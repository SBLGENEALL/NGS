"""Batch upload, staging, execution, and summaries for plasmid ONT runs."""

from __future__ import annotations

import csv
import io
import json
import os
import re
import shutil
import subprocess
import uuid
import zipfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path, PurePosixPath
from typing import BinaryIO, Callable, Sequence

from .fastq_qc import load_fastq_qc
from .results import (
    parse_consensus_fallback_variants,
    parse_flagstat,
    parse_vcf_variants,
    read_depth,
)
from .sequences import SequenceRecord, parse_single_sequence, sanitize_name, write_fasta


BARCODE_RE = re.compile(r"^barcode0*([0-9]+)$", re.IGNORECASE)
FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz", ".gz")
REFERENCE_SUFFIXES = (".fasta", ".fa", ".fna", ".txt")
GENERIC_READ_NAMES = {"read", "reads", "fastq", "fastq_pass", "pass", "fail"}


class BatchPreparationError(ValueError):
    pass


class BatchExecutionError(RuntimeError):
    pass


def build_sample_results_zip(sample_dir: Path, destination: Path) -> Path:
    """Create an on-demand ZIP containing every regular sample result file."""
    sample_dir = sample_dir.expanduser().resolve()
    if not sample_dir.is_dir():
        raise BatchPreparationError("분석 결과 폴더를 찾을 수 없습니다.")
    files = sorted(
        path
        for path in sample_dir.rglob("*")
        if path.is_file() and not path.is_symlink()
    )
    if not files:
        raise BatchPreparationError("다운로드할 분석 결과 파일이 없습니다.")

    destination = destination.expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".tmp")
    try:
        with zipfile.ZipFile(temporary, "w") as archive:
            archive.writestr(
                "README.txt",
                "ONT Plasmid Analyzer complete sample results\n"
                "Includes merged/filtered reads, BAM/index, depth, QC, raw/final "
                "consensus, VCF/index, and report files when generated.\n",
                compress_type=zipfile.ZIP_DEFLATED,
            )
            for path in files:
                relative = Path(sample_dir.name) / path.relative_to(sample_dir)
                lower = path.name.casefold()
                compression = (
                    zipfile.ZIP_STORED
                    if lower.endswith((".gz", ".bam", ".bai", ".csi", ".zip"))
                    else zipfile.ZIP_DEFLATED
                )
                archive.write(path, relative.as_posix(), compress_type=compression)
        os.replace(temporary, destination)
    finally:
        if temporary.exists():
            temporary.unlink()
    return destination


class ServerPathUpload:
    """Small UploadedFile-compatible wrapper around a readable server file."""

    def __init__(self, source_path: Path, name: str | None = None):
        resolved = source_path.expanduser().resolve()
        if not resolved.is_file():
            raise BatchPreparationError(f"서버 파일을 찾을 수 없습니다: {resolved}")
        self.source_path = resolved
        self.name = name or resolved.name
        self.size = resolved.stat().st_size

    def getvalue(self) -> bytes:
        return self.source_path.read_bytes()


def _server_files(
    value: str,
    suffixes: tuple[str, ...],
    description: str,
) -> list[ServerPathUpload]:
    raw = (value or "").strip()
    if not raw:
        raise BatchPreparationError(f"{description} 경로를 입력하세요.")
    root = Path(raw).expanduser().resolve()
    if not root.exists():
        raise BatchPreparationError(f"서버 경로를 찾을 수 없습니다: {root}")
    if root.is_file():
        paths = [root] if root.name.lower().endswith(suffixes) else []
        names = [root.name]
    elif root.is_dir():
        paths = sorted(
            (
                path
                for path in root.rglob("*")
                if path.is_file() and path.name.lower().endswith(suffixes)
            ),
            key=lambda path: natural_key(path.relative_to(root).as_posix()),
        )
        names = [path.relative_to(root).as_posix() for path in paths]
    else:
        paths = []
        names = []
    if not paths:
        expected = ", ".join(suffixes)
        raise BatchPreparationError(
            f"{root}에서 {description} 파일을 찾지 못했습니다. 지원 형식: {expected}"
        )
    return [ServerPathUpload(path, name) for path, name in zip(paths, names)]


def server_reference_uploads(value: str) -> list[ServerPathUpload]:
    """Load reference files from one server file or a server directory."""
    return _server_files(value, REFERENCE_SUFFIXES, "Reference")


def server_fastq_uploads(value: str) -> list[ServerPathUpload]:
    """Load passing FASTQs from a standard ONT run folder when available."""
    raw = (value or "").strip()
    if raw:
        root = Path(raw).expanduser().resolve()
        if root.is_dir() and root.name.casefold() != "fastq_pass":
            try:
                pass_dir = next(
                    (
                        child
                        for child in root.iterdir()
                        if child.is_dir() and child.name.casefold() == "fastq_pass"
                    ),
                    None,
                )
            except OSError:
                pass_dir = None
            if pass_dir is not None:
                return _server_files(str(pass_dir), FASTQ_SUFFIXES, "ONT FASTQ")
    return _server_files(value, FASTQ_SUFFIXES, "ONT FASTQ")


@dataclass(frozen=True)
class BatchSettings:
    experiment_name: str
    threads: int = 8
    parallel_jobs: int = 16
    min_read_length: int = 500
    min_read_quality: int = 10
    min_variant_quality: float = 20.0
    min_variant_depth: int = 10
    min_allele_fraction: float = 0.80
    circular: bool = True
    edge_margin: int = 50


@dataclass(frozen=True)
class BatchJob:
    job_dir: Path
    config_path: Path
    manifest_path: Path
    experiment_name: str
    sample_count: int


def natural_key(value: str) -> tuple[object, ...]:
    return tuple(
        int(part) if part.isdigit() else part.casefold()
        for part in re.split(r"(\d+)", value)
    )


def normalize_barcode(value: str) -> str | None:
    match = BARCODE_RE.match((value or "").strip())
    if not match:
        return None
    return f"barcode{int(match.group(1)):02d}"


def barcode_from_upload_name(value: str) -> str | None:
    normalized = value.replace("\\", "/")
    for part in PurePosixPath(normalized).parts:
        barcode = normalize_barcode(part)
        if barcode:
            return barcode
    match = re.search(r"barcode0*([0-9]+)", Path(normalized).name, re.IGNORECASE)
    return f"barcode{int(match.group(1)):02d}" if match else None


def is_fastq_name(value: str) -> bool:
    return value.lower().endswith(FASTQ_SUFFIXES)


def _strip_fastq_suffix(value: str) -> str:
    lowered = value.lower()
    for suffix in sorted(FASTQ_SUFFIXES, key=len, reverse=True):
        if lowered.endswith(suffix):
            return value[: -len(suffix)]
    return value


def sample_name_from_upload_name(value: str) -> str:
    """Infer an ONT sample ID from a barcode, folder name, or FASTQ filename."""
    barcode = barcode_from_upload_name(value)
    if barcode:
        return barcode
    normalized = value.replace("\\", "/")
    parts = [part for part in PurePosixPath(normalized).parts if part not in {"", "."}]
    for part in reversed(parts[:-1]):
        if part.casefold() not in GENERIC_READ_NAMES:
            return sanitize_name(part, "sample")
    filename = parts[-1] if parts else "sample"
    return sanitize_name(_strip_fastq_suffix(filename), "sample")


def uploaded_samples(files: Sequence[object]) -> dict[str, list[object]]:
    """Group uploaded FASTQs by barcode, ONT alias folder, or filename.

    If the browser removes all folder paths and sends repeated generic names such
    as ``reads.fastq``, each upload is kept as a separate sample instead of being
    silently merged.
    """
    grouped: dict[str, list[object]] = {}
    fallback_number = 1
    for uploaded in files:
        name = str(getattr(uploaded, "name", ""))
        if not is_fastq_name(name):
            continue
        sample_id = sample_name_from_upload_name(name)
        normalized = name.replace("\\", "/")
        has_folder = len(PurePosixPath(normalized).parts) > 1
        if not has_folder and sample_id.casefold() in GENERIC_READ_NAMES:
            while True:
                candidate = f"sample_{fallback_number:02d}"
                fallback_number += 1
                if candidate not in grouped:
                    sample_id = candidate
                    break
        grouped.setdefault(sample_id, []).append(uploaded)
    return dict(sorted(grouped.items(), key=lambda item: natural_key(item[0])))


def uploaded_barcodes(files: Sequence[object]) -> dict[str, list[object]]:
    grouped: dict[str, list[object]] = {}
    for uploaded in files:
        name = str(getattr(uploaded, "name", ""))
        if not is_fastq_name(name):
            continue
        barcode = barcode_from_upload_name(name)
        if barcode:
            grouped.setdefault(barcode, []).append(uploaded)
    return dict(sorted(grouped.items(), key=lambda item: natural_key(item[0])))


def parse_uploaded_reference(uploaded: object) -> SequenceRecord:
    name = str(getattr(uploaded, "name", "reference.fasta"))
    try:
        raw = uploaded.getvalue()  # type: ignore[attr-defined]
        text = raw.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise BatchPreparationError(f"{name}: UTF-8 FASTA/plain DNA text is required.") from exc
    except AttributeError as exc:
        raise BatchPreparationError(f"{name}: could not read uploaded reference.") from exc
    return parse_single_sequence(text, default_name=Path(name).stem)


def _safe_upload_path(value: str, fallback: str) -> Path:
    parts = [part for part in PurePosixPath(value.replace("\\", "/")).parts if part not in {"", "."}]
    if not parts or any(part == ".." for part in parts):
        return Path(fallback)
    return Path(*(sanitize_name(part, fallback, max_length=150) for part in parts))


def _copy_upload(uploaded: object, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    source_path = getattr(uploaded, "source_path", None)
    if isinstance(source_path, Path) and source_path.is_file():
        destination.symlink_to(source_path.resolve())
        return
    if hasattr(uploaded, "seek"):
        uploaded.seek(0)  # type: ignore[attr-defined]
    with destination.open("wb") as output:
        shutil.copyfileobj(uploaded, output, length=1024 * 1024)  # type: ignore[arg-type]


def _write_config(path: Path, job: BatchJob, settings: BatchSettings) -> None:
    staged = job.job_dir / "staged"
    values: list[tuple[str, object]] = [
        ("references_dir", str((staged / "references").resolve())),
        ("data_dir", str((staged / "data").resolve())),
        ("results_dir", str((job.job_dir / "results").resolve())),
        ("minimap2_preset", "map-ont"),
        ("threads", max(1, min(64, settings.threads))),
        ("parallel_jobs", max(1, min(64, settings.parallel_jobs))),
        ("min_read_length", max(0, settings.min_read_length)),
        ("min_read_quality", max(0, settings.min_read_quality)),
        ("variant_caller", "bcftools"),
        ("pilon_jar", ""),
        ("pilon_mem", "16G"),
    ]
    with path.open("w", encoding="utf-8") as handle:
        handle.write("# Auto-generated by the ONT batch UI.\n")
        for key, value in values:
            rendered = json.dumps(value) if isinstance(value, str) else str(value)
            handle.write(f"{key}: {rendered}\n")


def prepare_batch_job(
    run_root: Path,
    references: Sequence[object],
    reads: Sequence[object],
    mappings: Sequence[dict[str, object]],
    settings: BatchSettings,
) -> BatchJob:
    if not mappings:
        raise BatchPreparationError("No reference-barcode mapping was provided.")
    experiment = sanitize_name(settings.experiment_name, "ONT_batch")
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    job_dir = run_root.resolve() / f"{timestamp}_{experiment}_{uuid.uuid4().hex[:8]}"
    reference_dir = job_dir / "staged" / "references" / experiment
    data_dir = job_dir / "staged" / "data" / experiment
    upload_dir = job_dir / "uploads" / "reads"
    results_dir = job_dir / "results"
    for directory in (reference_dir, data_dir, upload_dir, results_dir):
        directory.mkdir(parents=True, exist_ok=True)

    reference_by_file = {Path(str(getattr(item, "name", ""))).name: item for item in references}
    if len(reference_by_file) != len(references):
        raise BatchPreparationError(
            "Reference file names must be unique, even when they came from different folders."
        )
    legacy_grouped_reads = uploaded_barcodes(reads)
    used_legacy_barcodes: set[str] = set()
    used_samples: set[str] = set()
    manifest_samples: list[dict[str, object]] = []
    for row_index, row in enumerate(mappings, 1):
        reference_file = Path(str(row.get("reference", ""))).name
        uploaded_reference = reference_by_file.get(reference_file)
        if uploaded_reference is None:
            raise BatchPreparationError(f"Reference upload not found: {reference_file}")
        record = parse_uploaded_reference(uploaded_reference)
        raw_samples = row.get("samples")
        sample_entries: list[tuple[str, list[object]]] = []
        if isinstance(raw_samples, list) and raw_samples:
            for raw_sample in raw_samples:
                if not isinstance(raw_sample, dict):
                    raise BatchPreparationError(f"{reference_file}: invalid sample assignment.")
                sample_id = sanitize_name(str(raw_sample.get("name", "")), "sample")
                sample_files = raw_sample.get("files", [])
                if not isinstance(sample_files, list) or not sample_files:
                    raise BatchPreparationError(
                        f"{reference_file}: no FASTQ was uploaded for {sample_id}."
                    )
                sample_entries.append((sample_id, sample_files))
        else:
            raw_barcodes = row.get("barcodes", [])
            if not isinstance(raw_barcodes, list) or not raw_barcodes:
                raise BatchPreparationError(f"{reference_file}: assign at least one sample.")
            for raw_barcode in raw_barcodes:
                barcode = normalize_barcode(str(raw_barcode))
                if not barcode:
                    raise BatchPreparationError(
                        f"{reference_file}: invalid barcode '{raw_barcode}'."
                    )
                if barcode in used_legacy_barcodes:
                    raise BatchPreparationError(f"{barcode} is assigned more than once.")
                barcode_files = legacy_grouped_reads.get(barcode, [])
                if not barcode_files:
                    raise BatchPreparationError(f"No uploaded FASTQ was detected for {barcode}.")
                used_legacy_barcodes.add(barcode)
                sample_entries.append((barcode, barcode_files))

        seen_ids: set[str] = set()
        for replicate, (sample_id, sample_files) in enumerate(sample_entries, 1):
            if sample_id in seen_ids:
                raise BatchPreparationError(
                    f"{reference_file}: sample name '{sample_id}' is duplicated."
                )
            seen_ids.add(sample_id)
            base = sanitize_name(Path(reference_file).stem, f"reference_{row_index}", max_length=52)
            sample_name = sanitize_name(f"{base}__S{replicate}_{sample_id}", max_length=80)
            if sample_name in used_samples:
                raise BatchPreparationError(f"Duplicate output name: {sample_name}")
            used_samples.add(sample_name)
            write_fasta(reference_dir / f"{sample_name}.fasta", SequenceRecord(sample_name, record.sequence))
            sample_data_dir = data_dir / sample_name
            sample_data_dir.mkdir(parents=True, exist_ok=True)
            for index, uploaded in enumerate(sample_files, 1):
                original = str(getattr(uploaded, "name", f"reads_{index}.fastq.gz"))
                relative = _safe_upload_path(original, f"reads_{index}.fastq.gz")
                source = upload_dir / sample_name / f"{index:04d}_{relative.name}"
                _copy_upload(uploaded, source)
                destination = sample_data_dir / f"{index:04d}_{source.name}"
                destination.symlink_to(source.resolve())
            manifest_samples.append(
                {
                    "reference_file": reference_file,
                    "reference_name": Path(reference_file).stem,
                    "sample_id": sample_id,
                    "barcode": sample_id,
                    "replicate": replicate,
                    "sample_name": sample_name,
                    "reference_length": len(record.sequence),
                    "input_file_count": len(sample_files),
                }
            )

    config_path = job_dir / "batch_config.yaml"
    manifest_path = job_dir / "batch_manifest.json"
    provisional = BatchJob(job_dir, config_path, manifest_path, experiment, len(manifest_samples))
    _write_config(config_path, provisional, settings)
    manifest_path.write_text(
        json.dumps(
            {
                "experiment_name": experiment,
                "thresholds": {
                    "min_quality": settings.min_variant_quality,
                    "min_depth": settings.min_variant_depth,
                    "min_af": settings.min_allele_fraction,
                    "circular": settings.circular,
                    "edge_margin": settings.edge_margin if settings.circular else 0,
                },
                "pipeline_settings": {
                    "min_read_length": settings.min_read_length,
                    "min_read_quality": settings.min_read_quality,
                },
                "samples": manifest_samples,
            },
            ensure_ascii=False,
            indent=2,
        ),
        encoding="utf-8",
    )
    return provisional


def run_batch_job(
    repository_root: Path,
    job: BatchJob,
    on_line: Callable[[str], None] | None = None,
) -> Path:
    log_path = job.job_dir / "pipeline.log"
    command = ["bash", str(repository_root / "run_pipeline.sh"), str(job.config_path)]
    recent: list[str] = []
    with log_path.open("w", encoding="utf-8") as log_handle:
        process = subprocess.Popen(
            command,
            cwd=repository_root,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert process.stdout is not None
        for line in process.stdout:
            log_handle.write(line)
            log_handle.flush()
            clean = line.rstrip("\n")
            recent.append(clean)
            recent = recent[-40:]
            if on_line:
                on_line(clean)
        return_code = process.wait()
    if return_code != 0:
        raise BatchExecutionError(
            f"Batch pipeline exited with code {return_code}.\n" + "\n".join(recent)
        )
    return log_path


def collect_batch_results(job: BatchJob) -> dict[str, object]:
    manifest = json.loads(job.manifest_path.read_text(encoding="utf-8"))
    thresholds = manifest.get("thresholds", {})
    result_roots = sorted(
        [path for path in (job.job_dir / "results").glob(f"{job.experiment_name}_*") if path.is_dir()],
        key=lambda path: path.stat().st_mtime,
    )
    result_root = result_roots[-1] if result_roots else None
    staged_reference_dir = job.job_dir / "staged" / "references" / job.experiment_name
    samples: list[dict[str, object]] = []
    for item in manifest["samples"]:
        sample_name = item["sample_name"]
        sample_dir = result_root / sample_name if result_root else Path("/__missing__")
        flagstat = sample_dir / f"{sample_name}.flagstat.txt"
        depth_file = sample_dir / f"{sample_name}.depth.txt"
        vcf = sample_dir / f"{sample_name}.vcf.gz"
        fastq_qc = load_fastq_qc(sample_dir / f"{sample_name}.fastq_qc.json")
        merged_qc = fastq_qc.get("merged")
        filtered_qc = fastq_qc.get("filtered")
        if not isinstance(merged_qc, dict):
            merged_qc = {}
        if not isinstance(filtered_qc, dict):
            filtered_qc = {}
        effective_qc = filtered_qc or merged_qc
        merged_reads = int(merged_qc.get("reads", 0) or 0)
        filtered_reads = int(filtered_qc.get("reads", 0) or 0) if filtered_qc else None
        retained_fraction = (
            filtered_reads / merged_reads
            if filtered_reads is not None and merged_reads > 0
            else None
        )
        read_n50 = int(effective_qc.get("n50", 0) or 0) or None
        mean_read_quality = (
            float(effective_qc.get("mean_read_quality", 0) or 0)
            if effective_qc
            else None
        )
        q10_read_fraction = (
            float(effective_qc.get("q10_read_fraction", 0) or 0)
            if effective_qc
            else None
        )
        q20_read_fraction = (
            float(effective_qc.get("q20_read_fraction", 0) or 0)
            if effective_qc
            else None
        )
        q30_read_fraction = (
            float(effective_qc.get("q30_read_fraction", 0) or 0)
            if effective_qc
            else None
        )
        if not flagstat.is_file() or not depth_file.is_file():
            samples.append(
                {
                    **item,
                    "status": "ERROR",
                    "merged_reads": merged_reads or None,
                    "filtered_reads": filtered_reads,
                    "retained_fraction": retained_fraction,
                    "read_n50": read_n50,
                    "mean_read_quality": mean_read_quality,
                    "q10_read_fraction": q10_read_fraction,
                    "q20_read_fraction": q20_read_fraction,
                    "q30_read_fraction": q30_read_fraction,
                    "message": (
                        "Mapping 결과가 생성되지 않았습니다. Pipeline log에서 해당 "
                        "sample의 [SKIP] 또는 error 메시지를 확인하세요."
                    ),
                }
            )
            continue
        mapping = parse_flagstat(flagstat)
        depth, _ = read_depth(depth_file)
        reference_path = staged_reference_dir / f"{sample_name}.fasta"
        reference_text = reference_path.read_text(encoding="utf-8")
        reference = parse_single_sequence(reference_text, default_name=sample_name)
        events = (
            parse_vcf_variants(
                vcf,
                reference_sequence=reference.sequence,
                min_quality=float(thresholds.get("min_quality", 20.0)),
                min_depth=int(thresholds.get("min_depth", 10)),
                min_allele_fraction=float(thresholds.get("min_af", 0.8)),
                edge_margin=int(thresholds.get("edge_margin", 50)),
                circular=bool(thresholds.get("circular", True)),
            )
            if vcf.is_file()
            else []
        )
        variant_source = "bcftools"
        if not events:
            events = parse_consensus_fallback_variants(
                sample_dir / f"{sample_name}.samtools.consensus.fasta",
                reference,
                circular=bool(thresholds.get("circular", True)),
            )
            if events:
                variant_source = "samtools consensus review"
        passed = sum(event.status == "PASS" for event in events)
        review = len(events) - passed
        low_qc = (
            mapping.mapping_rate is None
            or mapping.mapping_rate < 0.8
            or depth.coverage_1x < 0.95
            or review > 0
        )
        status = "REVIEW" if low_qc else ("VARIANT" if passed else "CLEAN")
        samples.append(
            {
                **item,
                "status": status,
                "total_reads": mapping.total_reads,
                "mapping_rate": mapping.mapping_rate,
                "mean_depth": depth.mean_depth,
                "coverage_1x": depth.coverage_1x,
                "coverage_10x": depth.coverage_10x,
                "merged_reads": merged_reads or None,
                "filtered_reads": filtered_reads,
                "retained_fraction": retained_fraction,
                "read_n50": read_n50,
                "mean_read_quality": mean_read_quality,
                "q10_read_fraction": q10_read_fraction,
                "q20_read_fraction": q20_read_fraction,
                "q30_read_fraction": q30_read_fraction,
                "variants": len(events),
                "pass_variants": passed,
                "review_variants": review,
                "snp": sum(event.kind == "SNP" for event in events),
                "insertion": sum(event.kind == "Insertion" for event in events),
                "deletion": sum(event.kind == "Deletion" for event in events),
                "variant_details": [event.as_row() for event in events],
                "variant_source": variant_source,
            }
        )

    columns = [
        "Reference", "Sample #", "Sample ID", "Output name", "Status", "Total reads",
        "Merged reads", "Filtered reads", "Read retention", "Read N50", "Mean read Q",
        "Q10 reads", "Q20 reads", "Q30 reads", "Mapping rate", "Mean depth",
        "Coverage 1x", "Coverage 10x", "Variants",
        "PASS variants", "REVIEW variants", "SNP", "Insertion", "Deletion",
    ]
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=columns)
    writer.writeheader()
    for sample in samples:
        writer.writerow(
            {
                "Reference": sample.get("reference_name"),
                "Sample #": sample.get("replicate"),
                "Sample ID": sample.get("sample_id", sample.get("barcode")),
                "Output name": sample.get("sample_name"),
                "Status": sample.get("status"),
                "Total reads": sample.get("total_reads"),
                "Merged reads": sample.get("merged_reads"),
                "Filtered reads": sample.get("filtered_reads"),
                "Read retention": sample.get("retained_fraction"),
                "Read N50": sample.get("read_n50"),
                "Mean read Q": sample.get("mean_read_quality"),
                "Q10 reads": sample.get("q10_read_fraction"),
                "Q20 reads": sample.get("q20_read_fraction"),
                "Q30 reads": sample.get("q30_read_fraction"),
                "Mapping rate": sample.get("mapping_rate"),
                "Mean depth": sample.get("mean_depth"),
                "Coverage 1x": sample.get("coverage_1x"),
                "Coverage 10x": sample.get("coverage_10x"),
                "Variants": sample.get("variants"),
                "PASS variants": sample.get("pass_variants"),
                "REVIEW variants": sample.get("review_variants"),
                "SNP": sample.get("snp"),
                "Insertion": sample.get("insertion"),
                "Deletion": sample.get("deletion"),
            }
        )
    groups: dict[str, list[dict[str, object]]] = {}
    for sample in samples:
        groups.setdefault(str(sample["reference_name"]), []).append(sample)
    return {
        "samples": samples,
        "groups": groups,
        "summary_csv": output.getvalue(),
        "result_root": result_root,
        "pipeline_settings": manifest.get("pipeline_settings", {}),
    }
