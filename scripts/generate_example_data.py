#!/usr/bin/env python3
"""Generate deterministic ONT Plasmid Analyzer example data."""

from __future__ import annotations

import csv
import random
import shutil
import zipfile
from pathlib import Path


REFERENCE_COUNT = 5
READS_PER_SAMPLE = 12
REFERENCE_LENGTH = 1400


def wrapped_fasta(name: str, sequence: str) -> str:
    lines = [sequence[index : index + 70] for index in range(0, len(sequence), 70)]
    return f">{name}\n" + "\n".join(lines) + "\n"


def alternate_base(base: str) -> str:
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[base]


def fastq(sample: str, sequence: str) -> str:
    quality = "I" * len(sequence)
    records: list[str] = []
    for read_index in range(1, READS_PER_SAMPLE + 1):
        records.extend(
            [f"@{sample}_read_{read_index:02d}", sequence, "+", quality]
        )
    return "\n".join(records) + "\n"


def build_dataset(root: Path) -> None:
    if root.exists():
        shutil.rmtree(root)
    references_dir = root / "references"
    samples_dir = root / "samples"
    references_dir.mkdir(parents=True)
    samples_dir.mkdir(parents=True)

    expected_rows: list[dict[str, object]] = []
    for reference_index in range(1, REFERENCE_COUNT + 1):
        rng = random.Random(20260902 + reference_index)
        reference = "".join(rng.choice("ACGT") for _ in range(REFERENCE_LENGTH))
        reference_name = f"demo_plasmid_{reference_index:02d}"
        (references_dir / f"{reference_name}.fasta").write_text(
            wrapped_fasta(reference_name, reference), encoding="utf-8"
        )

        prefix = f"P{reference_index:02d}"
        sample_variants: list[tuple[str, str, str, str]] = [
            (f"{prefix}_Control", reference, "None", "Reference match"),
        ]

        snp_position = 301
        snp_alt = alternate_base(reference[snp_position - 1])
        snp_sequence = (
            reference[: snp_position - 1] + snp_alt + reference[snp_position:]
        )
        sample_variants.append(
            (
                f"{prefix}_SNP",
                snp_sequence,
                "SNP",
                f"{snp_position}:{reference[snp_position - 1]}>{snp_alt}",
            )
        )

        indel_position = 701
        if reference_index % 2:
            indel_sample = f"{prefix}_Insertion"
            indel_sequence = reference[:indel_position] + "A" + reference[indel_position:]
            indel_type = "Insertion"
            indel_detail = f"A inserted after reference position {indel_position}"
        else:
            indel_sample = f"{prefix}_Deletion"
            deleted_base = reference[indel_position - 1]
            indel_sequence = reference[: indel_position - 1] + reference[indel_position:]
            indel_type = "Deletion"
            indel_detail = f"{deleted_base} deleted at reference position {indel_position}"
        sample_variants.append(
            (indel_sample, indel_sequence, indel_type, indel_detail)
        )

        for sample_name, sequence, variant_type, expected in sample_variants:
            sample_dir = samples_dir / reference_name / sample_name
            sample_dir.mkdir(parents=True)
            (sample_dir / "reads.fastq").write_text(
                fastq(sample_name, sequence), encoding="utf-8"
            )
            expected_rows.append(
                {
                    "Reference": f"{reference_name}.fasta",
                    "Sample ID": sample_name,
                    "Expected variant": variant_type,
                    "Expected detail": expected,
                }
            )

    with (root / "expected_variants.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(expected_rows[0]))
        writer.writeheader()
        writer.writerows(expected_rows)

    (root / "README.md").write_text(
        """# ONT Plasmid Analyzer example data

이 데이터는 UI와 pipeline 실행 확인을 위한 synthetic sequence입니다.

1. `references` 안의 FASTA 5개를 한 번에 업로드합니다.
2. 각 reference 영역에 같은 이름의 `samples/demo_plasmid_XX` 폴더를 드래그합니다.
3. Reference별로 Control, SNP, Insertion 또는 Deletion sample 3개가 감지되는지 확인합니다.
4. `Batch analysis 실행`을 누릅니다.
5. 결과를 `expected_variants.csv`와 비교합니다.

기본 Analysis settings에서 각 sample은 12 reads이므로 minimum depth 10 조건을 통과합니다.
실제 생물학적 데이터가 아닌 소프트웨어 테스트 전용 데이터입니다.
""",
        encoding="utf-8",
    )


def build_zip(dataset_root: Path, zip_path: Path) -> None:
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(dataset_root.rglob("*")):
            if path.is_file():
                archive.write(path, Path(dataset_root.name) / path.relative_to(dataset_root))


def main() -> None:
    repository_root = Path(__file__).resolve().parents[1]
    dataset_root = repository_root / "example_data" / "ONT_Plasmid_Analyzer_demo"
    build_dataset(dataset_root)
    build_zip(dataset_root, repository_root / "example_data" / "ONT_Plasmid_Analyzer_demo.zip")


if __name__ == "__main__":
    main()
