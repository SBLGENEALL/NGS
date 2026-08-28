#!/usr/bin/env python3
"""Local, offline Streamlit UI for the ONT reference-mapping pipeline."""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import pandas as pd
import streamlit as st

from ont_ui.compare import SequenceComparisonError, compare_sequences
from ont_ui.models import AlignmentSummary, PipelineRunResult, VariantEvent
from ont_ui.pipeline import (
    PipelineExecutionError,
    PipelinePreparationError,
    RawAnalysisSettings,
    find_fastq_files,
    missing_dependencies,
    prepare_job,
    run_job,
)
from ont_ui.results import parse_flagstat, parse_vcf_variants, read_depth, variants_csv
from ont_ui.sequences import (
    SequenceRecord,
    SequenceValidationError,
    parse_single_sequence,
    sanitize_name,
)


REPOSITORY_ROOT = Path(__file__).resolve().parent
UI_RUN_ROOT = REPOSITORY_ROOT / "ui_runs"


def _uploaded_text(uploaded) -> str:
    if uploaded is None:
        return ""
    try:
        return uploaded.getvalue().decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise SequenceValidationError(
            "The sequence file is not UTF-8 text. Export it as FASTA or plain text."
        ) from exc


def _sequence_from_inputs(uploaded, pasted: str, default_name: str) -> SequenceRecord:
    if uploaded is not None:
        default_name = Path(uploaded.name).name.split(".")[0] or default_name
        return parse_single_sequence(_uploaded_text(uploaded), default_name=default_name)
    return parse_single_sequence(pasted, default_name=default_name)


def _variant_counts(events: tuple[VariantEvent, ...] | list[VariantEvent]) -> dict[str, int]:
    return {
        kind: sum(event.kind == kind for event in events)
        for kind in ("SNP", "Insertion", "Deletion", "Complex")
    }


def _variants_frame(events: tuple[VariantEvent, ...] | list[VariantEvent]) -> pd.DataFrame:
    rows = [event.as_row() for event in events]
    columns = [
        "Position",
        "Type",
        "REF",
        "ALT",
        "Length (bp)",
        "QUAL",
        "Depth",
        "Allele fraction",
        "Status",
        "Warnings",
        "Reference context",
    ]
    return pd.DataFrame(rows, columns=columns)


def _render_quick_result(result: AlignmentSummary, reference: SequenceRecord) -> None:
    st.divider()
    st.subheader("Comparison result")
    count = _variant_counts(result.variants)
    metrics = st.columns(6)
    metrics[0].metric("Identity", f"{result.identity:.3%}")
    metrics[1].metric("Query coverage", f"{result.query_coverage:.2%}")
    metrics[2].metric("SNP", count["SNP"])
    metrics[3].metric("Insertion", count["Insertion"])
    metrics[4].metric("Deletion", count["Deletion"])
    metrics[5].metric("MAPQ", result.mapq)

    orientation = "Forward" if result.orientation == "Forward" else "Reverse complement"
    st.caption(
        f"Orientation: **{orientation}** · Reference start: **{result.reference_start:,}** · "
        f"Reference end: **{result.reference_end:,}** · Reference length: **{len(reference.sequence):,} bp**"
    )
    if result.query_coverage < 0.95:
        st.warning(
            "Less than 95% of the query aligned. Variants outside the aligned region cannot be reported."
        )

    if result.variants:
        st.dataframe(_variants_frame(result.variants), use_container_width=True, hide_index=True)
    else:
        st.success("No SNP, insertion, or deletion was found in the aligned sequence.")

    filename = f"{sanitize_name(reference.name)}_sequence_comparison.csv"
    st.download_button(
        "Download variant table (CSV)",
        data=variants_csv(result.variants).encode("utf-8-sig"),
        file_name=filename,
        mime="text/csv",
        use_container_width=False,
    )
    with st.expander("Alignment details"):
        st.code(result.raw_paf, language="text")
        st.caption(
            "SNP means a one-base substitution (point mutation). Insertions are positioned at the "
            "preceding reference base; deletions are positioned at the first deleted base."
        )


def _artifact_path(sample_dir: Path, sample_name: str, suffix: str) -> Path:
    return sample_dir / f"{sample_name}{suffix}"


def _download_file_button(label: str, path: Path, mime: str, key: str) -> None:
    if not path.is_file():
        return
    st.download_button(
        label,
        data=path.read_bytes(),
        file_name=path.name,
        mime=mime,
        key=key,
        use_container_width=True,
    )


def _render_raw_result(state: dict[str, object]) -> None:
    run = state["run"]
    reference = state["reference"]
    assert isinstance(run, PipelineRunResult)
    assert isinstance(reference, SequenceRecord)
    sample = run.sample_name
    sample_dir = run.sample_dir
    thresholds = state["thresholds"]
    assert isinstance(thresholds, dict)

    flagstat_path = _artifact_path(sample_dir, sample, ".flagstat.txt")
    depth_path = _artifact_path(sample_dir, sample, ".depth.txt")
    vcf_path = _artifact_path(sample_dir, sample, ".vcf.gz")
    mapping = parse_flagstat(flagstat_path)
    depth, chart_points = read_depth(depth_path)
    events: list[VariantEvent] = []
    if vcf_path.is_file():
        events = parse_vcf_variants(
            vcf_path,
            reference_sequence=reference.sequence,
            min_quality=float(thresholds["min_quality"]),
            min_depth=int(thresholds["min_depth"]),
            min_allele_fraction=float(thresholds["min_af"]),
            edge_margin=int(thresholds["edge_margin"]),
            circular=bool(thresholds["circular"]),
        )

    passed = sum(event.status == "PASS" for event in events)
    review = len(events) - passed
    count = _variant_counts(events)

    st.divider()
    st.subheader("ONT analysis result")
    first = st.columns(6)
    first[0].metric(
        "Mapping rate",
        f"{mapping.mapping_rate:.2%}" if mapping.mapping_rate is not None else "N/A",
    )
    first[1].metric("Mean depth", f"{depth.mean_depth:,.1f}×")
    first[2].metric("≥1× coverage", f"{depth.coverage_1x:.2%}")
    first[3].metric("PASS variants", passed)
    first[4].metric("Review variants", review)
    first[5].metric("Total reads", f"{mapping.total_reads:,}" if mapping.total_reads is not None else "N/A")

    second = st.columns(4)
    second[0].metric("SNP / point mutation", count["SNP"])
    second[1].metric("Insertion", count["Insertion"])
    second[2].metric("Deletion", count["Deletion"])
    second[3].metric("≥10× coverage", f"{depth.coverage_10x:.2%}")

    st.caption(
        "PASS requires the selected QUAL, depth, allele-fraction, and reference-edge thresholds. "
        "All calls remain visible; borderline calls are marked REVIEW instead of being discarded."
    )

    if chart_points:
        st.markdown("#### Depth across reference")
        depth_frame = pd.DataFrame(chart_points, columns=["Position", "Depth"]).set_index("Position")
        st.line_chart(depth_frame)

    st.markdown("#### Variant calls")
    if events:
        st.dataframe(_variants_frame(events), use_container_width=True, hide_index=True)
    elif vcf_path.is_file():
        st.success("No variant was called against the reference.")
    else:
        st.warning("No VCF was produced. Review the pipeline log and selected variant caller.")

    csv_data = variants_csv(events).encode("utf-8-sig")
    download_cols = st.columns(4)
    with download_cols[0]:
        st.download_button(
            "Variant table (CSV)",
            data=csv_data,
            file_name=f"{sample}_variants.csv",
            mime="text/csv",
            key=f"csv_{run.job_dir.name}",
            use_container_width=True,
        )
    with download_cols[1]:
        _download_file_button(
            "Original VCF.GZ", vcf_path, "application/gzip", f"vcf_{run.job_dir.name}"
        )
    with download_cols[2]:
        _download_file_button(
            "Consensus FASTA",
            _artifact_path(sample_dir, sample, ".consensus.fasta"),
            "text/plain",
            f"consensus_{run.job_dir.name}",
        )
    with download_cols[3]:
        _download_file_button(
            "Sample report",
            _artifact_path(sample_dir, sample, "_report.md"),
            "text/markdown",
            f"report_{run.job_dir.name}",
        )

    st.info(f"Complete result folder: `{sample_dir}`")
    st.caption("BAM/BAM.BAI and depth files remain in this folder to avoid loading large files into the browser.")
    with st.expander("Pipeline log"):
        st.code(run.log_path.read_text(encoding="utf-8", errors="replace"), language="text")


def _quick_compare_tab() -> None:
    st.subheader("Quick reference–sequence comparison")
    st.write(
        "Paste or upload one reference and one query/consensus sequence. The app automatically "
        "detects orientation and reports exact SNPs, insertions, and deletions."
    )
    left, right = st.columns(2)
    with left:
        st.markdown("#### Reference")
        quick_ref_file = st.file_uploader(
            "Reference FASTA or text", type=["fasta", "fa", "fna", "txt"], key="quick_ref_file"
        )
        quick_ref_text = st.text_area(
            "Or paste reference DNA",
            height=220,
            placeholder=">reference\nACGT...",
            key="quick_ref_text",
        )
    with right:
        st.markdown("#### Query / consensus")
        quick_query_file = st.file_uploader(
            "Query FASTA or text", type=["fasta", "fa", "fna", "txt"], key="quick_query_file"
        )
        quick_query_text = st.text_area(
            "Or paste query DNA",
            height=220,
            placeholder=">query\nACGT...",
            key="quick_query_text",
        )

    circular = st.checkbox(
        "Circular plasmid/vector",
        value=True,
        help="The reference is temporarily doubled for alignment, then coordinates are normalized back.",
        key="quick_circular",
    )
    if st.button("Compare sequences", type="primary", key="quick_run"):
        try:
            reference = _sequence_from_inputs(quick_ref_file, quick_ref_text, "reference")
            query = _sequence_from_inputs(quick_query_file, quick_query_text, "query")
            with st.spinner("Aligning sequences and extracting variants…"):
                result = compare_sequences(reference, query, circular=circular)
            st.session_state["quick_result"] = {
                "result": result,
                "reference": reference,
            }
        except (SequenceValidationError, SequenceComparisonError) as exc:
            st.error(str(exc))

    quick_state = st.session_state.get("quick_result")
    if quick_state:
        _render_quick_result(quick_state["result"], quick_state["reference"])


def _raw_ont_tab() -> None:
    st.subheader("Raw Oxford Nanopore read analysis")
    st.write(
        "Run the existing NanoFilt → minimap2 → samtools → bcftools pipeline from the browser. "
        "For large runs, selecting a server folder avoids uploading or duplicating FASTQ files."
    )

    identity_cols = st.columns(2)
    with identity_cols[0]:
        experiment_name = st.text_input("Experiment name", value="ONT_experiment")
    with identity_cols[1]:
        sample_name = st.text_input("Sample / vector name", value="sample01")

    reference_file = st.file_uploader(
        "Reference FASTA or text", type=["fasta", "fa", "fna", "txt"], key="raw_ref_file"
    )
    reference_text = st.text_area(
        "Or paste reference DNA",
        height=150,
        placeholder=">sample01\nACGT...",
        key="raw_ref_text",
    )

    reads_mode = st.radio(
        "FASTQ input",
        ["Server folder or file path (recommended)", "Upload FASTQ files"],
        horizontal=True,
    )
    reads_path = ""
    uploaded_reads = []
    if reads_mode.startswith("Server"):
        reads_path = st.text_input(
            "Path on this Linux server",
            placeholder="/data/user/MCET03/ONT_run/barcode01",
            help="Subfolders are searched recursively for FASTQ, FQ, and gzipped reads.",
        )
    else:
        uploaded_reads = st.file_uploader(
            "FASTQ / FASTQ.GZ",
            type=["fastq", "fq", "gz"],
            accept_multiple_files=True,
            key="raw_reads_upload",
        )
        st.caption("Browser upload is intended for smaller datasets; large ONT runs should use a server path.")

    with st.expander("Analysis settings", expanded=True):
        settings_cols = st.columns(4)
        threads = settings_cols[0].number_input(
            "CPU threads", min_value=1, max_value=max(1, os.cpu_count() or 1), value=min(8, os.cpu_count() or 1)
        )
        min_length = settings_cols[1].number_input(
            "Minimum read length", min_value=0, value=500, step=50, help="Use 300 if coverage is low."
        )
        min_read_quality = settings_cols[2].number_input(
            "Minimum mean Q", min_value=0, value=10, step=1, help="Use Q8 if coverage is low."
        )
        caller_label = settings_cols[3].selectbox(
            "Variant caller", ["bcftools (recommended)", "medaka", "pilon"]
        )
        caller = caller_label.split()[0]
        pilon_jar = ""
        pilon_mem = "16G"
        if caller == "pilon":
            pilon_cols = st.columns(2)
            pilon_jar = pilon_cols[0].text_input("Pilon JAR path")
            pilon_mem = pilon_cols[1].text_input("Pilon Java memory", value="16G")

        st.markdown("**Variant review thresholds**")
        threshold_cols = st.columns(5)
        min_variant_quality = threshold_cols[0].number_input("Minimum QUAL", min_value=0.0, value=20.0)
        min_variant_depth = threshold_cols[1].number_input("Minimum DP", min_value=0, value=10)
        min_af = threshold_cols[2].number_input(
            "Minimum allele fraction", min_value=0.0, max_value=1.0, value=0.80, step=0.05
        )
        circular = threshold_cols[3].checkbox("Circular reference", value=True)
        edge_margin = threshold_cols[4].number_input(
            "Edge margin (bp)", min_value=0, value=50, disabled=not circular
        )

    if st.button("Run ONT analysis", type="primary", key="raw_run"):
        try:
            reference = _sequence_from_inputs(reference_file, reference_text, sample_name or "sample")
            analysis_settings = RawAnalysisSettings(
                experiment_name=experiment_name,
                sample_name=sample_name or reference.name,
                threads=int(threads),
                min_read_length=int(min_length),
                min_read_quality=int(min_read_quality),
                variant_caller=caller,
                pilon_jar=pilon_jar,
                pilon_mem=pilon_mem,
            )
            missing = missing_dependencies(analysis_settings)
            if missing:
                raise PipelinePreparationError(
                    "Missing required program(s): " + ", ".join(missing) + ". Activate/update NGS_env first."
                )

            read_paths: list[Path] = []
            upload_streams = []
            if reads_mode.startswith("Server"):
                if not reads_path.strip():
                    raise PipelinePreparationError("Enter a server FASTQ file or folder path.")
                read_paths = find_fastq_files(Path(reads_path.strip()))
            else:
                if not uploaded_reads:
                    raise PipelinePreparationError("Upload at least one FASTQ/FASTQ.GZ file.")
                upload_streams = [(item.name, item) for item in uploaded_reads]

            job = prepare_job(
                UI_RUN_ROOT,
                reference,
                analysis_settings,
                read_paths=read_paths,
                uploaded_reads=upload_streams,
            )
            log_placeholder = st.empty()
            progress_placeholder = st.empty()
            log_lines: list[str] = []

            def show_line(line: str) -> None:
                log_lines.append(line)
                log_placeholder.code("\n".join(log_lines[-80:]), language="text")
                if "[" in line and "/6]" in line:
                    progress_placeholder.info(line)

            with st.spinner("ONT pipeline is running. Keep this browser tab open…"):
                run = run_job(REPOSITORY_ROOT, job, on_line=show_line)
            progress_placeholder.success("Pipeline completed.")
            st.session_state["raw_result"] = {
                "run": run,
                "reference": SequenceRecord(run.sample_name, reference.sequence),
                "thresholds": {
                    "min_quality": float(min_variant_quality),
                    "min_depth": int(min_variant_depth),
                    "min_af": float(min_af),
                    "circular": bool(circular),
                    "edge_margin": int(edge_margin) if circular else 0,
                },
            }
        except (
            SequenceValidationError,
            PipelinePreparationError,
            PipelineExecutionError,
        ) as exc:
            st.error(str(exc))

    raw_state = st.session_state.get("raw_result")
    if raw_state:
        _render_raw_result(raw_state)


def _sidebar() -> None:
    with st.sidebar:
        st.header("Environment")
        st.caption("All analysis runs locally. Sequence/read data is not sent outside this server.")
        for executable in ("minimap2", "samtools", "bcftools", "NanoFilt"):
            if shutil.which(executable):
                st.success(f"{executable}: ready")
            else:
                st.warning(f"{executable}: not found")
        st.divider()
        st.caption(f"UI result root\n`{UI_RUN_ROOT}`")


def main() -> None:
    st.set_page_config(
        page_title="ONT Variant Explorer",
        page_icon="🧬",
        layout="wide",
    )
    st.title("ONT Variant Explorer")
    st.caption("Offline plasmid/vector sequence comparison and Oxford Nanopore variant analysis")
    _sidebar()
    quick_tab, raw_tab = st.tabs(["Quick FASTA comparison", "Raw ONT analysis"])
    with quick_tab:
        _quick_compare_tab()
    with raw_tab:
        _raw_ont_tab()


if __name__ == "__main__":
    main()
