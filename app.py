#!/usr/bin/env python3
"""Local, offline Streamlit UI for the ONT reference-mapping pipeline."""

from __future__ import annotations

import os
import shutil
import base64
import html
import json
from pathlib import Path

import pandas as pd
import streamlit as st

from ont_ui.batch import (
    BatchExecutionError,
    BatchPreparationError,
    BatchSettings,
    collect_batch_results,
    natural_key,
    parse_uploaded_reference,
    prepare_batch_job,
    run_batch_job,
    uploaded_barcodes,
)
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


def _load_branding() -> dict[str, str]:
    branding = {
        "organization": "PLASMID SEQUENCING",
        "title": "ONT Plasmid Analyzer",
        "subtitle": "Reference-to-barcode batch mapping and variant review",
        "badge": "LOCAL RESEARCH TOOL",
        "distributed_by": "Jongin Baek",
        "primary_color": "#2446C8",
        "secondary_color": "#23398C",
    }
    config_path = REPOSITORY_ROOT / "branding.local.json"
    if config_path.is_file():
        try:
            values = json.loads(config_path.read_text(encoding="utf-8"))
            if isinstance(values, dict):
                for key in branding:
                    if isinstance(values.get(key), str) and values[key].strip():
                        branding[key] = values[key].strip()
        except (OSError, json.JSONDecodeError):
            pass
    return branding


def _brand_identity(branding: dict[str, str]) -> str:
    for filename, mime in (
        ("branding_logo.svg", "image/svg+xml"),
        ("branding_logo.png", "image/png"),
    ):
        path = REPOSITORY_ROOT / filename
        if path.is_file():
            encoded = base64.b64encode(path.read_bytes()).decode("ascii")
            return (
                f'<img src="data:{mime};base64,{encoded}" '
                'alt="Organization logo" class="brand-logo-image">'
            )
    return f'<div class="brand-wordmark">{html.escape(branding["organization"])}</div>'


def _sidebar_brand_logo() -> str:
    for filename, mime in (
        ("branding_sidebar_logo.svg", "image/svg+xml"),
        ("branding_sidebar_logo.png", "image/png"),
    ):
        path = REPOSITORY_ROOT / filename
        if path.is_file():
            encoded = base64.b64encode(path.read_bytes()).decode("ascii")
            return (
                f'<img src="data:{mime};base64,{encoded}" '
                'alt="Organization logo" class="sidebar-brand-logo">'
            )
    return ""


def _brand_header() -> None:
    branding = _load_branding()
    identity = _brand_identity(branding)
    markup = """
        <style>
        :root {
            --brand-primary: __PRIMARY__;
            --brand-secondary: __SECONDARY__;
            --ink: #17213A;
        }
        .stApp { background: linear-gradient(180deg, #F6F8FC 0, #FFFFFF 240px); }
        .brand-shell {
            display:flex; align-items:center; justify-content:space-between; gap:24px;
            padding:22px 28px; border-radius:18px; color:white; margin-bottom:20px;
            background:linear-gradient(120deg, var(--brand-primary) 0%, var(--brand-secondary) 62%, #5268D3 100%);
            box-shadow:0 12px 34px rgba(20,40,160,.20);
        }
        .brand-wordmark { font:700 19px/1 Arial,sans-serif; letter-spacing:.02em; }
        .brand-logo-image { display:block; max-width:240px; max-height:52px; object-fit:contain; }
        .brand-title { font:700 30px/1.15 "Segoe UI",Arial,sans-serif; margin-top:12px; }
        .brand-subtitle { opacity:.85; margin-top:6px; font-size:14px; }
        .brand-pill {
            border:1px solid rgba(255,255,255,.45); border-radius:999px;
            padding:8px 13px; white-space:nowrap; font-size:12px; font-weight:600;
        }
        [data-testid="stMetric"] {
            background:#FFFFFF; border:1px solid #E3E8F3; border-radius:13px; padding:12px;
        }
        div[data-testid="stFileUploader"] section {
            background:#FBFCFF; border:1.5px dashed #8A98D6; border-radius:14px;
        }
        .status-clean { color:#087A55; font-weight:700; }
        .status-variant { color:var(--brand-primary); font-weight:700; }
        .status-review { color:#B05D00; font-weight:700; }
        .status-error { color:#B42318; font-weight:700; }
        </style>
        <div class="brand-shell">
          <div>
            __IDENTITY__
            <div class="brand-title">__TITLE__</div>
            <div class="brand-subtitle">__SUBTITLE__</div>
          </div>
          <div class="brand-pill">__BADGE__</div>
        </div>
        """
    replacements = {
        "__PRIMARY__": html.escape(branding["primary_color"]),
        "__SECONDARY__": html.escape(branding["secondary_color"]),
        "__IDENTITY__": identity,
        "__TITLE__": html.escape(branding["title"]),
        "__SUBTITLE__": html.escape(branding["subtitle"]),
        "__BADGE__": html.escape(branding["badge"]),
    }
    for token, value in replacements.items():
        markup = markup.replace(token, value)
    st.markdown(
        markup,
        unsafe_allow_html=True,
    )


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


def _quick_compare_tab(settings: dict[str, object]) -> None:
    st.subheader("Quick reference–sequence comparison")
    st.write(
        "Paste or upload one reference and one query/consensus sequence. The app automatically "
        "detects orientation and reports exact SNPs, insertions, and deletions."
    )
    left, right = st.columns(2)
    with left:
        st.markdown("#### Reference")
        quick_ref_file = st.file_uploader(
            "Reference FASTA or text",
            type=["fasta", "fa", "fna", "txt"],
            key="quick_ref_file",
            help="Upload one reference plasmid sequence in FASTA or plain-text format.",
        )
        quick_ref_text = st.text_area(
            "Or paste reference DNA",
            height=220,
            placeholder=">reference\nACGT...",
            key="quick_ref_text",
            help="Use this field instead of uploading a reference file.",
        )
    with right:
        st.markdown("#### Query / consensus")
        quick_query_file = st.file_uploader(
            "Query FASTA or text",
            type=["fasta", "fa", "fna", "txt"],
            key="quick_query_file",
            help="Upload the assembled or consensus sequence to compare with the reference.",
        )
        quick_query_text = st.text_area(
            "Or paste query DNA",
            height=220,
            placeholder=">query\nACGT...",
            key="quick_query_text",
            help="Use this field instead of uploading a query file.",
        )

    circular = bool(settings["circular"])
    if st.button(
        "Compare sequences",
        type="primary",
        key="quick_run",
        help="Align the query against the reference and report SNPs, insertions, and deletions.",
    ):
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


def _default_batch_mapping(reference_files, barcode_names: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for index, uploaded in enumerate(reference_files):
        assigned = barcode_names[index * 3 : index * 3 + 3]
        rows.append(
            {
                "Order": index + 1,
                "Reference": Path(uploaded.name).name,
                "Barcode 1": assigned[0] if len(assigned) > 0 else "",
                "Barcode 2": assigned[1] if len(assigned) > 1 else "",
                "Barcode 3": assigned[2] if len(assigned) > 2 else "",
            }
        )
    return pd.DataFrame(rows)


def _result_status_html(status: str) -> str:
    css = {
        "CLEAN": "status-clean",
        "VARIANT": "status-variant",
        "REVIEW": "status-review",
        "ERROR": "status-error",
    }.get(status, "status-review")
    label = {
        "CLEAN": "● CLEAN",
        "VARIANT": "● VARIANT DETECTED",
        "REVIEW": "● REVIEW",
        "ERROR": "● ERROR",
    }.get(status, status)
    return f'<span class="{css}">{label}</span>'


def _render_batch_results(result: dict[str, object], job) -> None:
    st.divider()
    st.subheader("Batch analysis results")
    samples = result.get("samples", [])
    assert isinstance(samples, list)
    counts = {
        status: sum(sample.get("status") == status for sample in samples)
        for status in ("CLEAN", "VARIANT", "REVIEW", "ERROR")
    }
    metrics = st.columns(5)
    metrics[0].metric("Total samples", len(samples))
    metrics[1].metric("Clean", counts["CLEAN"])
    metrics[2].metric("Variant detected", counts["VARIANT"])
    metrics[3].metric("Review", counts["REVIEW"])
    metrics[4].metric("Error", counts["ERROR"])

    summary_rows = []
    for sample in samples:
        summary_rows.append(
            {
                "Reference": sample.get("reference_name"),
                "Replicate": sample.get("replicate"),
                "Barcode": sample.get("barcode"),
                "Status": sample.get("status"),
                "Mapping rate (%)": (
                    float(sample["mapping_rate"]) * 100
                    if sample.get("mapping_rate") is not None
                    else None
                ),
                "Mean depth": sample.get("mean_depth"),
                "Coverage ≥1× (%)": (
                    float(sample["coverage_1x"]) * 100
                    if sample.get("coverage_1x") is not None
                    else None
                ),
                "SNP": sample.get("snp", 0),
                "Insertion": sample.get("insertion", 0),
                "Deletion": sample.get("deletion", 0),
            }
        )
    if summary_rows:
        st.dataframe(
            pd.DataFrame(summary_rows),
            use_container_width=True,
            hide_index=True,
            column_config={
                "Mapping rate (%)": st.column_config.NumberColumn(format="%.2f%%"),
                "Coverage ≥1× (%)": st.column_config.NumberColumn(format="%.2f%%"),
                "Mean depth": st.column_config.NumberColumn(format="%.1f×"),
            },
        )
    st.download_button(
        "Download complete summary (CSV)",
        data=str(result.get("summary_csv", "")).encode("utf-8-sig"),
        file_name=f"{job.experiment_name}_ONT_summary.csv",
        mime="text/csv",
        type="primary",
    )
    st.caption(f"Complete BAM, VCF, consensus FASTA, and reports: `{result.get('result_root', '')}`")

    groups = result.get("groups", {})
    assert isinstance(groups, dict)
    st.markdown("#### Reference-by-reference review")
    for reference_name, group_samples in groups.items():
        assert isinstance(group_samples, list)
        labels = ", ".join(
            f"{sample.get('barcode')}: {sample.get('status')}" for sample in group_samples
        )
        with st.expander(f"{reference_name}  ·  {labels}"):
            columns = st.columns(max(1, len(group_samples)))
            for column, sample in zip(columns, group_samples):
                with column:
                    st.markdown(_result_status_html(str(sample.get("status"))), unsafe_allow_html=True)
                    st.write(f"**{sample.get('barcode')} · Replicate {sample.get('replicate')}**")
                    if sample.get("status") == "ERROR":
                        st.error(str(sample.get("message", "Analysis failed.")))
                        continue
                    mapping_rate = sample.get("mapping_rate")
                    st.metric("Mapping", f"{float(mapping_rate):.1%}" if mapping_rate is not None else "N/A")
                    st.metric("Mean depth", f"{float(sample.get('mean_depth', 0)):.1f}×")
                    st.write(
                        f"SNP **{sample.get('snp', 0)}** · INS **{sample.get('insertion', 0)}** · "
                        f"DEL **{sample.get('deletion', 0)}**"
                    )
                    details = sample.get("variant_details", [])
                    if details:
                        st.dataframe(pd.DataFrame(details), hide_index=True, use_container_width=True)


def _batch_ont_tab(settings: dict[str, object]) -> None:
    reference_count = int(settings["reference_count"])
    expected_sample_count = reference_count * 3
    st.subheader(f"{reference_count}-plasmid / {expected_sample_count}-barcode batch analysis")
    st.write(
        "Upload reference sequences and the completed ONT barcode directory. Barcodes are sorted "
        "numerically and assigned three at a time; review or edit the mapping before running."
    )
    st.info(
        "Recommended workflow: select the run folder that contains barcode13, barcode14, … folders. "
        "The folder hierarchy is used to identify each barcode automatically."
    )

    upload_columns = st.columns(2)
    with upload_columns[0]:
        st.markdown("#### 1 · Reference sequences")
        reference_uploads = st.file_uploader(
            f"Drag {reference_count} reference file(s)",
            type=["fasta", "fa", "fna", "txt"],
            accept_multiple_files=True,
            key="batch_references",
            help=(
                f"Upload exactly {reference_count} reference file(s). The file name becomes the "
                "plasmid/reference name in the final report."
            ),
        )
    with upload_columns[1]:
        st.markdown("#### 2 · ONT barcode directory")
        read_uploads = st.file_uploader(
            "Select or drag the directory containing barcode folders",
            type=["fastq", "fq", "gz"],
            accept_multiple_files="directory",
            key="batch_reads",
            help="Select a parent folder that contains barcode01, barcode02, ... subfolders.",
        )

    reference_uploads = list(reference_uploads or [])
    read_uploads = list(read_uploads or [])
    grouped_reads = uploaded_barcodes(read_uploads)
    barcode_names = list(grouped_reads)
    reference_signature = tuple((item.name, item.size) for item in reference_uploads)
    read_signature = tuple((item.name, item.size) for item in read_uploads)
    batch_signature = (reference_signature, read_signature)

    validation_errors: list[str] = []
    reference_rows: list[dict[str, object]] = []
    for item in reference_uploads:
        try:
            record = parse_uploaded_reference(item)
            reference_rows.append(
                {"Reference file": Path(item.name).name, "Length (bp)": len(record.sequence)}
            )
        except (BatchPreparationError, SequenceValidationError) as exc:
            validation_errors.append(str(exc))

    upload_metrics = st.columns(4)
    upload_metrics[0].metric("References", f"{len(reference_uploads)} / {reference_count}")
    upload_metrics[1].metric("Detected barcodes", f"{len(barcode_names)} / {expected_sample_count}")
    upload_metrics[2].metric("Planned samples", expected_sample_count)
    upload_metrics[3].metric("Uploaded FASTQ files", len(read_uploads))

    if reference_uploads and len(reference_uploads) != reference_count:
        st.warning(
            f"The sidebar is set to {reference_count} reference(s), but "
            f"{len(reference_uploads)} file(s) were uploaded."
        )

    if reference_rows:
        with st.expander("Reference file check"):
            st.dataframe(pd.DataFrame(reference_rows), hide_index=True, use_container_width=True)
    for error in validation_errors:
        st.error(error)
    if read_uploads and not grouped_reads:
        st.error(
            "No barcode folder name was detected. Select the parent directory that contains "
            "barcode01/barcode02/... folders rather than selecting FASTQ files individually."
        )

    if st.session_state.get("batch_signature") != batch_signature:
        st.session_state["batch_signature"] = batch_signature
        st.session_state["batch_mapping"] = _default_batch_mapping(reference_uploads, barcode_names)

    if reference_uploads and barcode_names:
        st.markdown("#### 3 · Review reference ↔ barcode mapping")
        st.caption(
            "References are shown in upload order and barcodes are assigned numerically in groups of three. "
            "Change Order or any barcode cell if the experimental order is different."
        )
        mapping_frame = st.session_state.get("batch_mapping")
        if not isinstance(mapping_frame, pd.DataFrame):
            mapping_frame = _default_batch_mapping(reference_uploads, barcode_names)
        edited = st.data_editor(
            mapping_frame,
            hide_index=True,
            use_container_width=True,
            num_rows="fixed",
            key="batch_mapping_editor",
            column_config={
                "Order": st.column_config.NumberColumn(
                    min_value=1,
                    step=1,
                    required=True,
                    help="Analysis and report order for the reference.",
                ),
                "Reference": st.column_config.TextColumn(
                    disabled=True, help="Reference file name used in the result report."
                ),
                "Barcode 1": st.column_config.SelectboxColumn(
                    options=barcode_names, required=True, help="First ONT replicate for this plasmid."
                ),
                "Barcode 2": st.column_config.SelectboxColumn(
                    options=barcode_names, required=True, help="Second ONT replicate for this plasmid."
                ),
                "Barcode 3": st.column_config.SelectboxColumn(
                    options=barcode_names, required=True, help="Third ONT replicate for this plasmid."
                ),
            },
        )
        st.session_state["batch_mapping"] = edited

        sorted_mapping = edited.sort_values("Order", kind="stable")
        assigned = [
            str(value)
            for column in ("Barcode 1", "Barcode 2", "Barcode 3")
            for value in sorted_mapping[column].tolist()
            if str(value) not in {"", "nan", "None"}
        ]
        invalid_assignments = sorted(
            {value for value in assigned if value not in barcode_names}, key=natural_key
        )
        duplicates = sorted({value for value in assigned if assigned.count(value) > 1}, key=natural_key)
        missing = sorted(set(barcode_names) - set(assigned), key=natural_key)
        expected_assignment_count = expected_sample_count
        if invalid_assignments or len(assigned) != expected_assignment_count:
            st.error(
                f"Every reference needs three valid barcodes ({expected_assignment_count} total assignments)."
            )
        if duplicates:
            st.error("Duplicate assignments: " + ", ".join(duplicates))
        if missing:
            st.warning("Uploaded but not assigned: " + ", ".join(missing))
        if len(barcode_names) != expected_sample_count:
            st.warning(
                f"The selected 3-per-reference layout expects {expected_sample_count} barcodes, "
                f"but {len(barcode_names)} were detected. You can still edit and run the mapping."
            )

        can_run = (
            not validation_errors
            and not duplicates
            and not invalid_assignments
            and len(reference_uploads) == reference_count
            and len(assigned) == expected_assignment_count
        )
        if st.button(
            f"Run batch analysis ({len(assigned)} samples)",
            type="primary",
            disabled=not can_run,
            key="batch_run",
            use_container_width=True,
            help="Run all mapped barcode samples with the Analysis settings selected in the sidebar.",
        ):
            try:
                mappings = []
                for _, row in sorted_mapping.iterrows():
                    mappings.append(
                        {
                            "reference": str(row["Reference"]),
                            "barcodes": [
                                str(row["Barcode 1"]),
                                str(row["Barcode 2"]),
                                str(row["Barcode 3"]),
                            ],
                        }
                    )
                missing_tools = [
                    executable
                    for executable in ("bash", "minimap2", "samtools", "bcftools")
                    if not shutil.which(executable)
                ]
                if missing_tools:
                    raise BatchPreparationError(
                        "Missing analysis tool(s): " + ", ".join(missing_tools)
                    )
                batch_settings = BatchSettings(
                    experiment_name=str(settings["experiment_name"]),
                    threads=int(settings["threads"]),
                    parallel_jobs=int(settings["parallel_jobs"]),
                    min_read_length=int(settings["min_length"]),
                    min_read_quality=int(settings["min_quality"]),
                )
                with st.spinner("Staging uploaded references and FASTQ files on the server…"):
                    job = prepare_batch_job(
                        UI_RUN_ROOT,
                        reference_uploads,
                        read_uploads,
                        mappings,
                        batch_settings,
                    )
                log_box = st.empty()
                progress_box = st.empty()
                log_lines: list[str] = []
                completed = 0

                def show_line(line: str) -> None:
                    nonlocal completed
                    log_lines.append(line)
                    log_box.code("\n".join(log_lines[-60:]), language="text")
                    if " Done: " in line:
                        completed += 1
                    progress_box.progress(
                        min(1.0, completed / max(1, job.sample_count)),
                        text=f"Completed {completed} / {job.sample_count} samples",
                    )

                with st.spinner("Batch analysis is running. Keep this browser tab open…"):
                    run_batch_job(REPOSITORY_ROOT, job, on_line=show_line)
                    result = collect_batch_results(job)
                progress_box.success(f"Completed {job.sample_count} samples.")
                st.session_state["batch_result"] = {"job": job, "result": result}
            except (
                BatchPreparationError,
                BatchExecutionError,
                SequenceValidationError,
                OSError,
            ) as exc:
                st.error(str(exc))

    batch_state = st.session_state.get("batch_result")
    if batch_state:
        _render_batch_results(batch_state["result"], batch_state["job"])


def _raw_ont_tab(settings: dict[str, object]) -> None:
    st.subheader("Raw Oxford Nanopore read analysis")
    st.write(
        "Run the existing NanoFilt → minimap2 → samtools → bcftools pipeline from the browser. "
        "For large runs, selecting a server folder avoids uploading or duplicating FASTQ files."
    )

    identity_cols = st.columns(2)
    with identity_cols[0]:
        experiment_name = st.text_input(
            "Experiment name",
            value="ONT_experiment",
            help="Name of the output experiment folder created under ui_runs.",
        )
    with identity_cols[1]:
        sample_name = st.text_input(
            "Sample / vector name",
            value="sample01",
            help="Short identifier used as the output file prefix.",
        )

    reference_file = st.file_uploader(
        "Reference FASTA or text",
        type=["fasta", "fa", "fna", "txt"],
        key="raw_ref_file",
        help="Upload one reference plasmid sequence in FASTA or plain-text format.",
    )
    reference_text = st.text_area(
        "Or paste reference DNA",
        height=150,
        placeholder=">sample01\nACGT...",
        key="raw_ref_text",
        help="Use this field instead of uploading a reference file.",
    )

    reads_mode = st.radio(
        "FASTQ input",
        ["Server folder or file path (recommended)", "Upload FASTQ files"],
        horizontal=True,
        help="Server paths avoid transferring large read files through the browser.",
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
            help="Upload one or more FASTQ/FASTQ.GZ files for this sample.",
        )
        st.caption("Browser upload is intended for smaller datasets; large ONT runs should use a server path.")

    threads = int(settings["threads"])
    min_length = int(settings["min_length"])
    min_read_quality = int(settings["min_read_quality"])
    caller = str(settings["caller"])
    pilon_jar = str(settings.get("pilon_jar", ""))
    pilon_mem = str(settings.get("pilon_mem", "16G"))
    min_variant_quality = float(settings["min_variant_quality"])
    min_variant_depth = int(settings["min_variant_depth"])
    min_af = float(settings["min_af"])
    circular = bool(settings["circular"])
    edge_margin = int(settings["edge_margin"])

    if st.button(
        "Run ONT analysis",
        type="primary",
        key="raw_run",
        help="Run the pipeline using the Analysis settings selected in the sidebar.",
    ):
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


def _sidebar() -> tuple[str, dict[str, object]]:
    branding = _load_branding()
    logo = _sidebar_brand_logo()
    with st.sidebar:
        st.markdown(
            """
            <style>
            .sidebar-brand-footer {
                width:100%; box-sizing:border-box; padding:12px 14px 10px;
                margin-top:24px;
                background:rgba(255,255,255,.96); border:1px solid #E3E8F3;
                border-radius:12px; box-shadow:0 6px 20px rgba(23,33,58,.08);
            }
            .sidebar-brand-logo {
                display:block; width:100%; max-width:210px; max-height:58px;
                margin:0 auto 8px; object-fit:contain;
            }
            .sidebar-brand-credit {
                color:#667085; font:500 11px/1.35 "Segoe UI",Arial,sans-serif;
                text-align:center;
            }
            </style>
            """,
            unsafe_allow_html=True,
        )
        st.markdown("### ONT Plasmid Analyzer")
        st.caption("Local ONT plasmid variant analysis")

        mode = st.radio(
            "Analysis menu",
            ["Batch plasmid analysis", "Quick sequence comparison", "Single-sample ONT analysis"],
            key="analysis_mode",
            help=(
                "Batch maps three barcode replicates to each plasmid. Quick comparison aligns two "
                "assembled sequences. Single-sample runs the complete ONT read pipeline."
            ),
        )

        settings: dict[str, object] = {}
        if mode == "Batch plasmid analysis":
            st.markdown("#### Run setup")
            reference_count = st.number_input(
                "Number of references",
                min_value=1,
                max_value=32,
                value=1,
                step=1,
                key="batch_reference_count",
                help=(
                    "Select how many different plasmids will be analyzed in this run. "
                    "Each reference is matched to three barcode samples."
                ),
            )
            st.caption(f"{int(reference_count)} references × 3 replicates = {int(reference_count) * 3} samples")
            st.markdown("#### Analysis settings")
            with st.expander("Batch pipeline", expanded=True):
                settings["experiment_name"] = st.text_input(
                    "Experiment name",
                    value="ONT_plasmid_batch",
                    key="batch_experiment_name",
                    help="Name of the output folder created under ui_runs.",
                )
                settings["threads"] = st.number_input(
                    "Threads per sample",
                    min_value=1,
                    max_value=64,
                    value=8,
                    key="batch_threads",
                    help="CPU threads assigned to each barcode sample. Eight is a safe default.",
                )
                settings["parallel_jobs"] = st.number_input(
                    "Parallel samples",
                    min_value=1,
                    max_value=96,
                    value=16,
                    key="batch_parallel_jobs",
                    help=(
                        "Maximum number of barcode samples processed at the same time. Reduce this "
                        "on a laptop or a low-memory system."
                    ),
                )
                settings["min_length"] = st.number_input(
                    "Minimum read length",
                    min_value=0,
                    value=500,
                    step=50,
                    key="batch_min_length",
                    help="Reads shorter than this length are removed before mapping. Use 300 if coverage is low.",
                )
                settings["min_quality"] = st.number_input(
                    "Minimum mean Q",
                    min_value=0,
                    value=10,
                    step=1,
                    key="batch_min_quality",
                    help="Reads below this mean Phred quality score are removed. Use Q8 if coverage is low.",
                )
            settings["reference_count"] = int(reference_count)
        elif mode == "Quick sequence comparison":
            st.markdown("#### Analysis settings")
            with st.expander("Sequence alignment", expanded=True):
                settings["circular"] = st.checkbox(
                    "Circular plasmid/vector",
                    value=True,
                    key="quick_circular",
                    help=(
                        "Temporarily doubles the reference during alignment so variants spanning the "
                        "plasmid origin can be detected, then restores the original coordinates."
                    ),
                )
        else:
            st.markdown("#### Analysis settings")
            with st.expander("Read filtering and caller", expanded=True):
                settings["threads"] = st.number_input(
                    "CPU threads",
                    min_value=1,
                    max_value=max(1, os.cpu_count() or 1),
                    value=min(8, os.cpu_count() or 1),
                    key="raw_threads",
                    help="CPU threads used for this single sample.",
                )
                settings["min_length"] = st.number_input(
                    "Minimum read length",
                    min_value=0,
                    value=500,
                    step=50,
                    key="raw_min_length",
                    help="Reads shorter than this value are removed. Use 300 if coverage is low.",
                )
                settings["min_read_quality"] = st.number_input(
                    "Minimum mean Q",
                    min_value=0,
                    value=10,
                    step=1,
                    key="raw_min_quality",
                    help="Reads below this mean Phred score are removed. Use Q8 if coverage is low.",
                )
                caller_label = st.selectbox(
                    "Variant caller",
                    ["bcftools (recommended)", "medaka", "pilon"],
                    key="raw_caller",
                    help="bcftools is the default lightweight caller. Medaka and Pilon require separate installations.",
                )
                settings["caller"] = caller_label.split()[0]
                if settings["caller"] == "pilon":
                    settings["pilon_jar"] = st.text_input(
                        "Pilon JAR path",
                        key="raw_pilon_jar",
                        help="Full Linux path to the installed pilon.jar file.",
                    )
                    settings["pilon_mem"] = st.text_input(
                        "Pilon Java memory",
                        value="16G",
                        key="raw_pilon_mem",
                        help="Maximum Java heap memory allocated to Pilon.",
                    )
            with st.expander("Variant review thresholds", expanded=False):
                settings["min_variant_quality"] = st.number_input(
                    "Minimum QUAL",
                    min_value=0.0,
                    value=20.0,
                    key="raw_variant_quality",
                    help="Calls below this variant quality score are marked REVIEW.",
                )
                settings["min_variant_depth"] = st.number_input(
                    "Minimum DP",
                    min_value=0,
                    value=10,
                    key="raw_variant_depth",
                    help="Calls supported by fewer reads than this depth are marked REVIEW.",
                )
                settings["min_af"] = st.number_input(
                    "Minimum allele fraction",
                    min_value=0.0,
                    max_value=1.0,
                    value=0.80,
                    step=0.05,
                    key="raw_variant_af",
                    help="Minimum fraction of reads supporting the alternative allele for PASS status.",
                )
                settings["circular"] = st.checkbox(
                    "Circular reference",
                    value=True,
                    key="raw_circular",
                    help="Enable for plasmids so origin-spanning alignments are handled as circular.",
                )
                settings["edge_margin"] = st.number_input(
                    "Edge margin (bp)",
                    min_value=0,
                    value=50,
                    disabled=not bool(settings["circular"]),
                    key="raw_edge_margin",
                    help="Calls this close to a linear reference edge are marked for review.",
                )

        st.divider()
        with st.expander("System status", expanded=False):
            tool_states = [
                f"{'✓' if shutil.which(executable) else '✕'} {executable}"
                for executable in ("minimap2", "samtools", "bcftools", "NanoFilt")
            ]
            st.caption(" · ".join(tool_states))
            st.caption(f"Results: `{UI_RUN_ROOT}`")
        st.markdown(
            (
                '<div class="sidebar-brand-footer">'
                f'{logo}<div class="sidebar-brand-credit">'
                f'Distributed by {html.escape(branding["distributed_by"])}</div></div>'
            ),
            unsafe_allow_html=True,
        )
    return mode, settings


def main() -> None:
    st.set_page_config(
        page_title="ONT Plasmid Analyzer",
        page_icon="🧬",
        layout="wide",
    )
    mode, settings = _sidebar()
    _brand_header()
    if mode == "Batch plasmid analysis":
        _batch_ont_tab(settings)
    elif mode == "Quick sequence comparison":
        _quick_compare_tab(settings)
    else:
        _raw_ont_tab(settings)


if __name__ == "__main__":
    main()
