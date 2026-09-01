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
    uploaded_samples,
)
from ont_ui.compare import SequenceComparisonError, compare_sequences
from ont_ui.demo import build_demo_batch_zip
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

DEFAULT_ANALYSIS_SETTINGS: dict[str, dict[str, object]] = {
    "Batch analysis": {
        "experiment_name": "ONT_plasmid_batch",
        "threads": 8,
        "parallel_jobs": 16,
        "min_length": 500,
        "min_quality": 10,
        "min_variant_quality": 20.0,
        "min_variant_depth": 10,
        "min_af": 0.80,
        "circular": True,
        "edge_margin": 50,
    },
    "Quick sequence comparison": {"circular": True},
    "Single-sample ONT analysis": {
        "threads": min(8, os.cpu_count() or 1),
        "min_length": 500,
        "min_read_quality": 10,
        "caller": "bcftools",
        "pilon_jar": "",
        "pilon_mem": "16G",
        "min_variant_quality": 20.0,
        "min_variant_depth": 10,
        "min_af": 0.80,
        "circular": True,
        "edge_margin": 50,
    },
}


def _analysis_settings(mode: str) -> dict[str, object]:
    all_settings = st.session_state.setdefault("saved_analysis_settings", {})
    if not isinstance(all_settings, dict):
        all_settings = {}
        st.session_state["saved_analysis_settings"] = all_settings
    current = all_settings.get(mode)
    if not isinstance(current, dict):
        current = dict(DEFAULT_ANALYSIS_SETTINGS[mode])
        all_settings[mode] = current
    return current


def _close_settings_page() -> None:
    st.session_state["show_analysis_settings"] = False


def _open_settings_page() -> None:
    st.session_state["show_analysis_settings"] = True


def _load_branding() -> dict[str, str]:
    branding = {
        "organization": "PLASMID SEQUENCING",
        "title": "ONT Plasmid Analyzer",
        "subtitle": "Reference와 ONT sample 매칭 및 variant 분석",
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
        html { color-scheme:light; }
        .stApp { background: linear-gradient(180deg, #F6F8FC 0, #FFFFFF 240px); }
        #MainMenu, [data-testid="stToolbar"] { visibility:hidden; }
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
    st.subheader("Quick sequence comparison")
    st.write(
        "완성된 두 DNA sequence를 직접 비교해 SNP, insertion, deletion을 확인합니다."
    )
    with st.expander("Query에는 무엇을 넣나요?", expanded=True):
        st.write(
            "**Query**에는 ONT read 분석 후 만들어진 **consensus FASTA**, assembly 결과, "
            "또는 확인하려는 plasmid의 완성 서열을 넣습니다. Raw FASTQ 파일을 넣는 곳은 "
            "아닙니다. Raw ONT read부터 분석하려면 `Batch analysis`를 사용하세요."
        )
        st.caption(
            "예: Reference = 설계한 pDNA FASTA · Query = 해당 sample에서 얻은 consensus FASTA"
        )
    left, right = st.columns(2)
    with left:
        st.markdown("#### Reference")
        quick_ref_file = st.file_uploader(
            "설계 Reference FASTA or text",
            type=["fasta", "fa", "fna", "txt"],
            key="quick_ref_file",
            help="비교 기준이 되는 원래 plasmid/vector 설계 서열입니다.",
        )
        quick_ref_text = st.text_area(
            "Or paste reference DNA",
            height=220,
            placeholder=">reference\nACGT...",
            key="quick_ref_text",
            help="파일 대신 reference sequence를 직접 붙여넣을 수 있습니다.",
        )
    with right:
        st.markdown("#### Query / consensus")
        quick_query_file = st.file_uploader(
            "Query / consensus FASTA or text",
            type=["fasta", "fa", "fna", "txt"],
            key="quick_query_file",
            help="ONT 분석으로 얻은 consensus 또는 assembled plasmid 전체 서열입니다. Raw FASTQ는 사용할 수 없습니다.",
        )
        quick_query_text = st.text_area(
            "Or paste query DNA",
            height=220,
            placeholder=">query\nACGT...",
            key="quick_query_text",
            help="파일 대신 query sequence를 직접 붙여넣을 수 있습니다.",
        )

    circular = bool(settings["circular"])
    if st.button(
        "Sequence 비교",
        type="primary",
        key="quick_run",
        help="Query를 reference에 alignment하여 SNP, insertion, deletion을 표시합니다.",
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
                "Sample #": sample.get("replicate"),
                "Sample ID": sample.get("sample_id", sample.get("barcode")),
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
            f"{sample.get('sample_id', sample.get('barcode'))}: {sample.get('status')}"
            for sample in group_samples
        )
        with st.expander(f"{reference_name}  ·  {labels}"):
            for start in range(0, len(group_samples), 3):
                sample_row = group_samples[start : start + 3]
                columns = st.columns(len(sample_row))
                for column, sample in zip(columns, sample_row):
                    with column:
                        st.markdown(
                            _result_status_html(str(sample.get("status"))),
                            unsafe_allow_html=True,
                        )
                        st.write(
                            f"**{sample.get('sample_id', sample.get('barcode'))} · "
                            f"Sample {sample.get('replicate')}**"
                        )
                        if sample.get("status") == "ERROR":
                            st.error(str(sample.get("message", "Analysis failed.")))
                            continue
                        mapping_rate = sample.get("mapping_rate")
                        st.metric(
                            "Mapping",
                            f"{float(mapping_rate):.1%}"
                            if mapping_rate is not None
                            else "N/A",
                        )
                        st.metric("Mean depth", f"{float(sample.get('mean_depth', 0)):.1f}×")
                        st.write(
                            f"SNP **{sample.get('snp', 0)}** · INS **{sample.get('insertion', 0)}** · "
                            f"DEL **{sample.get('deletion', 0)}**"
                        )
                        details = sample.get("variant_details", [])
                        if details:
                            st.dataframe(
                                pd.DataFrame(details), hide_index=True, use_container_width=True
                            )
    log_path = job.job_dir / "pipeline.log"
    if log_path.is_file():
        with st.expander("Analysis log · 필요할 때만 열기", expanded=False):
            st.code(
                log_path.read_text(encoding="utf-8", errors="replace"),
                language="text",
            )


def _batch_ont_tab(settings: dict[str, object]) -> None:
    st.subheader("Batch analysis")
    st.write(
        "Reference를 먼저 올리면 파일명별 ONT sample upload 영역이 자동으로 생성됩니다."
    )

    with st.expander("예제 데이터로 테스트하기 · 5 references / 15 samples"):
        st.write(
            "ZIP을 내려받아 압축을 풉니다. `references`의 FASTA 5개를 먼저 올린 뒤, "
            "각 reference 영역에 `demo_reads/<reference 이름>/`의 FASTQ 파일 또는 "
            "폴더를 그대로 드래그하세요."
        )
        st.download_button(
            "예제 데이터 다운로드 (ZIP)",
            data=build_demo_batch_zip(),
            file_name="ONT_plasmid_demo_5ref_15samples.zip",
            mime="application/zip",
            key="download_batch_demo",
            help="Software 테스트용 synthetic data이며 biological control로 사용할 수 없습니다.",
        )

    st.markdown("#### 1 · Reference")
    reference_uploads = st.file_uploader(
        "이번 run에 사용할 reference 파일을 모두 올리세요",
        type=["fasta", "fa", "fna", "txt"],
        accept_multiple_files=True,
        key="batch_references",
        help="업로드한 reference 개수와 파일명에 따라 ONT sample 영역이 생성됩니다.",
    )
    reference_uploads = sorted(
        list(reference_uploads or []),
        key=lambda item: natural_key(
            Path(str(getattr(item, "name", "reference"))).name
        ),
    )
    reference_names = [Path(item.name).name for item in reference_uploads]

    validation_errors: list[str] = []
    reference_rows: list[dict[str, object]] = []
    if len(set(reference_names)) != len(reference_names):
        validation_errors.append("Reference 파일명은 중복될 수 없습니다.")
    for item in reference_uploads:
        try:
            record = parse_uploaded_reference(item)
            reference_rows.append(
                {"Reference file": Path(item.name).name, "Length (bp)": len(record.sequence)}
            )
        except (BatchPreparationError, SequenceValidationError) as exc:
            validation_errors.append(str(exc))

    if reference_uploads:
        with st.expander("Reference file check", expanded=True):
            st.dataframe(pd.DataFrame(reference_rows), hide_index=True, use_container_width=True)
    for error in validation_errors:
        st.error(error)

    all_reads: list[object] = []
    mappings: list[dict[str, object]] = []
    assignment_rows: list[dict[str, object]] = []

    if reference_uploads:
        st.markdown("#### 2 · Reference별 ONT sample")
        st.caption(
            "각 영역에 해당 FASTQ 파일 또는 sample 폴더를 그대로 드래그하세요. "
            "barcode 번호와 사용자 지정 ONT sample name을 모두 인식합니다."
        )

    for index, reference_file in enumerate(reference_names, 1):
        reference_label = Path(reference_file).stem
        with st.expander(f"{index} · {reference_label}", expanded=True):
            uploaded = st.file_uploader(
                f"{reference_label}의 ONT FASTQ / sample folder",
                type=["fastq", "fq", "gz"],
                accept_multiple_files="directory",
                key=f"batch_reads_{index}_{sanitize_name(reference_label)}",
                help=(
                    "FASTQ 파일 또는 FASTQ가 들어 있는 폴더를 드래그합니다. Sample ID는 "
                    "barcode 번호, 폴더명 또는 FASTQ 파일명에서 자동으로 가져옵니다."
                ),
            )
            uploaded_files = list(uploaded or [])
            grouped = uploaded_samples(uploaded_files)

            if grouped:
                detected_ids = list(grouped)
                with st.expander("Sample ID 확인/수정", expanded=False):
                    st.caption(
                        "자동 인식된 이름이 실제 ONT sample name과 다르면 Sample ID를 수정하세요. "
                        "여러 FASTQ chunk를 하나로 묶으려면 같은 Sample ID를 입력합니다."
                    )
                    editor_rows = pd.DataFrame(
                        [
                            {
                                "Detected input": sample_id,
                                "Sample ID": sample_id,
                                "FASTQ files": len(grouped[sample_id]),
                            }
                            for sample_id in detected_ids
                        ]
                    )
                    editor_key_suffix = abs(
                        hash(tuple(str(getattr(item, "name", "")) for item in uploaded_files))
                    )
                    edited_rows = st.data_editor(
                        editor_rows,
                        hide_index=True,
                        use_container_width=True,
                        disabled=["Detected input", "FASTQ files"],
                        key=(
                            f"sample_ids_{index}_{sanitize_name(reference_label)}_"
                            f"{editor_key_suffix}"
                        ),
                        column_config={
                            "Sample ID": st.column_config.TextColumn(
                                help="결과에 표시할 ONT sample name입니다.", required=True
                            )
                        },
                    )
                renamed_groups: dict[str, list[object]] = {}
                for row in edited_rows.to_dict("records"):
                    detected_id = str(row["Detected input"])
                    sample_id = sanitize_name(str(row["Sample ID"]), detected_id)
                    renamed_groups.setdefault(sample_id, []).extend(grouped[detected_id])
                grouped = dict(
                    sorted(renamed_groups.items(), key=lambda item: natural_key(item[0]))
                )
                sample_ids = list(grouped)
                st.success(
                    f"{len(sample_ids)} samples 감지: " + ", ".join(sample_ids)
                )
                mappings.append(
                    {
                        "reference": reference_file,
                        "samples": [
                            {"name": sample_id, "files": grouped[sample_id]}
                            for sample_id in sample_ids
                        ],
                    }
                )
                all_reads.extend(uploaded_files)
                assignment_rows.append(
                    {
                        "Reference": reference_file,
                        "Samples": len(sample_ids),
                        "Sample ID": ", ".join(sample_ids),
                        "FASTQ files": len(uploaded_files),
                    }
                )
            elif uploaded_files:
                st.error("분석 가능한 FASTQ/FASTQ.GZ 파일을 찾지 못했습니다.")
            else:
                st.caption("아직 ONT sample을 올리지 않았습니다.")

    if assignment_rows:
        st.markdown("#### 3 · 분석 전 최종 확인")
        st.dataframe(
            pd.DataFrame(assignment_rows),
            hide_index=True,
            use_container_width=True,
        )

    total_samples = sum(int(row["Samples"]) for row in assignment_rows)
    metrics = st.columns(3)
    metrics[0].metric("Reference", len(reference_uploads))
    metrics[1].metric("감지된 samples", total_samples)
    metrics[2].metric("FASTQ files", len(all_reads))

    every_reference_has_reads = (
        bool(reference_uploads) and len(mappings) == len(reference_uploads)
    )
    if reference_uploads and not every_reference_has_reads:
        st.warning("모든 reference 영역에 한 개 이상의 ONT sample을 올리세요.")

    can_run = (
        every_reference_has_reads
        and not validation_errors
        and total_samples > 0
    )
    if st.button(
        f"Batch analysis 실행 ({total_samples} samples)",
        type="primary",
        disabled=not can_run,
        key="batch_run",
        use_container_width=True,
        help="각 reference 영역에 연결된 ONT sample을 분석합니다.",
    ):
        try:
            missing_tools = [
                executable
                for executable in ("bash", "minimap2", "samtools", "bcftools")
                if not shutil.which(executable)
            ]
            if missing_tools:
                raise BatchPreparationError(
                    "필수 분석 프로그램이 없습니다: " + ", ".join(missing_tools)
                )
            batch_settings = BatchSettings(
                experiment_name=str(settings["experiment_name"]),
                threads=int(settings["threads"]),
                parallel_jobs=int(settings["parallel_jobs"]),
                min_read_length=int(settings["min_length"]),
                min_read_quality=int(settings["min_quality"]),
                min_variant_quality=float(settings.get("min_variant_quality", 20.0)),
                min_variant_depth=int(settings.get("min_variant_depth", 10)),
                min_allele_fraction=float(settings.get("min_af", 0.80)),
                circular=bool(settings.get("circular", True)),
                edge_margin=int(settings.get("edge_margin", 50)),
            )
            with st.spinner("Reference와 FASTQ를 분석 폴더에 준비하고 있습니다…"):
                job = prepare_batch_job(
                    UI_RUN_ROOT,
                    reference_uploads,
                    all_reads,
                    mappings,
                    batch_settings,
                )
            progress_box = st.empty()
            with st.expander("Analysis log · 실행 중 상세 command", expanded=False):
                log_box = st.empty()
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
                    text=f"완료 {completed} / {job.sample_count} samples",
                )

            with st.spinner("Batch analysis를 실행하고 있습니다. 브라우저를 닫지 마세요…"):
                run_batch_job(REPOSITORY_ROOT, job, on_line=show_line)
                result = collect_batch_results(job)
            progress_box.success(f"{job.sample_count} samples 분석을 완료했습니다.")
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


def _settings_page(mode: str) -> None:
    st.subheader("Analysis settings")
    st.caption("각 `?`에 커서를 올리면 설정 설명을 볼 수 있습니다.")
    batch_tab, quick_tab = st.tabs(["Batch settings", "Quick comparison settings"])

    with batch_tab:
        current = _analysis_settings("Batch analysis")
        with st.form("batch_settings_form"):
            updated: dict[str, object] = {}
            pipeline_col, threshold_col = st.columns(2)
            with pipeline_col:
                st.markdown("#### Pipeline")
                updated["experiment_name"] = st.text_input(
                    "Experiment name",
                    value=str(current["experiment_name"]),
                    help="ui_runs 아래에 생성되는 결과 폴더명입니다.",
                )
                updated["threads"] = st.number_input(
                    "Threads per sample",
                    min_value=1,
                    max_value=64,
                    value=int(current["threads"]),
                    help="ONT sample 하나에 사용할 CPU threads입니다. 기본값은 8입니다.",
                )
                updated["parallel_jobs"] = st.number_input(
                    "Parallel samples",
                    min_value=1,
                    max_value=64,
                    value=int(current["parallel_jobs"]),
                    help="동시에 처리할 ONT sample의 최대 개수입니다.",
                )
                updated["min_length"] = st.number_input(
                    "Minimum read length",
                    min_value=0,
                    value=int(current["min_length"]),
                    step=50,
                    help="이 길이보다 짧은 read를 제외합니다. Coverage가 낮으면 300을 고려하세요.",
                )
                updated["min_quality"] = st.number_input(
                    "Minimum mean Q",
                    min_value=0,
                    value=int(current["min_quality"]),
                    step=1,
                    help="평균 Phred quality가 이 값보다 낮은 read를 제외합니다.",
                )
            with threshold_col:
                st.markdown("#### Variant review")
                updated["min_variant_quality"] = st.number_input(
                    "Minimum QUAL",
                    min_value=0.0,
                    value=float(current.get("min_variant_quality", 20.0)),
                    help="이 QUAL보다 낮은 variant는 REVIEW로 표시합니다.",
                )
                updated["min_variant_depth"] = st.number_input(
                    "Minimum DP",
                    min_value=0,
                    value=int(current.get("min_variant_depth", 10)),
                    help="이 read depth보다 낮은 variant는 REVIEW로 표시합니다.",
                )
                updated["min_af"] = st.number_input(
                    "Minimum allele fraction",
                    min_value=0.0,
                    max_value=1.0,
                    value=float(current.get("min_af", 0.80)),
                    step=0.05,
                    help="ALT allele을 지지하는 read 비율의 PASS 기준입니다.",
                )
                updated["circular"] = st.checkbox(
                    "Circular plasmid/vector",
                    value=bool(current.get("circular", True)),
                    help="Plasmid origin을 가로지르는 alignment와 variant 좌표를 처리합니다.",
                )
                updated["edge_margin"] = st.number_input(
                    "Edge margin (bp)",
                    min_value=0,
                    value=int(current.get("edge_margin", 50)),
                    help="Linear reference일 때 양 끝에서 REVIEW 처리할 범위입니다.",
                )
            save_batch = st.form_submit_button("Batch settings 저장", type="primary")
        if save_batch:
            st.session_state["saved_analysis_settings"]["Batch analysis"] = updated
            st.success("Batch settings를 저장했습니다.")

    with quick_tab:
        current_quick = _analysis_settings("Quick sequence comparison")
        with st.form("quick_settings_form"):
            quick_circular = st.checkbox(
                "Circular plasmid/vector",
                value=bool(current_quick["circular"]),
                help="Reference를 임시로 이어 붙여 plasmid origin을 통과하는 alignment를 허용합니다.",
            )
            save_quick = st.form_submit_button("Quick settings 저장", type="primary")
        if save_quick:
            st.session_state["saved_analysis_settings"]["Quick sequence comparison"] = {
                "circular": quick_circular
            }
            st.success("Quick comparison settings를 저장했습니다.")

    st.button("← 분석 화면으로", on_click=_close_settings_page)

def _sidebar() -> tuple[str, dict[str, object]]:
    branding = _load_branding()
    logo = _sidebar_brand_logo()
    fallback_logo = (
        f'<div class="sidebar-brand-wordmark">{html.escape(branding["organization"])}</div>'
    )
    with st.sidebar:
        st.markdown(
            """
            <style>
            [data-testid="stSidebar"] {
                background:
                    radial-gradient(circle at 78% 18%, rgba(212,180,255,.72), transparent 31%),
                    linear-gradient(155deg, #172EAE 0%, #4A38D0 53%, #A978EE 100%);
            }
            [data-testid="stSidebar"] [data-testid="stSidebarContent"] {
                padding-bottom:155px;
            }
            [data-testid="stSidebar"] h1,
            [data-testid="stSidebar"] h2,
            [data-testid="stSidebar"] h3,
            [data-testid="stSidebar"] h4,
            [data-testid="stSidebar"] p,
            [data-testid="stSidebar"] label,
            [data-testid="stSidebar"] span {
                color:#FFFFFF;
            }
            [data-testid="stSidebar"] [data-testid="stCaptionContainer"] p {
                color:rgba(255,255,255,.76);
            }
            [data-testid="stSidebar"] [data-testid="stNumberInput"] input {
                background:rgba(255,255,255,.96);
                color:#17213A;
            }
            [data-testid="stSidebar"] button {
                background:rgba(255,255,255,.14);
                color:#FFFFFF;
                border-color:rgba(255,255,255,.38);
            }
            [data-testid="stSidebar"] button:hover {
                background:rgba(255,255,255,.23);
                border-color:rgba(255,255,255,.72);
            }
            .sidebar-brand-footer {
                position:fixed; left:1.5rem; bottom:1rem; width:220px;
                box-sizing:border-box; z-index:999; padding:10px 12px 8px;
                background:rgba(255,255,255,.96); border-radius:14px;
                box-shadow:0 8px 24px rgba(16,25,90,.18);
            }
            .sidebar-brand-logo {
                display:block; width:100%; height:52px; object-fit:contain;
                object-position:center; margin:0 auto;
            }
            .sidebar-brand-wordmark {
                color:#244092; font:700 17px/1.2 Arial,sans-serif;
                letter-spacing:.04em; text-align:center; padding:18px 8px 10px;
            }
            .sidebar-brand-credit {
                color:#667085; font:500 11px/1.35 "Segoe UI",Arial,sans-serif;
                text-align:center; margin-top:5px;
            }
            </style>
            """,
            unsafe_allow_html=True,
        )
        st.markdown("### ONT Plasmid Analyzer")
        st.caption("ONT plasmid 변이 분석")
        mode = st.radio(
            "분석 메뉴",
            ["Batch analysis", "Quick sequence comparison"],
            key="analysis_mode",
            on_change=_close_settings_page,
            help=(
                "Batch analysis는 raw ONT sample을 분석합니다. Quick comparison은 완성된 "
                "두 sequence의 차이만 빠르게 확인합니다."
            ),
        )

        st.button(
            "⚙ Analysis settings",
            use_container_width=True,
            on_click=_open_settings_page,
            help="선택한 분석의 설정을 가운데 화면에서 변경합니다.",
        )
        st.markdown(
            (
                '<div class="sidebar-brand-footer">'
                f'{logo or fallback_logo}<div class="sidebar-brand-credit">'
                f'Distributed by {html.escape(branding["distributed_by"])}</div></div>'
            ),
            unsafe_allow_html=True,
        )

    return mode, dict(_analysis_settings(mode))

def main() -> None:
    st.set_page_config(
        page_title="ONT Plasmid Analyzer",
        page_icon="🧬",
        layout="wide",
    )
    mode, settings = _sidebar()
    _brand_header()
    if st.session_state.get("show_analysis_settings", False):
        _settings_page(mode)
    elif mode == "Batch analysis":
        _batch_ont_tab(settings)
    elif mode == "Quick sequence comparison":
        _quick_compare_tab(settings)
    else:
        _raw_ont_tab(settings)


if __name__ == "__main__":
    main()
