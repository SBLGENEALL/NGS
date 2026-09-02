#!/usr/bin/env python3
"""Local, offline Streamlit UI for the ONT reference-mapping pipeline."""

from __future__ import annotations

import shutil
import base64
import html
import json
import os
from pathlib import Path

import pandas as pd
import streamlit as st

try:
    import pyarrow  # noqa: F401
except (ImportError, ModuleNotFoundError):
    PYARROW_AVAILABLE = False
else:
    PYARROW_AVAILABLE = True

from ont_ui.batch import (
    BatchExecutionError,
    BatchPreparationError,
    BatchSettings,
    ServerPathUpload,
    collect_batch_results,
    natural_key,
    parse_uploaded_reference,
    prepare_batch_job,
    run_batch_job,
    server_fastq_uploads,
    server_reference_uploads,
    uploaded_samples,
)
from ont_ui.sequences import (
    SequenceValidationError,
    sanitize_name,
)


REPOSITORY_ROOT = Path(__file__).resolve().parent
UI_RUN_ROOT = REPOSITORY_ROOT / "ui_runs"
SERVER_DATA_ROOT = Path(
    os.environ.get("ONT_SERVER_ROOT", "/data")
).expanduser().resolve()
SERVER_BROWSER_START = SERVER_DATA_ROOT

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
}


def _display_table(frame: pd.DataFrame) -> None:
    """Render a table without making pyarrow mandatory in offline installs."""
    if PYARROW_AVAILABLE:
        try:
            st.dataframe(frame, use_container_width=True, hide_index=True)
            return
        except (ImportError, ModuleNotFoundError):
            pass
    def markdown_cell(value: object) -> str:
        if pd.isna(value):
            return ""
        return str(value).replace("\\", "\\\\").replace("|", "\\|").replace("\n", " ")

    columns = [markdown_cell(column) for column in frame.columns]
    rows = [
        "| " + " | ".join(markdown_cell(value) for value in row) + " |"
        for row in frame.itertuples(index=False, name=None)
    ]
    st.markdown(
        "\n".join(
            [
                "| " + " | ".join(columns) + " |",
                "| " + " | ".join("---" for _ in columns) + " |",
                *rows,
            ]
        )
    )


def _remember_server_uploads(
    state_key: str,
    uploads: list[ServerPathUpload],
) -> None:
    st.session_state[state_key] = [
        (str(item.source_path), item.name) for item in uploads
    ]


def _restore_server_uploads(state_key: str) -> list[ServerPathUpload]:
    stored = st.session_state.get(state_key, [])
    if not isinstance(stored, list):
        return []
    restored: list[ServerPathUpload] = []
    for item in stored:
        if not isinstance(item, (list, tuple)) or len(item) != 2:
            continue
        try:
            restored.append(ServerPathUpload(Path(str(item[0])), str(item[1])))
        except BatchPreparationError:
            continue
    return restored


def _clear_batch_assignments() -> None:
    for key in list(st.session_state):
        if str(key).startswith("assignment_"):
            st.session_state.pop(key, None)


def _inside_server_root(path: Path) -> bool:
    """Return True only for paths contained by the configured data root."""
    try:
        path.resolve().relative_to(SERVER_DATA_ROOT)
    except (OSError, ValueError):
        return False
    return True


def _select_browser_folder(state_key: str, folder: Path) -> None:
    """Persist a validated browser folder and clear any previous selection."""
    resolved = folder.resolve()
    if resolved.is_dir() and _inside_server_root(resolved):
        st.session_state[f"{state_key}_current"] = str(resolved)


def _server_folder_browser(
    title: str,
    state_key: str,
    file_suffixes: tuple[str, ...],
) -> Path | None:
    """Render a click-only server folder picker rooted at SERVER_DATA_ROOT."""
    if not SERVER_DATA_ROOT.is_dir():
        st.error("서버 data 폴더에 접근할 수 없습니다. 관리자에게 확인하세요.")
        return None

    current_raw = st.session_state.get(
        f"{state_key}_current", str(SERVER_BROWSER_START)
    )
    current = Path(str(current_raw)).expanduser().resolve()
    if not current.is_dir() or not _inside_server_root(current):
        current = SERVER_DATA_ROOT
        st.session_state[f"{state_key}_current"] = str(current)

    selected_raw = st.session_state.get(f"{state_key}_selected")
    selected = Path(str(selected_raw)).resolve() if selected_raw else None
    if selected is not None and (not selected.is_dir() or not _inside_server_root(selected)):
        selected = None
        st.session_state.pop(f"{state_key}_selected", None)

    with st.container(border=True):
        st.markdown(f"**{title}**")
        relative = current.relative_to(SERVER_DATA_ROOT)
        shown_path = SERVER_DATA_ROOT if str(relative) == "." else SERVER_DATA_ROOT / relative
        st.caption(f"현재 위치: {shown_path}")

        navigation = st.columns([1, 1, 3])
        with navigation[0]:
            st.button(
                "⬆ 상위",
                key=f"{state_key}_up",
                disabled=current == SERVER_DATA_ROOT,
                on_click=_select_browser_folder,
                args=(state_key, current.parent),
                use_container_width=True,
            )
        with navigation[1]:
            if st.button(
                "⌂ 시작 위치",
                key=f"{state_key}_root",
                disabled=current == SERVER_BROWSER_START,
                use_container_width=True,
            ):
                st.session_state[f"{state_key}_current"] = str(SERVER_BROWSER_START)
                st.rerun()

        try:
            child_folders = sorted(
                (item for item in current.iterdir() if item.is_dir()),
                key=lambda item: natural_key(item.name),
            )
            matching_files = sum(
                1
                for item in current.iterdir()
                if item.is_file() and item.name.lower().endswith(file_suffixes)
            )
        except PermissionError:
            st.error("이 폴더를 열 권한이 없습니다. 다른 폴더를 선택하세요.")
            child_folders = []
            matching_files = 0
        except OSError:
            st.error("폴더 목록을 읽을 수 없습니다.")
            child_folders = []
            matching_files = 0

        if child_folders:
            st.caption("폴더를 클릭해 이동하세요.")
            folder_columns = st.columns(3)
            for index, folder in enumerate(child_folders[:150]):
                with folder_columns[index % 3]:
                    st.button(
                        f"📁 {folder.name}",
                        key=f"{state_key}_dir_{index}_{folder.name}",
                        on_click=_select_browser_folder,
                        args=(state_key, folder),
                        use_container_width=True,
                    )
            if len(child_folders) > 150:
                st.warning("하위 폴더가 많아 이름순으로 150개만 표시합니다.")
        else:
            st.caption("표시할 하위 폴더가 없습니다.")

        st.caption(
            f"현재 폴더 바로 아래에서 지원 파일 {matching_files}개를 확인했습니다. "
            "선택하면 하위 폴더까지 함께 검색합니다."
        )
        choose, clear = st.columns([3, 1])
        with choose:
            if st.button(
                "✓ 이 폴더 사용",
                key=f"{state_key}_choose",
                type="primary",
                use_container_width=True,
            ):
                st.session_state[f"{state_key}_selected"] = str(current)
                st.session_state[f"{state_key}_reload"] = (
                    int(st.session_state.get(f"{state_key}_reload", 0)) + 1
                )
                st.rerun()
        with clear:
            if st.button(
                "선택 해제",
                key=f"{state_key}_clear",
                use_container_width=True,
            ):
                st.session_state.pop(f"{state_key}_selected", None)
                st.session_state[f"{state_key}_reload"] = (
                    int(st.session_state.get(f"{state_key}_reload", 0)) + 1
                )
                st.rerun()

        if selected is not None:
            st.success(f"선택됨: {selected.name or selected}")
    return selected


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
            data = path.read_bytes()
            valid = (
                filename.endswith(".png") and data.startswith(b"\x89PNG\r\n\x1a\n")
            ) or (
                filename.endswith(".svg") and b"<svg" in data[:2048].lower()
            )
            if not valid:
                continue
            encoded = base64.b64encode(data).decode("ascii")
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
        #MainMenu { visibility:hidden; }
        [data-testid="stSidebarCollapsedControl"],
        [data-testid="stSidebarCollapseButton"] {
            visibility:visible !important; display:flex !important; z-index:1002;
        }
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
        if PYARROW_AVAILABLE:
            st.line_chart(depth_frame)
        else:
            st.caption(
                "현재 offline 환경의 pyarrow가 호환되지 않아 depth chart는 생략합니다. "
                "Mean depth와 coverage 결과에는 영향이 없습니다."
            )

    st.markdown("#### Variant calls")
    if events:
        _display_table(_variants_frame(events))
    elif vcf_path.is_file():
        st.success("No variant was called against the reference.")
    else:
        st.warning("VCF가 생성되지 않았습니다. 입력 파일과 Analysis settings를 확인하세요.")

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

    st.caption("대용량 BAM과 depth 파일은 서버에 안전하게 보관됩니다.")


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
        summary_frame = pd.DataFrame(summary_rows)
        for column, suffix, decimals in (
            ("Mapping rate (%)", "%", 2),
            ("Coverage ≥1× (%)", "%", 2),
            ("Mean depth", "×", 1),
        ):
            summary_frame[column] = summary_frame[column].map(
                lambda value, s=suffix, d=decimals: ""
                if pd.isna(value)
                else f"{float(value):.{d}f}{s}"
            )
        _display_table(summary_frame)
    st.download_button(
        "Download complete summary (CSV)",
        data=str(result.get("summary_csv", "")).encode("utf-8-sig"),
        file_name=f"{job.experiment_name}_ONT_summary.csv",
        mime="text/csv",
        type="primary",
    )
    groups = result.get("groups", {})
    assert isinstance(groups, dict)
    show_details = st.checkbox(
        "Reference별 상세 결과 표시",
        value=False,
        key=f"show_batch_details_{job.job_dir.name}",
        help="필요할 때만 상세 결과 card와 variant table을 생성해 화면 속도를 유지합니다.",
    )
    if show_details:
        st.markdown("#### Reference-by-reference review")
    else:
        groups = {}
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
                            st.error("분석 결과를 생성하지 못했습니다. 입력 파일과 설정을 확인하세요.")
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
                            _display_table(pd.DataFrame(details))


def _batch_ont_tab(settings: dict[str, object]) -> None:
    st.title("Batch analysis")
    st.caption("폴더를 선택한 뒤 각 Reference에 분석할 ONT sample을 직접 지정하세요.")

    reference_is_selected = bool(st.session_state.get("reference_browser_selected"))
    with st.expander(
        "1 · Reference 폴더 선택",
        expanded=not reference_is_selected,
    ):
        reference_folder = _server_folder_browser(
            "서버 폴더",
            "reference_browser",
            (".fasta", ".fa", ".fna", ".txt"),
        )
    reference_reload = (
        str(reference_folder) if reference_folder else None,
        st.session_state.get("reference_browser_reload", 0),
    )
    if st.session_state.get("loaded_reference_folder") != reference_reload:
        st.session_state["loaded_reference_folder"] = reference_reload
        st.session_state.pop("loaded_server_references", None)
        st.session_state.pop("batch_result", None)
        _clear_batch_assignments()
        if reference_folder is not None:
            try:
                _remember_server_uploads(
                    "loaded_server_references",
                    list(server_reference_uploads(str(reference_folder))),
                )
            except BatchPreparationError as exc:
                st.error(str(exc))

    reference_uploads: list[object] = _restore_server_uploads(
        "loaded_server_references"
    )
    reference_uploads = sorted(
        reference_uploads,
        key=lambda item: natural_key(
            Path(str(getattr(item, "name", "reference"))).name
        ),
    )
    reference_names = [Path(item.name).name for item in reference_uploads]
    reference_signature = tuple(
        (Path(str(getattr(item, "name", ""))).name, int(getattr(item, "size", 0)))
        for item in reference_uploads
    )
    previous_signature = st.session_state.get("batch_reference_signature")
    if previous_signature != reference_signature:
        st.session_state["batch_reference_signature"] = reference_signature
        st.session_state.pop("batch_result", None)

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
        st.success(f"Reference {len(reference_uploads)}개 선택 완료")
        with st.expander(f"Reference 확인 · {len(reference_uploads)}개", expanded=False):
            for row in reference_rows:
                st.write(f"• {row['Reference file']}  ·  {row['Length (bp)']:,} bp")
    for error in validation_errors:
        st.error(error)

    ont_is_selected = bool(st.session_state.get("ont_browser_selected"))
    with st.expander(
        "2 · ONT 결과 폴더 선택",
        expanded=not ont_is_selected,
    ):
        st.caption("sample 폴더들이 들어 있는 상위 폴더를 선택하세요.")
        ont_folder = _server_folder_browser(
            "서버 폴더",
            "ont_browser",
            (".fastq", ".fq", ".fastq.gz", ".fq.gz", ".gz"),
        )
    ont_reload = (
        str(ont_folder) if ont_folder else None,
        st.session_state.get("ont_browser_reload", 0),
    )
    if st.session_state.get("loaded_ont_folder") != ont_reload:
        st.session_state["loaded_ont_folder"] = ont_reload
        st.session_state.pop("loaded_server_fastqs", None)
        st.session_state.pop("batch_result", None)
        _clear_batch_assignments()
        if ont_folder is not None:
            try:
                _remember_server_uploads(
                    "loaded_server_fastqs",
                    list(server_fastq_uploads(str(ont_folder))),
                )
            except BatchPreparationError as exc:
                st.error(str(exc))

    all_reads: list[object] = _restore_server_uploads("loaded_server_fastqs")
    grouped_reads = uploaded_samples(all_reads)
    detected_sample_ids = list(grouped_reads)
    if grouped_reads:
        st.success(f"ONT sample {len(grouped_reads)}개 선택 완료")

    mappings: list[dict[str, object]] = []
    selected_by_reference: dict[str, list[str]] = {}
    if reference_uploads and grouped_reads:
        st.markdown("### 3 · Reference별 ONT sample 선택")
        st.caption("이름이나 순서를 추정하지 않습니다. 분석할 조합을 직접 선택하세요.")
        assignment_columns = st.columns(2)
        for index, reference_file in enumerate(reference_names, 1):
            with assignment_columns[(index - 1) % 2]:
                reference_label = Path(reference_file).stem
                selected_ids = st.multiselect(
                    reference_label,
                    options=detected_sample_ids,
                    key=f"assignment_{index}_{sanitize_name(reference_file, str(index))}",
                    help=(
                        "이 Reference와 비교할 ONT sample을 직접 선택합니다. "
                        "한 sample을 여러 Reference에 선택해도 됩니다."
                    ),
                    placeholder="ONT sample 선택",
                )
                selected_by_reference[reference_file] = list(selected_ids)

    assignment_signature = tuple(
        (reference_file, tuple(sample_ids))
        for reference_file, sample_ids in selected_by_reference.items()
    )
    if st.session_state.get("batch_assignment_signature") != assignment_signature:
        st.session_state["batch_assignment_signature"] = assignment_signature
        st.session_state.pop("batch_result", None)

    for reference_file in reference_names:
        sample_ids = selected_by_reference.get(reference_file, [])
        if not sample_ids:
            continue
        mappings.append(
            {
                "reference": reference_file,
                "samples": [
                    {"name": sample_id, "files": grouped_reads[sample_id]}
                    for sample_id in sample_ids
                ],
            }
        )

    total_samples = sum(
        len(sample_ids) for sample_ids in selected_by_reference.values()
    )
    if total_samples:
        st.success(
            f"분석 조합 {total_samples}개 선택 완료 · Reference {len(mappings)}개"
        )

    can_run = (
        bool(mappings)
        and not validation_errors
        and total_samples > 0
    )
    if st.button(
        f"Batch analysis 실행 ({total_samples} samples)",
        type="primary",
        disabled=not can_run,
        key="batch_run",
        use_container_width=True,
        help="확인한 Reference와 ONT sample 배정대로 분석합니다.",
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
            progress_box.progress(0.0, text=f"완료 0 / {job.sample_count} samples")
            completed = 0

            def show_line(line: str) -> None:
                nonlocal completed
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
        except BatchExecutionError:
            st.error("분석 실행 중 오류가 발생했습니다. 입력 파일과 설정을 확인하세요.")
        except (BatchPreparationError, SequenceValidationError, OSError) as exc:
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
            progress_placeholder = st.empty()

            def show_line(line: str) -> None:
                if "[" in line and "/6]" in line:
                    progress_placeholder.info("분석을 진행하고 있습니다…")

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


def _settings_page() -> None:
    st.subheader("Analysis settings")
    st.caption("각 `?`에 커서를 올리면 설정 설명을 볼 수 있습니다.")
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

    st.button("← 분석 화면으로", on_click=_close_settings_page)


def _sidebar_tool_status() -> None:
    st.markdown("##### 분석 도구")
    required = ("minimap2", "samtools", "bcftools", "gzip", "bash")
    for tool in required:
        icon = "✅" if shutil.which(tool) else "❌"
        st.caption(f"{icon} {tool}")
    nanofilt_icon = "✅" if shutil.which("NanoFilt") else "⚪"
    st.caption(f"{nanofilt_icon} NanoFilt · optional QC")

def _sidebar() -> dict[str, object]:
    branding = _load_branding()
    logo = _sidebar_brand_logo()
    wordmark = html.escape(branding["organization"]).replace(" ", "<br>", 1)
    fallback_logo = f'<div class="sidebar-brand-wordmark">{wordmark}</div>'
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
                padding-bottom:210px;
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
            [data-testid="stSidebar"] [data-testid="stButton"] button[kind="primary"] {
                min-height:3.45rem; font-size:1.12rem; font-weight:750;
                background:rgba(255,255,255,.96); color:#2639A7;
                border:0; box-shadow:0 8px 22px rgba(17,26,100,.24);
            }
            [data-testid="stSidebar"] [data-testid="stButton"] button[kind="primary"] p {
                color:#2639A7;
            }
            [data-testid="stSidebar"] [data-testid="stButton"]:has(button[kind="secondary"]) {
                position:fixed; left:1.5rem; bottom:6.6rem; width:220px;
                z-index:1000;
            }
            .sidebar-brand-footer {
                position:fixed; left:1.5rem; bottom:1rem; width:220px;
                box-sizing:border-box; z-index:999; padding:4px 8px;
                background:transparent; border:0; box-shadow:none;
            }
            .sidebar-brand-logo {
                display:block; width:100%; height:52px; object-fit:contain;
                object-position:center; margin:0 auto;
            }
            .sidebar-brand-wordmark {
                color:#FFFFFF; font:800 18px/1.12 "Segoe UI",Arial,sans-serif;
                letter-spacing:.08em; text-align:left; padding:0;
                white-space:pre-line; text-shadow:0 2px 8px rgba(21,27,88,.25);
            }
            .sidebar-brand-credit {
                color:rgba(255,255,255,.78); font:500 11px/1.35 "Segoe UI",Arial,sans-serif;
                text-align:left; margin-top:7px;
            }
            </style>
            """,
            unsafe_allow_html=True,
        )
        st.markdown("### ONT Plasmid Analyzer")
        st.caption("ONT plasmid 변이 분석")

        st.button(
            "🧬  Batch analysis",
            type="primary",
            use_container_width=True,
            on_click=_close_settings_page,
            help="분석 화면으로 이동합니다.",
        )
        _sidebar_tool_status()

        st.button(
            "⚙ Analysis settings",
            type="secondary",
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

    return dict(_analysis_settings("Batch analysis"))

def main() -> None:
    st.set_page_config(
        page_title="ONT Plasmid Analyzer",
        page_icon="🧬",
        layout="wide",
    )
    settings = _sidebar()
    _brand_header()
    if st.session_state.get("show_analysis_settings", False):
        _settings_page()
    else:
        _batch_ont_tab(settings)


if __name__ == "__main__":
    main()
