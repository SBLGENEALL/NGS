"""Small append-only usage log for the local ONT analysis UI."""

from __future__ import annotations

import csv
import fcntl
from datetime import datetime, timedelta, timezone
from pathlib import Path


KST = timezone(timedelta(hours=9))
LOG_COLUMNS = (
    "timestamp_kst",
    "event",
    "project_name",
    "user_name",
    "reference_count",
    "analysis_count",
    "job_id",
)


def usage_timestamp() -> str:
    return datetime.now(KST).isoformat(timespec="seconds")


def _safe_csv_text(value: object, limit: int = 120) -> str:
    text = " ".join(str(value or "").split())[:limit]
    if text.startswith(("=", "+", "-", "@")):
        text = "'" + text
    return text


def append_usage_event(
    log_path: Path,
    *,
    event: str,
    project_name: str,
    user_name: str,
    reference_count: int,
    analysis_count: int,
    job_id: str = "",
) -> str:
    """Append one concurrency-safe CSV row and return its KST timestamp."""
    timestamp = usage_timestamp()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a+", encoding="utf-8", newline="") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            handle.seek(0, 2)
            is_empty = handle.tell() == 0
            writer = csv.DictWriter(handle, fieldnames=LOG_COLUMNS, lineterminator="\n")
            if is_empty:
                writer.writeheader()
            writer.writerow(
                {
                    "timestamp_kst": timestamp,
                    "event": _safe_csv_text(event, 30),
                    "project_name": _safe_csv_text(project_name),
                    "user_name": _safe_csv_text(user_name),
                    "reference_count": max(0, int(reference_count)),
                    "analysis_count": max(0, int(analysis_count)),
                    "job_id": _safe_csv_text(job_id, 160),
                }
            )
            handle.flush()
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    return timestamp
