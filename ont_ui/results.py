"""Parse ONT pipeline QC, depth, and VCF outputs for the local UI."""

from __future__ import annotations

import csv
import gzip
import io
import math
import re
from pathlib import Path
from typing import Iterable, TextIO

from .models import DepthMetrics, MappingMetrics, VariantEvent
from .sequences import max_homopolymer_run_near, sequence_context


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def parse_flagstat(path: Path) -> MappingMetrics:
    total = mapped = None
    mapping_rate = None
    if not path.is_file():
        return MappingMetrics()
    total_pattern = re.compile(r"^(\d+)\s+\+\s+\d+\s+in total")
    mapped_pattern = re.compile(r"^(\d+)\s+\+\s+\d+\s+mapped\s+\(([0-9.]+)%")
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if total is None:
                match = total_pattern.search(line)
                if match:
                    total = int(match.group(1))
            match = mapped_pattern.search(line)
            if match:
                mapped = int(match.group(1))
                mapping_rate = float(match.group(2)) / 100.0
    return MappingMetrics(total_reads=total, mapped_reads=mapped, mapping_rate=mapping_rate)


def read_depth(
    path: Path,
    max_chart_points: int = 5000,
) -> tuple[DepthMetrics, list[tuple[int, int]]]:
    depths: list[tuple[int, int]] = []
    total_depth = 0
    positions = covered_1x = covered_10x = 0
    min_depth: int | None = None
    max_depth = 0
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            try:
                position = int(cols[1])
                depth = int(cols[2])
            except ValueError:
                continue
            depths.append((position, depth))
            positions += 1
            total_depth += depth
            covered_1x += depth >= 1
            covered_10x += depth >= 10
            min_depth = depth if min_depth is None else min(min_depth, depth)
            max_depth = max(max_depth, depth)

    if not positions:
        metrics = DepthMetrics(0, 0.0, 0, 0, 0.0, 0.0)
        return metrics, []
    metrics = DepthMetrics(
        positions=positions,
        mean_depth=total_depth / positions,
        min_depth=min_depth or 0,
        max_depth=max_depth,
        coverage_1x=covered_1x / positions,
        coverage_10x=covered_10x / positions,
    )
    step = max(1, math.ceil(len(depths) / max_chart_points))
    chart_points = depths[::step]
    if depths[-1] not in chart_points:
        chart_points.append(depths[-1])
    return metrics, chart_points


def _parse_info(value: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for entry in value.split(";"):
        if "=" in entry:
            key, item = entry.split("=", 1)
            parsed[key] = item
        elif entry:
            parsed[entry] = "true"
    return parsed


def _number(value: str | None, cast):
    if value in {None, "", "."}:
        return None
    try:
        return cast(value)
    except (TypeError, ValueError):
        return None


def _classify_vcf(pos: int, ref: str, alt: str) -> tuple[int, str, str, str, int]:
    if len(ref) == 1 and len(alt) == 1:
        return pos, "SNP", ref, alt, 1
    if len(alt) > len(ref) and alt.startswith(ref):
        inserted = alt[len(ref) :]
        return pos + len(ref) - 1, "Insertion", "-", inserted, len(inserted)
    if len(ref) > len(alt) and ref.startswith(alt):
        deleted = ref[len(alt) :]
        return pos + len(alt), "Deletion", deleted, "-", len(deleted)
    if len(alt) > len(ref):
        return pos, "Insertion", ref, alt, len(alt) - len(ref)
    if len(ref) > len(alt):
        return pos, "Deletion", ref, alt, len(ref) - len(alt)
    return pos, "Complex", ref, alt, len(ref)


def parse_vcf_variants(
    path: Path,
    reference_sequence: str = "",
    min_quality: float = 20.0,
    min_depth: int = 10,
    min_allele_fraction: float = 0.80,
    edge_margin: int = 50,
    circular: bool = True,
    homopolymer_threshold: int = 4,
) -> list[VariantEvent]:
    """Read every ALT allele and mark calls that need review.

    Calls are not hidden by the thresholds. Instead, each call receives a
    `PASS` or `REVIEW` status and explicit warnings so users can inspect
    borderline ONT homopolymer/coverage events rather than silently losing
    them.
    """

    events: list[VariantEvent] = []
    ref_len = len(reference_sequence)
    with _open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            try:
                raw_pos = int(cols[1])
            except ValueError:
                continue
            ref = cols[3].upper()
            alts = [alt.upper() for alt in cols[4].split(",") if alt not in {"", "."}]
            quality = _number(cols[5], float)
            info = _parse_info(cols[7])

            sample_values: dict[str, str] = {}
            if len(cols) >= 10:
                keys = cols[8].split(":")
                values = cols[9].split(":")
                sample_values = dict(zip(keys, values))

            depth = _number(info.get("DP"), int)
            if depth is None:
                depth = _number(sample_values.get("DP"), int)

            info_af = info.get("AF", "").split(",") if info.get("AF") else []
            ad_values: list[int] = []
            if sample_values.get("AD") not in {None, "", "."}:
                try:
                    ad_values = [int(item) for item in sample_values["AD"].split(",")]
                except ValueError:
                    ad_values = []

            for alt_index, alt in enumerate(alts):
                allele_fraction = None
                if alt_index < len(info_af):
                    allele_fraction = _number(info_af[alt_index], float)
                if allele_fraction is None and len(ad_values) > alt_index + 1:
                    ad_total = sum(ad_values)
                    if ad_total > 0:
                        allele_fraction = ad_values[alt_index + 1] / ad_total

                pos, kind, display_ref, display_alt, length = _classify_vcf(
                    raw_pos, ref, alt
                )
                warnings: list[str] = []
                threshold_failed = False
                if quality is None:
                    warnings.append("QUAL unavailable")
                    threshold_failed = min_quality > 0
                elif quality < min_quality:
                    warnings.append(f"Low QUAL (<{min_quality:g})")
                    threshold_failed = True

                if depth is None:
                    warnings.append("Depth unavailable")
                    threshold_failed = min_depth > 0
                elif depth < min_depth:
                    warnings.append(f"Low depth (<{min_depth})")
                    threshold_failed = True

                if allele_fraction is None:
                    warnings.append("Allele fraction unavailable")
                    threshold_failed = min_allele_fraction > 0
                elif allele_fraction < min_allele_fraction:
                    warnings.append(f"Low allele fraction (<{min_allele_fraction:.2f})")
                    threshold_failed = True

                if circular and ref_len and edge_margin > 0:
                    if pos <= edge_margin or pos > ref_len - edge_margin:
                        warnings.append(f"Reference-edge region ({edge_margin} bp)")
                        threshold_failed = True

                if reference_sequence:
                    run = max_homopolymer_run_near(
                        reference_sequence, pos, circular=circular
                    )
                    if run >= homopolymer_threshold:
                        warnings.append(f"Nearby homopolymer ({run} bp)")

                context = (
                    sequence_context(reference_sequence, pos, circular=circular)
                    if reference_sequence
                    else ""
                )
                events.append(
                    VariantEvent(
                        position=pos,
                        kind=kind,
                        ref=display_ref,
                        alt=display_alt,
                        length=length,
                        context=context,
                        quality=quality,
                        depth=depth,
                        allele_fraction=allele_fraction,
                        status="REVIEW" if threshold_failed else "PASS",
                        warnings=tuple(warnings),
                    )
                )
    return events


def variants_csv(events: Iterable[VariantEvent]) -> str:
    rows = [event.as_row() for event in events]
    output = io.StringIO()
    fieldnames = [
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
    writer = csv.DictWriter(output, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()
