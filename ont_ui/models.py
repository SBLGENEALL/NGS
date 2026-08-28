"""Small data models shared by the UI and analysis backends."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass(frozen=True)
class VariantEvent:
    """A human-readable sequence difference in reference coordinates."""

    position: int
    kind: str
    ref: str
    alt: str
    length: int
    context: str = ""
    quality: Optional[float] = None
    depth: Optional[int] = None
    allele_fraction: Optional[float] = None
    status: str = "PASS"
    warnings: tuple[str, ...] = field(default_factory=tuple)

    def as_row(self) -> dict[str, object]:
        """Return stable, presentation-ready column names."""

        return {
            "Position": self.position,
            "Type": self.kind,
            "REF": self.ref,
            "ALT": self.alt,
            "Length (bp)": self.length,
            "QUAL": self.quality,
            "Depth": self.depth,
            "Allele fraction": self.allele_fraction,
            "Status": self.status,
            "Warnings": "; ".join(self.warnings),
            "Reference context": self.context,
        }


@dataclass(frozen=True)
class AlignmentSummary:
    orientation: str
    identity: float
    query_coverage: float
    mapq: int
    reference_start: int
    reference_end: int
    query_start: int
    query_end: int
    query_length: int
    variants: tuple[VariantEvent, ...]
    raw_paf: str = ""


@dataclass(frozen=True)
class MappingMetrics:
    total_reads: Optional[int] = None
    mapped_reads: Optional[int] = None
    mapping_rate: Optional[float] = None


@dataclass(frozen=True)
class DepthMetrics:
    positions: int
    mean_depth: float
    min_depth: int
    max_depth: int
    coverage_1x: float
    coverage_10x: float


@dataclass(frozen=True)
class PipelineRunResult:
    job_dir: Path
    result_dir: Path
    sample_dir: Path
    sample_name: str
    log_path: Path
