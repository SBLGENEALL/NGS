"""DNA input validation and small sequence utilities."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


IUPAC_DNA = frozenset("ACGTRYSWKMBDHVN")
_SAFE_NAME = re.compile(r"[^A-Za-z0-9._-]+")


class SequenceValidationError(ValueError):
    """Raised when pasted or uploaded sequence content is not valid DNA."""


@dataclass(frozen=True)
class SequenceRecord:
    name: str
    sequence: str


def sanitize_name(value: str, default: str = "sample", max_length: int = 80) -> str:
    """Make a user label safe to use as one path component."""

    cleaned = _SAFE_NAME.sub("_", (value or "").strip()).strip("._-")
    if not cleaned or cleaned in {".", ".."}:
        cleaned = default
    return cleaned[:max_length]


def parse_single_sequence(text: str, default_name: str = "sequence") -> SequenceRecord:
    """Parse one FASTA record or a plain DNA sequence.

    Spaces and line breaks are ignored. RNA `U` is accepted and normalized
    to DNA `T`. Multiple FASTA records are rejected because the UI compares
    exactly one reference with one query at a time.
    """

    if text is None:
        raise SequenceValidationError("No sequence was provided.")
    text = text.lstrip("\ufeff").strip()
    if not text:
        raise SequenceValidationError("No sequence was provided.")

    lines = text.splitlines()
    fasta_headers = [i for i, line in enumerate(lines) if line.lstrip().startswith(">")]
    name = default_name
    if fasta_headers:
        if fasta_headers[0] != 0:
            raise SequenceValidationError("FASTA content must start with a '>' header line.")
        if len(fasta_headers) != 1:
            raise SequenceValidationError(
                "Exactly one FASTA record is allowed in each input. "
                "Upload or paste one sequence at a time."
            )
        header = lines[0].strip()[1:].strip()
        if header:
            name = header.split()[0]
        sequence_text = "".join(lines[1:])
    else:
        sequence_text = "".join(lines)

    sequence = re.sub(r"\s+", "", sequence_text).upper().replace("U", "T")
    if not sequence:
        raise SequenceValidationError("The sequence contains no bases.")

    invalid = sorted(set(sequence) - IUPAC_DNA)
    if invalid:
        shown = " ".join(repr(char) for char in invalid[:10])
        raise SequenceValidationError(
            f"Unsupported character(s) in DNA sequence: {shown}. "
            "Use IUPAC DNA bases only."
        )

    return SequenceRecord(sanitize_name(name, default_name), sequence)


def write_fasta(path: Path, record: SequenceRecord, width: int = 80) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write(f">{record.name}\n")
        for start in range(0, len(record.sequence), width):
            handle.write(record.sequence[start : start + width] + "\n")


def normalize_position(position: int, reference_length: int, circular: bool) -> int:
    if reference_length <= 0:
        raise ValueError("Reference length must be positive.")
    if circular:
        return ((position - 1) % reference_length) + 1
    return max(1, min(position, reference_length))


def sequence_context(
    sequence: str,
    position: int,
    radius: int = 8,
    circular: bool = False,
) -> str:
    """Return context with the base at `position` enclosed in brackets."""

    if not sequence:
        return ""
    position = normalize_position(position, len(sequence), circular)
    center = position - 1
    if circular:
        offsets = range(-radius, radius + 1)
        chars = [sequence[(center + offset) % len(sequence)] for offset in offsets]
        chars[radius] = f"[{chars[radius]}]"
        return "".join(chars)

    start = max(0, center - radius)
    end = min(len(sequence), center + radius + 1)
    left = sequence[start:center]
    base = sequence[center]
    right = sequence[center + 1 : end]
    return f"{left}[{base}]{right}"


def max_homopolymer_run_near(
    sequence: str,
    position: int,
    radius: int = 8,
    circular: bool = False,
) -> int:
    """Return the longest identical-base run in a small variant window."""

    if not sequence:
        return 0
    position = normalize_position(position, len(sequence), circular)
    center = position - 1
    if circular:
        window = "".join(
            sequence[(center + offset) % len(sequence)]
            for offset in range(-radius, radius + 1)
        )
    else:
        window = sequence[max(0, center - radius) : min(len(sequence), center + radius + 1)]

    longest = current = 0
    previous = ""
    for base in window:
        if base == previous:
            current += 1
        else:
            previous = base
            current = 1
        longest = max(longest, current)
    return longest
