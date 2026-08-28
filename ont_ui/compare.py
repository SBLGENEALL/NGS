"""Reference-vs-query comparison using minimap2's base-level `cs` tag."""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from .models import AlignmentSummary, VariantEvent
from .sequences import (
    SequenceRecord,
    normalize_position,
    sequence_context,
    write_fasta,
)


class SequenceComparisonError(RuntimeError):
    """Raised when minimap2 cannot produce a usable comparison."""


_CS_TOKEN = re.compile(
    r":\d+|=[A-Za-z]+|\*[A-Za-z][A-Za-z]|\+[A-Za-z]+|-[A-Za-z]+|~[A-Za-z]{2}\d+[A-Za-z]{2}"
)


@dataclass(frozen=True)
class PafHit:
    query_name: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_length: int
    target_start: int
    target_end: int
    matches: int
    block_length: int
    mapq: int
    cs: str
    raw: str


def parse_paf_line(line: str) -> PafHit:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 12:
        raise SequenceComparisonError("minimap2 returned a malformed PAF alignment.")
    tags: dict[str, str] = {}
    for raw_tag in fields[12:]:
        pieces = raw_tag.split(":", 2)
        if len(pieces) == 3:
            tags[pieces[0]] = pieces[2]
    cs = tags.get("cs", "")
    if not cs:
        raise SequenceComparisonError("The minimap2 alignment is missing its cs tag.")
    return PafHit(
        query_name=fields[0],
        query_length=int(fields[1]),
        query_start=int(fields[2]),
        query_end=int(fields[3]),
        strand=fields[4],
        target_name=fields[5],
        target_length=int(fields[6]),
        target_start=int(fields[7]),
        target_end=int(fields[8]),
        matches=int(fields[9]),
        block_length=int(fields[10]),
        mapq=int(fields[11]),
        cs=cs,
        raw=line.rstrip("\n"),
    )


def _tokenize_cs(cs: str) -> list[str]:
    tokens: list[str] = []
    cursor = 0
    for match in _CS_TOKEN.finditer(cs):
        if match.start() != cursor:
            raise SequenceComparisonError(
                f"Unsupported minimap2 cs operation near: {cs[cursor:cursor + 20]}"
            )
        tokens.append(match.group(0))
        cursor = match.end()
    if cursor != len(cs):
        raise SequenceComparisonError(
            f"Unsupported minimap2 cs operation near: {cs[cursor:cursor + 20]}"
        )
    return tokens


def parse_cs_variants(
    cs: str,
    target_start: int,
    reference: str,
    circular: bool = False,
) -> tuple[VariantEvent, ...]:
    """Convert a minimap2 `cs` string into reference-coordinate events.

    `target_start` is the zero-based PAF target start. Insertions are
    reported at the reference base immediately before the inserted sequence;
    deletions are reported at the first deleted base.
    """

    ref_cursor = target_start
    ref_len = len(reference)
    events: list[VariantEvent] = []

    for token in _tokenize_cs(cs):
        operation = token[0]
        if operation == ":":
            ref_cursor += int(token[1:])
        elif operation == "=":
            ref_cursor += len(token) - 1
        elif operation == "*":
            pos = normalize_position(ref_cursor + 1, ref_len, circular)
            events.append(
                VariantEvent(
                    position=pos,
                    kind="SNP",
                    ref=token[1].upper(),
                    alt=token[2].upper(),
                    length=1,
                    context=sequence_context(reference, pos, circular=circular),
                )
            )
            ref_cursor += 1
        elif operation == "+":
            inserted = token[1:].upper()
            raw_anchor = ref_cursor if ref_cursor > 0 else 1
            pos = normalize_position(raw_anchor, ref_len, circular)
            events.append(
                VariantEvent(
                    position=pos,
                    kind="Insertion",
                    ref="-",
                    alt=inserted,
                    length=len(inserted),
                    context=sequence_context(reference, pos, circular=circular),
                )
            )
        elif operation == "-":
            deleted = token[1:].upper()
            pos = normalize_position(ref_cursor + 1, ref_len, circular)
            events.append(
                VariantEvent(
                    position=pos,
                    kind="Deletion",
                    ref=deleted,
                    alt="-",
                    length=len(deleted),
                    context=sequence_context(reference, pos, circular=circular),
                )
            )
            ref_cursor += len(deleted)
        elif operation == "~":
            match = re.fullmatch(r"~[A-Za-z]{2}(\d+)[A-Za-z]{2}", token)
            if not match:
                raise SequenceComparisonError(f"Could not parse cs splice operation: {token}")
            length = int(match.group(1))
            pos = normalize_position(ref_cursor + 1, ref_len, circular)
            events.append(
                VariantEvent(
                    position=pos,
                    kind="Deletion",
                    ref=f"<{length} bp>",
                    alt="-",
                    length=length,
                    context=sequence_context(reference, pos, circular=circular),
                    warnings=("Large skipped region",),
                )
            )
            ref_cursor += length
    return tuple(events)


def compare_sequences(
    reference: SequenceRecord,
    query: SequenceRecord,
    circular: bool = True,
    minimap2: str = "minimap2",
    preset: str = "asm5",
) -> AlignmentSummary:
    """Align a query sequence and return exact SNP/indel events.

    For circular references the reference is doubled before alignment, which
    lets a query cross the arbitrary FASTA start/end junction. Returned
    coordinates are normalized back to the original reference length.
    """

    executable = shutil.which(minimap2) if not Path(minimap2).is_file() else minimap2
    if not executable:
        raise SequenceComparisonError(
            "minimap2 was not found. Activate the NGS_env conda environment first."
        )

    target_sequence = reference.sequence * 2 if circular else reference.sequence
    target = SequenceRecord(reference.name, target_sequence)
    with tempfile.TemporaryDirectory(prefix="ont_sequence_compare_") as tmp:
        tmpdir = Path(tmp)
        ref_path = tmpdir / "reference.fasta"
        query_path = tmpdir / "query.fasta"
        write_fasta(ref_path, target)
        write_fasta(query_path, query)
        command = [
            str(executable),
            "-x",
            preset,
            "-c",
            "--cs=long",
            "--secondary=no",
            str(ref_path),
            str(query_path),
        ]
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
        )

    if completed.returncode != 0:
        detail = completed.stderr.strip() or "unknown minimap2 error"
        raise SequenceComparisonError(f"minimap2 failed: {detail}")

    hits: list[PafHit] = []
    for line in completed.stdout.splitlines():
        if line.strip():
            hits.append(parse_paf_line(line))
    if not hits:
        raise SequenceComparisonError(
            "No reliable alignment was found. Check that the sequences are related "
            "and contain enough unambiguous DNA."
        )

    hit = max(hits, key=lambda h: (h.query_end - h.query_start, h.matches, h.mapq))
    variants = parse_cs_variants(
        hit.cs,
        target_start=hit.target_start,
        reference=reference.sequence,
        circular=circular,
    )
    ref_len = len(reference.sequence)
    ref_start = normalize_position(hit.target_start + 1, ref_len, circular)
    ref_end = normalize_position(max(hit.target_start + 1, hit.target_end), ref_len, circular)
    identity = hit.matches / hit.block_length if hit.block_length else 0.0
    query_coverage = (
        (hit.query_end - hit.query_start) / hit.query_length if hit.query_length else 0.0
    )
    return AlignmentSummary(
        orientation="Forward" if hit.strand == "+" else "Reverse complement",
        identity=identity,
        query_coverage=query_coverage,
        mapq=hit.mapq,
        reference_start=ref_start,
        reference_end=ref_end,
        query_start=hit.query_start + 1,
        query_end=hit.query_end,
        query_length=hit.query_length,
        variants=variants,
        raw_paf=hit.raw,
    )
