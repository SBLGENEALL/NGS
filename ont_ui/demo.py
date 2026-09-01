"""Generate a small synthetic batch dataset for UI and pipeline smoke tests."""

from __future__ import annotations

import csv
import io
import random
import zipfile
from functools import lru_cache


REFERENCE_COUNT = 5
SAMPLES_PER_REFERENCE = 3
READS_PER_SAMPLE = 12
REFERENCE_LENGTH = 1400


def _wrapped_fasta(name: str, sequence: str) -> str:
    lines = [sequence[index : index + 70] for index in range(0, len(sequence), 70)]
    return f">{name}\n" + "\n".join(lines) + "\n"


def _alternate_base(base: str) -> str:
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[base]


def _sample_sequence(reference: str, sample_within_reference: int) -> tuple[str, str]:
    if sample_within_reference == 1:
        return reference, "Matches reference"
    if sample_within_reference == 2:
        position = 301
        sequence = (
            reference[: position - 1]
            + _alternate_base(reference[position - 1])
            + reference[position:]
        )
        return sequence, f"Synthetic SNP at reference position {position}"
    position = 701
    return (
        reference[:position] + "A" + reference[position:],
        f"Synthetic 1-bp insertion after reference position {position}",
    )


def _fastq(barcode: str, sequence: str) -> str:
    quality = "I" * len(sequence)
    records = []
    for read_index in range(1, READS_PER_SAMPLE + 1):
        records.extend(
            [
                f"@{barcode}_read_{read_index:02d}",
                sequence,
                "+",
                quality,
            ]
        )
    return "\n".join(records) + "\n"


@lru_cache(maxsize=1)
def build_demo_batch_zip() -> bytes:
    """Return a ZIP containing five references and fifteen valid barcode FASTQs."""
    output = io.BytesIO()
    mapping_buffer = io.StringIO()
    mapping_writer = csv.writer(mapping_buffer, lineterminator="\n")
    mapping_writer.writerow(["Barcode", "Reference", "Expected result"])

    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        barcode_number = 1
        for reference_index in range(1, REFERENCE_COUNT + 1):
            rng = random.Random(20260902 + reference_index)
            reference = "".join(rng.choice("ACGT") for _ in range(REFERENCE_LENGTH))
            reference_name = f"demo_plasmid_{reference_index:02d}"
            archive.writestr(
                f"references/{reference_name}.fasta",
                _wrapped_fasta(reference_name, reference),
            )
            for sample_index in range(1, SAMPLES_PER_REFERENCE + 1):
                barcode = f"barcode{barcode_number:02d}"
                sequence, expected = _sample_sequence(reference, sample_index)
                archive.writestr(
                    f"demo_reads/{barcode}/reads.fastq",
                    _fastq(barcode, sequence),
                )
                mapping_writer.writerow(
                    [barcode, f"{reference_name}.fasta", expected]
                )
                barcode_number += 1

        archive.writestr("expected_mapping.csv", mapping_buffer.getvalue())
        archive.writestr(
            "README.txt",
            "ONT Plasmid Analyzer synthetic demo\n"
            "===================================\n\n"
            "1. Set Number of samples to 15.\n"
            "2. Upload all five FASTA files in references/.\n"
            "3. Select demo_reads/ in the barcode directory uploader.\n"
            "4. Review the automatically balanced mapping against expected_mapping.csv.\n"
            "5. Run the batch analysis.\n\n"
            "Each reference has three example samples: exact, SNP, and 1-bp insertion.\n"
            "These synthetic sequences are for software testing only.\n",
        )
    return output.getvalue()
