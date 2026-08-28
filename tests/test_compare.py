import tempfile
import unittest
from pathlib import Path

from ont_ui.compare import compare_sequences, parse_cs_variants, parse_paf_line
from ont_ui.sequences import SequenceRecord


class ComparisonTests(unittest.TestCase):
    def test_cs_snp_insertion_and_deletion_positions(self):
        reference = "ACGTACGTACGTACGTACGT"
        events = parse_cs_variants(":3*ag:2+tt:1-cc:2", 0, reference, circular=False)
        self.assertEqual([event.kind for event in events], ["SNP", "Insertion", "Deletion"])
        self.assertEqual([event.position for event in events], [4, 6, 8])
        self.assertEqual(events[0].ref, "A")
        self.assertEqual(events[0].alt, "G")
        self.assertEqual(events[1].alt, "TT")
        self.assertEqual(events[2].ref, "CC")

    def test_cs_coordinates_are_normalized_for_circular_reference(self):
        events = parse_cs_variants(":3*ac", 8, "ACGTACGTAC", circular=True)
        self.assertEqual(events[0].position, 2)

    def test_parse_paf_line(self):
        line = "query\t8\t0\t8\t-\tref\t16\t6\t14\t7\t8\t60\tcs:Z::3*ag:4"
        hit = parse_paf_line(line)
        self.assertEqual(hit.strand, "-")
        self.assertEqual(hit.cs, ":3*ag:4")
        self.assertEqual(hit.mapq, 60)

    def test_compare_sequences_with_fake_minimap2(self):
        with tempfile.TemporaryDirectory() as tmp:
            fake = Path(tmp) / "minimap2"
            fake.write_text(
                "#!/usr/bin/env bash\n"
                "printf 'query\\t8\\t0\\t8\\t+\\treference\\t8\\t0\\t8\\t8\\t8\\t60\\tcs:Z:=ACGTACGT\\n'\n",
                encoding="utf-8",
            )
            fake.chmod(0o755)
            result = compare_sequences(
                SequenceRecord("reference", "ACGTACGT"),
                SequenceRecord("query", "ACGTACGT"),
                circular=False,
                minimap2=str(fake),
            )
        self.assertEqual(result.identity, 1.0)
        self.assertEqual(result.query_coverage, 1.0)
        self.assertEqual(result.variants, ())


if __name__ == "__main__":
    unittest.main()
