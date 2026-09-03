import gzip
import tempfile
import unittest
from pathlib import Path

from ont_ui.results import (
    parse_consensus_fallback_variants,
    parse_flagstat,
    parse_vcf_variants,
    read_depth,
    variants_csv,
)
from ont_ui.sequences import SequenceRecord


class ResultParsingTests(unittest.TestCase):
    def test_vcf_quality_depth_and_allele_fraction(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "calls.vcf.gz"
            content = (
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
                "ref\t5\t.\tA\tG\t40\tPASS\tDP=30\tGT:DP:AD\t1:30:2,28\n"
                "ref\t10\t.\tA\tATT\t15\tPASS\t.\tGT:DP:AD\t1:12:3,9\n"
                "ref\t20\t.\tATC\tA\t50\tPASS\tDP=25\tGT:DP:AD\t1:25:1,24\n"
            )
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write(content)
            events = parse_vcf_variants(
                path,
                reference_sequence="ACGT" * 10,
                min_quality=20,
                min_depth=10,
                min_allele_fraction=0.8,
                edge_margin=0,
                circular=False,
                homopolymer_threshold=10,
            )

        self.assertEqual([event.kind for event in events], ["SNP", "Insertion", "Deletion"])
        self.assertEqual(events[0].status, "PASS")
        self.assertAlmostEqual(events[0].allele_fraction, 28 / 30)
        self.assertEqual(events[1].position, 10)
        self.assertEqual(events[1].alt, "TT")
        self.assertEqual(events[1].status, "REVIEW")
        self.assertIn("Low QUAL (<20)", events[1].warnings)
        self.assertEqual(events[2].position, 21)
        self.assertEqual(events[2].ref, "TC")

    def test_flagstat_and_depth(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            flagstat = tmp_path / "sample.flagstat.txt"
            flagstat.write_text(
                "100 + 0 in total (QC-passed reads + QC-failed reads)\n"
                "80 + 0 mapped (80.00% : N/A)\n",
                encoding="utf-8",
            )
            mapping = parse_flagstat(flagstat)
            self.assertEqual(mapping.total_reads, 100)
            self.assertEqual(mapping.mapped_reads, 80)
            self.assertEqual(mapping.mapping_rate, 0.8)

            depth_path = tmp_path / "sample.depth.txt"
            depth_path.write_text(
                "ref\t1\t0\nref\t2\t5\nref\t3\t10\nref\t4\t20\n",
                encoding="utf-8",
            )
            depth, points = read_depth(depth_path, max_chart_points=2)
            self.assertEqual(depth.positions, 4)
            self.assertEqual(depth.mean_depth, 8.75)
            self.assertEqual(depth.coverage_1x, 0.75)
            self.assertEqual(depth.coverage_10x, 0.5)
            self.assertEqual(points[-1], (4, 20))

    def test_variant_csv_has_stable_headers(self):
        text = variants_csv([])
        self.assertTrue(text.startswith("Position,Type,REF,ALT"))

    def test_consensus_fallback_marks_differences_for_review(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            consensus = root / "consensus.fasta"
            consensus.write_text(">query\nACGGACGT\n", encoding="utf-8")
            fake = root / "minimap2"
            fake.write_text(
                "#!/usr/bin/env bash\n"
                "printf 'query\\t8\\t0\\t8\\t+\\treference\\t8\\t0\\t8\\t7\\t8\\t60\\tcs:Z::3*ta:4\\n'\n",
                encoding="utf-8",
            )
            fake.chmod(0o755)
            events = parse_consensus_fallback_variants(
                consensus,
                SequenceRecord("reference", "ACGTACGT"),
                circular=False,
                minimap2=str(fake),
            )

        self.assertEqual(len(events), 1)
        self.assertEqual(events[0].kind, "SNP")
        self.assertEqual(events[0].status, "REVIEW")
        self.assertIn("Consensus-only call", events[0].warnings[0])


if __name__ == "__main__":
    unittest.main()
