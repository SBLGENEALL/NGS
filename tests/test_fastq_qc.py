import gzip
import json
import tempfile
import unittest
from pathlib import Path

from ont_ui.fastq_qc import load_fastq_qc, summarize_fastq


class FastqQcTests(unittest.TestCase):
    def test_summarize_fastq_reports_read_level_qscore_fractions(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "reads.fastq.gz"
            records = (
                "@q40\nAAAA\n+\nIIII\n"
                "@q20\nAAAAAA\n+\n555555\n"
                "@q10\nAA\n+\n++\n"
            )
            path.write_bytes(gzip.compress(records.encode("utf-8")))

            result = summarize_fastq(path)

            self.assertEqual(result["reads"], 3)
            self.assertEqual(result["total_bases"], 12)
            self.assertEqual(result["n50"], 6)
            self.assertAlmostEqual(result["mean_read_quality"], 70 / 3)
            self.assertAlmostEqual(result["q10_read_fraction"], 1.0)
            self.assertAlmostEqual(result["q20_read_fraction"], 2 / 3)
            self.assertAlmostEqual(result["q30_read_fraction"], 1 / 3)

    def test_load_fastq_qc_returns_empty_dict_for_invalid_or_missing_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "qc.json"
            self.assertEqual(load_fastq_qc(path), {})
            path.write_text("not-json", encoding="utf-8")
            self.assertEqual(load_fastq_qc(path), {})
            path.write_text(json.dumps({"merged": {"reads": 2}}), encoding="utf-8")
            self.assertEqual(load_fastq_qc(path)["merged"]["reads"], 2)


if __name__ == "__main__":
    unittest.main()
