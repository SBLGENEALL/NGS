import csv
import gzip
import io
import unittest
import zipfile

from ont_ui.demo import build_demo_batch_zip


class DemoDatasetTests(unittest.TestCase):
    def test_demo_zip_has_realistic_ont_tree_and_fifteen_samples(self):
        with zipfile.ZipFile(io.BytesIO(build_demo_batch_zip())) as archive:
            names = archive.namelist()
            references = [name for name in names if name.startswith("references/")]
            fastqs = [
                name
                for name in names
                if name.startswith("demo_ont_run/fastq_pass/")
                and name.endswith(".fastq.gz")
            ]
            self.assertEqual(len(references), 5)
            self.assertEqual(len(fastqs), 15)
            self.assertIn("expected_mapping.csv", names)

            mapping = list(
                csv.DictReader(
                    io.StringIO(archive.read("expected_mapping.csv").decode("utf-8"))
                )
            )
            self.assertEqual(len(mapping), 15)
            self.assertEqual(mapping[0]["Sample"], "ONT_sample_01")
            self.assertEqual(mapping[-1]["Sample"], "ONT_sample_15")
            for folder in (
                "fastq_fail",
                "other_reports",
                "pod5_fail",
                "pod5_pass",
                "pod5_skip",
            ):
                self.assertIn(f"demo_ont_run/{folder}/README.txt", names)

    def test_demo_fastq_records_have_matching_sequence_and_quality_lengths(self):
        with zipfile.ZipFile(io.BytesIO(build_demo_batch_zip())) as archive:
            compressed = archive.read(
                "demo_ont_run/fastq_pass/ONT_sample_03/reads_0001.fastq.gz"
            )
            text = gzip.decompress(compressed).decode("utf-8")
            lines = text.splitlines()
            self.assertEqual(len(lines), 12 * 4)
            for index in range(0, len(lines), 4):
                self.assertTrue(lines[index].startswith("@ONT_sample_03_read_"))
                self.assertEqual(lines[index + 2], "+")
                self.assertEqual(len(lines[index + 1]), len(lines[index + 3]))


if __name__ == "__main__":
    unittest.main()
