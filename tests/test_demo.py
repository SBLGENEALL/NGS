import csv
import io
import unittest
import zipfile

from ont_ui.demo import build_demo_batch_zip


class DemoDatasetTests(unittest.TestCase):
    def test_demo_zip_has_five_references_and_fifteen_barcodes(self):
        with zipfile.ZipFile(io.BytesIO(build_demo_batch_zip())) as archive:
            names = archive.namelist()
            references = [name for name in names if name.startswith("references/")]
            fastqs = [
                name
                for name in names
                if name.startswith("demo_reads/") and name.endswith(".fastq")
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
            self.assertEqual(mapping[0]["Barcode"], "barcode01")
            self.assertEqual(mapping[-1]["Barcode"], "barcode15")

    def test_demo_fastq_records_have_matching_sequence_and_quality_lengths(self):
        with zipfile.ZipFile(io.BytesIO(build_demo_batch_zip())) as archive:
            text = archive.read(
                "demo_reads/demo_plasmid_01/barcode03.fastq"
            ).decode("utf-8")
            lines = text.splitlines()
            self.assertEqual(len(lines), 12 * 4)
            for index in range(0, len(lines), 4):
                self.assertTrue(lines[index].startswith("@barcode03_read_"))
                self.assertEqual(lines[index + 2], "+")
                self.assertEqual(len(lines[index + 1]), len(lines[index + 3]))


if __name__ == "__main__":
    unittest.main()
