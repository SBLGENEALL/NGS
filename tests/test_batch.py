import io
import json
import tempfile
import unittest
from pathlib import Path

from ont_ui.batch import (
    BatchPreparationError,
    BatchSettings,
    barcode_from_upload_name,
    prepare_batch_job,
    uploaded_barcodes,
)


class Upload(io.BytesIO):
    def __init__(self, name: str, data: bytes):
        super().__init__(data)
        self.name = name
        self.size = len(data)

    def getvalue(self):
        return super().getvalue()


class BatchTests(unittest.TestCase):
    def test_barcode_detection_from_directory_upload(self):
        self.assertEqual(
            barcode_from_upload_name("run/barcode013/fastq_pass/reads.fastq.gz"),
            "barcode13",
        )
        self.assertEqual(
            barcode_from_upload_name("barcode7_reads.fastq.gz"),
            "barcode07",
        )
        self.assertIsNone(barcode_from_upload_name("reads.fastq.gz"))

    def test_uploaded_reads_group_numerically(self):
        uploads = [
            Upload("run/barcode10/b.fastq.gz", b"x"),
            Upload("run/barcode2/a.fastq.gz", b"x"),
        ]
        self.assertEqual(list(uploaded_barcodes(uploads)), ["barcode02", "barcode10"])

    def test_prepare_three_replicates_for_one_reference(self):
        reference = Upload("C-5.pPB-161.fasta", b">ref\nACGTACGT\n")
        reads = [
            Upload(f"run/barcode{number}/reads.fastq", b"@r\nACGT\n+\nIIII\n")
            for number in (13, 14, 15)
        ]
        mappings = [
            {
                "reference": "C-5.pPB-161.fasta",
                "barcodes": ["barcode13", "barcode14", "barcode15"],
            }
        ]
        with tempfile.TemporaryDirectory() as tmp:
            job = prepare_batch_job(
                Path(tmp),
                [reference],
                reads,
                mappings,
                BatchSettings("test", threads=8, parallel_jobs=3),
            )
            manifest = json.loads(job.manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(job.sample_count, 3)
            self.assertEqual(manifest["thresholds"]["min_quality"], 20.0)
            self.assertEqual(manifest["thresholds"]["min_depth"], 10)
            self.assertEqual(manifest["thresholds"]["min_af"], 0.8)
            self.assertEqual(
                [sample["barcode"] for sample in manifest["samples"]],
                ["barcode13", "barcode14", "barcode15"],
            )
            reference_files = list(
                (job.job_dir / "staged" / "references" / "test").glob("*.fasta")
            )
            self.assertEqual(len(reference_files), 3)
            for sample in manifest["samples"]:
                links = list(
                    (job.job_dir / "staged" / "data" / "test" / sample["sample_name"]).iterdir()
                )
                self.assertEqual(len(links), 1)
                self.assertTrue(links[0].is_symlink())

    def test_duplicate_barcode_is_rejected(self):
        references = [
            Upload("a.fasta", b">a\nACGT\n"),
            Upload("b.fasta", b">b\nACGT\n"),
        ]
        reads = [Upload("run/barcode1/reads.fastq", b"@r\nA\n+\nI\n")]
        mappings = [
            {"reference": "a.fasta", "barcodes": ["barcode01"]},
            {"reference": "b.fasta", "barcodes": ["barcode01"]},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(BatchPreparationError):
                prepare_batch_job(
                    Path(tmp), references, reads, mappings, BatchSettings("test")
                )

    def test_variable_sample_counts_per_reference(self):
        references = [
            Upload("a.fasta", b">a\nACGTACGT\n"),
            Upload("b.fasta", b">b\nTGCATGCA\n"),
        ]
        reads = [
            Upload(f"run/barcode{number}/reads.fastq", b"@r\nACGT\n+\nIIII\n")
            for number in (1, 2, 3)
        ]
        mappings = [
            {"reference": "a.fasta", "barcodes": ["barcode01", "barcode02"]},
            {"reference": "b.fasta", "barcodes": ["barcode03"]},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            job = prepare_batch_job(
                Path(tmp), references, reads, mappings, BatchSettings("variable")
            )
            manifest = json.loads(job.manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(job.sample_count, 3)
            self.assertEqual(
                [sample["reference_file"] for sample in manifest["samples"]],
                ["a.fasta", "a.fasta", "b.fasta"],
            )


if __name__ == "__main__":
    unittest.main()
