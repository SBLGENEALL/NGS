import io
import json
import tempfile
import unittest
from pathlib import Path

from ont_ui.batch import (
    BatchPreparationError,
    BatchSettings,
    barcode_from_upload_name,
    natural_key,
    prepare_batch_job,
    server_fastq_uploads,
    server_reference_uploads,
    uploaded_barcodes,
    uploaded_samples,
)


class Upload(io.BytesIO):
    def __init__(self, name: str, data: bytes):
        super().__init__(data)
        self.name = name
        self.size = len(data)

    def getvalue(self):
        return super().getvalue()


class BatchTests(unittest.TestCase):
    def test_server_paths_load_references_and_group_sample_folders(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            references = root / "references"
            reads = root / "reads"
            references.mkdir()
            (references / "vector10.fasta").write_text(">v10\nACGT\n")
            (references / "vector2.fasta").write_text(">v2\nACGT\n")
            (reads / "Clone_B").mkdir(parents=True)
            (reads / "Clone_A").mkdir(parents=True)
            (reads / "Clone_B" / "reads.fastq").write_text("@r\nA\n+\nI\n")
            (reads / "Clone_A" / "reads.fastq.gz").write_bytes(b"test")

            reference_uploads = server_reference_uploads(str(references))
            fastq_uploads = server_fastq_uploads(str(reads))

            self.assertEqual(
                [Path(item.name).name for item in reference_uploads],
                ["vector2.fasta", "vector10.fasta"],
            )
            self.assertEqual(list(uploaded_samples(fastq_uploads)), ["Clone_A", "Clone_B"])

    def test_server_fastq_is_linked_without_copying_large_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reference_path = root / "vector.fasta"
            reads_dir = root / "Clone_A"
            reads_dir.mkdir()
            reference_path.write_text(">ref\nACGTACGT\n")
            reads_path = reads_dir / "reads.fastq"
            reads_path.write_text("@r\nACGT\n+\nIIII\n")
            references = server_reference_uploads(str(reference_path))
            reads = server_fastq_uploads(str(reads_dir))
            mappings = [
                {
                    "reference": "vector.fasta",
                    "samples": [{"name": "Clone_A", "files": reads}],
                }
            ]
            run_root = root / "runs"
            job = prepare_batch_job(
                run_root, references, reads, mappings, BatchSettings("server")
            )
            staged = next((job.job_dir / "staged" / "data" / "server").rglob("*.fastq"))
            self.assertTrue(staged.is_symlink())
            self.assertEqual(staged.resolve(), reads_path.resolve())

    def test_reference_names_use_natural_order(self):
        names = ["demo_plasmid_10.fasta", "demo_plasmid_02.fasta", "demo_plasmid_01.fasta"]
        self.assertEqual(
            sorted(names, key=natural_key),
            ["demo_plasmid_01.fasta", "demo_plasmid_02.fasta", "demo_plasmid_10.fasta"],
        )

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

    def test_uploaded_samples_accept_alias_folder_and_filename(self):
        uploads = [
            Upload("Clone_A/fastq_pass/chunk01.fastq.gz", b"x"),
            Upload("Clone_A/fastq_pass/chunk02.fastq.gz", b"x"),
            Upload("Custom_plasmid.fastq", b"x"),
        ]
        grouped = uploaded_samples(uploads)
        self.assertEqual(list(grouped), ["Clone_A", "Custom_plasmid"])
        self.assertEqual(len(grouped["Clone_A"]), 2)

    def test_generic_reads_are_not_silently_merged(self):
        grouped = uploaded_samples(
            [Upload("reads.fastq", b"x"), Upload("reads.fastq", b"y")]
        )
        self.assertEqual(list(grouped), ["sample_01", "sample_02"])

    def test_prepare_custom_named_sample(self):
        reference = Upload("vector.fasta", b">ref\nACGTACGT\n")
        reads = [Upload("Clone_A.fastq", b"@r\nACGT\n+\nIIII\n")]
        mappings = [
            {
                "reference": "vector.fasta",
                "samples": [{"name": "Clone_A", "files": reads}],
            }
        ]
        with tempfile.TemporaryDirectory() as tmp:
            job = prepare_batch_job(
                Path(tmp), [reference], reads, mappings, BatchSettings("custom")
            )
            manifest = json.loads(job.manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(manifest["samples"][0]["sample_id"], "Clone_A")
            self.assertIn("Clone_A", manifest["samples"][0]["sample_name"])

    def test_same_sample_can_be_compared_to_multiple_references(self):
        references = [
            Upload("reference_a.fasta", b">a\nACGTACGT\n"),
            Upload("reference_b.fasta", b">b\nACGTACGT\n"),
        ]
        reads = [Upload("Clone_A/reads.fastq", b"@r\nACGT\n+\nIIII\n")]
        mappings = [
            {
                "reference": reference.name,
                "samples": [{"name": "Clone_A", "files": reads}],
            }
            for reference in references
        ]
        with tempfile.TemporaryDirectory() as tmp:
            job = prepare_batch_job(
                Path(tmp), references, reads, mappings, BatchSettings("multi_ref")
            )
            manifest = json.loads(job.manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(job.sample_count, 2)
            self.assertEqual(
                [sample["reference_file"] for sample in manifest["samples"]],
                ["reference_a.fasta", "reference_b.fasta"],
            )
            self.assertEqual(
                [sample["sample_id"] for sample in manifest["samples"]],
                ["Clone_A", "Clone_A"],
            )

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
