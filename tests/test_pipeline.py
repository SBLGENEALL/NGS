import io
import tempfile
import unittest
from pathlib import Path

from ont_ui.pipeline import (
    PipelinePreparationError,
    RawAnalysisSettings,
    find_fastq_files,
    prepare_job,
)
from ont_ui.sequences import SequenceRecord


class PipelineStagingTests(unittest.TestCase):
    def test_find_fastq_files_recursively(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            nested = root / "barcode01"
            nested.mkdir()
            (nested / "reads.fastq.gz").touch()
            (nested / "ignore.txt").touch()
            files = find_fastq_files(root)
            self.assertEqual([path.name for path in files], ["reads.fastq.gz"])

    def test_prepare_job_uses_isolated_absolute_paths_and_symlinks(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            reads = root / "input.fastq.gz"
            reads.touch()
            job = prepare_job(
                root / "runs",
                SequenceRecord("ref", "ACGTACGT"),
                RawAnalysisSettings(experiment_name="test exp", sample_name="sample/01"),
                read_paths=[reads],
            )
            self.assertEqual(job.experiment_name, "test_exp")
            self.assertEqual(job.sample_name, "sample_01")
            staged = list((job.data_root / "test_exp" / "sample_01").iterdir())
            self.assertEqual(len(staged), 1)
            self.assertTrue(staged[0].is_symlink())
            config = job.config_path.read_text(encoding="utf-8")
            self.assertIn(str(job.references_root.resolve()), config)
            self.assertIn("min_read_length: 500", config)
            self.assertIn("min_read_quality: 10", config)

    def test_prepare_uploaded_read(self):
        with tempfile.TemporaryDirectory() as tmp:
            job = prepare_job(
                Path(tmp) / "runs",
                SequenceRecord("ref", "ACGT"),
                RawAnalysisSettings(experiment_name="exp", sample_name="sample"),
                uploaded_reads=[("reads.fastq", io.BytesIO(b"@r1\nACGT\n+\nIIII\n"))],
            )
            staged = list((job.data_root / "exp" / "sample").iterdir())
            self.assertEqual(len(staged), 1)
            self.assertEqual(staged[0].resolve().read_bytes(), b"@r1\nACGT\n+\nIIII\n")

    def test_prepare_job_requires_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(PipelinePreparationError):
                prepare_job(
                    Path(tmp),
                    SequenceRecord("ref", "ACGT"),
                    RawAnalysisSettings(experiment_name="exp", sample_name="sample"),
                )


if __name__ == "__main__":
    unittest.main()
