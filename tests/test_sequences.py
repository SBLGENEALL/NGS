import unittest

from ont_ui.sequences import (
    SequenceValidationError,
    max_homopolymer_run_near,
    parse_single_sequence,
    sanitize_name,
    sequence_context,
)


class SequenceTests(unittest.TestCase):
    def test_plain_sequence_is_normalized(self):
        record = parse_single_sequence("acgu\nNN", default_name="query")
        self.assertEqual(record.name, "query")
        self.assertEqual(record.sequence, "ACGTNN")

    def test_fasta_header_sets_name(self):
        record = parse_single_sequence(">vector one\nACGT\n")
        self.assertEqual(record.name, "vector")
        self.assertEqual(record.sequence, "ACGT")

    def test_multiple_records_are_rejected(self):
        with self.assertRaises(SequenceValidationError):
            parse_single_sequence(">one\nACGT\n>two\nTGCA\n")

    def test_invalid_character_is_rejected(self):
        with self.assertRaises(SequenceValidationError):
            parse_single_sequence("ACGT-Z")

    def test_safe_name_blocks_path_components(self):
        self.assertEqual(sanitize_name("../../my vector"), "my_vector")

    def test_circular_context_wraps(self):
        context = sequence_context("ACGT", 1, radius=2, circular=True)
        self.assertEqual(context, "GT[A]CG")

    def test_homopolymer_run(self):
        run = max_homopolymer_run_near("ACGAAAAATGC", 6, circular=False)
        self.assertEqual(run, 5)


if __name__ == "__main__":
    unittest.main()
