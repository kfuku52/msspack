import tempfile
import unittest
from pathlib import Path

from msspack.pipeline_actions import copy_input_fasta
from msspack.utils import MSSPackError


class PipelineActionTests(unittest.TestCase):
    def test_copy_input_fasta_rejects_duplicate_ids(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            source = base / "input.fa"
            source.write_text(">same\nATG\n>same\nTAA\n", encoding="utf-8")

            with self.assertRaisesRegex(MSSPackError, "Duplicate FASTA"):
                copy_input_fasta(
                    input_path=source,
                    output_path=base / "output.fa",
                    log_path=base / "copy.log",
                )
            self.assertFalse((base / "output.fa").exists())

    def test_copy_input_fasta_rejects_empty_records(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            source = base / "input.fa"
            source.write_text(">empty\n", encoding="utf-8")

            with self.assertRaisesRegex(MSSPackError, "sequence data"):
                copy_input_fasta(
                    input_path=source,
                    output_path=base / "output.fa",
                    log_path=base / "copy.log",
                )
            self.assertFalse((base / "output.fa").exists())

    def test_copy_input_fasta_rejects_a_file_without_records(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            source = base / "input.fa"
            source.write_text("\n", encoding="utf-8")

            with self.assertRaisesRegex(MSSPackError, "any sequence records"):
                copy_input_fasta(
                    input_path=source,
                    output_path=base / "output.fa",
                    log_path=base / "copy.log",
                )
            self.assertFalse((base / "output.fa").exists())


if __name__ == "__main__":
    unittest.main()
