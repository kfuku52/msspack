import tempfile
import unittest
from pathlib import Path

from msspack.pipeline_actions import copy_input_fasta
from msspack.utils import MSSPackError


class PipelineActionTests(unittest.TestCase):
    def test_invalid_fasta_is_rejected_without_publishing_output(self) -> None:
        cases = (
            (">same\nATG\n>same\nTAA\n", "Duplicate FASTA"),
            (">empty\n", "sequence data"),
            ("\n", "any sequence records"),
        )
        for content, message in cases:
            with self.subTest(content=content), tempfile.TemporaryDirectory() as tmp_dir:
                base = Path(tmp_dir)
                source = base / "input.fa"
                source.write_text(content, encoding="utf-8")
                with self.assertRaisesRegex(MSSPackError, message):
                    copy_input_fasta(
                        input_path=source,
                        output_path=base / "output.fa",
                        log_path=base / "copy.log",
                    )
                self.assertFalse((base / "output.fa").exists())


if __name__ == "__main__":
    unittest.main()
