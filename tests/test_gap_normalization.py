import tempfile
import unittest
from pathlib import Path

from msspack.gap_normalization import normalize_gap_lengths


class GapNormalizationTests(unittest.TestCase):
    def test_gapjust_updates_gff_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            fasta = base / "input.fa"
            gff = base / "input.gff3"
            out_fasta = base / "out.fa"
            out_gff = base / "out.gff3"
            log = base / "out.log"
            fasta.write_text(">chr1\nAAANNNNAAA\n", encoding="utf-8")
            gff.write_text(
                "chr1\tsrc\tgene\t8\t10\t.\t+\t.\tID=g1\n",
                encoding="utf-8",
            )

            normalize_gap_lengths(
                fasta_path=fasta,
                output_fasta_path=out_fasta,
                log_path=log,
                gap_len=2,
                gap_just_min=0,
                gap_just_max=10,
                input_gff_path=gff,
                output_gff_path=out_gff,
            )

            self.assertIn("AAANNAAA", out_fasta.read_text(encoding="utf-8"))
            self.assertIn("\t6\t8\t", out_gff.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
