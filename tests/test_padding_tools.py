import tempfile
import unittest
from pathlib import Path

from msspack.padding_tools import pad_fasta, write_spliced_cds_fasta


class PaddingToolsTests(unittest.TestCase):
    def test_write_spliced_cds_fasta_handles_plus_and_minus_strands(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            fasta = base / "input.fa"
            gff = base / "input.gff3"
            out = base / "out.fa"
            log = base / "out.log"
            fasta.write_text(
                ">chr1\nAAATTTCCCGGG\n>chr2\nATGCCCAAATTT\n",
                encoding="utf-8",
            )
            gff.write_text(
                "\n".join(
                    [
                        "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=tx1",
                        "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tID=tx1.cds1;Parent=tx1",
                        "chr1\tsrc\tCDS\t7\t9\t.\t+\t0\tID=tx1.cds2;Parent=tx1",
                        "chr2\tsrc\tmRNA\t1\t12\t.\t-\t.\tID=tx2",
                        "chr2\tsrc\tCDS\t1\t3\t.\t-\t0\tID=tx2.cds1;Parent=tx2",
                        "chr2\tsrc\tCDS\t10\t12\t.\t-\t0\tID=tx2.cds2;Parent=tx2",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            write_spliced_cds_fasta(
                fasta_path=fasta,
                gff_path=gff,
                output_path=out,
                log_path=log,
            )

            text = out.read_text(encoding="utf-8")
            self.assertIn(">tx1", text)
            self.assertIn("AAACCC", text)
            self.assertIn(">tx2", text)
            self.assertIn("AAACAT", text)

    def test_pad_fasta_reports_tail_padding(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            fasta = base / "input.fa"
            out = base / "out.fa"
            log = base / "out.log"
            fasta.write_text(">tx1\nATGA\n", encoding="utf-8")

            pad_fasta(
                input_path=fasta,
                output_path=out,
                log_path=log,
                genetic_code="1",
            )

            self.assertIn("tail_padding=2", log.read_text(encoding="utf-8"))
            self.assertIn("ATGANN", out.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
