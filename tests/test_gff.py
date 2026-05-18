import tempfile
import unittest
from pathlib import Path

from msspack.gff import sort_gff_file_precise


class GffSortTests(unittest.TestCase):
    def test_precise_sort_reorders_same_start_parent_before_child(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            source = base / "input.gff3"
            output = base / "output.gff3"
            source.write_text(
                "\n".join(
                    [
                        "##gff-version 3",
                        "chr2\tsrc\tgene\t5\t9\t.\t+\t.\tID=other",
                        "chr1\tsrc\texon\t1\t3\t.\t+\t.\tID=tx1.ex1;Parent=tx1",
                        "chr1\tsrc\tgene\t1\t20\t.\t+\t.\tID=g1",
                        "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tID=tx1.cds;Parent=tx1",
                        "chr1\tsrc\tmRNA\t1\t20\t.\t+\t.\tID=tx1;Parent=g1",
                        "chr1\tsrc\texon\t10\t12\t.\t+\t.\tID=tx1.ex2;Parent=tx1",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            sort_gff_file_precise(source, output)

            self.assertEqual(
                output.read_text(encoding="utf-8").splitlines(),
                [
                    "##gff-version 3",
                    "chr1\tsrc\tgene\t1\t20\t.\t+\t.\tID=g1",
                    "chr1\tsrc\tmRNA\t1\t20\t.\t+\t.\tID=tx1;Parent=g1",
                    "chr1\tsrc\texon\t1\t3\t.\t+\t.\tID=tx1.ex1;Parent=tx1",
                    "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tID=tx1.cds;Parent=tx1",
                    "chr1\tsrc\texon\t10\t12\t.\t+\t.\tID=tx1.ex2;Parent=tx1",
                    "chr2\tsrc\tgene\t5\t9\t.\t+\t.\tID=other",
                ],
            )

    def test_precise_sort_preserves_same_start_sibling_input_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            source = base / "input.gff3"
            output = base / "output.gff3"
            source.write_text(
                "\n".join(
                    [
                        "chr1\tsrc\tgene\t1\t20\t.\t+\t.\tID=g1",
                        "chr1\tsrc\tmRNA\t1\t20\t.\t+\t.\tID=tx1;Parent=g1",
                        "chr1\tsrc\texon\t1\t3\t.\t+\t.\tID=tx1.ex1;Parent=tx1",
                        "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tID=tx1.cds;Parent=tx1",
                        "chr1\tsrc\tCDS\t10\t12\t.\t+\t0\tID=tx1.cds;Parent=tx1",
                        "chr1\tsrc\texon\t10\t12\t.\t+\t.\tID=tx1.ex2;Parent=tx1",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            sort_gff_file_precise(source, output)

            self.assertEqual(
                output.read_text(encoding="utf-8").splitlines(),
                [
                    "chr1\tsrc\tgene\t1\t20\t.\t+\t.\tID=g1",
                    "chr1\tsrc\tmRNA\t1\t20\t.\t+\t.\tID=tx1;Parent=g1",
                    "chr1\tsrc\texon\t1\t3\t.\t+\t.\tID=tx1.ex1;Parent=tx1",
                    "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tID=tx1.cds;Parent=tx1",
                    "chr1\tsrc\tCDS\t10\t12\t.\t+\t0\tID=tx1.cds;Parent=tx1",
                    "chr1\tsrc\texon\t10\t12\t.\t+\t.\tID=tx1.ex2;Parent=tx1",
                ],
            )


if __name__ == "__main__":
    unittest.main()
