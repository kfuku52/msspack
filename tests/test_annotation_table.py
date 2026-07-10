import tempfile
import unittest
from pathlib import Path

from msspack.annotation_table import build_annotation_table


class AnnotationTableTests(unittest.TestCase):
    def test_uses_cds_product_when_mrna_has_no_product(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            gff = base / "input.gff3"
            output = base / "annotation.tsv"
            gff.write_text(
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=g1\n"
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n"
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=t1;product=real%20protein\n",
                encoding="utf-8",
            )

            build_annotation_table(
                gff_path=gff,
                output_path=output,
                locus_tag_prefix="X",
                log_path=base / "annotation.log",
            )

            self.assertIn("t1\treal protein\tX_g1", output.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
