import tempfile
import unittest
from pathlib import Path

from msspack.annotation_table import build_annotation_table


class AnnotationTableTests(unittest.TestCase):
    def test_root_transcripts_and_gene_direct_cds_have_annotation_rows(self) -> None:
        for parent_type in ("mRNA", "transcript", "gene"):
            with self.subTest(parent_type=parent_type), tempfile.TemporaryDirectory() as tmp:
                base = Path(tmp)
                gff = base / "input.gff3"
                # The CDS may precede its parent. Its product is still retained.
                gff.write_text(
                    "chr1\t.\tCDS\t1\t9\t.\t+\t0\tParent=t1;product=real%20protein\n"
                    f"chr1\t.\t{parent_type}\t1\t9\t.\t+\t.\tID=t1\n"
                )
                build_annotation_table(
                    gff_path=gff, output_path=base / "annotation.tsv",
                    locus_tag_prefix="X", log_path=base / "log",
                )
                self.assertEqual((base / "annotation.tsv").read_text(),
                                 "ID\tDescription\tLocus_tag\nt1\treal protein\tX_t1\n")

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
