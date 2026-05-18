import tempfile
import unittest
from pathlib import Path

from msspack.transcript_selection import select_one_mrna_per_gene

_GFF = """\
chr1	src	gene	100	400	.	+	.	ID=Gene1
chr1	src	mRNA	100	400	.	+	.	ID=Gene1-T1;Parent=Gene1
chr1	src	mRNA	100	400	.	+	.	ID=Gene1-T2;Parent=Gene1
chr1	src	exon	100	261	.	+	.	ID=Gene1-T1.exon1;Parent=Gene1-T1
chr1	src	CDS	120	221	.	+	0	ID=Gene1-T1.cds1;Parent=Gene1-T1
chr1	src	five_prime_UTR	100	119	.	+	.	ID=Gene1-T1.utr5;Parent=Gene1-T1
chr1	src	three_prime_UTR	222	261	.	+	.	ID=Gene1-T1.utr3;Parent=Gene1-T1
chr1	src	exon	100	200	.	+	.	ID=Gene1-T2.exon1;Parent=Gene1-T2
chr1	src	exon	320	360	.	+	.	ID=Gene1-T2.exon2;Parent=Gene1-T2
chr1	src	CDS	120	180	.	+	0	ID=Gene1-T2.cds1;Parent=Gene1-T2
chr1	src	CDS	320	360	.	+	0	ID=Gene1-T2.cds2;Parent=Gene1-T2
chr1	src	five_prime_UTR	100	119	.	+	.	ID=Gene1-T2.utr5;Parent=Gene1-T2
chr1	src	three_prime_UTR	181	200	.	+	.	ID=Gene1-T2.utr3a;Parent=Gene1-T2
chr1	src	three_prime_UTR	361	380	.	+	.	ID=Gene1-T2.utr3b;Parent=Gene1-T2
"""


class SelectOneMrnaTests(unittest.TestCase):
    def test_exact_tie_keeps_first_transcript(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            input_path = Path(tmp_dir) / "input.gff"
            output_path = Path(tmp_dir) / "output.gff"
            input_path.write_text(_GFF, encoding="utf-8")

            select_one_mrna_per_gene(
                input_gff_path=input_path,
                output_gff_path=output_path,
            )

            output = output_path.read_text(encoding="utf-8")
            self.assertIn("ID=Gene1-T1;Parent=Gene1", output)
            self.assertNotIn("ID=Gene1-T2;Parent=Gene1", output)

    def test_child_features_before_mrna_are_still_resolved(self) -> None:
        reordered = """\
chr1\tsrc\tgene\t100\t400\t.\t+\t.\tID=Gene1
chr1\tsrc\tCDS\t120\t221\t.\t+\t0\tID=Gene1-T1.cds1;Parent=Gene1-T1
chr1\tsrc\texon\t100\t261\t.\t+\t.\tID=Gene1-T1.exon1;Parent=Gene1-T1
chr1\tsrc\tmRNA\t100\t400\t.\t+\t.\tID=Gene1-T1;Parent=Gene1
"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            input_path = Path(tmp_dir) / "input.gff"
            output_path = Path(tmp_dir) / "output.gff"
            input_path.write_text(reordered, encoding="utf-8")

            select_one_mrna_per_gene(
                input_gff_path=input_path,
                output_gff_path=output_path,
            )

            output = output_path.read_text(encoding="utf-8")
            self.assertIn("ID=Gene1-T1;Parent=Gene1", output)
            self.assertIn("ID=Gene1-T1.cds1;Parent=Gene1-T1", output)
            self.assertIn("ID=Gene1-T1.exon1;Parent=Gene1-T1", output)

    def test_exact_tie_is_stable_when_mrna_input_order_flips(self) -> None:
        shuffled = _GFF.replace(
            "chr1\tsrc\tmRNA\t100\t400\t.\t+\t.\tID=Gene1-T1;Parent=Gene1\n"
            "chr1\tsrc\tmRNA\t100\t400\t.\t+\t.\tID=Gene1-T2;Parent=Gene1\n",
            "chr1\tsrc\tmRNA\t100\t400\t.\t+\t.\tID=Gene1-T2;Parent=Gene1\n"
            "chr1\tsrc\tmRNA\t100\t400\t.\t+\t.\tID=Gene1-T1;Parent=Gene1\n",
        )
        with tempfile.TemporaryDirectory() as tmp_dir:
            input_path = Path(tmp_dir) / "input.gff"
            output_path = Path(tmp_dir) / "output.gff"
            input_path.write_text(shuffled, encoding="utf-8")

            select_one_mrna_per_gene(
                input_gff_path=input_path,
                output_gff_path=output_path,
            )

            output = output_path.read_text(encoding="utf-8")
            self.assertIn("ID=Gene1-T1;Parent=Gene1", output)
            self.assertNotIn("ID=Gene1-T2;Parent=Gene1", output)


if __name__ == "__main__":
    unittest.main()
