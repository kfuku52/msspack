import csv
import tempfile
import unittest
from pathlib import Path

from msspack.submission_review import omit_redundant_mrna, write_submission_review


class SubmissionReviewTests(unittest.TestCase):
    def test_redundant_mrna_omission_preserves_transcript_evidence_and_partial_cds(self) -> None:
        source = "ctg\tsource\t1..900\n"
        mrna = "\tmRNA\tcomplement(join(1..40,61..99))\tlocus_tag\tTest_1\n"
        cds = "\tCDS\tcomplement(join(1..40,61..>99))\tlocus_tag\tTest_1\n"
        others = (
            "\tmRNA\t101..110\tlocus_tag\tTest_2\n"
            "\tmisc_feature\t101..110\tlocus_tag\tTest_2\n"
            "\tmRNA\t201..320\tlocus_tag\tTest_3\n"
            "\t\t\tnote\tRNA evidence\n"
            "\tCDS\t201..320\tlocus_tag\tTest_3\n"
            "\tmRNA\t401..520\tlocus_tag\tTest_4\n"
            "\tCDS\t404..520\tlocus_tag\tTest_4\n"
            "\tmRNA\t601..720\tlocus_tag\tTest_5\n"
            "\tmRNA\t610..720\tlocus_tag\tTest_5\n"
            "\tCDS\t601..720\tlocus_tag\tTest_5\n"
            "other\tsource\t1..900\n"
            "\tmRNA\tcomplement(join(1..40,61..99))\tlocus_tag\tTest_1\n"
        )
        result, count = omit_redundant_mrna(source + mrna + cds + others)
        self.assertEqual(count, 1)
        self.assertEqual(result, source + cds + others)

    def test_advisory_review_identifies_specific_loci_and_strands(self) -> None:
        ann = (
            "COMMON\tDATATYPE\t\ttype\tWGS\n"
            "ctg1\tsource\t1..1000\n"
            "\tCDS\tjoin(9,21..70)\tlocus_tag\tDemo_1\n"
            "\tCDS\tjoin(21..70,101..200)\tlocus_tag\tDemo_2\n"
            "\tCDS\tcomplement(join(21..70,301..400))\tlocus_tag\tDemo_3\n"
            "\tCDS\t501..>576\tlocus_tag\tDemo_4\n"
            "\ttRNA\t601..675\tlocus_tag\tDemo_Ile\n"
            "\t\t\tgene\ttrnI-CAU\n"
            "\t\t\tproduct\ttRNA-Met\n"
            "\ttRNA\t701..775\tlocus_tag\tDemo_fMet\n"
            "\t\t\tgene\ttrnfM-CAU\n"
            "\t\t\tproduct\ttRNA-Met\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path, review = Path(directory) / "test.ann", Path(directory) / "warnings.tsv"
            path.write_text(ann)
            self.assertEqual(write_submission_review(path, review), 6)
            with review.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(path.read_text(), ann)
        self.assertEqual(
            {(r["warning"], r["locus_tag"], r["location"]) for r in rows},
            {
                ("missing_genome_coverage", "", ""),
                ("short_cds", "Demo_1", "join(9,21..70)"),
                ("short_cds", "Demo_4", "501..>576"),
                ("single_base_cds_segment", "Demo_1", "9"),
                ("shared_cds_segment", "Demo_1,Demo_2", "21..70"),
                ("trna_gene_product_mismatch", "Demo_Ile", "601..675"),
            },
        )
        self.assertIn("partial CDS", next(r["detail"] for r in rows if r["locus_tag"] == "Demo_4"))
