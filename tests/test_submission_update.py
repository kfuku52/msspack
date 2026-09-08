import json
import tempfile
import unittest
from pathlib import Path

from msspack.cli import main
from msspack.submission_update import prepare_update
from msspack.utils import MSSPackError


class SubmissionUpdateTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.base = Path(self.tmp.name)
        self.ann = self.base / "input.ann.txt"
        self.fasta = self.base / "input.fasta"
        self.mapping = self.base / "mapping.tsv"
        self.output = self.base / "update"
        self.ann.write_text(
            "COMMON\tDBLINK\t\tproject\tPRJDB1\n"
            "scaffold1\tsource\t1..12\torganism\tTest species\n"
            "\t\t\tsubmitter_seqid\t@@[entry]@@\n"
            "\tmRNA\tjoin(1..6,10..12)\tlocus_tag\tTest_1\n"
            "\texon\t1..6\tlocus_tag\tTest_1\n"
            "\t\t\tnote\tremove this exon qualifier too\n"
            "\t5'UTR\t1..3\tlocus_tag\tTest_1\n"
            "\tCDS\tjoin(4..6,10..12)\tlocus_tag\tTest_1\n"
            "\t\t\tprotein_id\tGAB2236588.1\n"
            "\t\t\tproduct\thypothetical protein\n"
            "\texon\t7..9\tlocus_tag\tTest_1\n"
            "\texon\tcomplement(1..6)\tlocus_tag\tTest_1\n"
            "\t3'UTR\t10..12\tlocus_tag\tOther_1\n",
            encoding="utf-8",
        )
        self.fasta.write_text(">scaffold1 original header\nacgTNN\nACGTNN\n//\n", encoding="utf-8")
        self.mapping.write_text(
            "entry\taccession\tsubmitter_seqid\nscaffold1\tBAAGII010000001\tscaffold1\n",
            encoding="utf-8",
        )

    def run_update(self, **kwargs: bool) -> Path:
        return prepare_update(
            ann_path=self.ann, fasta_path=self.fasta, mapping_path=self.mapping,
            output_dir=self.output, **kwargs,
        )

    def test_preserves_bases_cds_protein_id_and_explicit_seqid(self) -> None:
        original = self.ann.read_text()
        self.run_update()
        ann = (self.output / "update.ann.txt").read_text()
        self.assertIn("BAAGII010000001\tsource", ann)
        self.assertIn("\tsubmitter_seqid\tscaffold1\n", ann)
        self.assertIn("\tprotein_id\tGAB2236588.1\n", ann)
        self.assertIn("\tCDS\tjoin(4..6,10..12)", ann)
        self.assertNotIn("remove this exon qualifier", ann)
        self.assertNotIn("\t5'UTR\t", ann)
        self.assertIn("\texon\t7..9", ann)
        self.assertIn("\texon\tcomplement(1..6)", ann)
        self.assertIn("\t3'UTR\t10..12\tlocus_tag\tOther_1", ann)
        self.assertEqual((self.output / "update.fasta").read_text(),
                         ">BAAGII010000001\nacgTNN\nACGTNN\n//\n")
        self.assertEqual(self.ann.read_text(), original)
        manifest = json.loads((self.output / "update-manifest.json").read_text())
        self.assertEqual(manifest["sequences"]["scaffold1"]["length"], 12)
        self.assertEqual(manifest["feature_counts"]["omitted_exon"], 1)

    def test_utr_opt_in(self) -> None:
        self.run_update(retain_utr_features=True)
        ann = (self.output / "update.ann.txt").read_text()
        self.assertIn("\t5'UTR\t1..3", ann)
        self.assertNotIn("\texon\t1..6", ann)

    def test_refuses_overwrite(self) -> None:
        self.run_update()
        with self.assertRaisesRegex(MSSPackError, "already exists"):
            self.run_update()

    def test_rejects_incomplete_or_duplicate_mapping(self) -> None:
        for extra in (
            "scaffold2\tBAAGII010000002\tscaffold2\n",
            "scaffold1\tBAAGII010000002\tscaffold2\n",
            "scaffold2\tBAAGII010000001\tscaffold2\n",
        ):
            with self.subTest(extra=extra):
                initial = self.mapping.read_text()
                self.mapping.write_text(initial + extra)
                with self.assertRaises(MSSPackError):
                    self.run_update()
                self.assertFalse(self.output.exists())
                self.mapping.write_text(initial)

    def test_rejects_mismatched_fasta_without_partial_outputs(self) -> None:
        self.fasta.write_text(">wrong\nACGT\n//\n")
        with self.assertRaisesRegex(MSSPackError, "FASTA entry"):
            self.run_update()
        self.assertFalse(self.output.exists())

    def test_rejects_conflicting_submitter_seqid(self) -> None:
        self.ann.write_text(self.ann.read_text().replace("@@[entry]@@", "different"))
        with self.assertRaisesRegex(MSSPackError, "conflicts"):
            self.run_update()
        self.assertFalse(self.output.exists())

    def test_cli(self) -> None:
        self.assertEqual(main([
            "prepare-update", "--ann", str(self.ann), "--fasta", str(self.fasta),
            "--mapping", str(self.mapping), "--output-dir", str(self.output),
        ]), 0)
        self.assertTrue((self.output / "update-manifest.json").is_file())

    def test_rejects_source_length_mismatch(self) -> None:
        self.ann.write_text(self.ann.read_text().replace("source\t1..12", "source\t1..13"))
        with self.assertRaisesRegex(MSSPackError, "lengths differ"):
            self.run_update()
        self.assertFalse(self.output.exists())

    def test_retains_fuzzy_and_unrepresented_locations(self) -> None:
        self.ann.write_text(self.ann.read_text() + "\texon\t<1..6\tlocus_tag\tTest_1\n")
        self.run_update()
        self.assertIn("\texon\t<1..6", (self.output / "update.ann.txt").read_text())

    def test_multiple_entries_and_duplicate_fasta(self) -> None:
        self.mapping.write_text(self.mapping.read_text()
                                + "scaffold2\tBAAGII010000002\tscaffold2\n")
        self.ann.write_text(self.ann.read_text()
                            + "scaffold2\tsource\t1..4\torganism\tTest species\n")
        self.fasta.write_text(self.fasta.read_text() + ">scaffold2\nACGT\n//\n")
        self.run_update()
        self.assertIn(">BAAGII010000002", (self.output / "update.fasta").read_text())
        self.assertIn("\tsubmitter_seqid\tscaffold2", (self.output / "update.ann.txt").read_text())
        self.output = self.base / "duplicate"
        self.fasta.write_text(self.fasta.read_text() + ">scaffold2\nACGT\n//\n")
        with self.assertRaisesRegex(MSSPackError, "duplicate FASTA"):
            self.run_update()
        self.assertFalse(self.output.exists())

    def test_accepts_feature_headers_without_qualifiers(self) -> None:
        self.ann.write_text(self.ann.read_text() + "\tmisc_feature\t7..9\n")
        self.run_update()
        self.assertIn("\tmisc_feature\t7..9\t\t\n",
                      (self.output / "update.ann.txt").read_text())
