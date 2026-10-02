import tempfile
import unittest
from pathlib import Path

from msspack.gff import child_ids, iter_gff_records
from msspack.pipeline_actions import copy_input_fasta, pad_locus_tags
from msspack.utils import MSSPackError


class PipelineActionTests(unittest.TestCase):
    def test_locus_padding_changes_identifiers_and_references_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            source = base / "input.gff"
            source.write_text(
                "##sequence-region Fix_1 1 9\n"
                "Fix_1\tFix_1\tgene\t1\t9\t.\t+\t.\t%49D=Fix%5F1;Name=Fix_1;Note=Fix_1\n"
                "Fix_1\tFix_1\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=Fix_1,Fix%5F2;Derives_from=Fix_1;product=Fix_1\n"
                "##FASTA\n>Fix_1\nATGAAATAA\n",
            )
            output = base / "output.gff"
            pad_locus_tags(source, output, "Fix", 6, log_path=base / "padding.log")
            records = list(iter_gff_records(output))
            self.assertEqual([row.seqid for row in records], ["Fix_1", "Fix_1"])
            self.assertEqual([row.source for row in records], ["Fix_1", "Fix_1"])
            self.assertEqual(records[0].attributes["ID"], "Fix_000001")
            self.assertEqual(child_ids(records[1].attributes["Parent"]), ["Fix_000001", "Fix_000002"])
            self.assertEqual(records[1].attributes["Derives_from"], "Fix_000001")
            self.assertEqual(records[0].attributes["Name"], "Fix_1")
            self.assertEqual(records[0].attributes["Note"], "Fix_1")
            self.assertEqual(records[1].attributes["product"], "Fix_1")
            self.assertTrue(output.read_text().startswith("##sequence-region Fix_1 1 9\n"))
            self.assertTrue(output.read_text().endswith("##FASTA\n>Fix_1\nATGAAATAA\n"))

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
