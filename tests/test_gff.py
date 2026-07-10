import tempfile
import unittest
from pathlib import Path

from msspack.gff import (
    iter_gff_records,
    parse_attributes,
    read_gff_document,
    sort_gff_file_precise,
    write_gff_document,
)
from msspack.gff_cleanup import fix_gff_semicolons_file, trim_gff_to_fasta_bounds
from msspack.utils import MSSPackError


class GffSortTests(unittest.TestCase):
    def test_attributes_reject_invalid_percent_escapes(self) -> None:
        with self.assertRaisesRegex(ValueError, "percent escape"):
            parse_attributes("ID=bad%ZZ")

    def test_dot_represents_an_empty_attribute_set(self) -> None:
        self.assertEqual(parse_attributes("."), {})

    def test_semicolon_repair_rejects_orphan_attribute_fragments(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            source = base / "input.gff3"
            source.write_text(
                "chr1\tsrc\tregion\t1\t3\t.\t+\t.\torphan-fragment\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(MSSPackError, "without a preceding key"):
                fix_gff_semicolons_file(
                    input_path=source,
                    output_path=base / "output.gff3",
                    log_path=base / "output.log",
                )

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

    def test_sort_preserves_embedded_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            source = base / "input.gff3"
            output = base / "output.gff3"
            source.write_text(
                "##gff-version 3\n"
                "chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=g1\n"
                "##FASTA\n>chr1\nATG\n",
                encoding="utf-8",
            )

            sort_gff_file_precise(source, output)

            self.assertIn("##FASTA\n>chr1\nATG\n", output.read_text(encoding="utf-8"))

    def test_gff_readers_stop_features_at_embedded_fasta(self) -> None:
        gff = (
            "##gff-version 3\n"
            "chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=g1\n"
            "##FASTA\n>chr1\nATG\n"
        )
        with tempfile.TemporaryDirectory() as tmp_dir:
            source = Path(tmp_dir) / "input.gff3"
            output = Path(tmp_dir) / "output.gff3"
            source.write_text(gff, encoding="utf-8")

            document = read_gff_document(source)
            self.assertEqual([record.attributes["ID"] for record in document.records], ["g1"])
            self.assertEqual(document.fasta_lines, ["##FASTA", ">chr1", "ATG"])
            self.assertEqual([record.attributes["ID"] for record in iter_gff_records(source)], ["g1"])

            write_gff_document(output, document)
            self.assertEqual(output.read_text(encoding="utf-8"), gff)

    def test_sort_rejects_malformed_feature_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            source = base / "input.gff3"
            source.write_text("chr1\tsrc\tgene\t1\n", encoding="utf-8")

            with self.assertRaises(MSSPackError):
                sort_gff_file_precise(source, base / "output.gff3")

    def test_trim_rejects_seqids_missing_from_fasta(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            fasta = base / "input.fa"
            gff = base / "input.gff3"
            fasta.write_text(">chr1\nATG\n", encoding="utf-8")
            gff.write_text(
                "chr2\tsrc\tgene\t1\t3\t.\t+\t.\tID=g1\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(MSSPackError, "chr2"):
                trim_gff_to_fasta_bounds(
                    gff_path=gff,
                    fasta_path=fasta,
                    output_path=base / "output.gff3",
                    log_path=base / "output.log",
                )

    def test_trim_preserves_embedded_fasta_without_parsing_it_as_features(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            fasta = base / "input.fa"
            gff = base / "input.gff3"
            output = base / "output.gff3"
            fasta.write_text(">chr1\nATG\n", encoding="utf-8")
            gff.write_text(
                "chr1\tsrc\tgene\t0\t4\t.\t+\t.\tID=g1\n"
                "##FASTA\n>chr1\nATG\n",
                encoding="utf-8",
            )

            trim_gff_to_fasta_bounds(
                gff_path=gff,
                fasta_path=fasta,
                output_path=output,
                log_path=base / "output.log",
            )

            self.assertEqual(
                output.read_text(encoding="utf-8"),
                "chr1\tsrc\tgene\t1\t3\t.\t+\t.\tID=g1\n"
                "##FASTA\n>chr1\nATG\n",
            )


if __name__ == "__main__":
    unittest.main()
