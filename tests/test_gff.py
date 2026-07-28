import tempfile
import unittest
from pathlib import Path

from msspack.coordinate_duplicates import read_coordinate_duplicate_map
from msspack.gff import (
    iter_gff_records,
    parse_attributes,
    read_gff_document,
    sort_gff_file_precise,
    write_gff_document,
)
from msspack.gff_cleanup import (
    drop_duplicate_coordinate_genes,
    fix_gff_semicolons_file,
    trim_gff_to_fasta_bounds,
)
from msspack.utils import MSSPackError


class GffSortTests(unittest.TestCase):
    def test_duplicate_coordinate_cleanup_records_kept_removed_pairs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            source = base / "input.gff3"
            output = base / "output.gff3"
            duplicate_map = base / "duplicates.tsv"
            source.write_text(
                "\n".join(
                    [
                        "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=kept",
                        "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=kept.t1;Parent=kept",
                        "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=removed1",
                        (
                            "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\t"
                            "ID=removed1.t1;Parent=removed1"
                        ),
                        "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=removed2",
                        "chr1\tsrc\tgene\t1\t100\t.\t-\t.\tID=opposite_strand",
                        "chr2\tsrc\tgene\t1\t100\t.\t+\t.\tID=other_seqid",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            result = drop_duplicate_coordinate_genes(
                input_path=source,
                output_path=output,
                log_path=base / "output.log",
                selection_policy="first",
                removed_gene_ids_path=base / "removed.txt",
                duplicate_map_path=duplicate_map,
                metrics_path=base / "metrics.json",
            )

            self.assertEqual(result["removed_gene_ids"], ["removed1", "removed2"])
            output_text = output.read_text(encoding="utf-8")
            self.assertIn("ID=kept\n", output_text)
            self.assertNotIn("ID=removed1", output_text)
            self.assertNotIn("ID=removed2", output_text)
            self.assertIn("ID=opposite_strand", output_text)
            map_lines = duplicate_map.read_text(encoding="utf-8").splitlines()
            self.assertEqual(len(map_lines), 3)
            self.assertIn("\tkept\tremoved1", map_lines[1])
            self.assertIn("\tkept\tremoved2", map_lines[2])

    def test_duplicate_coordinate_cleanup_prefers_longer_valid_cds(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            fasta = base / "input.fa"
            source = base / "input.gff3"
            output = base / "output.gff3"
            duplicate_map = base / "duplicates.tsv"
            fasta.write_text(">chr1\nATGAAACCCGGGTAA\n", encoding="utf-8")
            source.write_text(
                "\n".join(
                    [
                        "chr1\tsrc\tgene\t1\t15\t.\t+\t.\tID=short",
                        "chr1\tsrc\tmRNA\t1\t15\t.\t+\t.\tID=short.t1;Parent=short",
                        "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tParent=short.t1",
                        "chr1\tsrc\tCDS\t10\t15\t.\t+\t0\tParent=short.t1",
                        "chr1\tsrc\tgene\t1\t15\t.\t+\t.\tID=long",
                        "chr1\tsrc\tmRNA\t1\t15\t.\t+\t.\tID=long.t1;Parent=long",
                        "chr1\tsrc\tCDS\t1\t15\t.\t+\t0\tParent=long.t1",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            result = drop_duplicate_coordinate_genes(
                input_path=source,
                fasta_path=fasta,
                output_path=output,
                log_path=base / "output.log",
                duplicate_map_path=duplicate_map,
            )

            self.assertEqual(result["removed_gene_ids"], ["short"])
            pair = read_coordinate_duplicate_map(duplicate_map)[0]
            self.assertEqual(pair.kept_gene_id, "long")
            self.assertEqual(pair.removed_gene_id, "short")
            self.assertEqual(pair.kept_transcript_id, "long.t1")
            self.assertEqual(pair.removed_transcript_id, "short.t1")
            self.assertEqual(pair.kept_cds_length, 15)
            self.assertEqual(pair.removed_cds_length, 12)
            self.assertEqual(pair.kept_intron_count, 0)
            self.assertEqual(pair.removed_intron_count, 1)
            self.assertEqual(pair.selection_reason, "longer_cds")
            self.assertTrue(pair.low_confidence)
            output_text = output.read_text(encoding="utf-8")
            self.assertIn("ID=long", output_text)
            self.assertNotIn("ID=short", output_text)

    def test_duplicate_coordinate_cleanup_rejects_longer_cds_with_internal_stop(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            fasta = base / "input.fa"
            source = base / "input.gff3"
            output = base / "output.gff3"
            duplicate_map = base / "duplicates.tsv"
            fasta.write_text(">chr1\nATGTAACCCTAA\n", encoding="utf-8")
            source.write_text(
                "\n".join(
                    [
                        "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tID=long",
                        "chr1\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=long.t1;Parent=long",
                        "chr1\tsrc\tCDS\t1\t12\t.\t+\t0\tParent=long.t1",
                        "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tID=clean",
                        "chr1\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=clean.t1;Parent=clean",
                        "chr1\tsrc\tCDS\t1\t3\t.\t+\t0\tParent=clean.t1",
                        "chr1\tsrc\tCDS\t7\t12\t.\t+\t0\tParent=clean.t1",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            drop_duplicate_coordinate_genes(
                input_path=source,
                fasta_path=fasta,
                output_path=output,
                log_path=base / "output.log",
                duplicate_map_path=duplicate_map,
            )

            pair = read_coordinate_duplicate_map(duplicate_map)[0]
            self.assertEqual(pair.kept_gene_id, "clean")
            self.assertEqual(pair.removed_gene_id, "long")
            self.assertEqual(pair.kept_internal_stops, 0)
            self.assertEqual(pair.removed_internal_stops, 1)
            self.assertEqual(pair.kept_cds_length, 9)
            self.assertEqual(pair.removed_cds_length, 12)
            self.assertEqual(pair.selection_reason, "fewer_internal_stops")
            self.assertFalse(pair.low_confidence)

    def test_duplicate_coordinate_cleanup_uses_best_transcript_per_gene(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            fasta = base / "input.fa"
            source = base / "input.gff3"
            duplicate_map = base / "duplicates.tsv"
            fasta.write_text(">chr1\nATGAAACCCGGGTAA\n", encoding="utf-8")
            source.write_text(
                "\n".join(
                    [
                        "chr1\tsrc\tgene\t1\t15\t.\t+\t.\tID=multi",
                        "chr1\tsrc\tmRNA\t1\t15\t.\t+\t.\tID=multi.short;Parent=multi",
                        "chr1\tsrc\tCDS\t1\t6\t.\t+\t0\tParent=multi.short",
                        "chr1\tsrc\tCDS\t10\t15\t.\t+\t0\tParent=multi.short",
                        "chr1\tsrc\tmRNA\t1\t15\t.\t+\t.\tID=multi.long;Parent=multi",
                        "chr1\tsrc\tCDS\t1\t15\t.\t+\t0\tParent=multi.long",
                        "chr1\tsrc\tgene\t1\t15\t.\t+\t.\tID=partial",
                        "chr1\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=partial.t1;Parent=partial",
                        "chr1\tsrc\tCDS\t1\t12\t.\t+\t0\tParent=partial.t1",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            drop_duplicate_coordinate_genes(
                input_path=source,
                fasta_path=fasta,
                output_path=base / "output.gff3",
                log_path=base / "output.log",
                duplicate_map_path=duplicate_map,
            )

            pair = read_coordinate_duplicate_map(duplicate_map)[0]
            self.assertEqual(pair.kept_gene_id, "multi")
            self.assertEqual(pair.kept_transcript_id, "multi.long")
            self.assertEqual(pair.selection_reason, "more_complete_boundaries")

    def test_duplicate_coordinate_cleanup_translates_minus_strand_cds(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            fasta = base / "input.fa"
            source = base / "input.gff3"
            duplicate_map = base / "duplicates.tsv"
            fasta.write_text(">chr1\nTTAGGGTTTCAT\n", encoding="utf-8")
            source.write_text(
                "\n".join(
                    [
                        "chr1\tsrc\tgene\t1\t12\t.\t-\t.\tID=complete",
                        "chr1\tsrc\tmRNA\t1\t12\t.\t-\t.\tID=complete.t1;Parent=complete",
                        "chr1\tsrc\tCDS\t1\t12\t.\t-\t0\tParent=complete.t1",
                        "chr1\tsrc\tgene\t1\t12\t.\t-\t.\tID=partial",
                        "chr1\tsrc\tmRNA\t1\t12\t.\t-\t.\tID=partial.t1;Parent=partial",
                        "chr1\tsrc\tCDS\t4\t12\t.\t-\t0\tParent=partial.t1",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            drop_duplicate_coordinate_genes(
                input_path=source,
                fasta_path=fasta,
                output_path=base / "output.gff3",
                log_path=base / "output.log",
                duplicate_map_path=duplicate_map,
            )

            pair = read_coordinate_duplicate_map(duplicate_map)[0]
            self.assertEqual(pair.kept_gene_id, "complete")
            self.assertTrue(pair.kept_complete)
            self.assertFalse(pair.removed_complete)
            self.assertEqual(pair.selection_reason, "more_complete_boundaries")

    def test_duplicate_coordinate_cleanup_can_keep_all_collisions(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            source = base / "input.gff3"
            duplicate_map = base / "duplicates.tsv"
            source.write_text(
                "chr1\tsrc\tgene\t1\t10\t.\t+\t.\tID=first\n"
                "chr1\tsrc\tgene\t1\t10\t.\t+\t.\tID=second\n",
                encoding="utf-8",
            )

            result = drop_duplicate_coordinate_genes(
                input_path=source,
                output_path=base / "output.gff3",
                log_path=base / "output.log",
                selection_policy="keep_all",
                duplicate_map_path=duplicate_map,
            )

            self.assertEqual(result["removed_gene_ids"], [])
            self.assertEqual(result["coordinate_collision_groups"], 1)
            self.assertEqual(read_coordinate_duplicate_map(duplicate_map), ())
            self.assertEqual(
                base.joinpath("output.gff3").read_text(encoding="utf-8"),
                source.read_text(encoding="utf-8"),
            )

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
