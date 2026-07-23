import tempfile
import unittest
from pathlib import Path

from msspack.gff_adjustments import apply_padding_to_gff, fix_gff_to_inframe


class GffAdjustmentTests(unittest.TestCase):
    def test_fix_gff_to_inframe_synchronizes_explicit_terminal_codons(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            input_gff = base / "input.gff3"
            output_gff = base / "output.gff3"
            input_gff.write_text(
                "\n".join(
                    [
                        "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tID=g1",
                        "chr1\tsrc\tmRNA\t1\t8\t.\t+\t.\tID=t1;Parent=g1",
                        "chr1\tsrc\texon\t1\t8\t.\t+\t.\tID=e1;Parent=t1",
                        "chr1\tsrc\tCDS\t1\t8\t.\t+\t1\tID=c1;Parent=t1",
                        "chr1\tsrc\tstart_codon\t1\t3\t.\t+\t0\tID=s1;Parent=t1",
                        "chr1\tsrc\tstop_codon\t6\t8\t.\t+\t0\tID=p1;Parent=t1",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            summary = fix_gff_to_inframe(
                input_path=input_gff,
                output_path=output_gff,
                log_path=base / "summary.log",
            )

            text = output_gff.read_text(encoding="utf-8")
            self.assertIn("gene\t2\t7", text)
            self.assertIn("mRNA\t2\t7", text)
            self.assertIn("exon\t2\t7", text)
            self.assertIn("CDS\t2\t7\t.\t+\t0", text)
            self.assertIn("start_codon\t2\t4\t.\t+\t0", text)
            self.assertIn("stop_codon\t5\t7\t.\t+\t0", text)
            self.assertGreaterEqual(summary["synchronized_features"], 2)

    def test_inframe_removes_exhausted_terminal_segment_and_stale_intron(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            input_gff = base / "input.gff3"
            output_gff = base / "output.gff3"
            input_gff.write_text(
                "\n".join(
                    [
                        "chr1\tsrc\tgene\t1\t14\t.\t+\t.\tID=g1",
                        "chr1\tsrc\tmRNA\t1\t14\t.\t+\t.\tID=t1;Parent=g1",
                        "chr1\tsrc\texon\t1\t1\t.\t+\t.\tID=e1;Parent=t1",
                        "chr1\tsrc\texon\t10\t14\t.\t+\t.\tID=e2;Parent=t1",
                        "chr1\tsrc\tCDS\t1\t1\t.\t+\t2\tID=c1;Parent=t1",
                        "chr1\tsrc\tCDS\t10\t14\t.\t+\t0\tID=c2;Parent=t1",
                        "chr1\tsrc\tintron\t2\t9\t.\t+\t.\tID=i1;Parent=t1",
                        "chr1\tsrc\tstart_codon\t1\t1\t.\t+\t0\tID=s1;Parent=t1",
                        "chr1\tsrc\tstart_codon\t10\t11\t.\t+\t0\tID=s2;Parent=t1",
                        "chr1\tsrc\tstop_codon\t12\t14\t.\t+\t0\tID=p1;Parent=t1",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            fix_gff_to_inframe(
                input_path=input_gff,
                output_path=output_gff,
                log_path=base / "summary.log",
            )

            text = output_gff.read_text(encoding="utf-8")
            self.assertIn("gene\t10\t12", text)
            self.assertIn("mRNA\t10\t12", text)
            self.assertIn("CDS\t10\t12", text)
            self.assertNotIn("ID=c1", text)
            self.assertNotIn("ID=e1", text)
            self.assertNotIn("\tintron\t", text)
            self.assertEqual(text.count("\tstart_codon\t"), 1)
            self.assertIn("start_codon\t10\t12", text)
            self.assertIn("stop_codon\t10\t12", text)

    def test_inframe_synchronizes_terminal_codons_on_minus_strand(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            input_gff = base / "input.gff3"
            output_gff = base / "output.gff3"
            input_gff.write_text(
                "\n".join(
                    [
                        "chr1\tsrc\tgene\t1\t8\t.\t-\t.\tID=g1",
                        "chr1\tsrc\tmRNA\t1\t8\t.\t-\t.\tID=t1;Parent=g1",
                        "chr1\tsrc\texon\t1\t8\t.\t-\t.\tID=e1;Parent=t1",
                        "chr1\tsrc\tCDS\t1\t8\t.\t-\t1\tID=c1;Parent=t1",
                        "chr1\tsrc\tstart_codon\t6\t8\t.\t-\t0\tID=s1;Parent=t1",
                        "chr1\tsrc\tstop_codon\t1\t3\t.\t-\t0\tID=p1;Parent=t1",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            fix_gff_to_inframe(
                input_path=input_gff,
                output_path=output_gff,
                log_path=base / "summary.log",
            )

            text = output_gff.read_text(encoding="utf-8")
            self.assertIn("CDS\t2\t7\t.\t-\t0", text)
            self.assertIn("start_codon\t5\t7\t.\t-\t0", text)
            self.assertIn("stop_codon\t2\t4\t.\t-\t0", text)

    def test_fix_gff_to_inframe_updates_boundaries_and_phase(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            input_gff = base / "input.gff3"
            output_gff = base / "output.gff3"
            log_path = base / "summary.log"
            input_gff.write_text(
                "\n".join(
                    [
                        "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tID=g1",
                        "chr1\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=t1;Parent=g1",
                        "chr1\tsrc\texon\t1\t5\t.\t+\t.\tID=t1.ex1;Parent=t1",
                        "chr1\tsrc\texon\t10\t12\t.\t+\t.\tID=t1.ex2;Parent=t1",
                        "chr1\tsrc\tCDS\t1\t5\t.\t+\t1\tID=t1.cds;Parent=t1",
                        "chr1\tsrc\tCDS\t10\t12\t.\t+\t0\tID=t1.cds;Parent=t1",
                        "",
                    ]
                ),
                encoding="utf-8",
            )

            summary = fix_gff_to_inframe(
                input_path=input_gff,
                output_path=output_gff,
                log_path=log_path,
            )

            self.assertEqual(summary["updated_gene_models"], 1)
            text = output_gff.read_text(encoding="utf-8")
            self.assertIn("gene\t2\t11", text)
            self.assertIn("mRNA\t2\t11", text)
            self.assertIn("exon\t2\t5", text)
            self.assertIn("exon\t10\t11", text)
            self.assertIn("CDS\t2\t5\t.\t+\t0", text)
            self.assertIn("CDS\t10\t11\t.\t+\t0", text)

    def test_apply_padding_to_gff_updates_gene_and_lists_updated_gene(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            input_gff = base / "input.gff3"
            padding_log = base / "padding.log"
            output_gff = base / "output.gff3"
            with_stops = base / "with_stops.txt"
            updated = base / "updated.txt"
            input_gff.write_text(
                "\n".join(
                    [
                        "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=g1",
                        "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=g1-RA;Parent=g1",
                        "chr1\tsrc\texon\t1\t9\t.\t+\t.\tID=g1-RA.ex1;Parent=g1-RA",
                        "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=g1-RA.cds;Parent=g1-RA",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            padding_log.write_text(
                "g1, original_seqlen=9, head_padding=2, tail_padding=3, original_num_stop=1, new_num_stop=0\n",
                encoding="utf-8",
            )

            summary = apply_padding_to_gff(
                gff_path=input_gff,
                padding_log_path=padding_log,
                output_path=output_gff,
                genes_with_stops_path=with_stops,
                updated_genes_path=updated,
            )

            self.assertEqual(summary["updated_genes"], ["g1"])
            self.assertIn("g1", updated.read_text(encoding="utf-8"))
            self.assertEqual(with_stops.read_text(encoding="utf-8"), "")
            text = output_gff.read_text(encoding="utf-8")
            self.assertIn("gene\t2\t9", text)
            self.assertIn("mRNA\t2\t9", text)
            self.assertIn("exon\t2\t9", text)
            self.assertIn("CDS\t2\t9\t.\t+\t0", text)

    def test_apply_padding_to_gff_synchronizes_explicit_terminal_codons(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            input_gff = base / "input.gff3"
            input_gff.write_text(
                "\n".join(
                    [
                        "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=g1",
                        "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1",
                        "chr1\tsrc\texon\t1\t9\t.\t+\t.\tID=e1;Parent=t1",
                        "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=t1",
                        "chr1\tsrc\tstart_codon\t1\t3\t.\t+\t0\tID=s1;Parent=t1",
                        "chr1\tsrc\tstop_codon\t7\t9\t.\t+\t0\tID=p1;Parent=t1",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            (base / "padding.log").write_text(
                "t1, original_seqlen=9, head_padding=2, tail_padding=3, "
                "original_num_stop=1, new_num_stop=0\n",
                encoding="utf-8",
            )

            summary = apply_padding_to_gff(
                gff_path=input_gff,
                padding_log_path=base / "padding.log",
                output_path=base / "output.gff3",
                genes_with_stops_path=base / "stops.txt",
                updated_genes_path=base / "updated.txt",
            )

            text = (base / "output.gff3").read_text(encoding="utf-8")
            self.assertIn("start_codon\t2\t4\t.\t+\t0", text)
            self.assertIn("stop_codon\t7\t9\t.\t+\t0", text)
            self.assertGreaterEqual(summary["synchronized_features"], 1)

    def test_padding_clips_utrs_and_preserves_standalone_non_cds_features(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            input_gff = base / "input.gff3"
            input_gff.write_text(
                "\n".join(
                    [
                        "chr1\tsrc\tgene\t1\t12\t.\t+\t.\tID=g1",
                        "chr1\tsrc\tmRNA\t1\t12\t.\t+\t.\tID=t1;Parent=g1",
                        "chr1\tsrc\texon\t1\t12\t.\t+\t.\tID=e1;Parent=t1",
                        "chr1\tsrc\tCDS\t3\t10\t.\t+\t0\tID=c1;Parent=t1",
                        "chr1\tsrc\tfive_prime_UTR\t1\t2\t.\t+\t.\tID=u5;Parent=t1",
                        "chr1\tsrc\tthree_prime_UTR\t11\t12\t.\t+\t.\tID=u3;Parent=t1",
                        "chr1\tsrc\trepeat_region\t20\t30\t.\t+\t.\tID=rep1",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            (base / "padding.log").write_text(
                "t1, original_seqlen=8, head_padding=2, tail_padding=3, "
                "original_num_stop=1, new_num_stop=0\n",
                encoding="utf-8",
            )

            apply_padding_to_gff(
                gff_path=input_gff,
                padding_log_path=base / "padding.log",
                output_path=base / "output.gff3",
                genes_with_stops_path=base / "stops.txt",
                updated_genes_path=base / "updated.txt",
            )

            text = (base / "output.gff3").read_text(encoding="utf-8")
            self.assertIn("five_prime_UTR\t2\t2", text)
            self.assertIn("three_prime_UTR\t11\t12", text)
            self.assertIn("repeat_region\t20\t30", text)

    def test_padding_changes_only_the_target_transcript(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            input_gff = base / "input.gff3"
            input_gff.write_text(
                "\n".join(
                    [
                        "chr1\tsrc\tgene\t1\t28\t.\t+\t.\tID=g1",
                        "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1",
                        "chr1\tsrc\texon\t1\t9\t.\t+\t.\tID=e1;Parent=t1",
                        "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=t1",
                        "chr1\tsrc\tmRNA\t20\t28\t.\t+\t.\tID=t2;Parent=g1",
                        "chr1\tsrc\texon\t20\t28\t.\t+\t.\tID=e2;Parent=t2",
                        "chr1\tsrc\tCDS\t20\t28\t.\t+\t0\tID=c2;Parent=t2",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            (base / "padding.log").write_text(
                "t1, original_seqlen=9, head_padding=2, tail_padding=3, "
                "original_num_stop=1, new_num_stop=0\n",
                encoding="utf-8",
            )

            apply_padding_to_gff(
                gff_path=input_gff,
                padding_log_path=base / "padding.log",
                output_path=base / "output.gff3",
                genes_with_stops_path=base / "stops.txt",
                updated_genes_path=base / "updated.txt",
            )

            text = (base / "output.gff3").read_text(encoding="utf-8")
            self.assertIn("mRNA\t2\t9\t.\t+\t.\tID=t1", text)
            self.assertIn("CDS\t2\t9\t.\t+\t0\tID=c1", text)
            self.assertIn("mRNA\t20\t28\t.\t+\t.\tID=t2", text)
            self.assertIn("CDS\t20\t28\t.\t+\t0\tID=c2", text)

    def test_apply_padding_to_gff_keeps_gene_with_stops_unmodified(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            input_gff = base / "input.gff3"
            padding_log = base / "padding.log"
            output_gff = base / "output.gff3"
            with_stops = base / "with_stops.txt"
            updated = base / "updated.txt"
            original = "\n".join(
                    [
                        "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=g1",
                        "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=g1-RA;Parent=g1",
                        "chr1\tsrc\texon\t1\t9\t.\t+\t.\tID=g1-RA.ex1;Parent=g1-RA",
                        "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=g1-RA.cds;Parent=g1-RA",
                        "",
                    ]
                )
            input_gff.write_text(original, encoding="utf-8")
            padding_log.write_text(
                "g1, original_seqlen=9, head_padding=2, tail_padding=3, original_num_stop=1, new_num_stop=1\n",
                encoding="utf-8",
            )

            summary = apply_padding_to_gff(
                gff_path=input_gff,
                padding_log_path=padding_log,
                output_path=output_gff,
                genes_with_stops_path=with_stops,
                updated_genes_path=updated,
            )

            self.assertEqual(summary["genes_with_stops"], ["g1"])
            self.assertEqual(output_gff.read_text(encoding="utf-8"), original)
            self.assertIn("g1", with_stops.read_text(encoding="utf-8"))
            self.assertEqual(updated.read_text(encoding="utf-8"), "")

    def test_apply_padding_to_gff_resolves_transcript_ids_without_ra_suffix_rules(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            input_gff = base / "input.gff3"
            padding_log = base / "padding.log"
            output_gff = base / "output.gff3"
            with_stops = base / "with_stops.txt"
            updated = base / "updated.txt"
            input_gff.write_text(
                "\n".join(
                    [
                        "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=g1",
                        "chr1\tsrc\ttranscript\t1\t9\t.\t+\t.\tID=tx1;Parent=g1",
                        "chr1\tsrc\texon\t1\t9\t.\t+\t.\tID=tx1.ex1;Parent=tx1",
                        "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=tx1.cds;Parent=tx1",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            padding_log.write_text(
                "tx1, original_seqlen=9, head_padding=2, tail_padding=3, original_num_stop=1, new_num_stop=0\n",
                encoding="utf-8",
            )

            summary = apply_padding_to_gff(
                gff_path=input_gff,
                padding_log_path=padding_log,
                output_path=output_gff,
                genes_with_stops_path=with_stops,
                updated_genes_path=updated,
            )

            self.assertEqual(summary["updated_genes"], ["g1"])
            self.assertIn("g1", updated.read_text(encoding="utf-8"))
            text = output_gff.read_text(encoding="utf-8")
            self.assertIn("gene\t2\t9", text)
            self.assertIn("transcript\t2\t9", text)

    def test_partial_transcript_does_not_lose_an_extra_codon(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            input_gff = base / "input.gff3"
            padding_log = base / "padding.log"
            output_gff = base / "output.gff3"
            input_gff.write_text(
                "\n".join(
                    [
                        "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tID=g1",
                        "chr1\tsrc\ttranscript\t1\t8\t.\t+\t.\tID=t1;Parent=g1",
                        "chr1\tsrc\texon\t1\t8\t.\t+\t.\tID=e1;Parent=t1",
                        "chr1\tsrc\tCDS\t1\t8\t.\t+\t0\tID=c1;Parent=t1",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            padding_log.write_text(
                "t1, original_seqlen=8, head_padding=0, tail_padding=1, "
                "original_num_stop=0, new_num_stop=0\n",
                encoding="utf-8",
            )

            apply_padding_to_gff(
                gff_path=input_gff,
                padding_log_path=padding_log,
                output_path=output_gff,
                genes_with_stops_path=base / "stops.txt",
                updated_genes_path=base / "updated.txt",
            )

            text = output_gff.read_text(encoding="utf-8")
            self.assertIn("gene\t1\t6", text)
            self.assertIn("CDS\t1\t6", text)

    def test_padding_keeps_gene_without_explicit_exons(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            input_gff = base / "input.gff3"
            padding_log = base / "padding.log"
            output_gff = base / "output.gff3"
            input_gff.write_text(
                "\n".join(
                    [
                        "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=g1",
                        "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1",
                        "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=t1",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            padding_log.write_text(
                "t1, original_seqlen=9, head_padding=1, tail_padding=2, "
                "original_num_stop=1, new_num_stop=0\n",
                encoding="utf-8",
            )

            summary = apply_padding_to_gff(
                gff_path=input_gff,
                padding_log_path=padding_log,
                output_path=output_gff,
                genes_with_stops_path=base / "stops.txt",
                updated_genes_path=base / "updated.txt",
            )

            self.assertEqual(summary["updated_genes"], ["g1"])
            self.assertIn("gene", output_gff.read_text(encoding="utf-8"))

    def test_padding_keeps_embedded_fasta_after_feature_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            input_gff = base / "input.gff3"
            output_gff = base / "output.gff3"
            input_gff.write_text(
                "##gff-version 3\n"
                "chr1\tsrc\tgene\t1\t9\t.\t+\t.\tID=g1\n"
                "chr1\tsrc\tmRNA\t1\t9\t.\t+\t.\tID=t1;Parent=g1\n"
                "chr1\tsrc\tCDS\t1\t9\t.\t+\t0\tID=c1;Parent=t1\n"
                "##FASTA\n>chr1\nATGAAATAA\n",
                encoding="utf-8",
            )
            (base / "padding.log").write_text("", encoding="utf-8")

            apply_padding_to_gff(
                gff_path=input_gff,
                padding_log_path=base / "padding.log",
                output_path=output_gff,
                genes_with_stops_path=base / "stops.txt",
                updated_genes_path=base / "updated.txt",
            )

            output = output_gff.read_text(encoding="utf-8")
            self.assertLess(output.index("ID=c1"), output.index("##FASTA"))
            self.assertTrue(output.endswith("##FASTA\n>chr1\nATGAAATAA\n"))


if __name__ == "__main__":
    unittest.main()
