import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

from msspack.mss_converter.core import (
    ConversionOptions,
    _get_start_codons,
    build_gff_feature_indexes,
    convert_gff_to_mss,
    detect_gap_regions,
)


class MssConverterTests(unittest.TestCase):
    def test_start_codons_follow_selected_genetic_code(self) -> None:
        self.assertIn("ATA", _get_start_codons("2"))

    def test_detect_gap_regions_finds_contiguous_n_runs(self) -> None:
        record = SimpleNamespace(seq="AANNNTAANNNNN")
        out, gaps = detect_gap_regions(
            record,
            linkage_evidence="paired-ends",
            min_assembly_gap_size=3,
            gap_estimated_length="unknown",
        )
        self.assertEqual(gaps, [(3, 5), (9, 13)])
        self.assertIn("assembly_gap\t3..5", out)
        self.assertIn("assembly_gap\t9..13", out)

    def test_build_gff_feature_indexes_normalizes_parent_and_phase(self) -> None:
        rows = [
            {
                "seq_id": "chr1",
                "type": "gene",
                "ID": "Gene1",
                "Parent": None,
                "start": 10,
                "end": 50,
                "strand": "+",
                "phase": ".",
                "Name": None,
                "Type": None,
                "anticodon": None,
            },
            {
                "seq_id": "chr1",
                "type": "mRNA",
                "ID": "Gene1-T1",
                "Parent": "Gene1",
                "start": 10,
                "end": 50,
                "strand": "+",
                "phase": ".",
                "Name": "Gene1-T1",
                "Type": None,
                "anticodon": None,
            },
            {
                "seq_id": "chr1",
                "type": "CDS",
                "ID": "Gene1-T1.cds1",
                "Parent": "Gene1-T1",
                "start": 20,
                "end": 40,
                "strand": "+",
                "phase": "2",
                "Name": None,
                "Type": None,
                "anticodon": None,
            },
        ]

        gene_lookup, parent_lookup = build_gff_feature_indexes(rows)

        self.assertEqual(gene_lookup["chr1"][0].id, "Gene1")
        self.assertEqual(parent_lookup["Gene1"][0].id, "Gene1-T1")
        self.assertEqual(parent_lookup["Gene1-T1"][0].phase, 2)
        self.assertEqual(parent_lookup["Gene1-T1"][0].parent, "Gene1-T1")

    def test_build_gff_feature_indexes_splits_multiple_parents(self) -> None:
        shared = {
            "seq_id": "chr1",
            "type": "CDS",
            "ID": "shared",
            "Parent": "tx1,tx2",
            "start": 1,
            "end": 9,
            "strand": "+",
            "phase": "0",
        }

        _, parent_lookup = build_gff_feature_indexes([shared])

        self.assertEqual(parent_lookup["tx1"][0].id, "shared")
        self.assertEqual(parent_lookup["tx2"][0].id, "shared")

    def test_converter_renders_mixed_gff_feature_types_without_silent_loss(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            fasta = base / "genome.fa"
            gff = base / "annotation.gff3"
            annotation = base / "annotation.tsv"
            output = base / "output.ann.txt"
            fasta.write_text(">chr1\n" + "A" * 400 + "\n", encoding="utf-8")
            gff.write_text(
                "\n".join(
                    [
                        "##gff-version 3",
                        "chr1\tsrc\tgene\t1\t90\t.\t+\t.\tID=g1",
                        "chr1\tsrc\tmRNA\t1\t90\t.\t+\t.\tID=tx1;Parent=g1",
                        "chr1\tsrc\texon\t1\t90\t.\t+\t.\tID=ex1;Parent=tx1",
                        "chr1\tsrc\tfive_prime_UTR\t1\t9\t.\t+\t.\tID=utr1;Parent=tx1",
                        "chr1\tsrc\tCDS\t10\t90\t.\t+\t0\tID=cds1;Parent=tx1",
                        "chr1\tsrc\tgene\t100\t130\t.\t+\t.\tID=g2",
                        "chr1\tsrc\tmiRNA\t100\t130\t.\t+\t.\tID=mir1;Parent=g2;Name=mir-test",
                        "chr1\tsrc\texon\t100\t130\t.\t+\t.\tID=mir-ex;Parent=mir1",
                        "chr1\tsrc\tgene\t140\t170\t.\t-\t.\tID=g3",
                        "chr1\tsrc\trRNA\t140\t170\t.\t-\t.\tID=rrna1;Parent=g3;Type=16S",
                        "chr1\tsrc\tgene\t180\t210\t.\t+\t.\tID=g4",
                        "chr1\tsrc\ttRNA\t180\t210\t.\t+\t.\tID=trna1;Parent=g4;Name=tRNA-Lys",
                        "chr1\tsrc\tpseudogene\t220\t250\t.\t+\t.\tID=pg1;Name=abcP;pseudogene=processed",
                        "chr1\tsrc\trepeat_region\t260\t280\t.\t+\t.\tID=rep1;rpt_type=tandem",
                        "chr1\tsrc\tpromoter\t281\t290\t.\t+\t.\tID=prom1",
                        "chr1\tsrc\ttransposable_element\t291\t320\t.\t+\t.\tID=te1;Name=TestTE",
                        "chr1\tsrc\tnew_SO_feature\t321\t330\t.\t+\t.\tID=new1;Note=novel region",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            annotation.write_text(
                "ID\tDescription\n"
                "tx1\tATP synthase subunit alpha\n",
                encoding="utf-8",
            )

            summary = convert_gff_to_mss(
                ConversionOptions(
                    fasta_path=fasta,
                    gff_path=gff,
                    annotation_path=annotation,
                    output_path=output,
                    locus_tag_prefix="Mix",
                    organism_name="Test organism",
                )
            )

            text = output.read_text(encoding="utf-8")
            self.assertIn("\tmRNA\t1..90\tlocus_tag\tMix000000100", text)
            self.assertIn("\texon\t1..90\tlocus_tag\tMix000000100", text)
            self.assertIn("\t5'UTR\t1..9\tlocus_tag\tMix000000100", text)
            self.assertIn("\tCDS\t10..90\tlocus_tag\tMix000000100", text)
            self.assertIn("\tncRNA\t100..130\tlocus_tag\tMix000000200", text)
            self.assertIn("\t\t\tncRNA_class\tmiRNA", text)
            self.assertIn("\trRNA\tcomplement(140..170)\tproduct\t16S rRNA", text)
            self.assertIn("\ttRNA\t180..210\tproduct\ttRNA-Lys", text)
            self.assertIn("\tmisc_feature\t220..250\tlocus_tag\tMix000000500", text)
            self.assertIn("\t\t\tpseudogene\tprocessed", text)
            self.assertIn("\trepeat_region\t260..280\trpt_type\ttandem", text)
            self.assertIn("\tregulatory\t281..290\tregulatory_class\tpromoter", text)
            self.assertIn(
                "\tmobile_element\t291..320\tmobile_element_type\tother:TestTE",
                text,
            )
            self.assertIn(
                "\tmisc_feature\t321..330\tnote\tnovel region; GFF3 feature type: new_SO_feature",
                text,
            )
            self.assertEqual(summary.overall_counts["unknown_features"], 1)

    def test_suppresses_transcript_structure_when_exons_equal_cds(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            fasta = base / "genome.fa"
            gff = base / "annotation.gff3"
            annotation = base / "annotation.tsv"
            output = base / "output.ann.txt"
            fasta.write_text(">chr1\n" + "A" * 120 + "\n", encoding="utf-8")
            gff.write_text(
                "\n".join(
                    [
                        "##gff-version 3",
                        "chr1\tsrc\tgene\t1\t90\t.\t+\t.\tID=g1",
                        "chr1\tsrc\tmRNA\t1\t90\t.\t+\t.\tID=tx1;Parent=g1",
                        "chr1\tsrc\texon\t1\t30\t.\t+\t.\tID=ex1;Parent=tx1",
                        "chr1\tsrc\tCDS\t1\t30\t.\t+\t0\tID=cds1;Parent=tx1",
                        "chr1\tsrc\tintron\t31\t60\t.\t+\t.\tID=int1;Parent=tx1",
                        "chr1\tsrc\texon\t61\t90\t.\t+\t.\tID=ex2;Parent=tx1",
                        "chr1\tsrc\tCDS\t61\t90\t.\t+\t0\tID=cds2;Parent=tx1",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            annotation.write_text(
                "ID\tDescription\n"
                "tx1\ttest protein\n",
                encoding="utf-8",
            )

            convert_gff_to_mss(
                ConversionOptions(
                    fasta_path=fasta,
                    gff_path=gff,
                    annotation_path=annotation,
                    output_path=output,
                    locus_tag_prefix="Red",
                    organism_name="Test organism",
                )
            )

            text = output.read_text(encoding="utf-8")
            self.assertNotIn("\tmRNA\t", text)
            self.assertNotIn("\texon\t", text)
            self.assertNotIn("\tintron\t", text)
            self.assertIn("\tCDS\tjoin(1..30,61..90)\tlocus_tag\tRed000000100", text)

    def test_retains_mrna_without_cds(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            fasta = base / "genome.fa"
            gff = base / "annotation.gff3"
            annotation = base / "annotation.tsv"
            output = base / "output.ann.txt"
            fasta.write_text(">chr1\n" + "A" * 120 + "\n", encoding="utf-8")
            gff.write_text(
                "\n".join(
                    [
                        "##gff-version 3",
                        "chr1\tsrc\tgene\t1\t90\t.\t+\t.\tID=g1",
                        "chr1\tsrc\tmRNA\t1\t90\t.\t+\t.\tID=tx1;Parent=g1",
                        "chr1\tsrc\texon\t1\t30\t.\t+\t.\tID=ex1;Parent=tx1",
                        "chr1\tsrc\tintron\t31\t60\t.\t+\t.\tID=int1;Parent=tx1",
                        "chr1\tsrc\texon\t61\t90\t.\t+\t.\tID=ex2;Parent=tx1",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            annotation.write_text("ID\tDescription\n", encoding="utf-8")

            convert_gff_to_mss(
                ConversionOptions(
                    fasta_path=fasta,
                    gff_path=gff,
                    annotation_path=annotation,
                    output_path=output,
                    locus_tag_prefix="Rna",
                    organism_name="Test organism",
                )
            )

            text = output.read_text(encoding="utf-8")
            self.assertIn(
                "\tmRNA\tjoin(1..30,61..90)\tlocus_tag\tRna000000100",
                text,
            )
            self.assertIn("\texon\t1..30\tlocus_tag\tRna000000100", text)
            self.assertIn("\tintron\t31..60\tlocus_tag\tRna000000100", text)
            self.assertNotIn("\tCDS\t", text)

    def test_retains_mrna_for_alternative_isoforms(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            fasta = base / "genome.fa"
            gff = base / "annotation.gff3"
            annotation = base / "annotation.tsv"
            output = base / "output.ann.txt"
            fasta.write_text(">chr1\n" + "A" * 120 + "\n", encoding="utf-8")
            gff.write_text(
                "\n".join(
                    [
                        "##gff-version 3",
                        "chr1\tsrc\tgene\t1\t120\t.\t+\t.\tID=g1",
                        "chr1\tsrc\tmRNA\t1\t90\t.\t+\t.\tID=tx1;Parent=g1",
                        "chr1\tsrc\texon\t1\t90\t.\t+\t.\tID=ex1;Parent=tx1",
                        "chr1\tsrc\tCDS\t1\t90\t.\t+\t0\tID=cds1;Parent=tx1",
                        "chr1\tsrc\tmRNA\t31\t120\t.\t+\t.\tID=tx2;Parent=g1",
                        "chr1\tsrc\texon\t31\t120\t.\t+\t.\tID=ex2;Parent=tx2",
                        "chr1\tsrc\tCDS\t31\t120\t.\t+\t0\tID=cds2;Parent=tx2",
                        "",
                    ]
                ),
                encoding="utf-8",
            )
            annotation.write_text(
                "ID\tDescription\n"
                "tx1\tisoform 1\n"
                "tx2\tisoform 2\n",
                encoding="utf-8",
            )

            convert_gff_to_mss(
                ConversionOptions(
                    fasta_path=fasta,
                    gff_path=gff,
                    annotation_path=annotation,
                    output_path=output,
                    locus_tag_prefix="Alt",
                    organism_name="Test organism",
                )
            )

            text = output.read_text(encoding="utf-8")
            self.assertEqual(text.count("\tmRNA\t"), 2)
            self.assertIn("\tmRNA\t1..90\tlocus_tag\tAlt000000100", text)
            self.assertIn("\tmRNA\t31..120\tlocus_tag\tAlt000000100", text)
            self.assertEqual(text.count("\tCDS\t"), 2)


if __name__ == "__main__":
    unittest.main()
