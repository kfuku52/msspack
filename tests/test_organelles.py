from __future__ import annotations

import json
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition, CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from msspack.config import ConfigError, load_config
from msspack.organelles import prepare_organelles
from msspack.utils import MSSPackError
from msspack.workflow import run_all


class OrganelleTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        source = Path(__file__).parent / "fixtures" / "minimal_pack"
        shutil.copytree(source, self.root, dirs_exist_ok=True)
        self.config = self.root / "config.toml"
        self.base = self.config.read_text()

    def genbank(self, name: str = "annotated.gb", *, topology: str = "circular") -> Path:
        record = SeqRecord(Seq("ATG" + "GCT" * 8 + "TAA" + "A" * 170), id="ctg1")
        record.annotations.update(molecule_type="DNA", topology=topology)
        location = SimpleLocation(0, 30, strand=1)
        record.features = [
            SeqFeature(location, type="gene", qualifiers={"gene": ["abc"]}),
            SeqFeature(location, type="CDS", qualifiers={"gene": ["abc"], "product": ["test protein"],
                                                        "transl_table": ["1"], "translation": ["MAAAAAAAA"]}),
        ]
        path = self.root / name
        SeqIO.write([record], path, "genbank")
        return path

    def configure(self, text: str) -> None:
        self.config.write_text(self.base + "\n" + text)

    def test_genbank_preserves_sequences_features_and_separates_misc(self) -> None:
        self.genbank()
        self.configure('''[[organelles]]
name = "mt"
organelle = "mitochondrion"
genbank = "annotated.gb"
completeness = { ctg1 = "complete" }
''')
        config = load_config(self.config)
        nuclear = config.output_dir / "final" / "preserved.txt"
        nuclear.parent.mkdir(parents=True)
        nuclear.write_text("completed nuclear results")
        result = prepare_organelles(self.config, validate=False)
        annotation = result.ann_path.read_text()
        self.assertNotIn("\tDATATYPE\t", annotation)
        self.assertNotIn("\tgene\t1..30\t", annotation)
        self.assertNotIn("\tsubmitter_seqid\t", annotation)
        self.assertNotIn("\tkeyword\tWGS", annotation)
        self.assertIn("mt_ctg1\tTOPOLOGY\t\tcircular\t\n", annotation)
        self.assertIn("\t\t\torganelle\tmitochondrion\n", annotation)
        self.assertIn("\tCDS\t1..30\t", annotation)
        self.assertNotIn("\ttranslation\t", annotation)
        self.assertEqual(annotation.count("locus_tag\t"), 1)
        self.assertEqual(nuclear.read_text(), "completed nuclear results")
        self.assertEqual(result.final.parent, config.output_dir / "organelles")
        self.assertIn("ATGGCTGCT", result.fasta_path.read_text())
        summary = json.loads(result.summary_path.read_text())
        self.assertEqual(summary["records"][0]["feature_counts"], {"CDS": 1})
        self.assertEqual(summary["records"][0]["input_feature_counts"], {"gene": 1, "CDS": 1})
        self.assertFalse(summary["nuclear_analysis_executed"])
        self.assertEqual(summary["validation_status"], "not_run")
        audit = json.loads((result.final / "feature-audit.json").read_text())
        self.assertEqual(audit[1]["input_qualifiers"]["translation"], ["MAAAAAAAA"])
        self.assertIn("<svg", result.report_path.read_text())
        self.assertIn("Pipeline gene flow", result.report_path.read_text())
        self.assertTrue((result.final / "organelle-mt-gene-flow.sankey.pdf").is_file())
        flows = json.loads((result.final / "organelle-gene-flows.json").read_text())
        self.assertEqual(flows[0]["total"], 1)
        self.assertTrue((result.final / flows[0]["svg"]).is_file())

    def test_same_sequence_names_in_different_organelles_are_namespaced(self) -> None:
        self.genbank()
        self.genbank("linear.gb", topology="linear")
        self.configure('''[[organelles]]
name = "mt"
organelle = "mitochondrion"
genbank = "annotated.gb"
[[organelles]]
name = "pt"
organelle = "plastid:chloroplast"
genbank = "linear.gb"
''')
        result = prepare_organelles(self.config, validate=False)
        annotation = result.ann_path.read_text()
        self.assertEqual(annotation.count("\tTOPOLOGY\t"), 1)
        self.assertIn("pt_ctg1\tsource\t", annotation)
        self.assertEqual(len(json.loads(result.summary_path.read_text())["records"]), 2)

    def test_source_only_is_explicitly_reported_and_topology_is_not_guessed(self) -> None:
        (self.root / "org.fa").write_text(">ctg1 circular=true\naCGTacgtN\n")
        self.configure('''[[organelles]]
name = "mt"
organelle = "mitochondrion"
fasta = "org.fa"
''')
        with self.assertRaisesRegex(MSSPackError, "Specify topology"):
            prepare_organelles(self.config, validate=False)
        with self.config.open("a") as handle:
            handle.write('topology = { ctg1 = "linear" }\n')
        result = prepare_organelles(self.config, validate=False)
        self.assertIn("aCGTacgtN", result.fasta_path.read_text())
        warnings = json.loads(result.summary_path.read_text())["warnings"]
        self.assertTrue(any("source-only" in value for value in warnings))
        self.assertTrue(any("completeness" in value for value in warnings))

    def test_failed_validation_keeps_previous_published_pair(self) -> None:
        self.genbank()
        self.configure('''[[organelles]]
name = "mt"
organelle = "mitochondrion"
genbank = "annotated.gb"
''')
        first = prepare_organelles(self.config, validate=False)
        generation = first.final.resolve()
        original = first.ann_path.read_bytes()
        with patch("msspack.organelles.run_validation", side_effect=MSSPackError("validation failed")):
            with self.assertRaisesRegex(MSSPackError, "validation failed"):
                prepare_organelles(self.config)
        self.assertEqual(first.final.resolve(), generation)
        self.assertEqual(first.ann_path.read_bytes(), original)

    def test_gff_joined_negative_cds_uses_phase_without_moving_coordinates(self) -> None:
        (self.root / "org.fa").write_text(">ctg1\n" + "A" * 100 + "\n")
        (self.root / "org.gff").write_text(
            "ctg1\t.\tgene\t1\t100\t.\t-\t.\tID=g1;gene=abc\n"
            "ctg1\t.\tmRNA\t1\t100\t.\t-\t.\tID=t1;Parent=g1\n"
            "ctg1\t.\tCDS\t1\t8\t.\t-\t1\tID=c1;Parent=t1;product=test\n"
            "ctg1\t.\tCDS\t30\t40\t.\t-\t0\tID=c2;Parent=t1;product=test\n")
        self.configure('''[[organelles]]
name = "mt"
organelle = "mitochondrion"
fasta = "org.fa"
gff = "org.gff"
genetic_code = "1"
topology = { ctg1 = "linear" }
''')
        result = prepare_organelles(self.config, validate=False)
        self.assertIn("\tCDS\tcomplement(join(1..8,30..40))\t", result.ann_path.read_text())
        self.assertIn("codon_start\t1\n", result.ann_path.read_text())

    def test_circular_join_rna_editing_and_partial_locations_are_preserved(self) -> None:
        path = self.genbank()
        record = SeqIO.read(path, "genbank")
        record.features = [SeqFeature(
            CompoundLocation([SimpleLocation(180, 200, strand=1),
                              SimpleLocation(BeforePosition(0), 10, strand=1)]),
            type="CDS", qualifiers={"exception": ["RNA editing"], "trans_splicing": [""],
                                    "product": ["test protein"], "transl_table": ["1"],
                                    "translation": ["MKKKKKKKKK"]})]
        SeqIO.write(record, path, "genbank")
        self.configure('''[[organelles]]
name = "mt"
organelle = "mitochondrion"
genbank = "annotated.gb"
''')
        result = prepare_organelles(self.config, validate=False)
        self.assertIn("\tCDS\tjoin(181..200,<1..10)\t", result.ann_path.read_text())
        self.assertIn("exception\tRNA editing\n", result.ann_path.read_text())
        self.assertIn("trans_splicing\t\n", result.ann_path.read_text())
        self.assertIn("translation\tMKKKKKKKKK\n", result.ann_path.read_text())

    def test_invalid_configuration_and_stale_sequence_metadata_are_rejected(self) -> None:
        self.genbank()
        for text, message in [
            ('name="../unsafe"\norganelle="mitochondrion"\ngenbank="annotated.gb"', "organelle.name"),
            ('name="mt"\norganelle="wrong"\ngenbank="annotated.gb"', "organelle.organelle"),
            ('name="mt"\norganelle="mitochondrion"\ngenbank="annotated.gb"\nfasta="org.fa"', "either fasta"),
            ('name="mt"\norganelle="mitochondrion"\ngenbank="annotated.gb"\ngenetic_code="9999"', "genetic code"),
        ]:
            with self.subTest(text=text):
                self.configure("[[organelles]]\n" + text)
                with self.assertRaisesRegex(ConfigError, message):
                    load_config(self.config)
        self.configure('''[[organelles]]
name = "mt"
organelle = "mitochondrion"
genbank = "annotated.gb"
topology = { absent = "linear" }
''')
        with self.assertRaisesRegex(MSSPackError, "Unknown sequence ID"):
            prepare_organelles(self.config, validate=False)

    def test_run_keeps_nuclear_pair_compatible_and_links_organelle_report(self) -> None:
        self.genbank()
        self.configure('''[[organelles]]
name = "mt"
organelle = "mitochondrion"
genbank = "annotated.gb"
''')
        result = run_all(self.config, run_busco=False, validate=False)
        self.assertEqual(result.pipeline.ann_path.read_text(), (self.root / "expected.ann.txt").read_text())
        self.assertEqual(result.pipeline.fasta_path.read_text(), (self.root / "expected.fasta").read_text())
        self.assertIsNotNone(result.report)
        assert result.report is not None
        self.assertIn("../organelles/final/index.html", result.report.index_html.read_text())
        manifest = json.loads(result.pipeline.manifest_path.read_text())
        self.assertTrue(manifest["organelles"]["enabled"])

    def test_duplicate_gene_names_at_distinct_loci_get_distinct_tags(self) -> None:
        path = self.genbank()
        record = SeqIO.read(path, "genbank")
        loc = SimpleLocation(100, 130, strand=1)
        record.features.extend([
            SeqFeature(loc, type="gene", qualifiers={"gene": ["abc"]}),
            SeqFeature(loc, type="CDS", qualifiers={"gene": ["abc"], "product": ["test protein"],
                                                   "transl_table": ["1"]}),
        ])
        SeqIO.write(record, path, "genbank")
        self.configure('''[[organelles]]
name = "pt"
organelle = "plastid:chloroplast"
genbank = "annotated.gb"
''')
        result = prepare_organelles(self.config, validate=False)
        tags = [line.split("\t")[4] for line in result.ann_path.read_text().splitlines()
                if line.split("\t")[3] == "locus_tag"]
        self.assertEqual(len(tags), 2)
        self.assertEqual(len(set(tags)), 2)

    def test_nuclear_locus_tag_collision_prevents_publication(self) -> None:
        path = self.genbank()
        record = SeqIO.read(path, "genbank")
        for feature in record.features:
            feature.qualifiers["locus_tag"] = ["Fix_shared"]
        SeqIO.write(record, path, "genbank")
        self.configure('''[[organelles]]
name = "mt"
organelle = "mitochondrion"
genbank = "annotated.gb"
''')
        nuclear = self.root / "nuclear.ann"
        nuclear.write_text("\tCDS\t1..30\tlocus_tag\tFix_shared\n")
        with self.assertRaisesRegex(MSSPackError, "Nuclear/organelle locus_tag collision"):
            prepare_organelles(self.config, validate=False, nuclear_annotation=nuclear)
        self.assertFalse((load_config(self.config).output_dir / "organelles/final").exists())


if __name__ == "__main__":
    unittest.main()
