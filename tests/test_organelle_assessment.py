import copy
import csv
import hashlib
import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

from msspack.config import load_config
from msspack.config_models import OrganelleConfig
from msspack.organelle_assessment import render_assessment_html, write_organelle_assessment
from msspack.organelle_plots import write_organelle_gene_flows
from msspack.utils import MSSPackError


class OrganelleAssessmentTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        source = self.root / "annotation.gb"
        source.write_text("fixture external annotation")
        self.item = OrganelleConfig("pt", "plastid", genbank=str(source))
        self.config = replace(load_config(Path(__file__).parent / "fixtures/minimal_pack/config.toml"),
                              organelles=[self.item])
        self.records = [{"input_id": "c1", "sequence_sha256": "a" * 64, "feature_counts": {"CDS": 2}}]
        self.audit = [{"organelle": "pt", "source_record": "c1", "type": "CDS", "location": "1..90",
                       "output_qualifiers": {"gene": ["rpl2"], "product": ["ribosomal protein"], "locus_tag": [tag]}}
                      for tag in ["L1", "L2"]]
        self.payload = dict(schema_version=1, name="pt", organelle="plastid", method="External fixture assessment",
                            input_sha256=[hashlib.sha256(source.read_bytes()).hexdigest()], sequence_sha256={"c1": "a" * 64},
                            functional_classification={"assignments": [dict(entry="pt_c1", locus_tag=tag,
                                category="Translation", evidence="External supplied function") for tag in ["L1", "L2"]]},
                            gene_content={"scope": "supplied_annotations", "marker_set": dict(id="fixture", version="1",
                                definition="Single marker fixture, copied twice in an IR", references=["REF.1"]),
                                "markers": [dict(gene="rpl2", status="complete", matches=[dict(entry="pt_c1", locus_tag=tag)
                                      for tag in ["L1", "L2"]], evidence="Two full copies of the same reference marker")]})

    def run_payload(self, payload):
        p = self.root / "assessment.json"
        p.write_text(json.dumps(payload))
        self.item.assessment_json = str(p)
        return write_organelle_assessment(self.config, self.item, self.records, self.audit, self.root)

    def test_expected_ir_copies_do_not_inflate_gene_completeness(self):
        result = self.run_payload(self.payload)
        self.assertEqual(result["complete_percent"], 100)
        self.assertEqual(result["marker_count"], 1)
        self.assertEqual(result["functional_counts"], {"Translation": 2})
        self.assertTrue((self.root / result["functional_pdf"]).read_bytes().startswith(b"%PDF"))
        self.assertIn("1/1", render_assessment_html(result, self.root))

    def test_missing_assessment_is_unavailable_not_zero_percent(self):
        result = write_organelle_assessment(self.config, self.item, self.records, self.audit, self.root)
        self.assertIsNone(result["complete_percent"])
        self.assertIsNone(result["marker_count"])
        self.assertEqual(result["functional_counts"], {"Not classified": 2})
        self.assertIn("NOT ASSESSED", (self.root / result["gene_content_svg"]).read_text())

    def test_unassessed_marker_does_not_yield_a_percent(self):
        payload = copy.deepcopy(self.payload)
        payload["gene_content"]["markers"].append(dict(gene="rbcL", status="not_evaluated", matches=[], evidence="No assessment"))
        result = self.run_payload(payload)
        self.assertIsNone(result["complete_percent"])
        self.assertIsNone(result["detected_percent"])

    def test_stale_or_contradictory_external_results_are_rejected(self):
        changes = [lambda p: p.update(input_sha256=["0" * 64]),
                   lambda p: p.update(sequence_sha256={"c1": "0" * 64}),
                   lambda p: p["gene_content"]["markers"].append(copy.deepcopy(p["gene_content"]["markers"][0])),
                   lambda p: p["gene_content"]["markers"][0].update(status="not_detected"),
                   lambda p: p["functional_classification"]["assignments"][0].update(locus_tag="unknown")]
        for change in changes:
            with self.subTest(change=change):
                payload = copy.deepcopy(self.payload)
                change(payload)
                with self.assertRaises(MSSPackError):
                    self.run_payload(payload)
        for row in self.audit:
            row["location"] = "<1..90"
        with self.assertRaisesRegex(MSSPackError, "only partial"):
            self.run_payload(self.payload)

    def test_external_labels_are_escaped_in_svg_and_report(self):
        payload = copy.deepcopy(self.payload)
        payload["functional_classification"]["assignments"][0]["category"] = "<script>bad</script>"
        result = self.run_payload(payload)
        self.assertNotIn("<script>", (self.root / result["functional_svg"]).read_text())
        self.assertNotIn("<script>", render_assessment_html(result, self.root))

    def swiss_payload(self):
        payload = copy.deepcopy(self.payload)
        payload["swissprot"] = dict(database=dict(name="Swiss-Prot", version="fixture.1", fasta_sha256="b" * 64),
            method="External alignment of supplied CDS", thresholds=dict(evalue=1e-10), results=[
                dict(entry="pt_c1", locus_tag="L1", status="hit", accession="P12345", evidence="Accepted external match"),
                dict(entry="pt_c1", locus_tag="L2", status="no_hit", evidence="Searched, no match meets thresholds")])
        return payload

    def test_swissprot_requires_per_cds_evidence_and_database_provenance(self):
        payload = self.swiss_payload()
        changes = [lambda p: p["swissprot"].pop("database"),
                   lambda p: p["swissprot"]["database"].update(fasta_sha256="bad"),
                   lambda p: p["swissprot"].pop("thresholds"),
                   lambda p: p["swissprot"]["results"].pop(),
                   lambda p: p["swissprot"]["results"].append(p["swissprot"]["results"][0]),
                   lambda p: p["swissprot"]["results"][0].update(locus_tag="unknown"),
                   lambda p: p["swissprot"]["results"][0].pop("accession"),
                   lambda p: p["swissprot"]["results"][1].update(accession="P12345")]
        for change in changes:
            with self.subTest(change=change):
                bad = copy.deepcopy(payload)
                change(bad)
                with self.assertRaises(MSSPackError):
                    self.run_payload(bad)
        result = self.run_payload(payload)
        self.assertEqual(result["swissprot"]["state_counts"], {"hit": 1, "no_hit": 1, "not_evaluated": 0})
        self.assertIn("P12345", (self.root / result["swissprot_tsv"]).read_text())

    def test_swissprot_sankey_conserves_cds_rna_and_invalidates_changed_evidence(self):
        payload = self.swiss_payload()
        self.run_payload(payload)
        records = [{**self.records[0], "entry": "pt_c1", "organelle": "plastid",
                    "input_feature_counts": {"CDS": 2, "tRNA": 1}, "feature_counts": {"CDS": 2, "tRNA": 1}}]
        flow = write_organelle_gene_flows(self.config, {"records": records}, self.root, audit=self.audit)[0]
        self.assertEqual(flow["total"], 3)
        text = (self.root / flow["svg"]).read_text()
        self.assertIn("Swiss-Prot", text)
        self.assertIn("No Swiss-Prot hit", text)
        self.assertIn("applicable", text)
        tsv = (self.root / flow["tsv"]).read_text()
        self.assertIn("org_swissprot_hit", tsv)
        self.assertIn("org_annotation_RNA", tsv)
        with (self.root / flow["tsv"]).open(newline="") as stream:
            edges = list(csv.DictReader(stream, delimiter="\t"))
        self.assertTrue(all(None not in edge and None not in edge.values() for edge in edges))
        self.assertEqual(sum(int(edge["count"]) for edge in edges if edge["source"] == "CDS"), 2)
        self.assertEqual(sum(int(edge["count"]) for edge in edges if edge["target"] == "org_output"), 3)
        payload["swissprot"]["results"][1]["status"] = "not_evaluated"
        self.run_payload(payload)
        flow = write_organelle_gene_flows(self.config, {"records": records}, self.root, audit=self.audit)[0]
        self.assertEqual(flow["swissprot_counts"], {"hit": 1, "no_hit": 0, "not_evaluated": 1})
        self.assertIn("not assessed", (self.root / flow["svg"]).read_text())
        self.assertNotIn("No Swiss-Prot hit", (self.root / flow["svg"]).read_text())

    def test_existing_products_do_not_imply_swissprot_hits(self):
        result = self.run_payload(self.payload)
        self.assertEqual(result["supplied_product_count"], 2)
        self.assertEqual(result["swissprot"]["state_counts"], {"hit": 0, "no_hit": 0, "not_evaluated": 2})


if __name__ == "__main__":
    unittest.main()
