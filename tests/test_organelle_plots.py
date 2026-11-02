import hashlib
import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

from msspack.config import load_config
from msspack.config_models import OrganelleConfig
from msspack.organelle_plots import prepare_organelle_gene_flows, write_organelle_gene_flows
from msspack.utils import MSSPackError


class OrganellePlotTests(unittest.TestCase):
    def config(self, root):
        fixture = Path(__file__).parent / "fixtures/minimal_pack/config.toml"
        items = [OrganelleConfig("mt", "mitochondrion", genbank=str(root / "mt.gb")),
                 OrganelleConfig("pt", "plastid:chloroplast", genbank=str(root / "pt.gb"))]
        config = load_config(fixture)
        return replace(config, project=replace(config.project, output_dir=str(root)), organelles=items)

    def summary(self):
        return {"records": [
            {"entry": "mt_c1", "input_id": "c1", "organelle": "mitochondrion",
             "input_feature_counts": {"gene": 4, "CDS": 2, "tRNA": 1, "rRNA": 1, "exon": 8},
             "feature_counts": {"CDS": 2, "tRNA": 1, "rRNA": 1, "exon": 8}},
            {"entry": "mt_l2", "input_id": "l2", "organelle": "mitochondrion",
             "input_feature_counts": {"gene": 1, "CDS": 1}, "feature_counts": {"CDS": 1}},
            {"entry": "pt_c1", "input_id": "c1", "organelle": "plastid:chloroplast",
             "input_feature_counts": {"gene": 5, "CDS": 4, "rRNA": 1, "repeat_region": 2},
             "feature_counts": {"CDS": 4, "rRNA": 1, "repeat_region": 2}},
        ]}

    def test_groups_records_counts_features_once_and_reuses_unchanged_charts(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            config = self.config(root)
            summary = self.summary()
            final = root / "organelles/final"
            final.mkdir(parents=True)
            summary["inputs"] = []
            for item in config.organelles:
                path = Path(item.genbank)
                path.write_text("supplied external annotation " + item.name)
                summary["inputs"].append({"organelle": item.name, "path": str(path),
                                          "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
            (final / "summary.json").write_text(json.dumps(summary))
            flows = prepare_organelle_gene_flows(config)
            self.assertEqual([row["name"] for row in flows], ["mt", "pt"])
            self.assertEqual(flows[0]["input_counts"], {"CDS": 3, "tRNA": 1, "rRNA": 1})
            self.assertEqual(flows[0]["total"], 5)
            self.assertEqual(flows[0]["entries"], ["mt_c1", "mt_l2"])
            svg = Path(flows[0]["svg"])
            self.assertIn("CDS/RNA feature counts", svg.read_text())
            self.assertNotIn("BUSCO", svg.read_text())
            self.assertTrue(Path(flows[1]["pdf"]).read_bytes().startswith(b"%PDF"))
            stamp = svg.stat().st_mtime_ns
            prepare_organelle_gene_flows(config)
            self.assertEqual(svg.stat().st_mtime_ns, stamp)
            Path(config.organelles[0].genbank).write_text("changed input")
            with self.assertRaisesRegex(MSSPackError, "inputs changed"):
                prepare_organelle_gene_flows(config)

    def test_source_only_has_no_invented_gene_flow(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            config = self.config(root)
            config.organelles = config.organelles[:1]
            summary = {"records": [{"entry": "mt_c1", "input_id": "c1", "organelle": "mitochondrion",
                                    "input_feature_counts": {}, "feature_counts": {}}]}
            flow = write_organelle_gene_flows(config, summary, root)[0]
            self.assertTrue(flow["source_only"])
            self.assertEqual(flow["total"], 0)
            text = (root / flow["svg"]).read_text()
            self.assertIn("No CDS/RNA annotation supplied", text)
            self.assertNotIn("Retained", text)
            self.assertFalse(Path(flow["svg"]).is_absolute())

    def test_conversion_mismatch_or_missing_group_is_not_hidden(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            config = self.config(root)
            summary = self.summary()
            summary["records"][0]["feature_counts"]["CDS"] = 1
            with self.assertRaisesRegex(MSSPackError, "counts differ"):
                write_organelle_gene_flows(config, summary, root)
            with self.assertRaisesRegex(MSSPackError, "No prepared organelle records"):
                write_organelle_gene_flows(config, {"records": []}, root)
            with self.assertRaisesRegex(MSSPackError, "results are missing"):
                prepare_organelle_gene_flows(config)


if __name__ == "__main__":
    unittest.main()
