import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from msspack.annotation_taxonomy import (
    RELATION_WEIGHTS,
    TaxonomyRecord,
    resolve_annotation_taxonomy,
    taxonomy_relation,
)
from msspack.config_models import FunctionalAnnotationTaxonomyConfig


def _node(
    tax_id: int,
    name: str,
    rank: str,
    lineage: list[int],
) -> dict[str, object]:
    return {
        "query": [name],
        "taxonomy": {
            "tax_id": tax_id,
            "organism_name": name,
            "rank": rank,
            "lineage": lineage,
        },
    }


class AnnotationTaxonomyTests(unittest.TestCase):
    def test_taxonomy_relation_distinguishes_close_and_distant_hits(self) -> None:
        target = {
            "target": {
                "tax_id": 45172,
                "lineage": [
                    {"tax_id": 2759, "rank": "DOMAIN"},
                    {"tax_id": 33090, "rank": "KINGDOM"},
                    {"tax_id": 3398, "rank": "CLASS"},
                    {"tax_id": 4439, "rank": "GENUS"},
                    {"tax_id": 45172, "rank": "SPECIES"},
                ],
            }
        }
        same_genus = TaxonomyRecord(
            tax_id=999001,
            scientific_name="Plumbago testensis",
            rank="SPECIES",
            lineage_ids=(1, 2759, 33090, 3398, 4439),
        )
        animal = TaxonomyRecord(
            tax_id=9606,
            scientific_name="Homo sapiens",
            rank="SPECIES",
            lineage_ids=(1, 2759, 33208, 40674, 9605),
        )
        bacterium = TaxonomyRecord(
            tax_id=562,
            scientific_name="Escherichia coli",
            rank="SPECIES",
            lineage_ids=(1, 2, 1224, 561),
        )

        self.assertEqual(
            taxonomy_relation(target_context=target, subject=same_genus),
            ("same_genus", RELATION_WEIGHTS["same_genus"]),
        )
        self.assertEqual(
            taxonomy_relation(target_context=target, subject=animal),
            ("cross_kingdom", RELATION_WEIGHTS["cross_kingdom"]),
        )
        self.assertEqual(
            taxonomy_relation(target_context=target, subject=bacterium),
            ("cross_domain", RELATION_WEIGHTS["cross_domain"]),
        )

    def test_resolves_scientific_name_and_crosschecks_busco_lineage(self) -> None:
        records = {
            1: _node(1, "root", "NO_RANK", []),
            2759: _node(2759, "Eukaryota", "DOMAIN", [1]),
            33090: _node(33090, "Viridiplantae", "KINGDOM", [1, 2759]),
            71240: _node(
                71240,
                "eudicotyledons",
                "CLADE",
                [1, 2759, 33090],
            ),
            4439: _node(
                4439,
                "Plumbago",
                "GENUS",
                [1, 2759, 33090, 71240],
            ),
            45172: _node(
                45172,
                "Plumbago auriculata",
                "SPECIES",
                [1, 2759, 33090, 71240, 4439],
            ),
        }

        def fake_request(query: str) -> dict[str, object]:
            if query == "Plumbago auriculata":
                return {"taxonomy_nodes": [records[45172]]}
            if query == "eudicotyledons":
                return {"taxonomy_nodes": [records[71240]]}
            ids = [int(value) for value in query.split(",")]
            return {"taxonomy_nodes": [records[value] for value in ids]}

        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            output = base / "taxonomy.json"
            with patch(
                "msspack.annotation_taxonomy._request_taxonomy_payload",
                side_effect=fake_request,
            ):
                resolve_annotation_taxonomy(
                    scientific_name="Plumbago auriculata",
                    configured_busco_lineage="eudicotyledons_odb12",
                    busco_summary_paths=[],
                    output_path=output,
                    log_path=base / "taxonomy.log",
                    metrics_path=base / "taxonomy.metrics.json",
                    cache_dir=base / "cache",
                    config=FunctionalAnnotationTaxonomyConfig(),
                )

            payload = json.loads(output.read_text(encoding="utf-8"))
            self.assertEqual(payload["status"], "resolved")
            self.assertEqual(payload["source"], "scientific_name")
            self.assertEqual(payload["target"]["tax_id"], 45172)
            self.assertEqual(
                payload["busco_crosschecks"][0]["concordance"],
                "consistent",
            )


if __name__ == "__main__":
    unittest.main()
