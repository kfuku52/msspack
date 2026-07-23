import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from msspack.report import run_html_report


def _write_minimal_config(base: Path) -> Path:
    config_path = base / "demo.toml"
    (base / "input.fa").write_text(">chr1\nACGT\n", encoding="utf-8")
    (base / "input.gff").write_text("##gff-version 3\n", encoding="utf-8")
    config_path.write_text(
        "\n".join(
            [
                "[project]",
                'name = "Demo"',
                'output_dir = "build/Demo"',
                "",
                "[inputs]",
                'fasta = "input.fa"',
                'gff = "input.gff"',
                "",
                "[sample]",
                'locus_tag = "Demo"',
                "locus_tag_digits = 6",
                'scientific_name = "Demo demo"',
                "",
                "[submission]",
                'datatype = "CON"',
                'hold_date = "20270401"',
                'bioproject = "PRJDB000001"',
                'biosample = "SAMD000001"',
                "",
                "[submitter]",
                'ab_name = ["A. Author"]',
                'contact = "A. Author"',
                'institute = "Demo Institute"',
                'department = "Demo Department"',
                'country = "Japan"',
                'state = "Tokyo"',
                'city = "Tokyo"',
                'street = "1 Demo Street"',
                'zip = "100-0001"',
                'phone = "+81-3-0000-0000"',
                'email = "demo@example.org"',
                "",
                "[reference]",
                'title = "Demo title"',
                'ab_name = ["A. Author"]',
                "year = 2026",
                "",
                "[st_comment]",
                'assembly_method = "demo"',
                'assembly_name = "demo-v1"',
                'genome_coverage = "10x"',
                'sequencing_technology = "ONT"',
                "",
                "[pipeline]",
                "validate_with_parser = false",
                "validate_with_transchecker = false",
                "",
                "[busco]",
                "run_cds = true",
                "run_genome = false",
                'lineage_dataset = "embryophyta_odb12"',
                "auto_lineage = false",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return config_path


class ReportTests(unittest.TestCase):
    def test_summary_cards_escape_project_html(self) -> None:
        from msspack.report import _render_summary_cards

        html = _render_summary_cards(
            Path("/tmp/report"),
            {
                "project": "<img src=x onerror=alert(1)>",
                "outputs": {},
                "stage_summary": {},
                "runtime": {},
            },
        )

        self.assertNotIn("<img src=x", html)
        self.assertIn("&lt;img src=x", html)

    def test_run_html_report_writes_index_and_updates_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            config_path = _write_minimal_config(base)
            output_root = base / "build" / "Demo"
            logs_dir = output_root / "logs"
            plots_dir = output_root / "plots"
            busco_dir = output_root / "busco" / "cds"
            final_dir = output_root / "final"
            logs_dir.mkdir(parents=True)
            plots_dir.mkdir(parents=True)
            busco_dir.mkdir(parents=True)
            final_dir.mkdir(parents=True)

            metrics_path = logs_dir / "00.copy-input-fasta.metrics.json"
            metrics_path.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "step": "copy-input-fasta",
                        "status": "completed",
                        "count_unit": "sequences",
                        "input_total": 1,
                        "changed_total": 0,
                        "output_total": 1,
                        "details": {"source_file": str(base / "input.fa")},
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n",
                encoding="utf-8",
            )
            (logs_dir / "00.copy-input-fasta.log").write_text("$ demo\n", encoding="utf-8")

            for relative in (
                "pipeline-flow-summary.json",
                "pipeline-flow-summary.tsv",
                "pipeline-gene-flow.tsv",
                "pipeline-gene-flow.sankey.svg",
                "pipeline-gene-flow.sankey.pdf",
                "pipeline-event-counts.tsv",
                "pipeline-event-counts.svg",
                "pipeline-event-counts.pdf",
                "pipeline-gene-overlap.tsv",
                "pipeline-gene-overlap.svg",
                "pipeline-gene-overlap.pdf",
                "functional-annotation-name-consistency.tsv",
                "functional-annotation-name-consistency.svg",
                "functional-annotation-name-consistency.pdf",
                "functional-annotation-source-consistency.tsv",
                "functional-annotation-source-consistency.svg",
                "functional-annotation-source-consistency.pdf",
            ):
                (plots_dir / relative).write_text("demo\n", encoding="utf-8")

            for relative in (
                "comparison.tsv",
                "comparison.svg",
                "comparison.pdf",
                "input.summary.json",
                "processed.summary.json",
                "comparison.json",
            ):
                (busco_dir / relative).write_text("demo\n", encoding="utf-8")

            ann_path = final_dir / "Demo.ann.txt"
            fasta_path = final_dir / "Demo.fasta"
            evidence_path = final_dir / "functional-annotation.tsv"
            name_standardization_path = (
                final_dir / "functional-annotation-name-standardization.tsv"
            )
            taxonomy_path = final_dir / "functional-annotation-taxonomy.json"
            taxonomy_crosscheck_path = busco_dir / "taxonomy-crosscheck.json"
            ann_path.write_text("COMMON\n", encoding="utf-8")
            fasta_path.write_text(">chr1\nACGT\n//\n", encoding="utf-8")
            evidence_path.write_text("ID\tassigned_product\n", encoding="utf-8")
            name_standardization_path.write_text(
                "kind\tcode\tcount\n",
                encoding="utf-8",
            )
            taxonomy_path.write_text('{"status": "resolved"}\n', encoding="utf-8")
            taxonomy_crosscheck_path.write_text(
                '{"busco_crosschecks": []}\n',
                encoding="utf-8",
            )

            manifest_path = output_root / "build-manifest.json"
            manifest_path.write_text(
                json.dumps(
                    {
                        "project": "Demo",
                        "status": "completed",
                        "started_at": "2026-04-14T00:00:00+00:00",
                        "finished_at": "2026-04-14T00:01:00+00:00",
                        "duration_seconds": 60.0,
                        "runtime": {"python_version": "3.12.0"},
                        "outputs": {
                            "annotation": {"path": str(ann_path)},
                            "fasta": {"path": str(fasta_path)},
                            "functional_annotation_evidence": {"path": str(evidence_path)},
                            "functional_annotation_name_standardization": {
                                "path": str(name_standardization_path)
                            },
                            "functional_annotation_taxonomy": {
                                "path": str(taxonomy_path)
                            },
                        },
                        "stage_summary": {"count": 18, "ran": 18, "reused": 0},
                        "validation": {
                            "enabled": False,
                            "outputs": {},
                        },
                        "plots": {
                            "pipeline": {
                                "summary_json": str(plots_dir / "pipeline-flow-summary.json"),
                                "summary_tsv": str(plots_dir / "pipeline-flow-summary.tsv"),
                                "gene_flow_tsv": str(plots_dir / "pipeline-gene-flow.tsv"),
                                "gene_flow_svg": str(plots_dir / "pipeline-gene-flow.sankey.svg"),
                                "gene_flow_pdf": str(plots_dir / "pipeline-gene-flow.sankey.pdf"),
                                "event_counts_tsv": str(plots_dir / "pipeline-event-counts.tsv"),
                                "event_counts_svg": str(plots_dir / "pipeline-event-counts.svg"),
                                "event_counts_pdf": str(plots_dir / "pipeline-event-counts.pdf"),
                                "overlap_tsv": str(plots_dir / "pipeline-gene-overlap.tsv"),
                                "overlap_svg": str(plots_dir / "pipeline-gene-overlap.svg"),
                                "overlap_pdf": str(plots_dir / "pipeline-gene-overlap.pdf"),
                                "annotation_consistency": {
                                    "name_consistency_tsv": str(
                                        plots_dir / "functional-annotation-name-consistency.tsv"
                                    ),
                                    "name_consistency_svg": str(
                                        plots_dir / "functional-annotation-name-consistency.svg"
                                    ),
                                    "name_consistency_pdf": str(
                                        plots_dir / "functional-annotation-name-consistency.pdf"
                                    ),
                                    "source_consistency_tsv": str(
                                        plots_dir / "functional-annotation-source-consistency.tsv"
                                    ),
                                    "source_consistency_svg": str(
                                        plots_dir / "functional-annotation-source-consistency.svg"
                                    ),
                                    "source_consistency_pdf": str(
                                        plots_dir / "functional-annotation-source-consistency.pdf"
                                    ),
                                },
                            }
                        },
                        "busco": {
                            "enabled": True,
                            "taxonomy_crosscheck": str(taxonomy_crosscheck_path),
                            "comparisons": {
                                "cds": {
                                    "comparison_tsv": str(busco_dir / "comparison.tsv"),
                                    "comparison_svg": str(busco_dir / "comparison.svg"),
                                    "comparison_pdf": str(busco_dir / "comparison.pdf"),
                                    "lineage_dataset": "embryophyta_odb12",
                                    "mode": "transcriptome",
                                    "metrics": {
                                        "input": {"complete": 95.0},
                                        "processed": {"complete": 96.0},
                                        "delta": {"complete": 1.0},
                                    },
                                    "counts": {"input": {"complete": 950}, "processed": {"complete": 960}},
                                }
                            },
                        },
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n",
                encoding="utf-8",
            )

            with patch("msspack.report.run_pipeline_plots", return_value=None) as mocked:
                artifacts = run_html_report(config_path)

            mocked.assert_called_once_with(config_path, force=False)
            self.assertTrue(artifacts.index_html.exists())
            html = artifacts.index_html.read_text(encoding="utf-8")
            self.assertIn("Demo packaging report", html)
            self.assertIn("Pipeline gene flow", html)
            self.assertIn("BUSCO cds", html)
            self.assertIn("00.copy-input-fasta", html)
            self.assertIn("functional annotation evidence", html)
            self.assertIn("functional annotation taxonomy", html)
            self.assertIn("product-name standardization", html)
            self.assertIn("BUSCO taxonomy cross-check", html)
            self.assertIn("Functional annotation name consistency", html)
            self.assertIn("Name review rate by evidence source", html)
            updated_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(updated_manifest["report"]["index_html"], str(artifacts.index_html))


if __name__ == "__main__":
    unittest.main()
