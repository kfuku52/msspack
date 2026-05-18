import json
import tempfile
import unittest
from pathlib import Path

from msspack.pipeline_plots import (
    parse_pipeline_plot_metrics,
    run_pipeline_plots,
    summarize_pipeline_plots,
)


def _write_step_log(
    path: Path,
    *,
    step: str,
    unit: str,
    input_total: int,
    changed_total: int,
    output_total: int,
    details: list[str],
) -> None:
    text = "\n".join(
        [
            f"$ msspack internal {step}",
            "",
            "Summary table:",
            "| Metric | Value |",
            "| --- | --- |",
            f"| Step | {step} |",
            "| Status | completed |",
            "| Started at | 2026-04-09 10:00:00 |",
            "| Ended at | 2026-04-09 10:00:01 |",
            "| Duration seconds | 1.000 |",
            f"| Count unit | {unit} |",
            f"| Input total | {input_total:,} |",
            f"| Changed total | {changed_total:,} |",
            f"| Output total | {output_total:,} |",
            "",
            f"Step: {step}",
            "Status: completed",
            "Started at: 2026-04-09 10:00:00",
            "Ended at: 2026-04-09 10:00:01",
            "Duration seconds: 1.000",
            f"Count unit: {unit}",
            f"Input total: {input_total:,}",
            f"Changed total: {changed_total:,}",
            f"Output total: {output_total:,}",
            "",
            *details,
            "",
        ]
    )
    path.write_text(text, encoding="utf-8")


def _write_step_metrics(
    path: Path,
    *,
    step: str,
    unit: str,
    input_total: int,
    changed_total: int,
    output_total: int,
    details: dict[str, object],
) -> None:
    path.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "step": step,
                "status": "completed",
                "count_unit": unit,
                "input_total": input_total,
                "changed_total": changed_total,
                "output_total": output_total,
                "details": details,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def _write_id_file(path: Path, identifiers: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(identifiers) + ("\n" if identifiers else ""), encoding="utf-8")


class PipelinePlotTests(unittest.TestCase):
    def test_parse_pipeline_plot_metrics_extracts_stage_counts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            log_dir = Path(tmp_dir) / "logs"
            log_dir.mkdir()
            _write_step_log(
                log_dir / "06.drop-duplicate-coordinate-gene.log",
                step="drop-duplicate-coordinate-gene",
                unit="genes",
                input_total=100,
                changed_total=5,
                output_total=95,
                details=["Output feature count: 250"],
            )
            _write_step_log(
                log_dir / "07.select-one-mrna.log",
                step="select-one-mrna",
                unit="genes",
                input_total=95,
                changed_total=20,
                output_total=95,
                details=["Removed mRNAs: 23"],
            )
            _write_step_log(
                log_dir / "09.update-gff-to-inframe.log",
                step="update-gff-to-inframe",
                unit="genes",
                input_total=95,
                changed_total=3,
                output_total=95,
                details=["Number of unchanged gene models: 92"],
            )
            _write_step_log(
                log_dir / "11.update-gff-with-padding.log",
                step="update-gff-with-padding",
                unit="genes",
                input_total=95,
                changed_total=3,
                output_total=95,
                details=[
                    "Number of genes with stops (new_num_stop>0): 1",
                    "Number of updated genes (new_num_stop=0): 3",
                ],
            )
            _write_step_metrics(
                log_dir / "11.update-gff-with-padding.metrics.json",
                step="update-gff-with-padding",
                unit="genes",
                input_total=95,
                changed_total=4,
                output_total=95,
                details={
                    "updated_genes": 3,
                    "genes_with_stops": 1,
                },
            )
            _write_step_log(
                log_dir / "16.mss-cds-to-misc.log",
                step="mss-cds-to-misc",
                unit="feature blocks",
                input_total=110,
                changed_total=1,
                output_total=110,
                details=[
                    "Number of genes read from /tmp/genes-with-stops.txt: 1",
                    "Total number of CDS in input: 95",
                    "Total number of CDS in output: 94",
                    "Total number of misc_feature in output: 1",
                ],
            )

            metrics = parse_pipeline_plot_metrics(log_dir)

        self.assertEqual(metrics.initial_genes, 100)
        self.assertEqual(metrics.duplicate_removed_genes, 5)
        self.assertEqual(metrics.transcript_changed_genes, 20)
        self.assertEqual(metrics.transcript_unchanged_genes, 75)
        self.assertEqual(metrics.removed_mrnas, 23)
        self.assertEqual(metrics.inframe_updated_genes, 3)
        self.assertEqual(metrics.padding_updated_genes, 3)
        self.assertEqual(metrics.genes_with_stops, 1)
        self.assertEqual(metrics.converted_to_misc_genes, 1)
        self.assertEqual(metrics.final_cds_genes, 94)
        self.assertEqual(metrics.misc_feature_output, 1)

    def test_run_pipeline_plots_writes_outputs_and_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            config_path = base / "demo.toml"
            output_root = base / "build" / "Demo"
            log_dir = output_root / "logs"
            log_dir.mkdir(parents=True)
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
            manifest_path = output_root / "build-manifest.json"
            intermediate_dir = output_root / "intermediate"
            intermediate_dir.mkdir(parents=True)
            manifest_path.write_text('{"status":"completed"}\n', encoding="utf-8")
            _write_step_log(
                log_dir / "06.drop-duplicate-coordinate-gene.log",
                step="drop-duplicate-coordinate-gene",
                unit="genes",
                input_total=100,
                changed_total=5,
                output_total=95,
                details=["Output feature count: 250"],
            )
            _write_step_metrics(
                log_dir / "06.drop-duplicate-coordinate-gene.metrics.json",
                step="drop-duplicate-coordinate-gene",
                unit="genes",
                input_total=100,
                changed_total=5,
                output_total=95,
                details={
                    "output_feature_count": 250,
                    "removed_gene_ids_path": str(
                        log_dir / "06.drop-duplicate-coordinate-gene.changed-gene-ids.txt"
                    ),
                },
            )
            _write_step_log(
                log_dir / "07.select-one-mrna.log",
                step="select-one-mrna",
                unit="genes",
                input_total=95,
                changed_total=20,
                output_total=95,
                details=["Removed mRNAs: 23"],
            )
            _write_step_metrics(
                log_dir / "07.select-one-mrna.metrics.json",
                step="select-one-mrna",
                unit="genes",
                input_total=95,
                changed_total=20,
                output_total=95,
                details={
                    "removed_mrnas": 23,
                    "changed_gene_ids_path": str(log_dir / "07.select-one-mrna.changed-gene-ids.txt"),
                },
            )
            _write_step_log(
                log_dir / "09.update-gff-to-inframe.log",
                step="update-gff-to-inframe",
                unit="genes",
                input_total=95,
                changed_total=3,
                output_total=95,
                details=["Number of unchanged gene models: 92"],
            )
            _write_step_metrics(
                log_dir / "09.update-gff-to-inframe.metrics.json",
                step="update-gff-to-inframe",
                unit="genes",
                input_total=95,
                changed_total=3,
                output_total=95,
                details={
                    "unchanged_gene_models": 92,
                    "removed_features": 0,
                    "updated_gene_ids_path": str(log_dir / "09.update-gff-to-inframe.changed-gene-ids.txt"),
                },
            )
            _write_step_log(
                log_dir / "11.update-gff-with-padding.log",
                step="update-gff-with-padding",
                unit="genes",
                input_total=95,
                changed_total=3,
                output_total=95,
                details=[
                    "Number of genes with stops (new_num_stop>0): 1",
                    "Number of updated genes (new_num_stop=0): 3",
                ],
            )
            _write_step_metrics(
                log_dir / "11.update-gff-with-padding.metrics.json",
                step="update-gff-with-padding",
                unit="genes",
                input_total=95,
                changed_total=4,
                output_total=95,
                details={
                    "updated_genes": 3,
                    "genes_with_stops": 1,
                    "updated_genes_path": str(intermediate_dir / "11.gff.updated-genes.txt"),
                    "genes_with_stops_path": str(intermediate_dir / "11.gff.genes-with-stops.txt"),
                },
            )
            _write_step_log(
                log_dir / "16.mss-cds-to-misc.log",
                step="mss-cds-to-misc",
                unit="feature blocks",
                input_total=110,
                changed_total=1,
                output_total=110,
                details=[
                    "Number of genes read from /tmp/genes-with-stops.txt: 1",
                    "Total number of CDS in input: 95",
                    "Total number of CDS in output: 94",
                    "Total number of misc_feature in output: 1",
                ],
            )
            _write_step_metrics(
                log_dir / "16.mss-cds-to-misc.metrics.json",
                step="mss-cds-to-misc",
                unit="feature blocks",
                input_total=110,
                changed_total=1,
                output_total=110,
                details={
                    "genes_read": 1,
                    "cds_input": 95,
                    "cds_output": 94,
                    "misc_feature_output": 1,
                    "converted_gene_ids_path": str(log_dir / "16.mss-cds-to-misc.changed-gene-ids.txt"),
                },
            )
            _write_id_file(
                log_dir / "06.drop-duplicate-coordinate-gene.changed-gene-ids.txt",
                ["g1", "g2", "g3", "g4", "g5"],
            )
            _write_id_file(
                log_dir / "07.select-one-mrna.changed-gene-ids.txt",
                ["g4", "g5", "g6"],
            )
            _write_id_file(
                log_dir / "09.update-gff-to-inframe.changed-gene-ids.txt",
                ["g5", "g6", "g8"],
            )
            _write_id_file(intermediate_dir / "11.gff.updated-genes.txt", ["g6", "g8", "g9"])
            _write_id_file(intermediate_dir / "11.gff.genes-with-stops.txt", ["g10"])
            _write_id_file(
                log_dir / "16.mss-cds-to-misc.changed-gene-ids.txt",
                ["g10"],
            )

            artifacts = run_pipeline_plots(config_path)
            summary_lines = summarize_pipeline_plots(artifacts)
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))

            self.assertTrue(artifacts.summary_json.exists())
            self.assertTrue(artifacts.summary_tsv.exists())
            self.assertTrue(artifacts.gene_flow_tsv.exists())
            self.assertTrue(artifacts.gene_flow_svg.exists())
            self.assertTrue(artifacts.gene_flow_pdf.exists())
            self.assertTrue(artifacts.event_counts_tsv.exists())
            self.assertTrue(artifacts.event_counts_svg.exists())
            self.assertTrue(artifacts.event_counts_pdf.exists())
            self.assertTrue(artifacts.overlap_tsv.exists())
            self.assertTrue(artifacts.overlap_svg.exists())
            self.assertTrue(artifacts.overlap_pdf.exists())
            svg_texts = [
                artifacts.gene_flow_svg.read_text(encoding="utf-8"),
                artifacts.event_counts_svg.read_text(encoding="utf-8"),
                artifacts.overlap_svg.read_text(encoding="utf-8"),
            ]
            for svg_text in svg_texts:
                self.assertIn("font-size:8pt", svg_text)
                self.assertNotRegex(svg_text, r"font-size:\d+px")
            self.assertIn(
                "Stage-wise pipeline gene flow",
                svg_texts[0],
            )
            self.assertIn(
                "Pipeline event counts",
                svg_texts[1],
            )
            self.assertIn(
                "Changed-gene overlap",
                svg_texts[2],
            )
            payload = json.loads(artifacts.summary_json.read_text(encoding="utf-8"))
            self.assertEqual(payload["metrics"]["converted_to_misc_genes"], 1)
            self.assertEqual(payload["metrics"]["final_cds_genes"], 94)
            self.assertEqual(payload["overlap"]["row_count"], 7)
            self.assertIn("plots", manifest)
            self.assertIn("pipeline", manifest["plots"])
            self.assertEqual(
                manifest["plots"]["pipeline"]["gene_flow_pdf"],
                str(artifacts.gene_flow_pdf),
            )
            self.assertEqual(
                manifest["plots"]["pipeline"]["overlap_pdf"],
                str(artifacts.overlap_pdf),
            )
            self.assertEqual(len(summary_lines), 1)
            self.assertIn("dedup_removed=5", summary_lines[0])
            self.assertIn("misc_feature_genes=1", summary_lines[0])


if __name__ == "__main__":
    unittest.main()
