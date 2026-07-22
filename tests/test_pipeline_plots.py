import json
import tempfile
import unittest
from pathlib import Path

from msspack.chart_primitives import pdf_helvetica_text_width
from msspack.pipeline_plot_models import PipelineGeneSet, PipelinePlotMetrics
from msspack.pipeline_plot_render import (
    SANKEY_HEIGHT,
    _sankey_layout,
    build_sankey,
    load_sankey_busco_summaries,
)
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


def _write_busco_comparison(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "comparison_name": "cds",
                "summaries": [
                    {
                        "label": "input",
                        "lineage_dataset": "embryophyta_odb12",
                        "counts": {
                            "single_copy": 70,
                            "duplicated": 20,
                            "fragmented": 5,
                            "missing": 5,
                            "total_buscos": 100,
                        },
                    },
                    {
                        "label": "processed",
                        "lineage_dataset": "embryophyta_odb12",
                        "counts": {
                            "single_copy": 85,
                            "duplicated": 5,
                            "fragmented": 5,
                            "missing": 5,
                            "total_buscos": 100,
                        },
                    },
                ],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


def _sankey_gene_sets(groups: dict[str, list[str]]) -> tuple[PipelineGeneSet, ...]:
    keys = (
        "duplicate_removed_genes",
        "transcript_changed_genes",
        "inframe_updated_genes",
        "padding_updated_genes",
        "genes_with_stops",
        "converted_to_misc_genes",
    )
    return tuple(
        PipelineGeneSet(
            key=key,
            label=key,
            color="#000000",
            path=Path(f"/tmp/{key}.txt"),
            gene_ids=tuple(groups.get(key, [])),
        )
        for key in keys
    )


class PipelinePlotTests(unittest.TestCase):
    def test_pdf_helvetica_text_width_uses_real_character_widths(self) -> None:
        self.assertAlmostEqual(
            pdf_helvetica_text_width("Input", size=8, bold=True),
            19.552,
        )

    def test_sankey_layout_keeps_large_gene_flows_inside_plot_area(self) -> None:
        metrics = PipelinePlotMetrics(
            initial_genes=22_005,
            duplicate_removed_genes=4,
            genes_after_dedup=22_001,
            transcript_changed_genes=4_220,
            transcript_unchanged_genes=17_781,
            removed_mrnas=4_220,
            genes_after_single_mrna=22_001,
            inframe_updated_genes=3,
            inframe_unchanged_genes=21_998,
            genes_after_inframe=22_001,
            padding_updated_genes=0,
            genes_with_stops=0,
            padding_unchanged_genes=22_001,
            genes_after_padding=22_001,
            converted_to_misc_genes=0,
            final_cds_genes=22_001,
            total_cds_input=22_001,
            total_cds_output=22_001,
            misc_feature_output=0,
            sources={},
        )
        gene_sets = _sankey_gene_sets(
            {
                "duplicate_removed_genes": [f"d{index}" for index in range(4)],
                "transcript_changed_genes": [f"t{index}" for index in range(4_220)],
                "inframe_updated_genes": [f"f{index}" for index in range(3)],
            }
        )
        stage_labels, nodes, links = build_sankey(metrics, gene_sets)

        laid_out_nodes, laid_out_links, meta = _sankey_layout(stage_labels, nodes, links)

        self.assertEqual(len(stage_labels), 6)
        plot_top = meta["top"]
        plot_bottom = SANKEY_HEIGHT - 24.0
        for node_id in ("after_transcript", "after_inframe", "after_padding"):
            self.assertNotIn(node_id, laid_out_nodes)
        for node in laid_out_nodes.values():
            self.assertGreaterEqual(node.y, plot_top - 1e-6)
            self.assertLessEqual(node.y + node.height, plot_bottom + 1e-6)
        for link in laid_out_links:
            self.assertGreaterEqual(link.source_y, link.source.y - 1e-6)
            self.assertLessEqual(link.source_y + link.height, link.source.y + link.source.height + 1e-6)
            self.assertGreaterEqual(link.target_y, link.target.y - 1e-6)
            self.assertLessEqual(link.target_y + link.height, link.target.y + link.target.height + 1e-6)

    def test_load_sankey_busco_summaries_maps_results_to_measured_stages(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_root = Path(tmp_dir)
            _write_busco_comparison(output_root / "busco" / "cds" / "comparison.json")

            summaries = load_sankey_busco_summaries(output_root)

        self.assertEqual([summary.label for summary in summaries], ["Input CDS", "Adjusted CDS"])
        self.assertEqual([summary.stage for summary in summaries], [0, 4])
        self.assertEqual(summaries[0].complete_pct, 90.0)

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
                ["d1", "d2", "d3", "d4", "d5"],
            )
            _write_id_file(
                log_dir / "07.select-one-mrna.changed-gene-ids.txt",
                [f"g{index}" for index in range(1, 21)],
            )
            _write_id_file(
                log_dir / "09.update-gff-to-inframe.changed-gene-ids.txt",
                ["g18", "g19", "g21"],
            )
            _write_id_file(intermediate_dir / "11.gff.updated-genes.txt", ["g19", "g21", "g22"])
            _write_id_file(intermediate_dir / "11.gff.genes-with-stops.txt", ["g23"])
            _write_id_file(
                log_dir / "16.mss-cds-to-misc.changed-gene-ids.txt",
                ["g23"],
            )
            _write_busco_comparison(output_root / "busco" / "cds" / "comparison.json")

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
            gene_flow_tsv = artifacts.gene_flow_tsv.read_text(encoding="utf-8")
            self.assertNotIn("after_transcript", gene_flow_tsv)
            self.assertNotIn("after_inframe", gene_flow_tsv)
            self.assertNotIn("after_padding", gene_flow_tsv)
            self.assertIn("transcript_changed\tinframe_updated", gene_flow_tsv)
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
            self.assertIn(">mRNA selection</tspan>", svg_texts[0])
            self.assertNotIn(">One mRNA</tspan>", svg_texts[0])
            self.assertIn(">Frame</tspan>", svg_texts[0])
            self.assertIn(">correction</tspan>", svg_texts[0])
            self.assertIn(">Input GFF</tspan>", svg_texts[0])
            self.assertIn(">Coordinate duplicate</tspan>", svg_texts[0])
            self.assertIn(">removal</tspan>", svg_texts[0])
            self.assertIn(">CDS boundary</tspan>", svg_texts[0])
            self.assertIn(">adjustment</tspan>", svg_texts[0])
            self.assertIn(">Output ann.txt</tspan>", svg_texts[0])
            self.assertEqual(svg_texts[0].count('class="stage"'), 6)
            self.assertIn(">Already one</tspan>", svg_texts[0])
            self.assertIn(">Reduced to</tspan>", svg_texts[0])
            self.assertNotIn(">Selected<", svg_texts[0])
            self.assertNotIn(">Framed<", svg_texts[0])
            self.assertNotIn(">Pruned<", svg_texts[0])
            self.assertIn('width="7.2in"', svg_texts[0])
            self.assertIn('viewBox="0 0 518.40 420.00"', svg_texts[0])
            self.assertIn("BUSCO results (CDS; embryophyta_odb12; n=100)", svg_texts[0])
            self.assertIn("Single-copy 70.0%", svg_texts[0])
            self.assertIn(">Input CDS</text>", svg_texts[0])
            self.assertIn(">Adjusted CDS</text>", svg_texts[0])
            gene_flow_pdf = artifacts.gene_flow_pdf.read_bytes().decode("latin-1")
            self.assertIn("/MediaBox [0 0 518.40 420.00]", gene_flow_pdf)
            self.assertIn("(mRNA selection)", gene_flow_pdf)
            self.assertIn("(Coordinate duplicate)", gene_flow_pdf)
            self.assertIn("(CDS boundary)", gene_flow_pdf)
            self.assertIn("(Output ann.txt)", gene_flow_pdf)
            self.assertIn("(Adjusted CDS)", gene_flow_pdf)
            self.assertIn("(Already one)", gene_flow_pdf)
            self.assertIn("(Reduced to)", gene_flow_pdf)
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
