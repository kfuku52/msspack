import json
import tempfile
import unittest
from pathlib import Path

from msspack.pipeline_plot_data import (
    load_functional_annotation_summary,
    load_pipeline_validation_summary,
)
from msspack.pipeline_plot_models import (
    AnnotationConsistencyGroup,
    AnnotationConsistencySummary,
    FunctionalAnnotationGroup,
    FunctionalAnnotationSummary,
    PipelineGeneSet,
    PipelinePlotMetrics,
    SankeyBuscoSummary,
    SankeyLink,
    SankeyNode,
)
from msspack.pipeline_plot_render import (
    SANKEY_HEIGHT,
    _sankey_layout,
    build_sankey,
    load_sankey_busco_summaries,
    write_sankey_pdf,
    write_sankey_svg,
)
from msspack.pipeline_plots import (
    parse_pipeline_plot_metrics,
)
from msspack.validation import ValidationCheckResult, ValidationSummary


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


def _write_busco_comparison(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "comparison_name": "cds",
                "summaries": [
                    {
                        "label": "input",
                        "input_sequence_count": 95,
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
                        "input_sequence_count": 95,
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
    def test_load_pipeline_validation_summary_supports_explicit_validate_output(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_root = Path(tmp_dir)
            summary_path = (
                output_root / "final" / "validation" / "ddbj-validation-summary.json"
            )
            summary_path.parent.mkdir(parents=True)
            summary_path.write_text(
                json.dumps(
                    ValidationSummary(
                        path=summary_path,
                        status="passed",
                        checks=tuple(
                            ValidationCheckResult(
                                component=component,
                                label=label,
                                status="passed",
                                version=None,
                                log_path=None,
                                output_paths={},
                            )
                            for component, label in (
                                ("parser", "Parser"),
                                ("transchecker", "transChecker"),
                            )
                        ),
                    ).to_dict()
                ),
                encoding="utf-8",
            )

            summary = load_pipeline_validation_summary(output_root)

        self.assertIsNotNone(summary)
        assert summary is not None
        self.assertEqual(summary.path, summary_path)
        self.assertEqual(summary.status, "passed")

    def test_sankey_renders_ddbj_validation_band_in_svg_and_pdf(self) -> None:
        stage_labels = ["Input", "Final feature fate"]
        nodes = [
            SankeyNode("start", "Input genes", 0, 2, "#334155"),
            SankeyNode("final_cds", "Final gene models", 1, 2, "#059669"),
        ]
        links = [SankeyLink("start", "final_cds", 2, "#059669")]
        summary = ValidationSummary(
            path=Path("/tmp/ddbj-validation-summary.json"),
            status="passed",
            checks=(
                ValidationCheckResult(
                    component="parser",
                    label="Parser",
                    status="passed",
                    version="6.80",
                    log_path=Path("/tmp/parser.log"),
                    output_paths={},
                ),
                ValidationCheckResult(
                    component="transchecker",
                    label="transChecker",
                    status="passed",
                    version="2.26",
                    log_path=Path("/tmp/transchecker.log"),
                    output_paths={},
                    record_counts={"aa_fasta": 2, "nuc_fasta": 2},
                ),
            ),
        )
        with tempfile.TemporaryDirectory() as tmp_dir:
            svg_path = Path(tmp_dir) / "validation.svg"
            pdf_path = Path(tmp_dir) / "validation.pdf"
            write_sankey_svg(
                stage_labels,
                nodes,
                links,
                svg_path,
                validation_summary=summary,
            )
            write_sankey_pdf(
                stage_labels,
                nodes,
                links,
                pdf_path,
                validation_summary=summary,
            )
            svg_text = svg_path.read_text(encoding="utf-8")
            pdf_text = pdf_path.read_bytes().decode("latin-1")

        self.assertIn("DDBJ official validation (final MSS files)", svg_text)
        self.assertIn("Parser v6.80", svg_text)
        self.assertIn("transChecker v2.26", svg_text)
        self.assertEqual(svg_text.count(">PASS</text>"), 2)
        self.assertIn("Translated CDS: AA 2 / nucleotide 2", svg_text)
        self.assertEqual(pdf_text.count("/FontFile2"), 2)

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
            self.assertLessEqual(
                link.source_y + link.height, link.source.y + link.source.height + 1e-6
            )
            self.assertGreaterEqual(link.target_y, link.target.y - 1e-6)
            self.assertLessEqual(
                link.target_y + link.height, link.target.y + link.target.height + 1e-6
            )

    def test_sankey_includes_functional_annotation_outcomes(self) -> None:
        metrics = PipelinePlotMetrics(
            initial_genes=5,
            duplicate_removed_genes=1,
            genes_after_dedup=4,
            transcript_changed_genes=2,
            transcript_unchanged_genes=2,
            removed_mrnas=2,
            genes_after_single_mrna=4,
            inframe_updated_genes=1,
            inframe_unchanged_genes=3,
            genes_after_inframe=4,
            padding_updated_genes=1,
            genes_with_stops=1,
            padding_unchanged_genes=2,
            genes_after_padding=4,
            converted_to_misc_genes=1,
            final_cds_genes=3,
            total_cds_input=4,
            total_cds_output=3,
            misc_feature_output=1,
            sources={},
        )
        gene_sets = _sankey_gene_sets(
            {
                "duplicate_removed_genes": ["d1"],
                "transcript_changed_genes": ["g1", "g2"],
                "inframe_updated_genes": ["g2"],
                "padding_updated_genes": ["g3"],
                "genes_with_stops": ["g4"],
                "converted_to_misc_genes": ["g4"],
            }
        )
        summary = FunctionalAnnotationSummary(
            path=Path("/tmp/functional-annotation.tsv"),
            groups=(
                FunctionalAnnotationGroup(
                    key="swissprot",
                    source="swissprot",
                    label="Swiss-Prot assigned",
                    color="#16a34a",
                    locus_tags=("g1", "g3"),
                ),
                FunctionalAnnotationGroup(
                    key="pfam",
                    source="pfam",
                    label="Pfam fallback",
                    color="#7c3aed",
                    locus_tags=("g2",),
                ),
                FunctionalAnnotationGroup(
                    key="none",
                    source="none",
                    label="Unannotated",
                    color="#94a3b8",
                    locus_tags=("g4",),
                ),
            ),
        )

        stage_labels, nodes, links = build_sankey(metrics, gene_sets, summary)
        _laid_out_nodes, _laid_out_links, meta = _sankey_layout(stage_labels, nodes, links)

        self.assertEqual(stage_labels[-2:], ["Functional annotation", "Final feature fate"])
        self.assertEqual(len(stage_labels), 7)
        self.assertAlmostEqual(meta["width"], 7.2 * 72.0)
        node_counts = {node.id: node.count for node in nodes}
        self.assertEqual(node_counts["annotation_swissprot"], 2)
        self.assertEqual(node_counts["annotation_pfam"], 1)
        self.assertEqual(node_counts["annotation_none"], 1)
        link_counts = {(link.source, link.target): link.count for link in links}
        self.assertEqual(link_counts[("annotation_swissprot", "final_cds")], 2)
        self.assertEqual(link_counts[("annotation_pfam", "final_cds")], 1)
        self.assertEqual(link_counts[("annotation_none", "final_misc")], 1)

        consistency = AnnotationConsistencySummary(
            path=Path("/tmp/functional-annotation-consistency.tsv"),
            summary_path=Path("/tmp/functional-annotation-consistency-summary.tsv"),
            source_pair_path=Path("/tmp/functional-annotation-source-pairs.tsv"),
            groups=(
                AnnotationConsistencyGroup(
                    key="consistent",
                    label="Consistent",
                    color="#1d4ed8",
                    locus_tags=("g1", "g2"),
                ),
                AnnotationConsistencyGroup(
                    key="resolved",
                    label="Auto-resolved family variation",
                    color="#d97706",
                    locus_tags=("g3",),
                ),
                AnnotationConsistencyGroup(
                    key="unannotated",
                    label="Unannotated",
                    color="#cbd5e1",
                    locus_tags=("g4",),
                ),
            ),
        )
        stage_labels, nodes, links = build_sankey(
            metrics,
            gene_sets,
            summary,
            consistency,
        )
        _laid_out_nodes, _laid_out_links, meta = _sankey_layout(stage_labels, nodes, links)

        self.assertEqual(
            stage_labels[-2:],
            ["Functional annotation", "Final feature fate"],
        )
        self.assertAlmostEqual(meta["width"], 7.2 * 72.0)
        self.assertFalse(any(node.id.startswith("consistency_") for node in nodes))
        link_counts = {(link.source, link.target): link.count for link in links}
        self.assertEqual(link_counts[("annotation_swissprot", "final_cds")], 2)
        self.assertEqual(link_counts[("annotation_pfam", "final_cds")], 1)
        self.assertEqual(link_counts[("annotation_none", "final_misc")], 1)
        with tempfile.TemporaryDirectory() as tmp_dir:
            svg_path = Path(tmp_dir) / "consistency.svg"
            pdf_path = Path(tmp_dir) / "consistency.pdf"
            write_sankey_svg(
                stage_labels,
                nodes,
                links,
                svg_path,
                annotation_consistency=consistency,
            )
            write_sankey_pdf(
                stage_labels,
                nodes,
                links,
                pdf_path,
                annotation_consistency=consistency,
            )
            svg_text = svg_path.read_text(encoding="utf-8")
        self.assertIn("Name consistency (genes)", svg_text)
        self.assertIn("Close family peer threshold", svg_text)
        self.assertIn("&gt;=70% identity / &gt;=80% mutual coverage", svg_text)
        self.assertIn("Consistent 2 (50.0%)", svg_text)
        self.assertIn("Auto-resolved family variation 1 (25.0%)", svg_text)
        self.assertIn(">Adjusted</tspan>", svg_text)
        self.assertIn(">No adjustment</tspan>", svg_text)

        busco_summaries = (
            SankeyBuscoSummary(
                label="Input CDS",
                stage=0,
                lineage_dataset="embryophyta_odb12",
                input_sequences=4,
                total_buscos=100,
                single_copy=70,
                duplicated=20,
                fragmented=5,
                missing=5,
            ),
            SankeyBuscoSummary(
                label="Boundary-adjusted CDS",
                stage=4,
                lineage_dataset="embryophyta_odb12",
                input_sequences=4,
                total_buscos=100,
                single_copy=85,
                duplicated=5,
                fragmented=5,
                missing=5,
            ),
        )
        with tempfile.TemporaryDirectory() as tmp_dir:
            svg_path = Path(tmp_dir) / "combined.svg"
            pdf_path = Path(tmp_dir) / "combined.pdf"
            write_sankey_svg(
                stage_labels,
                nodes,
                links,
                svg_path,
                busco_summaries=busco_summaries,
                annotation_consistency=consistency,
            )
            write_sankey_pdf(
                stage_labels,
                nodes,
                links,
                pdf_path,
                busco_summaries=busco_summaries,
                annotation_consistency=consistency,
            )
            combined_svg = svg_path.read_text(encoding="utf-8")
        self.assertEqual(combined_svg.count('class="summary-pie-title"'), 3)
        self.assertIn("Input CDS BUSCO", combined_svg)
        self.assertIn("Boundary-adjusted CDS BUSCO", combined_svg)
        self.assertIn("BUSCO genes n=100", combined_svg)
        self.assertIn("CDS input n=4", combined_svg)
        self.assertIn("embryophyta_odb12", combined_svg)
        self.assertIn("Name consistency (n=4)", combined_svg)
        self.assertIn("Close family peer", combined_svg)
        self.assertIn("id&gt;=70%, cov&gt;=80%", combined_svg)
        self.assertIn("Consistent 50.0%", combined_svg)
        self.assertIn("Auto-resolved family", combined_svg)
        self.assertIn("variation 25.0%", combined_svg)

    def test_load_functional_annotation_summary_groups_evidence_sources(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_root = Path(tmp_dir)
            final_dir = output_root / "final"
            final_dir.mkdir()
            (final_dir / "functional-annotation.tsv").write_text(
                "ID\tLocus_tag\tsource\ng1.t1\tg1\tswissprot\ng2.t1\tg2\tpfam\ng3.t1\tg3\tnone\n",
                encoding="utf-8",
            )

            summary = load_functional_annotation_summary(output_root)

        self.assertIsNotNone(summary)
        assert summary is not None
        self.assertEqual(summary.total, 3)
        self.assertEqual(
            {group.key: group.count for group in summary.groups},
            {"swissprot": 1, "pfam": 1, "none": 1},
        )

    def test_load_functional_annotation_summary_maps_custom_locus_tags_to_gene_ids(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_root = Path(tmp_dir)
            final_dir = output_root / "final"
            intermediate_dir = output_root / "intermediate"
            final_dir.mkdir()
            intermediate_dir.mkdir()
            (intermediate_dir / "12.gff.final-sorted.gff").write_text(
                "\n".join(
                    [
                        "##gff-version 3",
                        "chr1\t.\tgene\t1\t90\t.\t+\t.\tID=g1",
                        "chr1\t.\tmRNA\t1\t90\t.\t+\t.\tID=g1.t1;Parent=g1",
                        "chr1\t.\tgene\t101\t190\t.\t+\t.\tID=g2",
                        "chr1\t.\tmRNA\t101\t190\t.\t+\t.\tID=g2.t1;Parent=g2",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            (final_dir / "functional-annotation.tsv").write_text(
                "ID\tLocus_tag\tsource\n"
                "g1.t1\tDemo_g1\tswissprot\n"
                "g2.t1\tDemo_g2\tpfam\n",
                encoding="utf-8",
            )

            summary = load_functional_annotation_summary(output_root)

        self.assertIsNotNone(summary)
        assert summary is not None
        self.assertEqual(
            {group.key: group.locus_tags for group in summary.groups},
            {"swissprot": ("g1",), "pfam": ("g2",)},
        )

    def test_load_sankey_busco_summaries_maps_results_to_measured_stages(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_root = Path(tmp_dir)
            _write_busco_comparison(output_root / "busco" / "cds" / "comparison.json")

            summaries = load_sankey_busco_summaries(output_root)

        self.assertEqual(
            [summary.label for summary in summaries],
            ["Input CDS", "Boundary-adjusted CDS"],
        )
        self.assertEqual([summary.stage for summary in summaries], [0, 4])
        self.assertEqual([summary.input_sequences for summary in summaries], [95, 95])
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

if __name__ == "__main__":
    unittest.main()
