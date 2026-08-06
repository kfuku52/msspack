import json
import tempfile
import unittest
from pathlib import Path

from msspack.chart_primitives import pdf_helvetica_text_width
from msspack.coordinate_duplicates import (
    CoordinateDuplicatePair,
    write_coordinate_duplicate_map,
)
from msspack.pipeline_plot_data import load_functional_annotation_summary
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
    EVENT_COUNTS_SUBTITLE,
    SANKEY_HEIGHT,
    _LaidOutNode,
    _sankey_label_y,
    _sankey_layout,
    build_sankey,
    load_sankey_busco_summaries,
    write_sankey_pdf,
    write_sankey_svg,
)
from msspack.pipeline_plots import (
    parse_pipeline_plot_metrics,
    run_pipeline_plots,
    summarize_pipeline_plots,
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


def _write_id_file(path: Path, identifiers: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(identifiers) + ("\n" if identifiers else ""), encoding="utf-8")


def _write_coordinate_duplicate_fixture(
    gff_path: Path,
    map_path: Path,
    *,
    removed_count: int,
) -> None:
    lines = [
        "##gff-version 3",
        "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=k1",
        "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=k1.t1;Parent=k1",
        "chr1\tsrc\texon\t1\t100\t.\t+\t.\tID=k1.ex1;Parent=k1.t1",
        "chr1\tsrc\tCDS\t1\t100\t.\t+\t0\tID=k1.cds1;Parent=k1.t1",
    ]
    pairs: list[CoordinateDuplicatePair] = []
    for index in range(1, removed_count + 1):
        lines.extend(
            [
                f"chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=d{index}",
                (
                    f"chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\t"
                    f"ID=d{index}.t1;Parent=d{index}"
                ),
                (
                    f"chr1\tsrc\texon\t1\t40\t.\t+\t.\t"
                    f"ID=d{index}.ex1;Parent=d{index}.t1"
                ),
                (
                    f"chr1\tsrc\tCDS\t1\t40\t.\t+\t0\t"
                    f"ID=d{index}.cds1;Parent=d{index}.t1"
                ),
                (
                    f"chr1\tsrc\texon\t60\t100\t.\t+\t.\t"
                    f"ID=d{index}.ex2;Parent=d{index}.t1"
                ),
                (
                    f"chr1\tsrc\tCDS\t60\t100\t.\t+\t0\t"
                    f"ID=d{index}.cds2;Parent=d{index}.t1"
                ),
            ]
        )
        pairs.append(
            CoordinateDuplicatePair(
                group_id="coordinate_duplicate_000001",
                seqid="chr1",
                start=1,
                end=100,
                strand="+",
                kept_gene_id="k1",
                removed_gene_id=f"d{index}",
            )
        )
    gff_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    write_coordinate_duplicate_map(map_path, pairs)


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
    def test_sankey_renders_ddbj_validation_band_in_svg_and_pdf(self) -> None:
        stage_labels = ["Input", "Final feature fate"]
        nodes = [
            SankeyNode("start", "Input genes", 0, 2, "#334155"),
            SankeyNode("final_cds", "Final CDS genes", 1, 2, "#059669"),
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

        self.assertIn('viewBox="0 0 518.40 368.00"', svg_text)
        self.assertIn("DDBJ official validation (final MSS files)", svg_text)
        self.assertIn("Parser v6.80", svg_text)
        self.assertIn("transChecker v2.26", svg_text)
        self.assertEqual(svg_text.count(">PASS</text>"), 2)
        self.assertIn("Translated CDS: AA 2 / nucleotide 2", svg_text)
        self.assertRegex(
            pdf_text,
            r"/MediaBox \[\s*0\s+0\s+518\.4(?:0)?\s+368(?:\.0+)?\s*\]",
        )
        self.assertEqual(pdf_text.count("/FontFile2"), 2)

    def test_large_no_adjustment_label_is_offset_from_the_center(self) -> None:
        node = _LaidOutNode(
            node=SankeyNode(
                "padding_unchanged",
                "No adjustment",
                4,
                100,
                "#22c55e",
            ),
            x=100.0,
            y=80.0,
            width=9.0,
            height=100.0,
        )

        self.assertEqual(_sankey_label_y(node, total_stages=6), 115.0)

    def test_event_counts_subtitle_fits_the_chart_width(self) -> None:
        self.assertLessEqual(
            pdf_helvetica_text_width(
                EVENT_COUNTS_SUBTITLE,
                size=8,
            ),
            7.2 * 72.0 - 32.0,
        )

    def test_pdf_helvetica_text_width_uses_real_character_widths(self) -> None:
        self.assertAlmostEqual(
            pdf_helvetica_text_width("Input", size=8, bold=True),
            23.91796875,
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
            self.assertLessEqual(
                link.source_y + link.height, link.source.y + link.source.height + 1e-6
            )
            self.assertGreaterEqual(link.target_y, link.target.y - 1e-6)
            self.assertLessEqual(
                link.target_y + link.height, link.target.y + link.target.height + 1e-6
            )

    def test_sankey_orders_annotation_nodes_by_priority(self) -> None:
        nodes = [
            SankeyNode("annotation_none", "Unannotated", 0, 1, "#000000"),
            SankeyNode("annotation_existing", "Existing", 0, 1, "#000000"),
            SankeyNode("annotation_cdd", "CDD", 0, 1, "#000000"),
            SankeyNode("annotation_pfam", "Pfam", 0, 1, "#000000"),
            SankeyNode("annotation_uniref90", "UniRef90", 0, 1, "#000000"),
            SankeyNode("annotation_close_reference", "Close reference", 0, 1, "#000000"),
            SankeyNode("annotation_swissprot", "Swiss-Prot", 0, 1, "#000000"),
        ]

        laid_out_nodes, _links, _meta = _sankey_layout(
            ["Functional annotation"],
            nodes,
            [],
        )

        annotation_order = sorted(
            (node for node in nodes if node.stage == 0),
            key=lambda node: laid_out_nodes[node.id].y,
        )
        self.assertEqual(
            [node.id for node in annotation_order],
            [
                "annotation_swissprot",
                "annotation_close_reference",
                "annotation_uniref90",
                "annotation_pfam",
                "annotation_cdd",
                "annotation_existing",
                "annotation_none",
            ],
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
            pdf_text = pdf_path.read_bytes().decode("latin-1")
        self.assertIn('width="7.2in"', svg_text)
        self.assertIn('viewBox="0 0 518.40 408.00"', svg_text)
        self.assertIn("Name consistency (genes)", svg_text)
        self.assertIn("Close family peer threshold", svg_text)
        self.assertIn("&gt;=70% identity / &gt;=80% mutual coverage", svg_text)
        self.assertIn("Consistent 2 (50.0%)", svg_text)
        self.assertIn("Auto-resolved family variation 1 (25.0%)", svg_text)
        self.assertIn(">Adjusted</tspan>", svg_text)
        self.assertIn(">No adjustment</tspan>", svg_text)
        self.assertRegex(
            pdf_text,
            r"/MediaBox \[\s*0\s+0\s+518\.4(?:0)?\s+408(?:\.0+)?\s*\]",
        )

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
            combined_pdf = pdf_path.read_bytes().decode("latin-1")
        self.assertIn('viewBox="0 0 518.40 432.00"', combined_svg)
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
        self.assertIn('fill-opacity="0.72"', combined_svg)
        self.assertEqual(combined_svg.count('stroke-dasharray="3 2"'), 3)
        self.assertEqual(combined_svg.count('height="132.00" rx="4"'), 3)
        self.assertRegex(
            combined_pdf,
            r"/MediaBox \[\s*0\s+0\s+518\.4(?:0)?\s+432(?:\.0+)?\s*\]",
        )

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
            duplicate_map_path = (
                log_dir / "06.drop-duplicate-coordinate-gene.duplicate-map.tsv"
            )
            _write_coordinate_duplicate_fixture(
                intermediate_dir / "05.gff.coordinates-trimmed.gff",
                duplicate_map_path,
                removed_count=5,
            )
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
                    "duplicate_map_path": str(duplicate_map_path),
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
                    "changed_gene_ids_path": str(
                        log_dir / "07.select-one-mrna.changed-gene-ids.txt"
                    ),
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
                    "updated_gene_ids_path": str(
                        log_dir / "09.update-gff-to-inframe.changed-gene-ids.txt"
                    ),
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
                    "converted_gene_ids_path": str(
                        log_dir / "16.mss-cds-to-misc.changed-gene-ids.txt"
                    ),
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
            plots_dir = output_root / "plots"
            plots_dir.mkdir()
            legacy_overlap_paths = tuple(
                plots_dir / f"pipeline-gene-overlap.{suffix}"
                for suffix in ("tsv", "svg", "pdf")
            )
            for legacy_path in legacy_overlap_paths:
                legacy_path.write_text("legacy\n", encoding="utf-8")

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
            self.assertTrue(artifacts.coordinate_duplicates_tsv.exists())
            self.assertTrue(artifacts.coordinate_duplicates_svg.exists())
            self.assertTrue(artifacts.coordinate_duplicates_pdf.exists())
            self.assertTrue(all(not path.exists() for path in legacy_overlap_paths))
            gene_flow_tsv = artifacts.gene_flow_tsv.read_text(encoding="utf-8")
            self.assertNotIn("after_transcript", gene_flow_tsv)
            self.assertNotIn("after_inframe", gene_flow_tsv)
            self.assertNotIn("after_padding", gene_flow_tsv)
            self.assertIn("transcript_changed\tinframe_updated", gene_flow_tsv)
            svg_texts = [
                artifacts.gene_flow_svg.read_text(encoding="utf-8"),
                artifacts.event_counts_svg.read_text(encoding="utf-8"),
            ]
            for svg_text in svg_texts:
                self.assertIn("font-size:8pt", svg_text)
                self.assertNotRegex(svg_text, r"font-size:\d+px")
            self.assertIn(
                "Stage-wise pipeline gene flow",
                svg_texts[0],
            )
            self.assertIn(">mRNA</tspan>", svg_texts[0])
            self.assertIn(">selection</tspan>", svg_texts[0])
            self.assertNotIn(">One mRNA</tspan>", svg_texts[0])
            self.assertIn(">Frame</tspan>", svg_texts[0])
            self.assertIn(">correction</tspan>", svg_texts[0])
            self.assertIn(">Input</tspan>", svg_texts[0])
            self.assertIn(">GFF</tspan>", svg_texts[0])
            self.assertIn(">Coordinate</tspan>", svg_texts[0])
            self.assertIn(">duplicate</tspan>", svg_texts[0])
            self.assertIn(">removal</tspan>", svg_texts[0])
            self.assertIn(">CDS</tspan>", svg_texts[0])
            self.assertIn(">boundary</tspan>", svg_texts[0])
            self.assertIn(">adjustment</tspan>", svg_texts[0])
            self.assertIn(">Output</tspan>", svg_texts[0])
            self.assertIn(">ann.txt</tspan>", svg_texts[0])
            self.assertEqual(svg_texts[0].count('class="stage"'), 6)
            self.assertIn(">Already one</tspan>", svg_texts[0])
            self.assertIn(">mRNA per gene</tspan>", svg_texts[0])
            self.assertIn(">Reduced to</tspan>", svg_texts[0])
            self.assertNotIn(">Selected<", svg_texts[0])
            self.assertNotIn(">Framed<", svg_texts[0])
            self.assertNotIn(">Pruned<", svg_texts[0])
            self.assertIn('width="7.2in"', svg_texts[0])
            self.assertIn('viewBox="0 0 518.40 420.00"', svg_texts[0])
            self.assertIn("BUSCO results (CDS; embryophyta_odb12; n=100)", svg_texts[0])
            self.assertIn("Single-copy 70.0%", svg_texts[0])
            self.assertIn(">Input CDS (CDS input n=95)</text>", svg_texts[0])
            self.assertIn(">Boundary-adjusted CDS (CDS input n=95)</text>", svg_texts[0])
            gene_flow_pdf = artifacts.gene_flow_pdf.read_bytes().decode("latin-1")
            self.assertRegex(
                gene_flow_pdf,
                r"/MediaBox \[\s*0\s+0\s+518\.4(?:0)?\s+420(?:\.0+)?\s*\]",
            )
            self.assertIn("(mRNA)", gene_flow_pdf)
            self.assertIn("(selection)", gene_flow_pdf)
            self.assertIn("(Coordinate)", gene_flow_pdf)
            self.assertIn("(duplicate)", gene_flow_pdf)
            self.assertIn("(CDS)", gene_flow_pdf)
            self.assertIn("(boundary)", gene_flow_pdf)
            self.assertIn("(Output)", gene_flow_pdf)
            self.assertIn("(ann.txt)", gene_flow_pdf)
            self.assertIn(
                r"(Boundary-adjusted CDS \(CDS input n=95\))",
                gene_flow_pdf,
            )
            self.assertIn("(Already one)", gene_flow_pdf)
            self.assertIn("(mRNA per gene)", gene_flow_pdf)
            self.assertIn("(Reduced to)", gene_flow_pdf)
            self.assertIn(
                "Pipeline event counts",
                svg_texts[1],
            )
            self.assertIn("CDS boundary-adjusted genes", svg_texts[1])
            self.assertIn('width="7.2in"', svg_texts[1])
            self.assertIn('viewBox="0 0 518.40 298.00"', svg_texts[1])
            event_counts_pdf = artifacts.event_counts_pdf.read_bytes().decode("latin-1")
            self.assertRegex(
                event_counts_pdf,
                r"/MediaBox \[\s*0\s+0\s+518\.4(?:0)?\s+298(?:\.0+)?\s*\]",
            )
            duplicate_svg = artifacts.coordinate_duplicates_svg.read_text(
                encoding="utf-8"
            )
            self.assertIn("Coordinate duplicate gene models", duplicate_svg)
            self.assertIn('data-full-label="KEPT k1 / k1.t1"', duplicate_svg)
            self.assertIn('data-full-label="REMOVED d5 / d5.t1"', duplicate_svg)
            payload = json.loads(artifacts.summary_json.read_text(encoding="utf-8"))
            self.assertEqual(payload["metrics"]["converted_to_misc_genes"], 1)
            self.assertEqual(payload["metrics"]["final_cds_genes"], 94)
            self.assertNotIn("overlap", payload)
            self.assertIn("plots", manifest)
            self.assertIn("pipeline", manifest["plots"])
            self.assertEqual(
                manifest["plots"]["pipeline"]["gene_flow_pdf"],
                str(artifacts.gene_flow_pdf),
            )
            self.assertNotIn("overlap_pdf", manifest["plots"]["pipeline"])
            duplicate_manifest = manifest["plots"]["pipeline"]["coordinate_duplicates"]
            self.assertEqual(duplicate_manifest["limit"], 50)
            self.assertEqual(duplicate_manifest["total_removed_genes"], 5)
            self.assertEqual(duplicate_manifest["shown_removed_genes"], 5)
            self.assertFalse(duplicate_manifest["truncated"])
            self.assertEqual(len(summary_lines), 1)
            self.assertIn("dedup_removed=5", summary_lines[0])
            self.assertIn("misc_feature_genes=1", summary_lines[0])


if __name__ == "__main__":
    unittest.main()
