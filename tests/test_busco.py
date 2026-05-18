import json
import tempfile
import unittest
from pathlib import Path

from msspack.busco import (
    BuscoArtifacts,
    BuscoComparisonArtifacts,
    _build_busco_command,
    _discover_short_summary,
    _publish_busco_workspace,
    _update_busco_manifest,
    _write_comparison_svg,
    _write_summary_json,
    busco_workspace_root,
    cleanup_busco_cache,
    parse_short_summary,
    summarize_busco_artifacts,
)
from msspack.config import BuscoConfig


class BuscoTests(unittest.TestCase):
    def test_parse_short_summary_reads_percentages_counts_and_dataset(self) -> None:
        summary_text = """\
# BUSCO version is: 6.0.0
# The lineage dataset is: embryophyta_odb12 (Creation date: 2024-11-15, number of genomes: 161, number of BUSCOs: 425)
# BUSCO was run in mode: genome
C:98.6%[S:97.4%,D:1.2%],F:0.5%,M:0.9%,n:425
419 Complete BUSCOs (C)
414 Complete and single-copy BUSCOs (S)
5 Complete and duplicated BUSCOs (D)
2 Fragmented BUSCOs (F)
4 Missing BUSCOs (M)
"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            summary_path = Path(tmp_dir) / "short_summary.txt"
            summary_path.write_text(summary_text, encoding="utf-8")

            summary = parse_short_summary(
                summary_path,
                label="input",
                input_fasta=Path("/tmp/input.fa"),
                raw_output_dir=Path("/tmp/raw"),
                selection_strategy="auto-lineage",
            )

        self.assertEqual(summary.lineage_dataset, "embryophyta_odb12")
        self.assertEqual(summary.mode, "genome")
        self.assertEqual(summary.busco_version, "6.0.0")
        self.assertEqual(summary.total_buscos, 425)
        self.assertEqual(summary.complete_count, 419)
        self.assertEqual(summary.single_copy_count, 414)
        self.assertAlmostEqual(summary.complete_pct, 98.6)
        self.assertEqual(summary.selection_strategy, "auto-lineage")

    def test_parse_short_summary_accepts_transcriptome_mode(self) -> None:
        summary_text = """\
# BUSCO version is: 6.0.0
# The lineage dataset is: embryophyta_odb12
# BUSCO was run in mode: transcriptome
C:95.0%[S:90.0%,D:5.0%],F:2.0%,M:3.0%,n:100
"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            summary_path = Path(tmp_dir) / "short_summary.txt"
            summary_path.write_text(summary_text, encoding="utf-8")
            summary = parse_short_summary(
                summary_path,
                label="cds",
                input_fasta=Path("/tmp/cds.fa"),
                raw_output_dir=Path("/tmp/raw"),
            )

        self.assertEqual(summary.mode, "transcriptome")

    def test_build_busco_command_supports_auto_lineage_scope(self) -> None:
        command = _build_busco_command(
            busco=BuscoConfig(
                command="busco",
                mode="genome",
                lineage_dataset="",
                auto_lineage=True,
                auto_lineage_scope="euk",
                threads=8,
                opt_out_run_stats=True,
            ),
            input_fasta=Path("/tmp/input.fa"),
            out_name="run",
            out_path=Path("/tmp/out"),
            use_auto_lineage=True,
        )

        self.assertIn("--auto-lineage-euk", command)
        self.assertIn("--opt-out-run-stats", command)
        self.assertIn("8", command)

    def test_write_comparison_svg_includes_labels_and_metrics(self) -> None:
        summary_text = """\
# BUSCO version is: 6.0.0
# The lineage dataset is: embryophyta_odb12
# BUSCO was run in mode: genome
C:98.6%[S:97.4%,D:1.2%],F:0.5%,M:0.9%,n:425
"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            short_summary = base / "short_summary.txt"
            short_summary.write_text(summary_text, encoding="utf-8")
            input_summary = parse_short_summary(
                short_summary,
                label="input",
                input_fasta=Path("/tmp/input.fa"),
                raw_output_dir=base / "input",
            )
            processed_summary = parse_short_summary(
                short_summary,
                label="processed",
                input_fasta=Path("/tmp/processed.fa"),
                raw_output_dir=base / "processed",
                selection_strategy="auto-lineage-from-input",
            )
            svg_path = base / "comparison.svg"
            _write_comparison_svg([input_summary, processed_summary], svg_path, comparison_name="genome")
            svg = svg_path.read_text(encoding="utf-8")

        self.assertIn("BUSCO comparison: genome", svg)
        self.assertIn("input", svg)
        self.assertIn("processed", svg)
        self.assertIn("embryophyta_odb12", svg)
        self.assertIn("font-size:8pt", svg)
        self.assertNotRegex(svg, r"font-size:\d+px")

    def test_summarize_busco_artifacts_reports_metric_changes(self) -> None:
        summary_text_input = """\
# BUSCO version is: 6.0.0
# The lineage dataset is: eudicotyledons_odb12
# BUSCO was run in mode: transcriptome
C:95.3%[S:69.9%,D:25.4%],F:1.7%,M:3.0%,n:2805
"""
        summary_text_processed = """\
# BUSCO version is: 6.0.0
# The lineage dataset is: eudicotyledons_odb12
# BUSCO was run in mode: transcriptome
C:95.2%[S:89.4%,D:5.8%],F:1.8%,M:3.0%,n:2805
"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            (base / "input.short_summary.txt").write_text(summary_text_input, encoding="utf-8")
            (base / "processed.short_summary.txt").write_text(summary_text_processed, encoding="utf-8")
            input_summary = parse_short_summary(
                base / "input.short_summary.txt",
                label="input",
                input_fasta=Path("/tmp/input.cds.fa"),
                raw_output_dir=base / "input",
            )
            processed_summary = parse_short_summary(
                base / "processed.short_summary.txt",
                label="processed",
                input_fasta=Path("/tmp/processed.cds.fa"),
                raw_output_dir=base / "processed",
                selection_strategy="auto-lineage-from-cds-input",
            )
            comparison = BuscoComparisonArtifacts(
                root=base / "cds",
                logs_dir=base / "cds" / "logs",
                raw_root=base / "cds" / "raw",
                input_fasta=base / "cds" / "input.cds.fasta",
                processed_fasta=base / "cds" / "processed.cds.fasta",
                input_summary_json=base / "cds" / "input.summary.json",
                processed_summary_json=base / "cds" / "processed.summary.json",
                comparison_json=base / "cds" / "comparison.json",
                comparison_tsv=base / "cds" / "comparison.tsv",
                comparison_svg=base / "cds" / "comparison.svg",
                comparison_pdf=base / "cds" / "comparison.pdf",
            )
            comparison.root.mkdir(parents=True)
            _write_summary_json(input_summary, comparison.input_summary_json)
            _write_summary_json(processed_summary, comparison.processed_summary_json)
            lines = summarize_busco_artifacts(
                BuscoArtifacts(root=base / "busco", cds=comparison, genome=None)
            )

        self.assertEqual(len(lines), 1)
        self.assertIn("cds", lines[0])
        self.assertIn("lineage=eudicotyledons_odb12", lines[0])
        self.assertIn("single_copy=69.9->89.4", lines[0])

    def test_update_busco_manifest_records_comparison_outputs(self) -> None:
        summary_text = """\
# BUSCO version is: 6.0.0
# The lineage dataset is: eudicotyledons_odb12
# BUSCO was run in mode: transcriptome
C:95.3%[S:69.9%,D:25.4%],F:1.7%,M:3.0%,n:2805
"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            summary_path = base / "summary.txt"
            summary_path.write_text(summary_text, encoding="utf-8")
            summary = parse_short_summary(
                summary_path,
                label="input",
                input_fasta=Path("/tmp/input.cds.fa"),
                raw_output_dir=base / "raw",
            )
            comparison = BuscoComparisonArtifacts(
                root=base / "cds",
                logs_dir=base / "cds" / "logs",
                raw_root=base / "cds" / "raw",
                input_fasta=base / "cds" / "input.cds.fasta",
                processed_fasta=base / "cds" / "processed.cds.fasta",
                input_summary_json=base / "cds" / "input.summary.json",
                processed_summary_json=base / "cds" / "processed.summary.json",
                comparison_json=base / "cds" / "comparison.json",
                comparison_tsv=base / "cds" / "comparison.tsv",
                comparison_svg=base / "cds" / "comparison.svg",
                comparison_pdf=base / "cds" / "comparison.pdf",
            )
            comparison.logs_dir.mkdir(parents=True)
            _write_summary_json(summary, comparison.input_summary_json)
            _write_summary_json(summary, comparison.processed_summary_json)
            manifest_path = base / "build-manifest.json"
            manifest_path.write_text("{\"status\": \"completed\"}\n", encoding="utf-8")
            _update_busco_manifest(
                manifest_path,
                busco=BuscoConfig(
                    command="/usr/bin/busco",
                    run_cds=True,
                    run_genome=False,
                    lineage_dataset="eudicotyledons_odb12",
                    auto_lineage=False,
                ),
                artifacts=BuscoArtifacts(root=base / "busco", cds=comparison, genome=None),
            )
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))

        self.assertIn("busco", payload)
        self.assertTrue(payload["busco"]["run_cds"])
        self.assertFalse(payload["busco"]["run_genome"])
        self.assertIn("cds", payload["busco"]["comparisons"])
        self.assertIn("comparison_pdf", payload["busco"]["comparisons"]["cds"])

    def test_cleanup_busco_cache_removes_workspace_root(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            workspace_root = Path(tmp_dir) / "busco-work"
            workspace_root.mkdir()
            (workspace_root / "stale-run").mkdir()

            from unittest.mock import patch

            with patch("msspack.busco.busco_workspace_root", return_value=workspace_root):
                cleaned = cleanup_busco_cache()

        self.assertEqual(cleaned, workspace_root)
        self.assertFalse(workspace_root.exists())

    def test_busco_workspace_root_uses_msspack_cache_dir(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            cache_root = Path(tmp_dir) / "cache"
            from unittest.mock import patch

            with patch("msspack.busco.default_cache_dir", return_value=cache_root):
                self.assertEqual(busco_workspace_root(), cache_root / "busco-work")

    def test_publish_busco_workspace_keeps_summary_discoverable_from_final_raw_dir(self) -> None:
        summary_text = """\
# BUSCO version is: 6.0.0
# The lineage dataset is: embryophyta_odb12
# BUSCO was run in mode: transcriptome
C:95.0%[S:90.0%,D:5.0%],F:2.0%,M:3.0%,n:100
"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            work_dir = base / "work"
            raw_dir = base / "raw"
            run_dir = work_dir / "run"
            run_dir.mkdir(parents=True)
            short_summary = run_dir / "short_summary.txt"
            short_summary.write_text(summary_text, encoding="utf-8")

            _publish_busco_workspace(work_dir, raw_dir)

            discovered = _discover_short_summary(raw_dir)

        self.assertEqual(discovered, raw_dir / "run" / "short_summary.txt")
