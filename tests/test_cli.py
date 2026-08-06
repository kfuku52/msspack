import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from msspack import cli
from msspack.busco import BuscoArtifacts, BuscoComparisonArtifacts
from msspack.cli import main
from msspack.databases import DatabaseStatus
from msspack.doctor import Check
from msspack.pipeline import PipelineOutputs
from msspack.pipeline_plots import PipelinePlotArtifacts
from msspack.report import ReportArtifacts
from msspack.workflow import RunArtifacts


class CliTests(unittest.TestCase):
    def test_main_version_prints_package_version(self) -> None:
        stdout = io.StringIO()

        with self.assertRaises(SystemExit) as raised:
            with contextlib.redirect_stdout(stdout):
                main(["--version"])

        self.assertEqual(raised.exception.code, 0)
        self.assertEqual(stdout.getvalue().strip(), "msspack 0.7.4")

    def test_main_init_writes_config(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            output = Path(tmp_dir) / "msspack.toml"
            stdout = io.StringIO()

            with contextlib.redirect_stdout(stdout):
                exit_code = main(["init", str(output)])

            self.assertEqual(exit_code, 0)
            self.assertTrue(output.exists())
            self.assertEqual(stdout.getvalue().strip(), str(output.resolve()))
            self.assertEqual(output.read_text(encoding="utf-8"), cli._example_config_text())

    def test_main_init_refuses_to_overwrite_without_force(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            output = Path(tmp_dir) / "msspack.toml"
            output.write_text("keep me\n", encoding="utf-8")

            with contextlib.redirect_stderr(io.StringIO()):
                exit_code = main(["init", str(output)])

            self.assertEqual(exit_code, 1)
            self.assertEqual(output.read_text(encoding="utf-8"), "keep me\n")

    def test_main_init_force_overwrites_existing_config(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            output = Path(tmp_dir) / "msspack.toml"
            output.write_text("replace me\n", encoding="utf-8")

            with contextlib.redirect_stdout(io.StringIO()):
                exit_code = main(["init", str(output), "--force"])

            self.assertEqual(exit_code, 0)
            self.assertEqual(output.read_text(encoding="utf-8"), cli._example_config_text())

    def test_main_demo_writes_bundled_dataset(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            output = Path(tmp_dir) / "demo"
            stdout = io.StringIO()

            with contextlib.redirect_stdout(stdout):
                exit_code = main(["demo", "--output", str(output)])

            self.assertEqual(exit_code, 0)
            self.assertTrue((output / "genome.fa").is_file())
            self.assertTrue((output / "annotation.gff3").is_file())
            self.assertEqual(
                stdout.getvalue().splitlines(),
                [
                    str(output.resolve()),
                    str(output.resolve() / "config.toml"),
                    str(output.resolve() / "config.functional.toml"),
                ],
            )

    def test_main_doctor_returns_nonzero_for_required_failure(self) -> None:
        stdout = io.StringIO()
        checks = [Check("required", False, "missing", required=True)]
        with patch("msspack.cli.run_doctor", return_value=checks), patch(
            "msspack.cli.render_doctor_report", return_value="[MISSING] required: missing"
        ), contextlib.redirect_stdout(stdout):
            exit_code = main(["doctor"])

        self.assertEqual(exit_code, 1)
        self.assertIn("[MISSING]", stdout.getvalue())

    def test_main_pack_dispatches_to_pipeline(self) -> None:
        outputs = PipelineOutputs(
            root=Path("/tmp/build"),
            intermediate=Path("/tmp/build/intermediate"),
            logs=Path("/tmp/build/logs"),
            final=Path("/tmp/build/final"),
            ann_path=Path("/tmp/build/final/sample.ann.txt"),
            fasta_path=Path("/tmp/build/final/sample.fasta"),
            manifest_path=Path("/tmp/build/build-manifest.json"),
        )
        stdout = io.StringIO()

        with patch("msspack.cli.run_pipeline", return_value=outputs) as mocked:
            with contextlib.redirect_stdout(stdout):
                exit_code = main(["pack", "--config", "/tmp/demo.toml", "--no-validate"])

        self.assertEqual(exit_code, 0)
        mocked.assert_called_once_with("/tmp/demo.toml", validate=False)
        self.assertEqual(
            stdout.getvalue().splitlines(),
            [str(outputs.ann_path), str(outputs.fasta_path)],
        )

    def test_main_run_dispatches_full_workflow(self) -> None:
        outputs = PipelineOutputs(
            root=Path("/tmp/build"),
            intermediate=Path("/tmp/build/intermediate"),
            logs=Path("/tmp/build/logs"),
            final=Path("/tmp/build/final"),
            ann_path=Path("/tmp/build/final/sample.ann.txt"),
            fasta_path=Path("/tmp/build/final/sample.fasta"),
            manifest_path=Path("/tmp/build/build-manifest.json"),
        )
        plots = PipelinePlotArtifacts(
            root=Path("/tmp/build/plots"),
            summary_json=Path("/tmp/build/plots/summary.json"),
            summary_tsv=Path("/tmp/build/plots/summary.tsv"),
            gene_flow_tsv=Path("/tmp/build/plots/flow.tsv"),
            gene_flow_svg=Path("/tmp/build/plots/flow.svg"),
            gene_flow_pdf=Path("/tmp/build/plots/flow.pdf"),
            event_counts_tsv=Path("/tmp/build/plots/events.tsv"),
            event_counts_svg=Path("/tmp/build/plots/events.svg"),
            event_counts_pdf=Path("/tmp/build/plots/events.pdf"),
            coordinate_duplicates_tsv=Path("/tmp/build/plots/duplicates.tsv"),
            coordinate_duplicates_svg=Path("/tmp/build/plots/duplicates.svg"),
            coordinate_duplicates_pdf=Path("/tmp/build/plots/duplicates.pdf"),
            name_consistency_tsv=Path("/tmp/build/plots/names.tsv"),
            name_consistency_svg=Path("/tmp/build/plots/names.svg"),
            name_consistency_pdf=Path("/tmp/build/plots/names.pdf"),
            source_consistency_tsv=Path("/tmp/build/plots/sources.tsv"),
            source_consistency_svg=Path("/tmp/build/plots/sources.svg"),
            source_consistency_pdf=Path("/tmp/build/plots/sources.pdf"),
        )
        status = DatabaseStatus(
            root=Path("/shared/msspack-db"),
            mode="shared",
            resources=(),
        )
        artifacts = RunArtifacts(
            pipeline=outputs,
            busco=None,
            plots=plots,
            report=None,
            database_status=status,
            duration_seconds=12.3456,
        )
        stdout = io.StringIO()

        with patch("msspack.cli.run_all", return_value=artifacts) as mocked:
            with contextlib.redirect_stdout(stdout):
                exit_code = main(
                    [
                        "run",
                        "--config",
                        "/tmp/demo.toml",
                        "--db-dir",
                        "/shared/msspack-db",
                        "--force-compute",
                        "--no-busco",
                        "--no-validate",
                        "--no-report",
                    ]
                )

        self.assertEqual(exit_code, 0)
        mocked.assert_called_once_with(
            "/tmp/demo.toml",
            database_dir="/shared/msspack-db",
            force_compute=True,
            run_busco=False,
            validate=False,
            write_report=False,
        )
        self.assertEqual(
            stdout.getvalue().splitlines(),
            [
                "Database root: /shared/msspack-db",
                "Database mode: shared",
                str(outputs.ann_path),
                str(outputs.fasta_path),
                str(plots.gene_flow_svg),
                "duration_seconds\t12.346",
            ],
        )

    def test_main_db_status_uses_database_override(self) -> None:
        status = DatabaseStatus(
            root=Path("/shared/msspack-db"),
            mode="shared",
            resources=(),
        )
        stdout = io.StringIO()

        with patch("msspack.cli.load_config", return_value=object()), patch(
            "msspack.cli.collect_database_status",
            return_value=status,
        ) as mocked:
            with contextlib.redirect_stdout(stdout):
                exit_code = main(
                    [
                        "db",
                        "status",
                        "--config",
                        "/tmp/demo.toml",
                        "--db-dir",
                        "/shared/msspack-db",
                    ]
                )

        self.assertEqual(exit_code, 0)
        mocked.assert_called_once()
        self.assertEqual(
            stdout.getvalue().splitlines(),
            [
                "Database root: /shared/msspack-db",
                "Database mode: shared",
            ],
        )

    def test_main_busco_dispatches_to_busco_module(self) -> None:
        artifacts = BuscoArtifacts(
            root=Path("/tmp/build/Demo/busco"),
            cds=BuscoComparisonArtifacts(
                root=Path("/tmp/build/Demo/busco/cds"),
                logs_dir=Path("/tmp/build/Demo/busco/cds/logs"),
                raw_root=Path("/tmp/build/Demo/busco/cds/raw"),
                input_fasta=Path("/tmp/build/Demo/busco/cds/input.cds.fasta"),
                processed_fasta=Path("/tmp/build/Demo/busco/cds/processed.cds.fasta"),
                input_summary_json=Path("/tmp/build/Demo/busco/cds/input.summary.json"),
                processed_summary_json=Path("/tmp/build/Demo/busco/cds/processed.summary.json"),
                comparison_json=Path("/tmp/build/Demo/busco/cds/comparison.json"),
                comparison_tsv=Path("/tmp/build/Demo/busco/cds/comparison.tsv"),
                comparison_svg=Path("/tmp/build/Demo/busco/cds/comparison.svg"),
                comparison_pdf=Path("/tmp/build/Demo/busco/cds/comparison.pdf"),
            ),
            genome=BuscoComparisonArtifacts(
                root=Path("/tmp/build/Demo/busco/genome"),
                logs_dir=Path("/tmp/build/Demo/busco/genome/logs"),
                raw_root=Path("/tmp/build/Demo/busco/genome/raw"),
                input_fasta=Path("/tmp/build/Demo/busco/genome/input.genome.fasta"),
                processed_fasta=Path("/tmp/build/Demo/busco/genome/processed.genome.fasta"),
                input_summary_json=Path("/tmp/build/Demo/busco/genome/input.summary.json"),
                processed_summary_json=Path(
                    "/tmp/build/Demo/busco/genome/processed.summary.json"
                ),
                comparison_json=Path("/tmp/build/Demo/busco/genome/comparison.json"),
                comparison_tsv=Path("/tmp/build/Demo/busco/genome/comparison.tsv"),
                comparison_svg=Path("/tmp/build/Demo/busco/genome/comparison.svg"),
                comparison_pdf=Path("/tmp/build/Demo/busco/genome/comparison.pdf"),
            ),
        )
        stdout = io.StringIO()

        with patch("msspack.cli.run_busco_comparison", return_value=artifacts) as mocked, patch(
            "msspack.cli.summarize_busco_artifacts",
            return_value=[
                "cds\tlineage=eudicotyledons_odb12\tcomplete=95.3->95.2",
                "genome\tlineage=eudicotyledons_odb12\tcomplete=98.0->98.0",
            ],
        ), patch("msspack.cli.cleanup_busco_cache", return_value=Path("/tmp/cache/busco-work")):
            with contextlib.redirect_stdout(stdout):
                exit_code = main(
                    [
                        "busco",
                        "--config",
                        "/tmp/demo.toml",
                        "--lineage-dataset",
                        "embryophyta_odb12",
                        "--threads",
                        "4",
                        "--force",
                        "--clean-cache",
                        "--genome",
                        "--auto-lineage-scope",
                        "euk",
                    ]
                )

        self.assertEqual(exit_code, 0)
        mocked.assert_called_once_with(
            "/tmp/demo.toml",
            lineage_dataset="embryophyta_odb12",
            threads=4,
            force=True,
            clean_cache=False,
            auto_lineage=None,
            auto_lineage_scope="euk",
            run_genome=True,
            run_cds=None,
        )
        self.assertEqual(
            stdout.getvalue().splitlines(),
            [
                "cleaned_staging_cache\t/tmp/cache/busco-work",
                str(artifacts.root),
                str(artifacts.cds.comparison_tsv),
                str(artifacts.cds.comparison_svg),
                str(artifacts.cds.comparison_pdf),
                str(artifacts.genome.comparison_tsv),
                str(artifacts.genome.comparison_svg),
                str(artifacts.genome.comparison_pdf),
                "cds\tlineage=eudicotyledons_odb12\tcomplete=95.3->95.2",
                "genome\tlineage=eudicotyledons_odb12\tcomplete=98.0->98.0",
            ],
        )

    def test_main_busco_prints_cds_only_when_genome_is_disabled(self) -> None:
        artifacts = BuscoArtifacts(
            root=Path("/tmp/build/Demo/busco"),
            cds=BuscoComparisonArtifacts(
                root=Path("/tmp/build/Demo/busco/cds"),
                logs_dir=Path("/tmp/build/Demo/busco/cds/logs"),
                raw_root=Path("/tmp/build/Demo/busco/cds/raw"),
                input_fasta=Path("/tmp/build/Demo/busco/cds/input.cds.fasta"),
                processed_fasta=Path("/tmp/build/Demo/busco/cds/processed.cds.fasta"),
                input_summary_json=Path("/tmp/build/Demo/busco/cds/input.summary.json"),
                processed_summary_json=Path("/tmp/build/Demo/busco/cds/processed.summary.json"),
                comparison_json=Path("/tmp/build/Demo/busco/cds/comparison.json"),
                comparison_tsv=Path("/tmp/build/Demo/busco/cds/comparison.tsv"),
                comparison_svg=Path("/tmp/build/Demo/busco/cds/comparison.svg"),
                comparison_pdf=Path("/tmp/build/Demo/busco/cds/comparison.pdf"),
            ),
            genome=None,
        )
        stdout = io.StringIO()

        with patch("msspack.cli.run_busco_comparison", return_value=artifacts), patch(
            "msspack.cli.summarize_busco_artifacts",
            return_value=["cds\tlineage=eudicotyledons_odb12\tcomplete=95.3->95.2"],
        ), patch("msspack.cli.cleanup_busco_cache", return_value=None):
            with contextlib.redirect_stdout(stdout):
                exit_code = main(["busco", "--config", "/tmp/demo.toml"])

        self.assertEqual(exit_code, 0)
        self.assertEqual(
            stdout.getvalue().splitlines(),
            [
                str(artifacts.root),
                str(artifacts.cds.comparison_tsv),
                str(artifacts.cds.comparison_svg),
                str(artifacts.cds.comparison_pdf),
                "cds\tlineage=eudicotyledons_odb12\tcomplete=95.3->95.2",
            ],
        )

    def test_main_plot_dispatches_to_pipeline_plots(self) -> None:
        artifacts = PipelinePlotArtifacts(
            root=Path("/tmp/build/Demo/plots"),
            summary_json=Path("/tmp/build/Demo/plots/pipeline-flow-summary.json"),
            summary_tsv=Path("/tmp/build/Demo/plots/pipeline-flow-summary.tsv"),
            gene_flow_tsv=Path("/tmp/build/Demo/plots/pipeline-gene-flow.tsv"),
            gene_flow_svg=Path("/tmp/build/Demo/plots/pipeline-gene-flow.sankey.svg"),
            gene_flow_pdf=Path("/tmp/build/Demo/plots/pipeline-gene-flow.sankey.pdf"),
            event_counts_tsv=Path("/tmp/build/Demo/plots/pipeline-event-counts.tsv"),
            event_counts_svg=Path("/tmp/build/Demo/plots/pipeline-event-counts.svg"),
            event_counts_pdf=Path("/tmp/build/Demo/plots/pipeline-event-counts.pdf"),
            coordinate_duplicates_tsv=Path(
                "/tmp/build/Demo/plots/coordinate-duplicate-gene-models.tsv"
            ),
            coordinate_duplicates_svg=Path(
                "/tmp/build/Demo/plots/coordinate-duplicate-gene-models.svg"
            ),
            coordinate_duplicates_pdf=Path(
                "/tmp/build/Demo/plots/coordinate-duplicate-gene-models.pdf"
            ),
            name_consistency_tsv=Path(
                "/tmp/build/Demo/plots/functional-annotation-name-consistency.tsv"
            ),
            name_consistency_svg=Path(
                "/tmp/build/Demo/plots/functional-annotation-name-consistency.svg"
            ),
            name_consistency_pdf=Path(
                "/tmp/build/Demo/plots/functional-annotation-name-consistency.pdf"
            ),
            source_consistency_tsv=Path(
                "/tmp/build/Demo/plots/functional-annotation-source-consistency.tsv"
            ),
            source_consistency_svg=Path(
                "/tmp/build/Demo/plots/functional-annotation-source-consistency.svg"
            ),
            source_consistency_pdf=Path(
                "/tmp/build/Demo/plots/functional-annotation-source-consistency.pdf"
            ),
        )
        stdout = io.StringIO()

        with patch("msspack.cli.run_pipeline_plots", return_value=artifacts) as mocked, patch(
            "msspack.cli.summarize_pipeline_plots",
            return_value=[
                "pipeline-flow\tgenes=100\tdedup_removed=5\ttranscript_changed=20\tremoved_mrnas=23"
            ],
        ):
            with contextlib.redirect_stdout(stdout):
                exit_code = main(["plot", "--config", "/tmp/demo.toml", "--force"])

        self.assertEqual(exit_code, 0)
        mocked.assert_called_once_with("/tmp/demo.toml", force=True)
        self.assertEqual(
            stdout.getvalue().splitlines(),
            [
                str(artifacts.root),
                str(artifacts.summary_tsv),
                str(artifacts.gene_flow_svg),
                str(artifacts.gene_flow_pdf),
                str(artifacts.event_counts_svg),
                str(artifacts.event_counts_pdf),
                str(artifacts.coordinate_duplicates_svg),
                str(artifacts.coordinate_duplicates_pdf),
                "pipeline-flow\tgenes=100\tdedup_removed=5\ttranscript_changed=20\tremoved_mrnas=23",
            ],
        )

    def test_main_report_dispatches_to_report_module(self) -> None:
        artifacts = ReportArtifacts(
            root=Path("/tmp/build/Demo/report"),
            index_html=Path("/tmp/build/Demo/report/index.html"),
        )
        stdout = io.StringIO()

        with patch("msspack.cli.run_html_report", return_value=artifacts) as mocked:
            with contextlib.redirect_stdout(stdout):
                exit_code = main(["report", "--config", "/tmp/demo.toml", "--force"])

        self.assertEqual(exit_code, 0)
        mocked.assert_called_once_with("/tmp/demo.toml", force=True)
        self.assertEqual(
            stdout.getvalue().splitlines(),
            [str(artifacts.root), str(artifacts.index_html)],
        )

    def test_main_validate_dispatches_to_validation_module(self) -> None:
        validation_outputs = {
            "parser_log": Path("/tmp/final/logs/parser.log"),
            "transchecker_log": Path("/tmp/final/logs/transchecker.log"),
            "aa_fasta": Path("/tmp/final/validation/transChecker.aa.fasta"),
            "nuc_fasta": Path("/tmp/final/validation/transChecker.nuc.fasta"),
        }
        stdout = io.StringIO()

        with patch("msspack.cli.validate_existing", return_value=validation_outputs) as mocked:
            with contextlib.redirect_stdout(stdout):
                exit_code = main(
                    [
                        "validate",
                        "--config",
                        "/tmp/demo.toml",
                        "--ann",
                        "/tmp/final/sample.ann.txt",
                        "--fasta",
                        "/tmp/final/sample.fasta",
                    ]
                )

        self.assertEqual(exit_code, 0)
        mocked.assert_called_once_with(
            config_file="/tmp/demo.toml",
            ann_path="/tmp/final/sample.ann.txt",
            fasta_path="/tmp/final/sample.fasta",
        )
        self.assertEqual(
            stdout.getvalue().splitlines(),
            [str(path) for path in validation_outputs.values()],
        )

    def test_tools_install_defaults_to_validation_components(self) -> None:
        stdout = io.StringIO()
        stderr = io.StringIO()

        with patch(
            "msspack.cli.install_component",
            side_effect=lambda component, **_kwargs: type(
                "Install", (), {"component": component, "version": "1.0", "root": Path("/tmp")}
            )(),
        ) as mocked:
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                exit_code = main(["tools", "install"])

        self.assertEqual(exit_code, 0)
        self.assertEqual(
            [call.args[0] for call in mocked.call_args_list],
            ["parser", "transchecker"],
        )
        self.assertIn("DDBJ agreement", stderr.getvalue())

    def test_main_internal_remove_trailing_ns_runs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            input_path = base / "input.fa"
            output_path = base / "output.fa"
            log_path = base / "remove.log"
            input_path.write_text(">seq1\nACGTNNN\n", encoding="utf-8")
            stdout = io.StringIO()

            with contextlib.redirect_stdout(stdout):
                exit_code = main(
                    [
                        "internal",
                        "remove-trailing-ns",
                        "--input",
                        str(input_path),
                        "--output",
                        str(output_path),
                        "--log",
                        str(log_path),
                    ]
                )

            self.assertEqual(exit_code, 0)
            self.assertEqual(output_path.read_text(encoding="utf-8"), ">seq1\nACGT\n")
            log_text = log_path.read_text(encoding="utf-8")
            self.assertIn("Summary table:", log_text)
            self.assertIn("| Metric | Value |", log_text)
            self.assertIn("| Step | remove-trailing-ns |", log_text)
            self.assertIn("| Input total | 1 |", log_text)
            self.assertIn("| Changed total | 1 |", log_text)
            self.assertIn("| Output total | 1 |", log_text)
            self.assertIn("Count unit: sequences", log_text)
            self.assertIn("Input total: 1", log_text)
            self.assertIn("Changed total: 1", log_text)
            self.assertIn("Output total: 1", log_text)
            self.assertIn("Started at:", log_text)
            self.assertIn("Ended at:", log_text)
            self.assertIn("Total Ns removed: 3", log_text)
            self.assertEqual(
                stdout.getvalue().splitlines(),
                [str(output_path.resolve()), str(log_path.resolve())],
            )

    def test_main_internal_select_one_mrna_runs(self) -> None:
        gff_text = """\
chr1\tsrc\tgene\t100\t400\t.\t+\t.\tID=Gene1
chr1\tsrc\tmRNA\t100\t400\t.\t+\t.\tID=Gene1-T1;Parent=Gene1
chr1\tsrc\tmRNA\t100\t400\t.\t+\t.\tID=Gene1-T2;Parent=Gene1
chr1\tsrc\texon\t100\t261\t.\t+\t.\tID=Gene1-T1.exon1;Parent=Gene1-T1
chr1\tsrc\tCDS\t120\t221\t.\t+\t0\tID=Gene1-T1.cds1;Parent=Gene1-T1
chr1\tsrc\tfive_prime_UTR\t100\t119\t.\t+\t.\tID=Gene1-T1.utr5;Parent=Gene1-T1
chr1\tsrc\tthree_prime_UTR\t222\t261\t.\t+\t.\tID=Gene1-T1.utr3;Parent=Gene1-T1
chr1\tsrc\texon\t100\t200\t.\t+\t.\tID=Gene1-T2.exon1;Parent=Gene1-T2
chr1\tsrc\texon\t320\t360\t.\t+\t.\tID=Gene1-T2.exon2;Parent=Gene1-T2
chr1\tsrc\tCDS\t120\t180\t.\t+\t0\tID=Gene1-T2.cds1;Parent=Gene1-T2
chr1\tsrc\tCDS\t320\t360\t.\t+\t0\tID=Gene1-T2.cds2;Parent=Gene1-T2
chr1\tsrc\tfive_prime_UTR\t100\t119\t.\t+\t.\tID=Gene1-T2.utr5;Parent=Gene1-T2
chr1\tsrc\tthree_prime_UTR\t181\t200\t.\t+\t.\tID=Gene1-T2.utr3a;Parent=Gene1-T2
chr1\tsrc\tthree_prime_UTR\t361\t380\t.\t+\t.\tID=Gene1-T2.utr3b;Parent=Gene1-T2
"""
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            input_path = base / "input.gff"
            output_path = base / "output.gff"
            log_path = base / "select.log"
            input_path.write_text(gff_text, encoding="utf-8")
            stdout = io.StringIO()

            with contextlib.redirect_stdout(stdout):
                exit_code = main(
                    [
                        "internal",
                        "select-one-mrna",
                        "--input",
                        str(input_path),
                        "--output",
                        str(output_path),
                        "--log",
                        str(log_path),
                    ]
                )

            self.assertEqual(exit_code, 0)
            output = output_path.read_text(encoding="utf-8")
            self.assertIn("ID=Gene1-T1;Parent=Gene1", output)
            self.assertNotIn("ID=Gene1-T2;Parent=Gene1", output)
            log_text = log_path.read_text(encoding="utf-8")
            self.assertIn("Summary table:", log_text)
            self.assertIn("| Metric | Value |", log_text)
            self.assertIn("| Step | select-one-mrna |", log_text)
            self.assertIn("| Input total | 1 |", log_text)
            self.assertIn("| Changed total | 1 |", log_text)
            self.assertIn("| Output total | 1 |", log_text)
            self.assertIn("Count unit: genes", log_text)
            self.assertIn("Input total: 1", log_text)
            self.assertIn("Changed total: 1", log_text)
            self.assertIn("Output total: 1", log_text)
            self.assertIn("Started at:", log_text)
            self.assertIn("Ended at:", log_text)
            self.assertIn("Removed mRNAs: 1", log_text)
            self.assertEqual(
                stdout.getvalue().splitlines(),
                [str(output_path.resolve()), str(log_path.resolve())],
            )
