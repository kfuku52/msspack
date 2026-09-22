import contextlib
import io
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from msspack import __version__
from msspack.busco import BuscoArtifacts
from msspack.cli import main
from msspack.databases import DatabaseStatus
from msspack.doctor import Check
from msspack.workflow import RunArtifacts


class CliTests(unittest.TestCase):
    def test_main_version_prints_package_version(self) -> None:
        stdout = io.StringIO()

        with self.assertRaises(SystemExit) as raised:
            with contextlib.redirect_stdout(stdout):
                main(["--version"])

        self.assertEqual(raised.exception.code, 0)
        self.assertEqual(stdout.getvalue().strip(), f"msspack {__version__}")

    def test_main_init_writes_config(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            output = Path(tmp_dir) / "msspack.toml"
            stdout = io.StringIO()

            with contextlib.redirect_stdout(stdout):
                exit_code = main(["init", str(output)])

            self.assertEqual(exit_code, 0)
            self.assertEqual(stdout.getvalue().strip(), str(output.resolve()))
            self.assertEqual(
                output.read_text(encoding="utf-8"),
                (
                    "# Replace the placeholder paths and submission metadata for your project.\n"
                    + (Path(__file__).parents[1] / "examples/msspack.example.toml").read_text()
                ),
            )

            expected = output.read_text(encoding="utf-8")
            output.write_text("keep me\n", encoding="utf-8")
            with contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(main(["init", str(output)]), 1)
            self.assertEqual(output.read_text(encoding="utf-8"), "keep me\n")
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(["init", str(output), "--force"]), 0)
            self.assertEqual(output.read_text(encoding="utf-8"), expected)

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
        with (
            patch("msspack.cli.run_doctor", return_value=checks),
            contextlib.redirect_stdout(stdout),
        ):
            exit_code = main(["doctor"])

        self.assertEqual(exit_code, 1)
        self.assertIn("[MISSING]", stdout.getvalue())

    def test_main_pack_dispatches_to_pipeline(self) -> None:
        outputs = SimpleNamespace(
            ann_path=Path("/tmp/build/final/sample.ann.txt"),
            fasta_path=Path("/tmp/build/final/sample.fasta"),
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
        outputs = SimpleNamespace(
            ann_path=Path("/tmp/build/final/sample.ann.txt"),
            fasta_path=Path("/tmp/build/final/sample.fasta"),
        )
        plots = SimpleNamespace(
            gene_flow_svg=Path("/tmp/build/plots/flow.svg"),
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
        fixture = Path(__file__).parent / "fixtures" / "minimal_pack" / "config.toml"
        stdout = io.StringIO()
        with tempfile.TemporaryDirectory() as tmp_dir, contextlib.redirect_stdout(stdout):
            exit_code = main(
                [
                    "db",
                    "status",
                    "--config",
                    str(fixture),
                    "--db-dir",
                    tmp_dir,
                ]
            )
            self.assertEqual(exit_code, 0)
            self.assertIn(f"Database root: {tmp_dir}", stdout.getvalue())
            self.assertIn("Database mode: shared", stdout.getvalue())

    def test_main_busco_dispatches_to_busco_module(self) -> None:
        artifacts = BuscoArtifacts(
            root=Path("/tmp/build/Demo/busco"),
            cds=SimpleNamespace(
                comparison_tsv=Path("/tmp/build/Demo/busco/cds/comparison.tsv"),
                comparison_svg=Path("/tmp/build/Demo/busco/cds/comparison.svg"),
                comparison_pdf=Path("/tmp/build/Demo/busco/cds/comparison.pdf"),
            ),
            genome=SimpleNamespace(
                comparison_tsv=Path("/tmp/build/Demo/busco/genome/comparison.tsv"),
                comparison_svg=Path("/tmp/build/Demo/busco/genome/comparison.svg"),
                comparison_pdf=Path("/tmp/build/Demo/busco/genome/comparison.pdf"),
            ),
        )
        stdout = io.StringIO()

        with (
            patch("msspack.cli.run_busco_comparison", return_value=artifacts) as mocked,
            patch(
                "msspack.cli.summarize_busco_artifacts",
                return_value=[
                    "cds\tlineage=eudicotyledons_odb12\tcomplete=95.3->95.2",
                    "genome\tlineage=eudicotyledons_odb12\tcomplete=98.0->98.0",
                ],
            ),
            patch("msspack.cli.cleanup_busco_cache", return_value=Path("/tmp/cache/busco-work")),
        ):
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
            cds=SimpleNamespace(
                comparison_tsv=Path("/tmp/build/Demo/busco/cds/comparison.tsv"),
                comparison_svg=Path("/tmp/build/Demo/busco/cds/comparison.svg"),
                comparison_pdf=Path("/tmp/build/Demo/busco/cds/comparison.pdf"),
            ),
            genome=None,
        )
        stdout = io.StringIO()

        with (
            patch("msspack.cli.run_busco_comparison", return_value=artifacts),
            patch(
                "msspack.cli.summarize_busco_artifacts",
                return_value=["cds\tlineage=eudicotyledons_odb12\tcomplete=95.3->95.2"],
            ),
            patch("msspack.cli.cleanup_busco_cache", return_value=None),
        ):
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
