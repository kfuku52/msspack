import os
import tempfile
import unittest
from pathlib import Path

from msspack.busco import run_busco_comparison
from msspack.ddbj_tools import install_component, run_parser, run_transchecker


def _validation_fixture(base: Path) -> tuple[Path, Path]:
    fixture = Path(__file__).parent / "fixtures" / "minimal_pack"
    annotation = base / "official.ann.txt"
    annotation.write_text(
        (fixture / "expected.ann.txt")
        .read_text(encoding="utf-8")
        .replace("source\t1..9", "source\t1..120")
        .replace("\tCDS\t1..9", "\tCDS\t1..120"),
        encoding="utf-8",
    )
    sequence = "ATG" + "AAA" * 38 + "TAA"
    fasta = base / "official.fasta"
    fasta.write_text(f">ctg1\n{sequence}\n//\n", encoding="utf-8")
    return annotation, fasta


@unittest.skipUnless(
    os.environ.get("MSSPACK_RUN_DDBJ_EXTERNAL") == "1",
    "set MSSPACK_RUN_DDBJ_EXTERNAL=1 to run official DDBJ tools",
)
class ExternalDdbjTests(unittest.TestCase):
    def test_official_parser_and_transchecker_accept_fixture(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            annotation, fasta = _validation_fixture(base)
            parser = install_component("parser", cache_dir=base / "cache")
            transchecker = install_component("transchecker", cache_dir=base / "cache")

            run_parser(
                parser,
                ann_path=annotation,
                fasta_path=fasta,
                heap="1G",
                java_cmd="java",
                log_path=base / "parser.log",
            )
            run_transchecker(
                transchecker,
                ann_path=annotation,
                fasta_path=fasta,
                aa_out=base / "transchecker.aa.fasta",
                nuc_out=base / "transchecker.nuc.fasta",
                heap="1G",
                java_cmd="java",
                log_path=base / "transchecker.log",
            )

            self.assertIn("jParser (Ver. 6.80) finished", (base / "parser.log").read_text())
            self.assertTrue((base / "transchecker.aa.fasta").exists())
            self.assertTrue((base / "transchecker.nuc.fasta").exists())


@unittest.skipUnless(
    os.environ.get("MSSPACK_RUN_BUSCO_EXTERNAL") == "1",
    "set MSSPACK_RUN_BUSCO_EXTERNAL=1 to run BUSCO",
)
class ExternalBuscoTests(unittest.TestCase):
    def test_busco_comparison_completes_with_real_binary(self) -> None:
        fixture = Path(__file__).parent / "fixtures" / "minimal_pack"
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            sequence = "ATG" + "AAA" * 38 + "TAA"
            (base / "input.fa").write_text(f">ctg1\n{sequence}\n", encoding="utf-8")
            (base / "input.gff3").write_text(
                "##gff-version 3\n"
                "ctg1\tsrc\tgene\t1\t120\t.\t+\t.\tID=gene1\n"
                "ctg1\tsrc\tmRNA\t1\t120\t.\t+\t.\tID=tx1;Parent=gene1;product=test%20protein\n"
                "ctg1\tsrc\texon\t1\t120\t.\t+\t.\tID=tx1.ex1;Parent=tx1\n"
                "ctg1\tsrc\tCDS\t1\t120\t.\t+\t0\tID=tx1.cds1;Parent=tx1\n",
                encoding="utf-8",
            )
            config_text = (fixture / "config.toml").read_text(encoding="utf-8")
            download_path = os.environ.get("MSSPACK_BUSCO_DOWNLOAD_PATH", "")
            config_text += (
                "\n[busco]\n"
                'command = "busco"\n'
                "run_cds = true\n"
                "run_genome = false\n"
                'cds_mode = "transcriptome"\n'
                f'lineage_dataset = "{os.environ.get("MSSPACK_BUSCO_LINEAGE", "bacteria_odb12")}"\n'
                "auto_lineage = false\n"
                "threads = 2\n"
                f'download_path = "{download_path}"\n'
                f'offline = {"true" if os.environ.get("MSSPACK_BUSCO_OFFLINE") == "1" else "false"}\n'
            )
            config_path = base / "config.toml"
            config_path.write_text(config_text, encoding="utf-8")

            artifacts = run_busco_comparison(config_path)

            self.assertIsNotNone(artifacts.cds)
            assert artifacts.cds is not None
            self.assertTrue(artifacts.cds.comparison_json.exists())
            self.assertTrue(artifacts.cds.comparison_svg.exists())


if __name__ == "__main__":
    unittest.main()
