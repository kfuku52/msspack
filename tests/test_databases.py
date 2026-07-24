import json
import tempfile
import unittest
from pathlib import Path

from msspack.config import load_config
from msspack.databases import (
    collect_database_status,
    render_database_status,
    write_database_manifest,
)
from msspack.demo import write_demo_dataset


class DatabaseStatusTests(unittest.TestCase):
    def test_reports_project_database_root_and_writes_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            demo_root = write_demo_dataset(Path(tmp_dir) / "demo")
            config = load_config(demo_root / "config.toml")
            status = collect_database_status(config)

            self.assertEqual(status.root, demo_root / "msspack_db")
            self.assertEqual(status.mode, "project")
            self.assertIn(f"Database root: {demo_root / 'msspack_db'}", render_database_status(status))

            manifest_path = write_database_manifest(config)
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["root"], str(demo_root / "msspack_db"))
            self.assertEqual(payload["mode"], "project")
            self.assertTrue(payload["resources"])

    def test_detects_prepared_annotation_databases(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            demo_root = write_demo_dataset(Path(tmp_dir) / "demo")
            config = load_config(demo_root / "config.toml")
            root = config.database_dir
            swissprot = root / "diamond" / "swissprot" / "swissprot.dmnd"
            swissprot.parent.mkdir(parents=True)
            swissprot.write_bytes(b"database")
            pfam = root / "pfam" / "objects" / "digest.hmm"
            pfam.parent.mkdir(parents=True)
            pfam.write_text("HMMER3/f\n", encoding="utf-8")
            for suffix in (".h3f", ".h3i", ".h3m", ".h3p"):
                Path(str(pfam) + suffix).write_bytes(b"index")

            status = collect_database_status(config)
            swissprot_resource = next(
                item for item in status.resources if item.name == "Swiss-Prot"
            )
            pfam_resource = next(item for item in status.resources if item.name == "Pfam")

            self.assertFalse(swissprot_resource.enabled)
            self.assertTrue(swissprot_resource.ready)
            self.assertFalse(pfam_resource.enabled)
            self.assertTrue(pfam_resource.ready)
            rendered = render_database_status(status)
            self.assertIn("Swiss-Prot: disabled (cached)", rendered)
            self.assertIn("Pfam: disabled (cached)", rendered)

    def test_requires_complete_busco_lineage_and_finds_versioned_cdd(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            demo_root = write_demo_dataset(Path(tmp_dir) / "demo")
            config = load_config(demo_root / "config.toml")
            config.busco.lineage_dataset = "eukaryota_odb12"
            lineage = config.busco_database_dir / "lineages" / "eukaryota_odb12"
            lineage.mkdir(parents=True)
            (lineage / "ancestral").write_text("partial\n", encoding="utf-8")

            cdd_prefix = (
                config.database_dir
                / "cdd"
                / "database-versions"
                / ("a" * 64)
                / "little_endian"
                / "Cdd"
            )
            cdd_prefix.parent.mkdir(parents=True)
            for suffix in (".aux", ".freq", ".loo", ".rps"):
                Path(str(cdd_prefix) + suffix).write_text(
                    f"{suffix}\n",
                    encoding="utf-8",
                )
            cdd_data = (
                config.database_dir
                / "cdd"
                / "data-versions"
                / ("b" * 64)
            )
            cdd_data.mkdir(parents=True)
            for name in (
                "cddid.tbl",
                "cdtrack.txt",
                "family_superfamily_links",
                "cddannot.dat",
                "cddannot_generic.dat",
                "bitscore_specific.txt",
            ):
                (cdd_data / name).write_text(f"{name}\n", encoding="utf-8")
            (cdd_data / "provenance.json").write_text("{}\n", encoding="utf-8")

            status = collect_database_status(config)
            busco = next(item for item in status.resources if item.name.startswith("BUSCO "))
            cdd = next(item for item in status.resources if item.name == "CDD")
            self.assertFalse(busco.ready)
            self.assertTrue(cdd.ready)

            (lineage / "dataset.cfg").write_text("[busco]\n", encoding="utf-8")
            status = collect_database_status(config)
            busco = next(item for item in status.resources if item.name.startswith("BUSCO "))
            self.assertTrue(busco.ready)
