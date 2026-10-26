import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from msspack.config import ConfigError, load_config


class ConfigTests(unittest.TestCase):
    def setUp(self) -> None:
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.base = Path(temporary.name)
        self.config_path = self.base / "msspack.toml"
        fixture = Path(__file__).parent / "fixtures" / "minimal_pack" / "config.toml"
        self.original = fixture.read_text(encoding="utf-8")
        self.config_path.write_text(self.original, encoding="utf-8")

    def test_load_config_populates_defaults_and_resolves_paths(self) -> None:
        self.config_path.write_text(
            self.original.replace("run_gapjust = false\n", "")
            .replace("locus_tag_digits = 6\n", "")
            .replace('datatype = "WGS"\n', ""),
            encoding="utf-8",
        )
        config = load_config(self.config_path)

        self.assertEqual(config.project.name, "Fixture")
        self.assertEqual(config.sample.locus_tag_digits, 6)
        self.assertEqual(config.submission.datatype, "WGS")
        self.assertEqual(config.submission.hold_date, "20271231")
        self.assertTrue(config.pipeline.run_gapjust)
        self.assertEqual(config.pipeline.coordinate_duplicate_policy, "longest_valid_cds")
        self.assertEqual(config.plots.coordinate_duplicate_limit, 50)
        self.assertEqual(config.tools.java, "java")
        self.assertEqual(config.busco.command, "busco")
        self.assertTrue(config.busco.auto_lineage)
        self.assertTrue(config.busco.run_cds)
        self.assertFalse(config.busco.run_genome)
        self.assertEqual(config.busco.cds_mode, "transcriptome")
        self.assertFalse(config.functional_annotation.enabled)
        self.assertTrue(config.functional_annotation.swissprot_enabled)
        self.assertTrue(config.functional_annotation.pfam_enabled)
        self.assertEqual(config.fasta_path, (self.base / "input.fa").resolve())
        self.assertEqual(config.gff_path, (self.base / "input.gff3").resolve())
        self.assertEqual(config.output_dir, (self.base / "build" / "Fixture").resolve())
        self.assertEqual(config.database_dir, (self.base / "msspack_db").resolve())
        self.assertEqual(config.busco_database_dir, (self.base / "msspack_db" / "busco").resolve())

        shared_root = self.base / "shared-databases"
        with patch.dict("os.environ", {"MSSPACK_DB_DIR": str(shared_root)}):
            self.assertEqual(config.database_dir, shared_root)
            self.assertEqual(config.busco_database_dir, shared_root / "busco")

    def test_load_config_preserves_explicit_settings(self) -> None:
        text = self.original.replace('locus_tag = "Fix"', 'locus_tag = "Test_123-demo.v1"')
        text = text.replace(
            'collection_date = "2026-01-02"',
            'collection_date = "2020-06-10/2020-10-14"',
        ).replace("run_gapjust = false", 'run_gapjust = false\ncoordinate_duplicate_policy = "first"')
        self.config_path.write_text(
            text + """
[plots]
coordinate_duplicate_limit = 125
[functional_annotation]
enabled = true
threads = 3
min_identity = 30
[functional_annotation.consistency]
enabled = true
harmonize_safe_equivalents = true
auto_resolve_conflicts = false
family_identity = 72.0
source_pair_min_pairs = 7
[functional_annotation.taxonomy]
enabled = true
target_taxon_id = 45172
resolve_scientific_name = false
offline = true
busco_crosscheck = false
strict = true
distant_specificity_identity = 55.0
""",
            encoding="utf-8",
        )
        config = load_config(self.config_path)

        self.assertEqual(config.sample.locus_tag, "Test_123-demo.v1")
        self.assertEqual(config.sample.collection_date, "2020-06-10/2020-10-14")
        self.assertEqual(config.pipeline.coordinate_duplicate_policy, "first")
        self.assertEqual(config.plots.coordinate_duplicate_limit, 125)
        annotation = config.functional_annotation
        self.assertTrue(annotation.enabled)
        self.assertEqual(annotation.threads, 3)
        self.assertEqual(annotation.min_identity, 30.0)
        self.assertTrue(annotation.consistency.enabled)
        self.assertTrue(annotation.consistency.harmonize_safe_equivalents)
        self.assertFalse(annotation.consistency.auto_resolve_conflicts)
        self.assertEqual(annotation.consistency.family_identity, 72.0)
        self.assertEqual(annotation.consistency.source_pair_min_pairs, 7)
        self.assertTrue(annotation.taxonomy.enabled)
        self.assertEqual(annotation.taxonomy.target_taxon_id, 45172)
        self.assertFalse(annotation.taxonomy.resolve_scientific_name)
        self.assertTrue(annotation.taxonomy.offline)
        self.assertFalse(annotation.taxonomy.busco_crosscheck)
        self.assertTrue(annotation.taxonomy.strict)
        self.assertEqual(annotation.taxonomy.distant_specificity_identity, 55.0)

    def test_load_config_allows_omitting_hold_date_for_immediate_release(self) -> None:
        self.config_path.write_text(
            self.original.replace('hold_date = "20271231"\n', ""), encoding="utf-8",
        )
        self.assertEqual(load_config(self.config_path).submission.hold_date, "")

    def test_load_config_rejects_invalid_settings(self) -> None:
        cases = (
            (
                '[project]\nname = "Fixture"\noutput_dir = "build/Fixture"',
                'project = "oops"',
                "project",
            ),
            ("", '\n[busco]\nauto_lineage = false\nlineage_dataset = ""\n', "lineage_dataset"),
            ('locus_tag = "Fix"', "locus_tag = '../victim'", "sample.locus_tag"),
            ('locus_tag = "Fix"', r"locus_tag = '..\victim'", "sample.locus_tag"),
            ('locus_tag = "Fix"', "locus_tag = 'C:victim'", "sample.locus_tag"),
            ('hold_date = "20271231"', 'hold_date = "2026-12-31"', "hold_date"),
            ('hold_date = "20271231"', 'hold_date = "20271399"', "hold_date"),
            ('collection_date = "2026-01-02"', 'collection_date = "2026/01/01"', "collection_date"),
            (
                'collection_date = "2026-01-02"',
                'collection_date = "2020-10-14/2020-06-10"',
                "collection_date",
            ),
            (
                'feature_with_gap = "misc_feature"',
                'feature_with_gap = "bad-value"',
                "feature_with_gap",
            ),
            (
                "run_gapjust = false",
                'run_gapjust = false\ncoordinate_duplicate_policy = "unsupported"',
                "coordinate_duplicate_policy",
            ),
            ("", "\n[plots]\ncoordinate_duplicate_limit = 0\n", "plots.coordinate_duplicate_limit"),
            ("", "\n[databases]\nlock_poll_seconds = 0\n", "databases.lock_poll_seconds"),
            ("", "\n[databases]\nlock_poll_seconds = nan\n", "databases.lock_poll_seconds"),
            ("", "\n[databases]\nlock_poll_seconds = inf\n", "databases.lock_poll_seconds"),
            (
                "", "\n[databases]\nlock_heartbeat_seconds = 60\nlock_stale_seconds = 60\n",
                "databases.lock_stale_seconds",
            ),
            ('genome_coverage = "10x"', 'genome_coverage = "tenfold"', "genome_coverage"),
            ("run_gapjust = false", 'run_gapjust = "false"', "pipeline.run_gapjust"),
            (
                "run_gapjust = false",
                "run_gapjust = false\nvalidate_with_parsr = false",
                "pipeline.validate_with_parsr",
            ),
            (
                "replace_product_patterns = []",
                'replace_product_patterns = ["("]',
                "replace_product_patterns",
            ),
            ('genetic_code = "1"', 'genetic_code = "999"', "genetic code"),
            (
                'title = "Fixture genome sequencing"',
                'title = """Fixture genome\nsequencing"""',
                "control characters",
            ),
            ("", "\n[functional_annotation]\nthreads = true\n", "functional_annotation.threads"),
            ("", "\n[functional_annotation]\nthreads = 3.0\n", "functional_annotation.threads"),
            ("", "\n[functional_annotation]\nmin_identity = false\n", "functional_annotation.min_identity"),
            ("", '\n[functional_annotation]\nmin_identity = "3"\n', "functional_annotation.min_identity"),
            ("", "\n[functional_annotation]\nenabled = 0\n", "functional_annotation.enabled"),
            ("", "\n[functional_annotation]\nevalue = nan\n", "functional_annotation.evalue"),
            ("", "\n[functional_annotation]\nevalue = inf\n", "functional_annotation.evalue"),
            ("", "\n[functional_annotation]\nevalue = -inf\n", "functional_annotation.evalue"),
            (
                "", "\n[functional_annotation.taxonomy]\ntarget_taxon_id = true\n",
                "functional_annotation.taxonomy.target_taxon_id",
            ),
            (
                "", "\n[functional_annotation.taxonomy]\ndistant_specificity_identity = nan\n",
                "functional_annotation.taxonomy.distant_specificity_identity",
            ),
            (
                "", "\n[functional_annotation.consistency]\nfamily_identity = false\n",
                "functional_annotation.consistency.family_identity",
            ),
            (
                "", "\n[functional_annotation.consistency]\nfamily_identity = inf\n",
                "functional_annotation.consistency.family_identity",
            ),
            (
                "", "\n[functional_annotation.taxonomy]\nunknown_taxonomy_option = true\n",
                "functional_annotation.taxonomy.unknown_taxonomy_option",
            ),
            (
                "", "\n[functional_annotation.consistency]\nunknown_threshold = 1\n",
                "functional_annotation.consistency.unknown_threshold",
            ),
        )
        for old, new, message in cases:
            with self.subTest(setting=message, value=new):
                text = self.original.replace(old, new) if old else self.original + new
                self.config_path.write_text(text, encoding="utf-8")
                with self.assertRaisesRegex(ConfigError, message):
                    load_config(self.config_path)


if __name__ == "__main__":
    unittest.main()
