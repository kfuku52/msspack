import tempfile
import unittest
from math import inf, nan
from pathlib import Path
from unittest.mock import patch

from msspack.config import ConfigError, DatabasesConfig, PipelineConfig, load_config
from msspack.config_loading import (
    _validate_raw_config,
    load_functional_annotation_config,
    load_pipeline_config,
)
from msspack.config_validation import (
    ensure_collection_date,
    validate_databases_config,
    validate_functional_annotation_config,
    validate_pipeline_config,
)


class ConfigTests(unittest.TestCase):
    def test_collection_date_accepts_an_iso_date_range(self) -> None:
        ensure_collection_date(
            "2020-06-10/2020-10-14",
            "sample.collection_date",
        )

    def test_collection_date_rejects_a_reversed_range(self) -> None:
        with self.assertRaises(ConfigError):
            ensure_collection_date(
                "2020-10-14/2020-06-10",
                "sample.collection_date",
            )

    def test_loads_and_validates_nested_annotation_consistency_config(self) -> None:
        raw = {
            "functional_annotation": {
                "enabled": True,
                "consistency": {
                    "enabled": True,
                    "harmonize_safe_equivalents": True,
                    "auto_resolve_conflicts": False,
                    "family_identity": 72.0,
                    "source_pair_min_pairs": 7,
                },
            }
        }

        _validate_raw_config(raw)
        config = load_functional_annotation_config(raw["functional_annotation"])
        validate_functional_annotation_config(config)

        self.assertTrue(config.consistency.enabled)
        self.assertTrue(config.consistency.harmonize_safe_equivalents)
        self.assertFalse(config.consistency.auto_resolve_conflicts)
        self.assertEqual(config.consistency.family_identity, 72.0)
        self.assertEqual(config.consistency.source_pair_min_pairs, 7)

    def test_loads_and_validates_annotation_taxonomy_config(self) -> None:
        raw = {
            "functional_annotation": {
                "taxonomy": {
                    "enabled": True,
                    "target_taxon_id": 45172,
                    "resolve_scientific_name": False,
                    "offline": True,
                    "busco_crosscheck": False,
                    "strict": True,
                    "distant_specificity_identity": 55.0,
                },
            }
        }

        _validate_raw_config(raw)
        config = load_functional_annotation_config(raw["functional_annotation"])
        validate_functional_annotation_config(config)

        self.assertEqual(config.taxonomy.target_taxon_id, 45172)
        self.assertFalse(config.taxonomy.resolve_scientific_name)
        self.assertTrue(config.taxonomy.offline)
        self.assertFalse(config.taxonomy.busco_crosscheck)
        self.assertTrue(config.taxonomy.strict)
        self.assertEqual(config.taxonomy.distant_specificity_identity, 55.0)

    def test_rejects_unknown_annotation_taxonomy_key(self) -> None:
        with self.assertRaises(ConfigError):
            _validate_raw_config(
                {
                    "functional_annotation": {
                        "taxonomy": {"unknown_taxonomy_option": True},
                    }
                }
            )

    def test_rejects_unknown_annotation_consistency_key(self) -> None:
        with self.assertRaises(ConfigError):
            _validate_raw_config(
                {
                    "functional_annotation": {
                        "consistency": {"unknown_threshold": 1},
                    }
                }
            )

    def test_load_config_populates_defaults_and_resolves_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            config_path = base / "msspack.toml"
            config_path.write_text(
                "\n".join(
                    [
                        "[project]",
                        'name = "Demo"',
                        "",
                        "[inputs]",
                        'fasta = "input.fa"',
                        'gff = "input.gff3"',
                        "",
                        "[sample]",
                        'locus_tag = "Demo"',
                        'scientific_name = "Demo species"',
                        "",
                        "[submission]",
                        'hold_date = "20261231"',
                        'bioproject = "PRJDB1"',
                        'biosample = "SAMD1"',
                        "",
                        "[submitter]",
                        'ab_name = ["Fukushima,K."]',
                        'contact = "Kenji Fukushima"',
                        'institute = "NIG"',
                        'department = "Lab"',
                        'country = "Japan"',
                        'state = "Shizuoka"',
                        'city = "Mishima"',
                        'street = "1111 Yata"',
                        'zip = "411-8540"',
                        'phone = "81-00-0000-0000"',
                        'email = "x@example.org"',
                        "",
                        "[reference]",
                        'title = "Demo sequencing"',
                        'ab_name = ["Fukushima,K."]',
                        "year = 2026",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            config = load_config(config_path)

            self.assertEqual(config.project.name, "Demo")
            self.assertEqual(config.sample.locus_tag_digits, 6)
            self.assertEqual(config.submission.datatype, "WGS")
            self.assertTrue(config.pipeline.run_gapjust)
            self.assertEqual(
                config.pipeline.coordinate_duplicate_policy,
                "longest_valid_cds",
            )
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
            self.assertEqual(config.fasta_path, (base / "input.fa").resolve())
            self.assertEqual(config.gff_path, (base / "input.gff3").resolve())
            self.assertEqual(config.output_dir, (base / "build" / "Demo").resolve())
            self.assertEqual(config.database_dir, (base / "msspack_db").resolve())
            self.assertEqual(
                config.busco_database_dir,
                (base / "msspack_db" / "busco").resolve(),
            )

            shared_root = base / "shared-databases"
            with patch.dict("os.environ", {"MSSPACK_DB_DIR": str(shared_root)}):
                self.assertEqual(config.database_dir, shared_root)
                self.assertEqual(
                    config.busco_database_dir,
                    shared_root / "busco",
                )

    def test_coordinate_duplicate_plot_limit_is_configurable(self) -> None:
        raw = {"plots": {"coordinate_duplicate_limit": 125}}

        _validate_raw_config(raw)

        from msspack.config_loading import load_plots_config
        from msspack.config_validation import validate_plots_config

        plots = load_plots_config(raw["plots"])
        validate_plots_config(plots)
        self.assertEqual(plots.coordinate_duplicate_limit, 125)

    def test_coordinate_duplicate_plot_limit_must_be_positive(self) -> None:
        from msspack.config import PlotsConfig
        from msspack.config_validation import validate_plots_config

        with self.assertRaises(ConfigError):
            validate_plots_config(PlotsConfig(coordinate_duplicate_limit=0))

    def test_coordinate_duplicate_selection_policy_is_configurable(self) -> None:
        raw = {"pipeline": {"coordinate_duplicate_policy": "first"}}

        _validate_raw_config(raw)
        pipeline = load_pipeline_config(raw["pipeline"])
        validate_pipeline_config(pipeline)

        self.assertEqual(pipeline.coordinate_duplicate_policy, "first")

    def test_coordinate_duplicate_selection_policy_rejects_unknown_value(
        self,
    ) -> None:
        with self.assertRaises(ConfigError):
            validate_pipeline_config(
                PipelineConfig(coordinate_duplicate_policy="unsupported")
            )

    def test_database_lock_settings_must_be_positive(self) -> None:
        with self.assertRaises(ConfigError):
            validate_databases_config(DatabasesConfig(lock_poll_seconds=0))

    def test_database_lock_settings_must_be_finite(self) -> None:
        for value in (nan, inf):
            with self.subTest(value=value), self.assertRaises(ConfigError):
                validate_databases_config(DatabasesConfig(lock_poll_seconds=value))

    def test_database_lock_stale_interval_must_exceed_heartbeat(self) -> None:
        with self.assertRaises(ConfigError):
            validate_databases_config(
                DatabasesConfig(
                    lock_heartbeat_seconds=60,
                    lock_stale_seconds=60,
                )
            )

    def test_load_config_rejects_non_table_section(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            config_path = Path(tmp_dir) / "msspack.toml"
            config_path.write_text(
                "\n".join(
                    [
                        'project = "oops"',
                        "",
                        "[inputs]",
                        'fasta = "input.fa"',
                        'gff = "input.gff3"',
                        "",
                        "[sample]",
                        'locus_tag = "Demo"',
                        'scientific_name = "Demo species"',
                        "",
                        "[submission]",
                        'hold_date = "20261231"',
                        'bioproject = "PRJDB1"',
                        'biosample = "SAMD1"',
                        "",
                        "[submitter]",
                        'ab_name = ["Fukushima,K."]',
                        'contact = "Kenji Fukushima"',
                        'institute = "NIG"',
                        'department = "Lab"',
                        'country = "Japan"',
                        'state = "Shizuoka"',
                        'city = "Mishima"',
                        'street = "1111 Yata"',
                        'zip = "411-8540"',
                        'phone = "81-00-0000-0000"',
                        'email = "x@example.org"',
                        "",
                        "[reference]",
                        'title = "Demo sequencing"',
                        'ab_name = ["Fukushima,K."]',
                        "year = 2026",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            with self.assertRaises(ConfigError):
                load_config(config_path)

    def test_load_config_rejects_busco_without_lineage_when_auto_disabled(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            config_path = Path(tmp_dir) / "msspack.toml"
            config_path.write_text(
                "\n".join(
                    [
                        "[project]",
                        'name = "Demo"',
                        "",
                        "[inputs]",
                        'fasta = "input.fa"',
                        'gff = "input.gff3"',
                        "",
                        "[sample]",
                        'locus_tag = "Demo"',
                        'scientific_name = "Demo species"',
                        "",
                        "[submission]",
                        'hold_date = "20261231"',
                        'bioproject = "PRJDB1"',
                        'biosample = "SAMD1"',
                        "",
                        "[submitter]",
                        'ab_name = ["Fukushima,K."]',
                        'contact = "Kenji Fukushima"',
                        'institute = "NIG"',
                        'department = "Lab"',
                        'country = "Japan"',
                        'state = "Shizuoka"',
                        'city = "Mishima"',
                        'street = "1111 Yata"',
                        'zip = "411-8540"',
                        'phone = "81-00-0000-0000"',
                        'email = "x@example.org"',
                        "",
                        "[reference]",
                        'title = "Demo sequencing"',
                        'ab_name = ["Fukushima,K."]',
                        "year = 2026",
                        "",
                        "[busco]",
                        "auto_lineage = false",
                        'lineage_dataset = ""',
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            with self.assertRaises(ConfigError):
                load_config(config_path)

    def test_load_config_rejects_busco_when_all_targets_are_disabled(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            config_path = Path(tmp_dir) / "msspack.toml"
            config_path.write_text(
                "\n".join(
                    [
                        "[project]",
                        'name = "Demo"',
                        "",
                        "[inputs]",
                        'fasta = "input.fa"',
                        'gff = "input.gff3"',
                        "",
                        "[sample]",
                        'locus_tag = "Demo"',
                        'scientific_name = "Demo species"',
                        "",
                        "[submission]",
                        'hold_date = "20261231"',
                        'bioproject = "PRJDB1"',
                        'biosample = "SAMD1"',
                        "",
                        "[submitter]",
                        'ab_name = ["Fukushima,K."]',
                        'contact = "Kenji Fukushima"',
                        'institute = "NIG"',
                        'department = "Lab"',
                        'country = "Japan"',
                        'state = "Shizuoka"',
                        'city = "Mishima"',
                        'street = "1111 Yata"',
                        'zip = "411-8540"',
                        'phone = "81-00-0000-0000"',
                        'email = "x@example.org"',
                        "",
                        "[reference]",
                        'title = "Demo sequencing"',
                        'ab_name = ["Fukushima,K."]',
                        "year = 2026",
                        "",
                        "[busco]",
                        "run_cds = false",
                        "run_genome = false",
                        'lineage_dataset = "embryophyta_odb12"',
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            with self.assertRaises(ConfigError):
                load_config(config_path)

    def test_load_config_rejects_bad_hold_date_format(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            config_path = Path(tmp_dir) / "msspack.toml"
            config_path.write_text(
                "\n".join(
                    [
                        "[project]",
                        'name = "Demo"',
                        "",
                        "[inputs]",
                        'fasta = "input.fa"',
                        'gff = "input.gff3"',
                        "",
                        "[sample]",
                        'locus_tag = "Demo"',
                        'scientific_name = "Demo species"',
                        "",
                        "[submission]",
                        'hold_date = "2026-12-31"',
                        'bioproject = "PRJDB1"',
                        'biosample = "SAMD1"',
                        "",
                        "[submitter]",
                        'ab_name = ["Fukushima,K."]',
                        'contact = "Kenji Fukushima"',
                        'institute = "NIG"',
                        'department = "Lab"',
                        'country = "Japan"',
                        'state = "Shizuoka"',
                        'city = "Mishima"',
                        'street = "1111 Yata"',
                        'zip = "411-8540"',
                        'phone = "81-00-0000-0000"',
                        'email = "x@example.org"',
                        "",
                        "[reference]",
                        'title = "Demo sequencing"',
                        'ab_name = ["Fukushima,K."]',
                        "year = 2026",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            with self.assertRaises(ConfigError):
                load_config(config_path)

    def test_load_config_rejects_invalid_collection_date(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            config_path = Path(tmp_dir) / "msspack.toml"
            config_path.write_text(
                "\n".join(
                    [
                        "[project]",
                        'name = "Demo"',
                        "",
                        "[inputs]",
                        'fasta = "input.fa"',
                        'gff = "input.gff3"',
                        "",
                        "[sample]",
                        'locus_tag = "Demo"',
                        'scientific_name = "Demo species"',
                        'collection_date = "2026/01/01"',
                        "",
                        "[submission]",
                        'hold_date = "20261231"',
                        'bioproject = "PRJDB1"',
                        'biosample = "SAMD1"',
                        "",
                        "[submitter]",
                        'ab_name = ["Fukushima,K."]',
                        'contact = "Kenji Fukushima"',
                        'institute = "NIG"',
                        'department = "Lab"',
                        'country = "Japan"',
                        'state = "Shizuoka"',
                        'city = "Mishima"',
                        'street = "1111 Yata"',
                        'zip = "411-8540"',
                        'phone = "81-00-0000-0000"',
                        'email = "x@example.org"',
                        "",
                        "[reference]",
                        'title = "Demo sequencing"',
                        'ab_name = ["Fukushima,K."]',
                        "year = 2026",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            with self.assertRaises(ConfigError):
                load_config(config_path)

    def test_load_config_rejects_invalid_pipeline_choice(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            config_path = Path(tmp_dir) / "msspack.toml"
            config_path.write_text(
                "\n".join(
                    [
                        "[project]",
                        'name = "Demo"',
                        "",
                        "[inputs]",
                        'fasta = "input.fa"',
                        'gff = "input.gff3"',
                        "",
                        "[sample]",
                        'locus_tag = "Demo"',
                        'scientific_name = "Demo species"',
                        "",
                        "[submission]",
                        'hold_date = "20261231"',
                        'bioproject = "PRJDB1"',
                        'biosample = "SAMD1"',
                        "",
                        "[submitter]",
                        'ab_name = ["Fukushima,K."]',
                        'contact = "Kenji Fukushima"',
                        'institute = "NIG"',
                        'department = "Lab"',
                        'country = "Japan"',
                        'state = "Shizuoka"',
                        'city = "Mishima"',
                        'street = "1111 Yata"',
                        'zip = "411-8540"',
                        'phone = "81-00-0000-0000"',
                        'email = "x@example.org"',
                        "",
                        "[reference]",
                        'title = "Demo sequencing"',
                        'ab_name = ["Fukushima,K."]',
                        "year = 2026",
                        "",
                        "[pipeline]",
                        'feature_with_gap = "bad-value"',
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            with self.assertRaises(ConfigError):
                load_config(config_path)

    def test_rejects_impossible_hold_date(self) -> None:
        self._assert_fixture_rejected(
            lambda text: text.replace('hold_date = "20271231"', 'hold_date = "20271399"'),
            "submission.hold_date",
        )

    def test_rejects_quoted_boolean(self) -> None:
        self._assert_fixture_rejected(
            lambda text: text.replace("run_gapjust = false", 'run_gapjust = "false"'),
            "pipeline.run_gapjust",
        )

    def test_rejects_unknown_key(self) -> None:
        self._assert_fixture_rejected(
            lambda text: text.replace(
                "run_gapjust = false",
                "run_gapjust = false\nvalidate_with_parsr = false",
            ),
            "pipeline.validate_with_parsr",
        )

    def test_rejects_invalid_product_regex(self) -> None:
        self._assert_fixture_rejected(
            lambda text: text.replace(
                "replace_product_patterns = []",
                'replace_product_patterns = ["("]',
            ),
            "replace_product_patterns",
        )

    def test_rejects_unknown_genetic_code(self) -> None:
        self._assert_fixture_rejected(
            lambda text: text.replace('genetic_code = "1"', 'genetic_code = "999"'),
            "genetic code",
        )

    def test_rejects_control_characters_in_optional_strings(self) -> None:
        self._assert_fixture_rejected(
            lambda text: text.replace(
                'title = "Fixture genome sequencing"',
                'title = """Fixture genome\nsequencing"""',
            ),
            "control characters",
        )

    def _assert_fixture_rejected(self, transform, message: str) -> None:
        fixture = Path(__file__).parent / "fixtures" / "minimal_pack" / "config.toml"
        with tempfile.TemporaryDirectory() as tmp_dir:
            config_path = Path(tmp_dir) / "msspack.toml"
            config_path.write_text(transform(fixture.read_text(encoding="utf-8")), encoding="utf-8")
            with self.assertRaisesRegex(ConfigError, message):
                load_config(config_path)


if __name__ == "__main__":
    unittest.main()
