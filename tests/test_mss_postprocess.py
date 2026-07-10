import tempfile
import unittest
from pathlib import Path

from msspack.mss_postprocess import convert_cds_features_to_misc


class MssPostprocessTests(unittest.TestCase):
    def test_convert_cds_features_to_misc_updates_target_block_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            mss_in = base / "input.mss"
            genes_in = base / "genes.txt"
            mss_out = base / "output.mss"
            changed_ids = base / "changed.txt"
            mss_in.write_text(
                "".join(
                    [
                        "chr1\tsource\t1..100\tff_definition\tExample\n",
                        "\tassembly_gap\t20..30\testimated_length\tunknown\n",
                        "\t\t\tgap_type\twithin scaffold\n",
                        "\tCDS\t1..9\tlocus_tag\tOrg_Gene1\n",
                        "\t\t\tnote\ttranscript_id:Gene1-T1\n",
                        "\t\t\tproduct\thypothetical protein\n",
                        "\t\t\ttransl_table\t1\n",
                        "\t\t\tcodon_start\t1\n",
                        "\tCDS\t40..48\tlocus_tag\tOrg_Gene2\n",
                        "\t\t\tnote\ttranscript_id:Gene2-T1\n",
                        "\t\t\tproduct\tknown protein\n",
                        "\t\t\ttransl_table\t1\n",
                        "\t\t\tcodon_start\t1\n",
                    ]
                ),
                encoding="utf-8",
            )
            genes_in.write_text("Gene1\n", encoding="utf-8")

            summary = convert_cds_features_to_misc(
                mss_input_path=mss_in,
                genes_input_path=genes_in,
                mss_output_path=mss_out,
                converted_gene_ids_path=changed_ids,
                locus_tag_prefix="Org",
            )

            text = mss_out.read_text(encoding="utf-8")
            self.assertIn("\tmisc_feature\t1..9\tlocus_tag\tOrg_Gene1\n", text)
            self.assertIn("\t\t\tnote\thypothetical protein\n", text)
            self.assertNotIn("transl_table\t1\n", text.split("\tCDS\t40..48", 1)[0])
            self.assertIn("\tCDS\t40..48\tlocus_tag\tOrg_Gene2\n", text)
            self.assertEqual(summary["edited_genes"], 1)
            self.assertEqual(summary["cds_input"], 2)
            self.assertEqual(summary["cds_output"], 1)
            self.assertEqual(summary["misc_feature_output"], 1)
            self.assertEqual(summary["converted_gene_ids"], ["Gene1"])
            self.assertEqual(changed_ids.read_text(encoding="utf-8"), "Gene1\n")

    def test_similar_numeric_suffix_does_not_convert_unrelated_gene(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            mss_in = base / "input.mss"
            genes_in = base / "genes.txt"
            mss_out = base / "output.mss"
            mss_in.write_text(
                "\tCDS\t1..9\tlocus_tag\tPrefix_bar_001\n"
                "\t\t\tproduct\tvalid protein\n"
                "\t\t\ttransl_table\t1\n"
                "\t\t\tcodon_start\t1\n",
                encoding="utf-8",
            )
            genes_in.write_text("foo_001\n", encoding="utf-8")

            summary = convert_cds_features_to_misc(
                mss_input_path=mss_in,
                genes_input_path=genes_in,
                mss_output_path=mss_out,
                locus_tag_prefix="Prefix",
            )

            self.assertEqual(summary["edited_genes"], 0)
            self.assertIn("\tCDS\t", mss_out.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
