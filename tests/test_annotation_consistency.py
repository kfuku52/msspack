from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from xml.etree import ElementTree

from msspack.annotation_consistency import (
    audit_annotation_consistency,
    classify_product_names,
    run_family_similarity_search,
)
from msspack.annotation_consistency_plots import (
    prepare_name_consistency_plot_data,
    prepare_source_consistency_plot_data,
    write_name_consistency_pdf,
    write_name_consistency_svg,
    write_source_consistency_pdf,
    write_source_consistency_svg,
)
from msspack.config_models import (
    FunctionalAnnotationConfig,
    FunctionalAnnotationConsistencyConfig,
)


class AnnotationConsistencyTests(unittest.TestCase):
    def test_classifies_safe_equivalents_granularity_and_specific_conflicts(self) -> None:
        exact = classify_product_names(
            "ATP synthase subunit alpha",
            "ATP-synthase subunit alpha",
        )
        canonical = classify_product_names(
            "BTB/POZ domain-containing protein",
            "BTB domain-containing protein",
        )
        compatible = classify_product_names(
            "protein kinase domain-containing protein",
            "serine threonine protein kinase",
            similarity_threshold=0.30,
        )
        conflict = classify_product_names(
            "cellulose synthase-like protein B4",
            "cellulose synthase-like protein H1",
        )

        self.assertEqual(exact.relation, "exact")
        self.assertEqual(canonical.relation, "canonical_equivalent")
        self.assertEqual(compatible.relation, "compatible_granularity")
        self.assertEqual(conflict.relation, "review")
        self.assertEqual(conflict.reason_code, "specific_modifier_conflict")

    def test_family_search_writes_one_broad_all_vs_all_table(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            proteins = base / "proteins.faa"
            proteins.write_text(">t1\nMPEPTIDE\n>t2\nMPEPTIDE\n", encoding="utf-8")

            def fake_run(command: list[str], *, log_path: Path) -> None:
                log_path.write_text("ok\n", encoding="utf-8")
                if "blastp" in command:
                    output = Path(command[command.index("--out") + 1])
                    output.write_text(
                        "t1\tt1\t100\t100\t100\t0\t20\t8\t8\n"
                        "t1\tt2\t95\t100\t100\t1e-20\t18\t8\t8\n",
                        encoding="utf-8",
                    )

            config = FunctionalAnnotationConfig(
                consistency=FunctionalAnnotationConsistencyConfig(enabled=True)
            )
            with (
                patch(
                    "msspack.annotation_consistency.run_command",
                    side_effect=fake_run,
                ),
                patch(
                    "msspack.annotation_consistency._diamond_version",
                    return_value="diamond test",
                ),
            ):
                run_family_similarity_search(
                    protein_fasta_path=proteins,
                    output_path=base / "similarity.tsv",
                    provenance_path=base / "provenance.json",
                    log_path=base / "search.log",
                    metrics_path=base / "search.json",
                    config=config,
                )

            lines = (base / "similarity.tsv").read_text(encoding="utf-8").splitlines()
            self.assertEqual(lines[0].split("\t")[0:2], ["qseqid", "sseqid"])
            self.assertEqual(len(lines), 3)
            self.assertIn("--max-target-seqs 0", (base / "search.log").read_text())

    def test_audit_partitions_genes_and_only_harmonizes_safe_equivalents(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            annotation = base / "annotation.tsv"
            evidence = base / "evidence.tsv"
            similarity = base / "similarity.tsv"
            annotation.write_text(
                "ID\tDescription\tLocus_tag\n"
                "t1\tATP synthase alpha subunit\tG1\n"
                "t2\tATP synthase subunit alpha\tG2\n"
                "t3\tcellulose synthase-like protein B4\tG3\n"
                "t4\tcellulose synthase-like protein H1\tG4\n"
                "t5\thypothetical protein\tG5\n"
                "t6\tprotein kinase domain-containing protein\tG6\n"
                "t7\tprotein kinase domain-containing protein\tG7\n",
                encoding="utf-8",
            )
            evidence.write_text(
                "ID\tLocus_tag\tassigned_product\tsource\tconfidence\tbitscore\tconsensus_score\n"
                "t1\tG1\tATP synthase alpha subunit\tswissprot\thigh\t500\t1\n"
                "t2\tG2\tATP synthase subunit alpha\tuniref90\tmoderate\t480\t0.9\n"
                "t3\tG3\tcellulose synthase-like protein B4\tswissprot\thigh\t400\t1\n"
                "t4\tG4\tcellulose synthase-like protein H1\tuniref90\tmoderate\t390\t0.9\n"
                "t5\tG5\thypothetical protein\tnone\tnone\t\t\n"
                "t6\tG6\tprotein kinase domain-containing protein\tpfam\tlow\t200\t0.7\n"
                "t7\tG7\tprotein kinase domain-containing protein\tuniref90\tmoderate\t250\t0.8\n",
                encoding="utf-8",
            )
            similarity.write_text(
                "qseqid\tsseqid\tpident\tqcovhsp\tscovhsp\tevalue\tbitscore\tqlen\tslen\n"
                "t1\tt2\t95\t100\t100\t1e-80\t500\t300\t300\n"
                "t3\tt4\t95\t100\t100\t1e-70\t450\t300\t300\n"
                "t2\tt3\t45\t65\t65\t1e-20\t100\t300\t300\n"
                "t6\tt7\t75\t90\t90\t1e-40\t300\t300\t300\n",
                encoding="utf-8",
            )

            output_paths = {
                "output_annotation_table_path": base / "consistent.tsv",
                "gene_output_path": base / "gene.tsv",
                "family_output_path": base / "family.tsv",
                "pair_output_path": base / "pair.tsv",
                "review_output_path": base / "review.tsv",
                "summary_output_path": base / "summary.tsv",
                "source_pair_output_path": base / "source.tsv",
                "log_path": base / "audit.log",
                "metrics_path": base / "audit.json",
            }
            config = FunctionalAnnotationConfig(
                consistency=FunctionalAnnotationConsistencyConfig(
                    enabled=True,
                    harmonize_safe_equivalents=True,
                )
            )
            audit_annotation_consistency(
                annotation_table_path=annotation,
                evidence_path=evidence,
                similarity_path=similarity,
                config=config,
                **output_paths,
            )

            with (base / "gene.tsv").open("r", encoding="utf-8", newline="") as handle:
                genes = {row["ID"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(genes["t1"]["name_consistency"], "consistent")
            self.assertEqual(genes["t3"]["name_consistency"], "resolved")
            self.assertEqual(genes["t3"]["conflict_severity"], "high")
            self.assertEqual(
                genes["t4"]["resolution_action"],
                "higher_priority_near_identical_product_propagated",
            )
            self.assertEqual(
                genes["t4"]["resolved_product"],
                "cellulose synthase-like protein B4",
            )
            self.assertEqual(genes["t5"]["name_consistency"], "unannotated")
            self.assertEqual(genes["t6"]["name_consistency"], "consistent")
            self.assertEqual(genes["t6"]["reason_code"], "close_family_names_compatible")
            self.assertEqual(len({row["name_consistency"] for row in genes.values()}), 3)

            with (base / "consistent.tsv").open(
                "r", encoding="utf-8", newline=""
            ) as handle:
                products = {
                    row["ID"]: row["Description"]
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(products["t2"], "ATP synthase alpha subunit")
            self.assertEqual(products["t4"], "cellulose synthase-like protein B4")

            with (base / "summary.tsv").open(
                "r", encoding="utf-8", newline=""
            ) as handle:
                summaries = {
                    row["tier"]: row for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(summaries["near_identical"]["evaluable_clusters"], "2")
            self.assertEqual(summaries["near_identical"]["review_clusters"], "1")
            self.assertEqual(summaries["near_identical"]["canonical_equivalent_clusters"], "1")

    def test_consistency_plots_write_chart_ready_tsv_svg_and_pdf(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_dir:
            base = Path(tmp_dir)
            summary = base / "summary.tsv"
            source = base / "source.tsv"
            summary.write_text(
                "tier\tthreshold\tevaluable_clusters\texact_clusters\tcanonical_equivalent_clusters\t"
                "compatible_granularity_clusters\treview_clusters\tstrict_equivalence_rate\t"
                "biological_compatibility_rate\n"
                "near_identical\tidentity>=90%; mutual coverage>=90%\t10\t6\t2\t1\t1\t0.8\t0.9\n"
                "family\tidentity>=70%; mutual coverage>=80%\t20\t10\t4\t3\t3\t0.7\t0.85\n"
                "broad\tidentity>=40%; mutual coverage>=60%\t40\t12\t8\t8\t12\t0.5\t0.7\n",
                encoding="utf-8",
            )
            source.write_text(
                "tier\tsource_a\tsource_b\tannotated_direct_pairs\treview_pairs\treview_rate\tmeets_minimum_pairs\n"
                "near_identical\tswissprot\tswissprot\t100\t2\t0.02\ttrue\n"
                "near_identical\tswissprot\tuniref90\t4\t2\t0.5\tfalse\n",
                encoding="utf-8",
            )
            name_rows = prepare_name_consistency_plot_data(summary, base / "name.tsv")
            source_rows = prepare_source_consistency_plot_data(source, base / "source-plot.tsv")
            write_name_consistency_svg(name_rows, base / "name.svg")
            write_name_consistency_pdf(name_rows, base / "name.pdf")
            write_source_consistency_svg(source_rows, base / "source.svg")
            write_source_consistency_pdf(source_rows, base / "source.pdf")

            ElementTree.parse(base / "name.svg")
            ElementTree.parse(base / "source.svg")
            self.assertTrue((base / "name.pdf").read_bytes().startswith(b"%PDF-1.4"))
            self.assertTrue((base / "source.pdf").read_bytes().startswith(b"%PDF-1.4"))
            name_svg = (base / "name.svg").read_text()
            source_svg = (base / "source.svg").read_text()
            self.assertIn("each bar totals 100%", name_svg)
            self.assertIn("id &gt;=90%", name_svg)
            self.assertIn("mutual cov &gt;=90%", name_svg)
            self.assertIn("id &gt;=70%", name_svg)
            self.assertIn("mutual cov &gt;=80%", name_svg)
            self.assertIn("id &gt;=40%", name_svg)
            self.assertIn("mutual cov &gt;=60%", name_svg)
            self.assertIn("grey cells", source_svg)
            for svg_text in (name_svg, source_svg):
                self.assertIn('width="3.6in"', svg_text)
                self.assertIn("font-size:8pt", svg_text)
                self.assertNotIn("font-size:7pt", svg_text)
                self.assertNotIn("font-size:14pt", svg_text)
            name_pdf = (base / "name.pdf").read_bytes().decode("latin-1")
            source_pdf = (base / "source.pdf").read_bytes().decode("latin-1")
            self.assertIn("/MediaBox [0 0 259.20 260.00]", name_pdf)
            self.assertIn("/MediaBox [0 0 259.20 186.00]", source_pdf)
            for pdf_text in (name_pdf, source_pdf):
                self.assertNotIn(" 7 Tf", pdf_text)
                self.assertNotIn(" 14 Tf", pdf_text)


if __name__ == "__main__":
    unittest.main()
