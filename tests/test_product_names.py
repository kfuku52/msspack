from __future__ import annotations

import unittest

from msspack.product_names import (
    ProductNameContext,
    audit_product_name,
    standardize_product_name,
)


class ProductNameTests(unittest.TestCase):
    def test_standardizes_case_localization_and_evolutionary_wording(self) -> None:
        localized = standardize_product_name(
            "Glyceraldehyde-3-phosphate dehydrogenase B, chloroplastic",
            source="swissprot",
        )
        homolog = standardize_product_name(
            "Protein unc-50 homolog",
            source="uniref90",
        )

        self.assertEqual(
            localized.product,
            "glyceraldehyde-3-phosphate dehydrogenase B",
        )
        self.assertIn("lowercased_initial", localized.actions)
        self.assertIn("removed_localization_suffix", localized.actions)
        self.assertEqual(homolog.product, "unc-50-like protein")
        self.assertIn("replaced_evolutionary_term", homolog.actions)
        self.assertEqual(
            standardize_product_name("PsbQ-like protein 3 chloroplastic").product,
            "PsbQ-like protein 3",
        )

    def test_keeps_acronyms_chemical_symbols_and_established_truncated_name(self) -> None:
        self.assertEqual(
            standardize_product_name("DNA repair protein RAD51").product,
            "DNA repair protein RAD51",
        )
        self.assertEqual(
            standardize_product_name("Fe(3+)/Cu(2+)-chelate reductase").product,
            "Fe(3+)/Cu(2+)-chelate reductase",
        )
        truncated = standardize_product_name("Group 2 truncated hemoglobin 3-2")
        self.assertEqual(truncated.product, "group 2 truncated hemoglobin 3-2")
        self.assertNotIn("partial_or_fragment_term", truncated.warnings)

    def test_reduces_cdd_explanations_and_decodes_html_entities(self) -> None:
        result = standardize_product_name(
            "KOW: an acronym for the authors&apos; surnames "
            "(Kyrpides, Ouzounis and Woese) domain-containing protein",
            source="cdd",
            evidence_id="cl00001",
        )

        self.assertEqual(result.product, "KOW domain-containing protein")
        self.assertIn("decoded_html_entity", result.actions)
        self.assertIn("removed_explanatory_text", result.actions)
        self.assertEqual(result.warnings, ())

    def test_removes_cog_category_and_normalizes_contains_description(self) -> None:
        result = standardize_product_name(
            "Lipoprotein NlpI, contains TPR repeats "
            "[Cell wall/membrane/envelope biogenesis] domain-containing protein",
            source="cdd",
            evidence_id="COG4785",
            context=ProductNameContext(domain="Eukaryota", kingdom="Viridiplantae"),
        )

        self.assertEqual(result.product, "TPR repeat-containing protein")
        self.assertIn("removed_cog_category", result.actions)
        self.assertIn("reduced_contains_description", result.actions)
        self.assertEqual(result.warnings, ())

    def test_generalizes_cross_domain_names_for_eukaryotes(self) -> None:
        result = standardize_product_name(
            "Type IV pilus assembly protein PilF/PilW "
            "[Cell motility, Extracellular structures] domain-containing protein",
            source="cdd",
            evidence_id="COG3063",
            context=ProductNameContext(domain="Eukaryota", kingdom="Viridiplantae"),
        )

        self.assertEqual(result.product, "PilF/PilW family protein")
        self.assertIn("generalized_cross_domain_name", result.actions)
        self.assertEqual(result.warnings, ())

    def test_uses_stable_model_identifier_for_uninformative_domain_label(self) -> None:
        result = standardize_product_name(
            "Plant protein 1589 of unknown function (A_thal_3526) "
            "domain-containing protein",
            source="pfam",
            evidence_id="PF12345.7",
            context=ProductNameContext(domain="Eukaryota", kingdom="Viridiplantae"),
        )

        self.assertEqual(result.product, "PF12345 domain-containing protein")
        self.assertIn("used_stable_model_identifier", result.actions)
        self.assertEqual(result.warnings, ())

    def test_replaces_uninformative_product_with_ddbj_recommended_fallback(self) -> None:
        result = standardize_product_name("Uncharacterized protein")

        self.assertEqual(result.product, "hypothetical protein")
        self.assertFalse(result.informative)
        self.assertIn("used_hypothetical_product", result.actions)

    def test_audit_reports_residual_nonstandard_terms(self) -> None:
        self.assertEqual(
            set(audit_product_name("probable enzyme homolog, mitochondrial")),
            {
                "evolutionary_term",
                "discouraged_uncertainty",
                "localization_suffix",
            },
        )


if __name__ == "__main__":
    unittest.main()
