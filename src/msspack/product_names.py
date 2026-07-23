from __future__ import annotations

import html
import re
import unicodedata
from dataclasses import dataclass

HYPOTHETICAL_PRODUCT = "hypothetical protein"

_UNINFORMATIVE_PRODUCTS = {
    "hypothetical protein",
    "uncharacterized protein",
    "unknown protein",
    "protein of unknown function",
    "unnamed protein product",
    "protein",
    "domain-containing protein",
}

_AMERICAN_SPELLING = (
    (re.compile(r"\buncharacterised\b", re.IGNORECASE), "uncharacterized"),
    (re.compile(r"\bcharacterisation\b", re.IGNORECASE), "characterization"),
    (re.compile(r"\bsulphur\b", re.IGNORECASE), "sulfur"),
    (re.compile(r"\bsulphate\b", re.IGNORECASE), "sulfate"),
    (re.compile(r"\bhaemoglobin\b", re.IGNORECASE), "hemoglobin"),
    (re.compile(r"\bhaem\b", re.IGNORECASE), "heme"),
    (re.compile(r"\btumour\b", re.IGNORECASE), "tumor"),
    (re.compile(r"\bfibre\b", re.IGNORECASE), "fiber"),
    (re.compile(r"\bcentre\b", re.IGNORECASE), "center"),
)

_COG_CATEGORY_WORDS = re.compile(
    r"\b(?:"
    r"biogenesis|cell cycle|cell motility|chromatin structure|coenzyme metabolism|"
    r"defense mechanisms|energy production|extracellular structures|"
    r"function prediction|inorganic ion transport|lipid metabolism|"
    r"nucleotide metabolism|posttranslational modification|replication|"
    r"secondary metabolites|signal transduction|transcription|translation|"
    r"transport and metabolism"
    r")\b",
    re.IGNORECASE,
)

_LOCALIZATION_SUFFIX = re.compile(
    r"(?:,\s*|\s+)(?:"
    r"amyloplastic|chloroplastic|chloroplast|chromoplastic|cytosolic|cytoplasmic|"
    r"endoplasmic reticulum|extracellular|glyoxysomal|Golgi|lysosomal|"
    r"mitochondrial|nuclear|peroxisomal|plasma membrane|plastidial|secreted|vacuolar"
    r")(?:/(?:"
    r"amyloplastic|chloroplastic|chloroplast|chromoplastic|cytosolic|cytoplasmic|"
    r"endoplasmic reticulum|extracellular|glyoxysomal|Golgi|lysosomal|"
    r"mitochondrial|nuclear|peroxisomal|plasma membrane|plastidial|secreted|vacuolar"
    r"))*(?:-like)?$",
    re.IGNORECASE,
)

_EVOLUTIONARY_TERM = re.compile(
    r"\b(?:homolog|homologue|ortholog|orthologue|paralog|paralogue)s?\b",
    re.IGNORECASE,
)
_DISCOURAGED_UNCERTAINTY = re.compile(
    r"\b(?:candidate|likely|possible|potential|predicted|probable)\b",
    re.IGNORECASE,
)
_PARTIAL_TERM = re.compile(r"\b(?:fragment|partial|truncated?)\b", re.IGNORECASE)
_DATABASE_IDENTIFIER = re.compile(
    r"\b(?:COG|FOG|KOG)\d+\b|"
    r"\bEC\s+\d+\.\d+\.\d+\.\d+\b|"
    r"\bGO:\d+\b",
    re.IGNORECASE,
)
_HTML_ENTITY = re.compile(r"&(?:amp|apos|gt|lt|quot);", re.IGNORECASE)
_EXPLANATORY_PHRASE = re.compile(
    r"\b(?:an acronym for|and related proteins|and similar proteins|found in)\b",
    re.IGNORECASE,
)

_CHEMICAL_SYMBOLS = {
    "Ag",
    "Al",
    "Ca",
    "Cd",
    "Cl",
    "Co",
    "Cr",
    "Cu",
    "Fe",
    "Hg",
    "K",
    "Li",
    "Mg",
    "Mn",
    "Mo",
    "Na",
    "Ni",
    "Pb",
    "Se",
    "Si",
    "Sn",
    "Zn",
}

_EUKARYOTE_BACTERIAL_GENERALIZATIONS = (
    (
        re.compile(
            r"^bacterial surface antigen\s*\(?(D15)\)?\s+domain-containing protein$",
            re.IGNORECASE,
        ),
        r"\1 domain-containing protein",
    ),
    (
        re.compile(
            r"^MSCRAMM family adhesin clumping factor\s+(ClfA)"
            r"\s+domain-containing protein$",
            re.IGNORECASE,
        ),
        r"\1-like domain-containing protein",
    ),
    (
        re.compile(
            r"^methyl-accepting chemotaxis protein\s*\(MCP\),?\s*"
            r"signaling domain-containing protein$",
            re.IGNORECASE,
        ),
        "MCP signaling domain-containing protein",
    ),
    (
        re.compile(
            r"^type IV pilus assembly protein\s+(PilF/PilW).*?domain-containing protein$",
            re.IGNORECASE,
        ),
        r"\1 family protein",
    ),
    (
        re.compile(
            r"^lipoprotein NlpI,?\s*contains TPR repeats.*?domain-containing protein$",
            re.IGNORECASE,
        ),
        "TPR repeat-containing protein",
    ),
    (
        re.compile(
            r"^pneumococcal surface protein\s+(PspC).*?domain-containing protein$",
            re.IGNORECASE,
        ),
        r"\1-like domain-containing protein",
    ),
    (
        re.compile(
            r"^Clostridium epsilon toxin ETX/Bacillus mosquitocidal toxin MTX2"
            r"\s+domain-containing protein$",
            re.IGNORECASE,
        ),
        "ETX/MTX2 domain-containing protein",
    ),
    (
        re.compile(
            r"^(?:bacterial\s+)?extracellular solute-binding protein"
            r"\s+domain-containing protein$",
            re.IGNORECASE,
        ),
        "solute-binding domain-containing protein",
    ),
    (
        re.compile(
            r"^bacterial transferase hexapeptide\s*\([^)]*\)"
            r"\s+domain-containing protein$",
            re.IGNORECASE,
        ),
        "hexapeptide repeat-containing protein",
    ),
    (
        re.compile(
            r"^bacterial\s+(PH|Ig-like)\s+domain-containing protein$",
            re.IGNORECASE,
        ),
        r"\1 domain-containing protein",
    ),
    (
        re.compile(
            r"^(RcsF)\s+lipoprotein\s+domain-containing protein$",
            re.IGNORECASE,
        ),
        r"\1-like domain-containing protein",
    ),
    (
        re.compile(
            r"^antigen I/II family LPXTG-anchored adhesin"
            r"\s+domain-containing protein$",
            re.IGNORECASE,
        ),
        "LPXTG-anchor domain-containing protein",
    ),
)

_RESIDUAL_AUDITS = (
    ("html_entity", _HTML_ENTITY),
    ("evolutionary_term", _EVOLUTIONARY_TERM),
    ("discouraged_uncertainty", _DISCOURAGED_UNCERTAINTY),
    ("partial_or_fragment_term", _PARTIAL_TERM),
    ("localization_suffix", _LOCALIZATION_SUFFIX),
    ("explanatory_phrase", _EXPLANATORY_PHRASE),
    ("forbidden_character", re.compile(r'[\x00-\x1f\x7f"\\]')),
)


@dataclass(frozen=True)
class ProductNameContext:
    domain: str = ""
    kingdom: str = ""

    @classmethod
    def from_taxonomy_context(cls, payload: dict[str, object]) -> ProductNameContext:
        target = payload.get("target")
        if not isinstance(target, dict):
            return cls()
        lineage = target.get("lineage")
        if not isinstance(lineage, list):
            return cls()
        ranks: dict[str, str] = {}
        for item in lineage:
            if not isinstance(item, dict):
                continue
            rank = str(item.get("rank", "")).upper()
            scientific_name = str(item.get("scientific_name", "")).strip()
            if rank and scientific_name:
                ranks[rank] = scientific_name
        return cls(
            domain=ranks.get("DOMAIN", ""),
            kingdom=ranks.get("KINGDOM", ""),
        )


@dataclass(frozen=True)
class ProductNameStandardization:
    product: str
    actions: tuple[str, ...]
    warnings: tuple[str, ...]

    @property
    def informative(self) -> bool:
        return self.product.casefold() not in _UNINFORMATIVE_PRODUCTS


def _replace(
    value: str,
    pattern: re.Pattern[str],
    replacement: str,
    action: str,
    actions: list[str],
) -> str:
    updated = pattern.sub(replacement, value)
    if updated != value:
        actions.append(action)
    return updated


def _strip_cog_categories(value: str, actions: list[str]) -> str:
    def replacement(match: re.Match[str]) -> str:
        content = match.group(0)[1:-1]
        return "" if _COG_CATEGORY_WORDS.search(content) else match.group(0)

    updated = re.sub(r"\[[^]]+]", replacement, value)
    if updated != value:
        actions.append("removed_cog_category")
    return updated


def _singular_family_label(value: str) -> str:
    if re.search(r"(?:ases|proteins)$", value, re.IGNORECASE):
        return value[:-1]
    return value


def _simplify_domain_description(value: str, actions: list[str]) -> str:
    match = re.fullmatch(
        r"([^:]{1,50}):\s*an acronym for\b.*?\s+domain-containing protein",
        value,
        flags=re.IGNORECASE,
    )
    if match:
        actions.append("removed_explanatory_text")
        value = f"{match.group(1).strip()} domain-containing protein"

    match = re.fullmatch(
        r"(?:[A-Z][a-z]+\s+[a-z-]+\s+)?"
        r"([A-Za-z][A-Za-z0-9_.-]*\d[A-Za-z0-9_.-]*)"
        r"\s+and similar proteins\s+domain-containing protein",
        value,
        flags=re.IGNORECASE,
    )
    if match:
        actions.extend(("removed_organism_name", "generalized_related_proteins"))
        value = f"{match.group(1)} family protein"

    match = re.fullmatch(
        r"(.+?\bdomain):\s*superfamily domain-containing protein",
        value,
        flags=re.IGNORECASE,
    )
    if match:
        actions.append("normalized_domain_architecture")
        value = f"{match.group(1).strip()}-containing protein"

    match = re.fullmatch(
        r"(.+?\b(?:domain|family))\s+found in\b.*?\s+domain-containing protein",
        value,
        flags=re.IGNORECASE,
    )
    if match:
        label = match.group(1).strip()
        actions.append("removed_explanatory_text")
        value = (
            f"{label}-containing protein"
            if label.casefold().endswith("domain")
            else f"{label} protein"
        )

    match = re.fullmatch(
        r"(.+?)\s+and (?:related|similar) proteins"
        r"(?:\s+with\b.*?)?\s+domain-containing protein",
        value,
        flags=re.IGNORECASE,
    )
    if match:
        label = _singular_family_label(match.group(1).strip())
        actions.append("generalized_related_proteins")
        value = f"{label} protein" if " family" in label.casefold() else f"{label} family protein"

    match = re.fullmatch(
        r"found in\s+(.+?)\s+protein family protein",
        value,
        flags=re.IGNORECASE,
    )
    if match:
        actions.append("removed_explanatory_text")
        value = f"{match.group(1).strip()} family protein"

    match = re.fullmatch(
        r".+?,\s*contains?\s+(.+?)\s+domain-containing protein",
        value,
        flags=re.IGNORECASE,
    )
    if match:
        label = match.group(1).strip()
        label = re.sub(r"\brepeats$", "repeat", label, flags=re.IGNORECASE)
        actions.append("reduced_contains_description")
        value = (
            f"{label}-containing protein"
            if re.search(r"\b(?:domain|repeat)$", label, flags=re.IGNORECASE)
            else f"{label} domain-containing protein"
        )

    value = _replace(
        value,
        re.compile(r"\s*\(or\s+[^)]+\)", re.IGNORECASE),
        "",
        "removed_alternative_name",
        actions,
    )
    value = _replace(
        value,
        re.compile(r"\brepeats-containing\b", re.IGNORECASE),
        "repeat-containing",
        "singularized_repeat_name",
        actions,
    )
    value = _replace(
        value,
        re.compile(r"\bhomeodomain domain-containing protein\b", re.IGNORECASE),
        "homeodomain-containing protein",
        "removed_repeated_domain_word",
        actions,
    )
    value = _replace(
        value,
        re.compile(r"\bdomain domain-containing protein\b", re.IGNORECASE),
        "domain-containing protein",
        "removed_repeated_domain_word",
        actions,
    )
    value = _replace(
        value,
        re.compile(r"\bsuperfamily domain-containing protein\b", re.IGNORECASE),
        "superfamily protein",
        "normalized_superfamily_name",
        actions,
    )
    return value


def _model_identifier_fallback(
    value: str,
    *,
    source: str,
    evidence_id: str,
    actions: list[str],
) -> str:
    if source.casefold() not in {"cdd", "pfam"} or not evidence_id:
        return value
    if not re.search(
        r"\b(?:"
        r"function unknown|of unknown function|plant protein \d+|"
        r"uncharacteri[sz]ed conserved protein|unknown protein"
        r")\b",
        value,
        re.IGNORECASE,
    ):
        return value
    stable_id = evidence_id.split(".", 1)[0]
    if not re.fullmatch(r"(?:PF\d{5}|cd\d+|cl\d+)", stable_id, flags=re.IGNORECASE):
        return value
    actions.append("used_stable_model_identifier")
    return f"{stable_id} domain-containing protein"


def _generalize_eukaryotic_bacterial_name(
    value: str,
    *,
    source: str,
    context: ProductNameContext,
    actions: list[str],
) -> str:
    if context.domain.casefold() != "eukaryota" or source.casefold() not in {
        "cdd",
        "pfam",
        "uniref90",
    }:
        return value
    for pattern, replacement in _EUKARYOTE_BACTERIAL_GENERALIZATIONS:
        updated = pattern.sub(replacement, value)
        if updated != value:
            actions.append("generalized_cross_domain_name")
            return updated
    if re.fullmatch(r"(?:membrane\s+)?lipoprotein", value, flags=re.IGNORECASE):
        actions.extend(("rejected_cross_domain_name", "used_hypothetical_product"))
        return HYPOTHETICAL_PRODUCT
    return value


def _remove_subject_organism(
    value: str,
    subject_organism: str,
    actions: list[str],
) -> str:
    if not subject_organism:
        return value
    pattern = re.compile(rf"\b{re.escape(subject_organism)}\b", re.IGNORECASE)
    updated = pattern.sub("", value)
    if updated != value:
        actions.append("removed_organism_name")
    return updated


def _normalize_evolutionary_wording(value: str, actions: list[str]) -> str:
    subfamily = re.sub(
        r"\b(DnaJ)\s+homolog(?:ue)?\s+subfamily\b",
        r"\1 subfamily",
        value,
        flags=re.IGNORECASE,
    )
    if subfamily != value:
        actions.append("replaced_evolutionary_term")
        value = subfamily
    leading = re.fullmatch(
        r"protein\s+(.+?)(?:-|\s+)"
        r"(?:homolog|homologue|ortholog|orthologue|paralog|paralogue)"
        r"(?:\s+([A-Za-z0-9.-]+))?",
        value,
        flags=re.IGNORECASE,
    )
    if leading:
        designator = f" {leading.group(2)}" if leading.group(2) else ""
        actions.append("replaced_evolutionary_term")
        return f"{leading.group(1).strip()}-like protein{designator}"
    match = re.fullmatch(
        r"(.+?)(?:-|\s+)"
        r"(?:homolog|homologue|ortholog|orthologue|paralog|paralogue)"
        r"(?:\s+protein)?(?:\s+([A-Za-z0-9.-]+))?",
        value,
        flags=re.IGNORECASE,
    )
    if match:
        designator = f" {match.group(2)}" if match.group(2) else ""
        actions.append("replaced_evolutionary_term")
        return f"{match.group(1).strip()}-like protein{designator}"
    match = re.fullmatch(
        r"(.+?)\s+(?:homolog|homologue|ortholog|orthologue|paralog|paralogue)"
        r"\s+isoform\s+([A-Za-z0-9.-]+)",
        value,
        flags=re.IGNORECASE,
    )
    if match:
        actions.append("replaced_evolutionary_term")
        return f"{match.group(1).strip()}-like protein isoform {match.group(2)}"
    return value


def _normalize_initial_case(value: str, actions: list[str]) -> str:
    match = re.match(r"^([A-Z][A-Za-z0-9]*)(.*)$", value)
    if match is None:
        return value
    token = match.group(1)
    if (
        token in _CHEMICAL_SYMBOLS
        or token.isupper()
        or any(character.isupper() for character in token[1:])
        or (len(token) <= 6 and any(character.isdigit() for character in token))
    ):
        return value
    updated = token[0].lower() + token[1:] + match.group(2)
    actions.append("lowercased_initial")
    return updated


def _qualified_database_identifier(value: str) -> bool:
    return bool(
        re.search(
            r"\b(?:PF\d{5}|cd\d+|cl\d+|smart\d+)\s+"
            r"(?:domain-containing|family)\s+protein\b",
            value,
            re.IGNORECASE,
        )
    )


def audit_product_name(
    value: str,
    *,
    subject_organism: str = "",
) -> tuple[str, ...]:
    warnings = [
        code
        for code, pattern in _RESIDUAL_AUDITS
        if pattern.search(value)
        and not (
            code == "partial_or_fragment_term"
            and re.search(r"\btruncated hemoglobin\b", value, re.IGNORECASE)
        )
    ]
    if _DATABASE_IDENTIFIER.search(value) and not _qualified_database_identifier(value):
        warnings.append("database_identifier")
    if any(
        _COG_CATEGORY_WORDS.search(match.group(0)[1:-1])
        for match in re.finditer(r"\[[^]]+]", value)
    ):
        warnings.append("cog_category")
    if subject_organism and re.search(
        rf"\b{re.escape(subject_organism)}\b",
        value,
        re.IGNORECASE,
    ):
        warnings.append("organism_name")
    if ";" in value:
        warnings.append("semicolon")
    if len(value) > 120:
        warnings.append("long_product_name")
    return tuple(dict.fromkeys(warnings))


def standardize_product_name(
    product: str,
    *,
    source: str = "",
    evidence_id: str = "",
    subject_organism: str = "",
    context: ProductNameContext | None = None,
) -> ProductNameStandardization:
    actions: list[str] = []
    value = unicodedata.normalize("NFKC", product)
    unescaped = html.unescape(value)
    if unescaped != value:
        actions.append("decoded_html_entity")
    value = unescaped
    value = _replace(
        value,
        re.compile(r"^LOW QUALITY PROTEIN:\s*", re.IGNORECASE),
        "",
        "removed_quality_prefix",
        actions,
    )
    value = _replace(
        value,
        re.compile(r"^LOW protein:\s*", re.IGNORECASE),
        "",
        "removed_quality_prefix",
        actions,
    )
    value = value.replace("\\", "/").replace('"', "")
    value = re.sub(r"[\x00-\x20\x7f]+", " ", value).strip().strip(".;")
    value = _strip_cog_categories(value, actions)
    value = _remove_subject_organism(value, subject_organism, actions)
    value = _simplify_domain_description(value, actions)
    value = _generalize_eukaryotic_bacterial_name(
        value,
        source=source,
        context=context or ProductNameContext(),
        actions=actions,
    )
    value = _model_identifier_fallback(
        value,
        source=source,
        evidence_id=evidence_id,
        actions=actions,
    )
    value = _replace(
        value,
        _LOCALIZATION_SUFFIX,
        "",
        "removed_localization_suffix",
        actions,
    )
    value = _replace(
        value,
        re.compile(r"\s+\b(?:COG|KOG|FOG)\d+\b$", re.IGNORECASE),
        "",
        "removed_database_identifier",
        actions,
    )
    value = _replace(
        value,
        re.compile(r"\s+(?:LOC\d+|At\d+g\d+|Os\d+g\d+)$", re.IGNORECASE),
        "",
        "removed_locus_identifier",
        actions,
    )
    value = _replace(
        value,
        re.compile(r"\s*\([A-Z]_[A-Za-z]+_\d+\)", re.IGNORECASE),
        "",
        "removed_species_locus_identifier",
        actions,
    )
    value = _replace(
        value,
        re.compile(r"\bhelicase conserved C-terminal domain", re.IGNORECASE),
        "helicase C-terminal domain",
        "removed_conserved_term",
        actions,
    )
    value = _replace(
        value,
        re.compile(
            r"\bplants and prokaryotes conserved\s*\(PCC\)\s+domain",
            re.IGNORECASE,
        ),
        "PCC domain",
        "removed_conserved_term",
        actions,
    )
    value = _replace(
        value,
        re.compile(r"\s+conserved(?=\s+domain-containing protein)", re.IGNORECASE),
        "",
        "removed_conserved_term",
        actions,
    )
    value = _replace(
        value,
        re.compile(r"^(?:fragment|partial)\s+", re.IGNORECASE),
        "",
        "removed_partial_term",
        actions,
    )
    if not re.search(r"\btruncated hemoglobin\b", value, re.IGNORECASE):
        value = _replace(
            value,
            re.compile(r"^truncated\s+", re.IGNORECASE),
            "",
            "removed_partial_term",
            actions,
        )
    uncertainty = re.match(
        r"^(candidate|likely|possible|potential|predicted|probable)\s+",
        value,
        flags=re.IGNORECASE,
    )
    if uncertainty:
        if "domain-containing protein" in value.casefold():
            value = value[uncertainty.end() :]
            actions.append("removed_redundant_uncertainty")
        else:
            value = "putative " + value[uncertainty.end() :]
            actions.append("normalized_uncertainty")
    value = _replace(
        value,
        re.compile(r"\bmesoderm development candidate\b", re.IGNORECASE),
        "mesoderm development protein",
        "removed_candidate_term",
        actions,
    )
    value = _replace(
        value,
        re.compile(r"\bgolgin candidate\b", re.IGNORECASE),
        "golgin protein",
        "removed_candidate_term",
        actions,
    )
    value = _replace(
        value,
        re.compile(
            r"^protein CANDIDATE G PROTEIN-COUPLED RECEPTOR\s*",
            re.IGNORECASE,
        ),
        "G protein-coupled receptor-like protein ",
        "generalized_phenotype_name",
        actions,
    )
    value = _replace(
        value,
        re.compile(r"^novel plant SNARE\s+", re.IGNORECASE),
        "SNARE protein ",
        "generalized_taxon_specific_name",
        actions,
    )
    value = _replace(
        value,
        re.compile(r"^Orf\d+\s+protein$", re.IGNORECASE),
        HYPOTHETICAL_PRODUCT,
        "used_hypothetical_product",
        actions,
    )
    value = _replace(
        value,
        re.compile(r"^pseudo histidine-containing phosphotransfer protein", re.IGNORECASE),
        "inactive histidine-containing phosphotransfer protein",
        "normalized_inactive_protein_name",
        actions,
    )
    value = _replace(
        value,
        re.compile(r"\bsecreted\b\s*", re.IGNORECASE),
        "",
        "removed_localization_term",
        actions,
    )
    value = _replace(
        value,
        re.compile(r"_HOMOLOG(?:UE)?\b", re.IGNORECASE),
        "-like protein",
        "replaced_evolutionary_term",
        actions,
    )
    value = _replace(
        value,
        re.compile(
            r"^protein HOMOLOG OF MAMMALIAN LYST-INTERACTING PROTEIN\s*",
            re.IGNORECASE,
        ),
        "LYST-interacting protein-like protein ",
        "generalized_phenotype_name",
        actions,
    )
    value = _replace(
        value,
        re.compile(
            r"^spen paralogue and orthologue\s+SPOC\s+",
            re.IGNORECASE,
        ),
        "SPOC ",
        "removed_acronym_expansion",
        actions,
    )
    value = _replace(
        value,
        re.compile(
            r"^(.+?)\s*;\s*solute binding domain-containing protein$",
            re.IGNORECASE,
        ),
        "solute-binding domain-containing protein",
        "removed_alternative_name",
        actions,
    )
    value = _replace(
        value,
        re.compile(
            r"^trifunctional UDP-glucose 4,6-dehydratase/"
            r"UDP-4-keto-6-deoxy-D-glucose 3,5-epimerase/"
            r"UDP-4-keto-L-rhamnose-reductase (RHM\d+)$",
            re.IGNORECASE,
        ),
        r"multifunctional protein \1",
        "shortened_multifunctional_name",
        actions,
    )
    value = _normalize_evolutionary_wording(value, actions)
    value = _replace(
        value,
        re.compile(r"(?<=\w);(?=\w)"),
        "-",
        "replaced_semicolon",
        actions,
    )
    for pattern, replacement in _AMERICAN_SPELLING:
        value = _replace(value, pattern, replacement, "american_spelling", actions)
    value = re.sub(r"\s+", " ", value).strip().strip(",;.")
    if value.casefold() in _UNINFORMATIVE_PRODUCTS or not value:
        if value != HYPOTHETICAL_PRODUCT:
            actions.append("used_hypothetical_product")
        value = HYPOTHETICAL_PRODUCT
    value = _normalize_initial_case(value, actions)
    warnings = audit_product_name(value, subject_organism=subject_organism)
    return ProductNameStandardization(
        product=value,
        actions=tuple(dict.fromkeys(actions)),
        warnings=warnings,
    )


__all__ = [
    "HYPOTHETICAL_PRODUCT",
    "ProductNameContext",
    "ProductNameStandardization",
    "audit_product_name",
    "standardize_product_name",
]
