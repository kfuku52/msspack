from __future__ import annotations

from collections.abc import Mapping

# GFF3 feature roles used by both coordinate adjustment and MSS rendering.
CODING_TRANSCRIPT_TYPES = frozenset({"mRNA", "transcript"})
EXON_TYPES = frozenset({"exon", "pseudogenic_exon"})
CDS_TYPES = frozenset({"CDS"})
START_CODON_TYPES = frozenset({"start_codon"})
STOP_CODON_TYPES = frozenset({"stop_codon"})
CODON_TYPES = START_CODON_TYPES | STOP_CODON_TYPES
INTRON_TYPES = frozenset({"intron"})
FIVE_PRIME_UTR_TYPES = frozenset({"five_prime_UTR", "5_prime_UTR", "5'UTR"})
THREE_PRIME_UTR_TYPES = frozenset({"three_prime_UTR", "3_prime_UTR", "3'UTR"})
UTR_TYPES = FIVE_PRIME_UTR_TYPES | THREE_PRIME_UTR_TYPES

# These rows describe the structure of a parent transcript.  They are kept and
# synchronized in GFF3, but start/stop codons are not INSDC feature keys.
TRANSCRIPT_STRUCTURE_TYPES = EXON_TYPES | CDS_TYPES | CODON_TYPES | INTRON_TYPES | UTR_TYPES

RNA_TYPE_TO_CLASS: Mapping[str, str] = {
    "miRNA": "miRNA",
    "microRNA": "miRNA",
    "siRNA": "siRNA",
    "scRNA": "scRNA",
    "snRNA": "snRNA",
    "snoRNA": "snoRNA",
    "antisense_RNA": "antisense_RNA",
    "lncRNA": "lncRNA",
    "lnc_RNA": "lncRNA",
    "long_noncoding_RNA": "lncRNA",
    "guide_RNA": "guide_RNA",
    "ribozyme": "ribozyme",
    "RNase_P_RNA": "RNase_P_RNA",
    "RNase_MRP_RNA": "RNase_MRP_RNA",
    "telomerase_RNA": "telomerase_RNA",
    "vault_RNA": "vault_RNA",
    "Y_RNA": "Y_RNA",
}

RNA_FEATURE_TYPES = frozenset(
    {
        "rRNA",
        "tRNA",
        "tmRNA",
        "ncRNA",
        "misc_RNA",
        "precursor_RNA",
        "primary_transcript",
        "pseudogenic_transcript",
        *RNA_TYPE_TO_CLASS,
    }
)

REGULATORY_TYPE_TO_CLASS: Mapping[str, str] = {
    "attenuator": "attenuator",
    "CAAT_signal": "CAAT_signal",
    "enhancer": "enhancer",
    "GC_signal": "GC_signal",
    "-35_signal": "minus_35_signal",
    "minus_35_signal": "minus_35_signal",
    "-10_signal": "minus_10_signal",
    "minus_10_signal": "minus_10_signal",
    "polyA_signal": "polyA_signal_sequence",
    "promoter": "promoter",
    "RBS": "ribosome_binding_site",
    "ribosome_entry_site": "ribosome_binding_site",
    "TATA_signal": "TATA_box",
    "TATA_box": "TATA_box",
    "terminator": "terminator",
    "response_element": "response_element",
    "insulator": "other",
    "silencer": "other",
    "regulatory_region": "other",
}

REPEAT_TYPE_ALIASES = frozenset(
    {
        "repeat",
        "repeat_unit",
        "tandem_repeat",
        "microsatellite",
        "minisatellite",
        "satellite_DNA",
        "low_complexity_region",
    }
)

MOBILE_ELEMENT_ALIASES = frozenset(
    {
        "transposable_element",
        "transposon",
        "retrotransposon",
        "LINE",
        "SINE",
        "MITE",
        "insertion_sequence",
        "integron",
    }
)

# INSDC feature keys accepted or recommended by DDBJ.  CDS is rendered through
# its transcript model so that joined locations, phases and products are handled
# consistently.
DIRECT_INSDC_FEATURE_TYPES = frozenset(
    {
        "C_region",
        "D-loop",
        "D_segment",
        "J_segment",
        "V_region",
        "V_segment",
        "centromere",
        "exon",
        "gap",
        "intron",
        "mat_peptide",
        "misc_binding",
        "misc_difference",
        "misc_feature",
        "misc_RNA",
        "misc_structure",
        "mobile_element",
        "modified_base",
        "mRNA",
        "ncRNA",
        "operon",
        "oriT",
        "polyA_site",
        "precursor_RNA",
        "prim_transcript",
        "primer_bind",
        "propeptide",
        "protein_bind",
        "regulatory",
        "repeat_region",
        "rep_origin",
        "rRNA",
        "sig_peptide",
        "stem_loop",
        "telomere",
        "tmRNA",
        "transit_peptide",
        "tRNA",
        "unsure",
        "variation",
        "3'UTR",
        "5'UTR",
    }
)


def canonical_insdc_feature(feature_type: str) -> str | None:
    """Map a GFF3/SO feature type to an INSDC feature key.

    ``None`` means that the row is structural metadata or a container and
    should not be emitted independently.
    """

    if feature_type in START_CODON_TYPES | STOP_CODON_TYPES | CDS_TYPES:
        return None
    if feature_type in {"source", "assembly_gap"}:
        # Both are derived directly from the submitted FASTA.
        return None
    if feature_type in FIVE_PRIME_UTR_TYPES:
        return "5'UTR"
    if feature_type in THREE_PRIME_UTR_TYPES:
        return "3'UTR"
    if feature_type == "transcript":
        return "mRNA"
    if feature_type == "primary_transcript":
        return "prim_transcript"
    if feature_type == "pseudogenic_transcript":
        return "misc_RNA"
    if feature_type == "pseudogenic_exon":
        return "exon"
    if feature_type in RNA_TYPE_TO_CLASS:
        return "ncRNA"
    if feature_type in REGULATORY_TYPE_TO_CLASS:
        return "regulatory"
    if feature_type in REPEAT_TYPE_ALIASES:
        return "repeat_region"
    if feature_type in MOBILE_ELEMENT_ALIASES:
        return "mobile_element"
    if feature_type == "pseudogene":
        # DDBJ's current MSS Parser forbids the INSDC gene key for new
        # submissions.  misc_feature accepts locus_tag, gene and pseudogene.
        return "misc_feature"
    if feature_type in DIRECT_INSDC_FEATURE_TYPES:
        return feature_type
    if feature_type == "gene":
        return None
    return "misc_feature"


def derived_ncRNA_class(feature_type: str, attributes: Mapping[str, str]) -> str:
    return (
        attributes.get("ncRNA_class")
        or attributes.get("ncrna_class")
        or RNA_TYPE_TO_CLASS.get(feature_type, "")
    )


def derived_regulatory_class(feature_type: str, attributes: Mapping[str, str]) -> str:
    return (
        attributes.get("regulatory_class")
        or attributes.get("regulatory_region")
        or REGULATORY_TYPE_TO_CLASS.get(feature_type, "")
    )


def is_pseudogenic(feature_type: str, attributes: Mapping[str, str]) -> bool:
    pseudo_value = attributes.get("pseudo", "").strip().lower()
    return (
        feature_type.startswith("pseudogen")
        or "pseudogene" in attributes
        or pseudo_value in {"1", "true", "yes"}
    )


def pseudogene_type(attributes: Mapping[str, str]) -> str:
    value = attributes.get("pseudogene", "").strip().lower()
    valid = {"processed", "unprocessed", "unitary", "allelic", "unknown"}
    return value if value in valid else "unknown"


def is_transcript_structure(feature_type: str) -> bool:
    return feature_type in TRANSCRIPT_STRUCTURE_TYPES
