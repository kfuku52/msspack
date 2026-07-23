from __future__ import annotations

import re
from collections import Counter
from dataclasses import replace

from ..fasta import reverse_complement
from ..gff_feature_policy import (
    CDS_TYPES,
    CODING_TRANSCRIPT_TYPES,
    CODON_TYPES,
    EXON_TYPES,
    RNA_FEATURE_TYPES,
    UTR_TYPES,
    canonical_insdc_feature,
    derived_ncRNA_class,
    derived_regulatory_class,
    is_pseudogenic,
    is_transcript_structure,
    pseudogene_type,
)
from .gaps import GapRegions, append_position
from .models import AnnotationEntry, FeatureRecord
from .render import (
    render_cds_feature,
    render_generic_feature,
    render_rrna_feature,
    render_trna_feature,
)

_CONTAINER_TYPES = frozenset({"chromosome", "contig", "scaffold", "supercontig", "region"})


def _attribute(feature: FeatureRecord, *keys: str) -> str:
    for key in keys:
        value = feature.attributes.get(key, "").strip()
        if value:
            return value
    return ""


def _clean_value(value: str) -> str:
    return value.replace("\t", " ").replace("\r", " ").replace("\n", " ").strip()


def _comma_values(value: str) -> list[str]:
    return [_clean_value(item) for item in value.split(",") if _clean_value(item)]


def _location_for_features(
    features: list[FeatureRecord],
    *,
    strand: str,
    gap_regions: GapRegions,
) -> tuple[str, bool]:
    if not features:
        return "", False
    count = 0
    position = ""
    joint_prefix = ""
    joint_suffix = ""
    out_gap_flag = False
    ordered = sorted(features, key=lambda feature: (feature.start, feature.end))
    for feature in ordered:
        count += 1
        position, joint_prefix, joint_suffix, feature_gap = append_position(
            feature,
            count,
            position,
            gap_regions,
            strand,
        )
        out_gap_flag = out_gap_flag or feature_gap
    strand_prefix = "complement(" if strand == "-" else ""
    strand_suffix = ")" if strand == "-" else ""
    return strand_prefix + joint_prefix + position + joint_suffix + strand_suffix, out_gap_flag


def _location_features(
    feature: FeatureRecord,
    parent_lookup: dict[str, list[FeatureRecord]],
) -> list[FeatureRecord]:
    children = parent_lookup.get(feature.id, []) if feature.id else []
    exons = [child for child in children if child.type in EXON_TYPES]
    if feature.type in CODING_TRANSCRIPT_TYPES and not exons:
        transcript_parts = [
            child for child in children if child.type in CDS_TYPES | UTR_TYPES
        ]
        if transcript_parts:
            return _merge_adjacent_features(transcript_parts)
    return exons or [feature]


def _merge_adjacent_features(features: list[FeatureRecord]) -> list[FeatureRecord]:
    """Merge overlapping or adjacent transcript parts into mature-RNA segments."""

    ordered = sorted(features, key=lambda feature: (feature.start, feature.end))
    if not ordered:
        return []
    merged = [ordered[0]]
    for feature in ordered[1:]:
        previous = merged[-1]
        if feature.start <= previous.end + 1:
            merged[-1] = replace(previous, end=max(previous.end, feature.end))
        else:
            merged.append(feature)
    return merged


def _covered_intervals(features: list[FeatureRecord]) -> tuple[tuple[int, int], ...]:
    return tuple(
        (feature.start, feature.end) for feature in _merge_adjacent_features(features)
    )


def _mrna_adds_submission_information(
    feature: FeatureRecord,
    parent_lookup: dict[str, list[FeatureRecord]],
    *,
    alternative_isoforms: bool,
) -> bool:
    """Return whether an mRNA conveys information beyond its CDS.

    A coding transcript whose explicit exon coverage is identical to its CDS
    coverage adds no mature-transcript information to a DDBJ feature table.
    Transcripts with UTR/non-coding exon sequence, alternative isoforms, or no
    CDS remain informative.
    """

    children = parent_lookup.get(feature.id, []) if feature.id else []
    cds_features = [child for child in children if child.type in CDS_TYPES]
    if not cds_features:
        return True
    if alternative_isoforms:
        return True
    if any(child.type in UTR_TYPES for child in children):
        return True
    exon_features = [child for child in children if child.type in EXON_TYPES]
    if not exon_features:
        # A bare GFF mRNA span includes introns and cannot safely describe the
        # mature transcript when no exon or UTR structure was supplied.
        return False
    return _covered_intervals(exon_features) != _covered_intervals(cds_features)


def _generic_qualifiers(
    feature: FeatureRecord,
    *,
    feature_key: str,
    locus_tag: str | None,
    product_override: str = "",
    pseudogene_override: str = "",
) -> list[tuple[str, str]]:
    qualifiers: list[tuple[str, str]] = []
    if locus_tag:
        qualifiers.append(("locus_tag", locus_tag))

    gene_name = _attribute(feature, "gene", "gene_name")
    if not gene_name and feature.type in {"gene", "pseudogene"}:
        gene_name = feature.name
    if gene_name:
        qualifiers.append(("gene", _clean_value(gene_name)))
    for synonym in _comma_values(_attribute(feature, "gene_synonym", "Alias")):
        qualifiers.append(("gene_synonym", synonym))

    product = product_override or _attribute(feature, "product", "Product")
    if not product and feature_key in {"rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA"}:
        product = feature.name or feature.rna_type
    if feature_key == "rRNA":
        product = product or "ribosomal RNA"
        if product.lower() != "ribosomal rna" and not product.lower().endswith("rrna"):
            product += " rRNA"
    elif feature_key == "tRNA":
        product = product or "tRNA"
    elif feature_key == "ncRNA":
        product = product or "non-coding RNA"
    # The INSDC definition permits /product on mRNA, but the current DDBJ
    # Parser usage matrix rejects it for new MSS submissions. Product names
    # remain on CDS and RNA product features.
    if product and feature_key != "mRNA":
        qualifiers.append(("product", _clean_value(product)))

    nc_rna_class = derived_ncRNA_class(feature.type, feature.attributes)
    if feature_key == "ncRNA":
        qualifiers.append(("ncRNA_class", _clean_value(nc_rna_class or "other")))

    regulatory_class = derived_regulatory_class(feature.type, feature.attributes)
    if feature_key == "regulatory":
        qualifiers.append(("regulatory_class", _clean_value(regulatory_class or "other")))

    if feature_key == "repeat_region":
        repeat_type = _attribute(feature, "rpt_type", "repeat_type")
        repeat_family = _attribute(feature, "rpt_family", "repeat_family")
        repeat_unit = _attribute(feature, "rpt_unit_seq", "repeat_unit")
        if repeat_type:
            qualifiers.append(("rpt_type", _clean_value(repeat_type)))
        if repeat_family:
            qualifiers.append(("rpt_family", _clean_value(repeat_family)))
        if repeat_unit:
            qualifiers.append(("rpt_unit_seq", _clean_value(repeat_unit)))

    if feature_key == "mobile_element":
        mobile_type = _attribute(feature, "mobile_element_type")
        if not mobile_type:
            label = feature.name or feature.type.replace("_", " ")
            mobile_type = f"other:{label}"
        qualifiers.append(("mobile_element_type", _clean_value(mobile_type)))

    if feature_key == "modified_base":
        modified_base = _attribute(feature, "mod_base")
        if modified_base:
            qualifiers.append(("mod_base", _clean_value(modified_base)))

    if feature_key == "gap":
        estimated_length = _attribute(feature, "estimated_length") or str(
            feature.end - feature.start + 1
        )
        qualifiers.append(("estimated_length", _clean_value(estimated_length)))

    if feature_key == "operon":
        operon_name = _attribute(feature, "operon") or feature.name or feature.id
        qualifiers.append(("operon", _clean_value(operon_name or "unnamed operon")))

    if feature_key in {"misc_binding", "protein_bind"}:
        bound_moiety = _attribute(feature, "bound_moiety") or feature.name or "unknown"
        qualifiers.append(("bound_moiety", _clean_value(bound_moiety)))

    anticodon = _attribute(feature, "anticodon")
    if feature_key == "tRNA" and anticodon:
        qualifiers.append(("anticodon", _clean_value(anticodon)))

    if pseudogene_override:
        qualifiers.append(("pseudogene", pseudogene_override))
    elif is_pseudogenic(feature.type, feature.attributes):
        qualifiers.append(("pseudogene", pseudogene_type(feature.attributes)))

    for db_xref in _comma_values(_attribute(feature, "Dbxref", "db_xref")):
        qualifiers.append(("db_xref", db_xref))

    note = _attribute(feature, "Note", "note", "description", "Description")
    if feature_key == "misc_feature":
        type_note = f"GFF3 feature type: {feature.type}"
        note = f"{note}; {type_note}" if note else type_note
    if note:
        qualifiers.append(("note", _clean_value(note)))
    if feature.id and feature_key in RNA_FEATURE_TYPES | {"mRNA"}:
        qualifiers.append(("note", f"transcript_id:{_clean_value(feature.id)}"))
    return qualifiers


def build_generic_feature_text(
    *,
    feature: FeatureRecord,
    parent_lookup: dict[str, list[FeatureRecord]],
    gap_regions: GapRegions,
    locus_tag: str | None = None,
    product_override: str = "",
    feature_key_override: str | None = None,
    pseudogene_override: str = "",
) -> tuple[str, bool]:
    feature_key = feature_key_override or canonical_insdc_feature(feature.type)
    if feature_key is None:
        return "", False
    if feature_key == "modified_base" and not _attribute(feature, "mod_base"):
        feature_key = "misc_feature"
    location, out_gap_flag = _location_for_features(
        _location_features(feature, parent_lookup),
        strand=feature.strand,
        gap_regions=gap_regions,
    )
    if not location:
        return "", False
    qualifiers = _generic_qualifiers(
        feature,
        feature_key=feature_key,
        locus_tag=locus_tag,
        product_override=product_override,
        pseudogene_override=pseudogene_override,
    )
    if out_gap_flag:
        qualifiers.append(("artificial_location", "low-quality sequence region"))
    return (
        render_generic_feature(
            feature_key=feature_key,
            position=location,
            qualifiers=qualifiers,
        ),
        feature_key == "misc_feature"
        and feature.type not in {"misc_feature", "gene", "pseudogene"},
    )


def _detect_incomplete_start(
    transcript_order_feature: FeatureRecord,
    count: int,
    strand_prefix: str,
    codon_start: int,
) -> tuple[bool, bool, int]:
    if count != 1:
        return False, False, codon_start
    codon_start = transcript_order_feature.phase + 1
    if codon_start == 1:
        return False, False, codon_start
    if strand_prefix:
        return False, True, codon_start
    return True, False, codon_start


def build_mrna_text(
    *,
    parent_lookup: dict[str, list[FeatureRecord]],
    rna_feature: FeatureRecord,
    locus_tag_prefix: str,
    locus_tag_counter: int,
    annotation_lookup: dict[str, AnnotationEntry],
    protein_lookup: dict[str, str] | None,
    gap_regions: GapRegions,
    contig_sequence: str,
    infer_boundary: bool,
    start_codons: set[str],
    stop_codons: set[str],
    feature_with_gap: str,
    minimum_intron_size_cutoff: int,
    transl_table: str,
    event_counts: Counter[str],
) -> str:
    count = 0
    incomplete_5 = False
    reverse_incomplete_5 = False
    incomplete_3 = False
    reverse_incomplete_3 = False
    artificial_location_flag = False
    codon_start = 1
    strand_prefix = ""
    strand_suffix = ""
    position = ""
    joint_prefix = ""
    joint_suffix = ""
    out_gap_flag = False

    strand = rna_feature.strand
    if strand == "-":
        strand_prefix = "complement("
        strand_suffix = ")"

    transcript_id = rna_feature.id
    annotation = annotation_lookup.get(
        transcript_id,
        AnnotationEntry(
            product_name=_attribute(rna_feature, "product", "Product") or "hypothetical protein",
            custom_locus_tag=None,
        ),
    )
    product_name = annotation.product_name
    custom_locus_tag = annotation.custom_locus_tag

    sub_features = parent_lookup.get(transcript_id)
    if not sub_features:
        return ""
    cds_features = [feature for feature in sub_features if feature.type == "CDS"]
    if not cds_features:
        return ""
    cds_by_position = sorted(cds_features, key=lambda feature: feature.start)
    cds_by_transcript = list(reversed(cds_by_position)) if strand == "-" else cds_by_position

    for genomic_feature, transcript_feature in zip(
        cds_by_position, cds_by_transcript, strict=True
    ):
        count += 1
        position, joint_prefix, joint_suffix, out_gap_flag = append_position(
            genomic_feature,
            count,
            position,
            gap_regions,
            strand,
        )
        incomplete_5_tmp, reverse_incomplete_5_tmp, codon_start = _detect_incomplete_start(
            transcript_feature,
            count,
            strand_prefix,
            codon_start,
        )
        incomplete_5 = incomplete_5 or incomplete_5_tmp
        reverse_incomplete_5 = reverse_incomplete_5 or reverse_incomplete_5_tmp

    if infer_boundary:
        spliced_parts = [contig_sequence[exon.start - 1 : exon.end] for exon in cds_by_position]
        spliced_cds = "".join(spliced_parts)
        if strand == "-":
            spliced_cds = reverse_complement(spliced_cds)
        start_codon_seq = spliced_cds[:3]
        stop_codon_seq = spliced_cds[-3:]
        if start_codon_seq not in start_codons:
            event_counts["start_codon_missing"] += 1
            if strand == "+":
                incomplete_5 = True
            else:
                reverse_incomplete_5 = True
        if stop_codon_seq not in stop_codons:
            event_counts["stop_codon_missing"] += 1
            if strand == "+":
                incomplete_3 = True
            else:
                reverse_incomplete_3 = True

    if incomplete_5 or reverse_incomplete_3:
        position = re.sub(r"^", "<", position)
    if reverse_incomplete_5 or incomplete_3:
        position = re.sub(r"([^.,]+$)", r">\1", position)

    joined_location = strand_prefix + joint_prefix + position + joint_suffix + strand_suffix

    intron_sizes: list[int] = []
    for end_start in position.split(".."):
        if "," in end_start:
            end_val, start_val = end_start.split(",")
            intron_size = int(start_val.strip("><")) - int(end_val.strip("><")) - 1
            intron_sizes.append(intron_size)
    out = render_cds_feature(
        joined_location=joined_location,
        locus_tag_prefix=locus_tag_prefix,
        locus_tag_counter=locus_tag_counter,
        mrna_id=transcript_id,
        product_name=product_name,
        custom_locus_tag=custom_locus_tag,
        transl_table=transl_table,
        codon_start=codon_start,
    )
    if intron_sizes and min(intron_sizes) < minimum_intron_size_cutoff:
        event_counts["small_introns"] += 1
        artificial_location_flag = True

    if out_gap_flag:
        if feature_with_gap == "asis":
            event_counts["gap_artificial_location"] += 1
            artificial_location_flag = True
        elif feature_with_gap == "misc_feature":
            event_counts["gap_misc_feature"] += 1
            out = re.sub(r"^\tCDS\t", "\tmisc_feature\t", out)
            out = re.sub(r"\n\t\t\tproduct\t.*\n", "\n", out)
            out = re.sub(r"\n\t\t\ttransl_table\t.*\n", "\n", out)
            out = re.sub(r"\n\t\t\tcodon_start\t.*\n", "\n", out)
            artificial_location_flag = False
    if artificial_location_flag:
        out += "\t\t\tartificial_location\tlow-quality sequence region\n"
    if protein_lookup:
        protein_id = protein_lookup.get(transcript_id)
        if protein_id:
            out += f"\t\t\tprotein_id\t{protein_id}\n"
    return out


def build_rrna_text(
    *,
    parent_lookup: dict[str, list[FeatureRecord]],
    rna_feature: FeatureRecord,
    locus_tag_prefix: str,
    locus_tag_counter: int,
    gap_regions: GapRegions,
) -> str:
    count = 0
    strand_prefix = ""
    strand_suffix = ""
    position = ""
    joint_prefix = ""
    joint_suffix = ""
    out_gap_flag = False
    if rna_feature.strand == "-":
        strand_prefix = "complement("
        strand_suffix = ")"

    sub_features = parent_lookup.get(rna_feature.id, [])
    location_features = [feature for feature in sub_features if feature.type == "exon"] or [
        rna_feature
    ]
    for sub_feature in location_features:
        count += 1
        position, joint_prefix, joint_suffix, out_gap_flag = append_position(
            sub_feature,
            count,
            position,
            gap_regions,
            rna_feature.strand,
        )
    joined_location = strand_prefix + joint_prefix + position + joint_suffix + strand_suffix
    rrna_name = (
        _attribute(rna_feature, "product", "Product")
        or rna_feature.rna_type
        or rna_feature.name
    )
    out = render_rrna_feature(rrna_name, joined_location)
    out += f"\t\t\tlocus_tag\t{locus_tag_prefix}{str(locus_tag_counter).zfill(9)}\n"
    out += f"\t\t\tnote\ttranscript_id:{rna_feature.id or rna_feature.name}\n"
    if out_gap_flag:
        out += "\t\t\tartificial_location\tlow-quality sequence region\n"
    return out


def build_trna_text(
    *,
    parent_lookup: dict[str, list[FeatureRecord]],
    rna_feature: FeatureRecord,
    locus_tag_prefix: str,
    locus_tag_counter: int,
    gap_regions: GapRegions,
) -> str:
    count = 0
    strand_prefix = ""
    strand_suffix = ""
    position = ""
    joint_prefix = ""
    joint_suffix = ""
    out_gap_flag = False
    if rna_feature.strand == "-":
        strand_prefix = "complement("
        strand_suffix = ")"

    sub_features = parent_lookup.get(rna_feature.id, [])
    location_features = [feature for feature in sub_features if feature.type == "exon"] or [
        rna_feature
    ]
    for sub_feature in location_features:
        count += 1
        position, joint_prefix, joint_suffix, out_gap_flag = append_position(
            sub_feature,
            count,
            position,
            gap_regions,
            rna_feature.strand,
        )
    joined_location = strand_prefix + joint_prefix + position + joint_suffix + strand_suffix
    out = render_trna_feature(
        position=joined_location,
        locus_tag_prefix=locus_tag_prefix,
        locus_tag_counter=locus_tag_counter,
        product=_attribute(rna_feature, "product", "Product") or rna_feature.name or "tRNA",
        anticodon=rna_feature.anticodon,
        note="",
    )
    if out_gap_flag:
        out += "\t\t\tartificial_location\tlow-quality sequence region\n"
    return out


def build_standalone_cds_text(
    *,
    cds_feature: FeatureRecord,
    locus_tag: str,
    annotation_lookup: dict[str, AnnotationEntry],
    gap_regions: GapRegions,
    transl_table: str,
) -> str:
    location, out_gap_flag = _location_for_features(
        [cds_feature],
        strand=cds_feature.strand,
        gap_regions=gap_regions,
    )
    annotation = annotation_lookup.get(cds_feature.id)
    product = (
        annotation.product_name
        if annotation is not None
        else _attribute(cds_feature, "product", "Product") or "hypothetical protein"
    )
    out = render_cds_feature(
        joined_location=location,
        locus_tag_prefix="",
        locus_tag_counter=0,
        mrna_id=cds_feature.parent or cds_feature.id,
        product_name=product,
        custom_locus_tag=locus_tag,
        transl_table=transl_table,
        codon_start=cds_feature.phase + 1,
    )
    if out_gap_flag:
        out += "\t\t\tartificial_location\tlow-quality sequence region\n"
    return out


def convert_contig_features(
    *,
    gene_lookup: dict[str, list[FeatureRecord]],
    parent_lookup: dict[str, list[FeatureRecord]],
    seq_lookup: dict[str, list[FeatureRecord]],
    contig_name: str,
    locus_tag_prefix: str,
    locus_tag_counter: int,
    annotation_lookup: dict[str, AnnotationEntry],
    protein_lookup: dict[str, str] | None,
    gap_regions: GapRegions,
    contig_sequence: str,
    infer_boundary: bool,
    start_codons: set[str],
    stop_codons: set[str],
    feature_with_gap: str,
    minimum_intron_size_cutoff: int,
    transl_table: str,
    event_counts: Counter[str],
) -> tuple[int, str]:
    genes = gene_lookup.get(contig_name, [])
    chunks: list[str] = []
    processed: set[int] = set()

    def current_locus_tag() -> str:
        return f"{locus_tag_prefix}{str(locus_tag_counter).zfill(9)}"

    def render_generic(
        feature: FeatureRecord,
        *,
        locus_tag: str | None,
        product_override: str = "",
        feature_key_override: str | None = None,
        pseudogene_override: str = "",
    ) -> None:
        text, used_misc_fallback = build_generic_feature_text(
            feature=feature,
            parent_lookup=parent_lookup,
            gap_regions=gap_regions,
            locus_tag=locus_tag,
            product_override=product_override,
            feature_key_override=feature_key_override,
            pseudogene_override=pseudogene_override,
        )
        if text:
            chunks.append(text)
            event_counts["non_cds_features"] += 1
            if used_misc_fallback:
                event_counts["unknown_features"] += 1

    def render_descendants(
        parent_id: str,
        *,
        locus_tag: str | None,
        emit_transcript_structure: bool = True,
    ) -> None:
        for descendant in parent_lookup.get(parent_id, []):
            if id(descendant) in processed:
                continue
            processed.add(id(descendant))
            if descendant.type in CDS_TYPES | CODON_TYPES:
                render_descendants(
                    descendant.id,
                    locus_tag=locus_tag,
                    emit_transcript_structure=emit_transcript_structure,
                )
                continue
            if (
                not emit_transcript_structure
                and is_transcript_structure(descendant.type)
            ):
                if descendant.id:
                    render_descendants(
                        descendant.id,
                        locus_tag=locus_tag,
                        emit_transcript_structure=False,
                    )
                continue
            render_generic(descendant, locus_tag=locus_tag)
            if descendant.id:
                render_descendants(
                    descendant.id,
                    locus_tag=locus_tag,
                    emit_transcript_structure=emit_transcript_structure,
                )

    for gene_feature in genes:
        processed.add(id(gene_feature))
        locus_tag_counter += 100
        locus_tag = current_locus_tag()
        gene_children = parent_lookup.get(gene_feature.id)
        if not gene_children:
            render_generic(
                gene_feature,
                locus_tag=locus_tag,
                feature_key_override="misc_feature",
            )
            continue
        gene_is_pseudogenic = is_pseudogenic(gene_feature.type, gene_feature.attributes)
        if gene_is_pseudogenic:
            render_generic(
                gene_feature,
                locus_tag=locus_tag,
                feature_key_override="misc_feature",
            )
        coding_transcripts = [
            child for child in gene_children if child.type in CODING_TRANSCRIPT_TYPES
        ]
        alternative_isoforms = len(coding_transcripts) > 1
        for child in gene_children:
            if id(child) in processed:
                continue
            processed.add(id(child))
            if child.type in CODING_TRANSCRIPT_TYPES:
                annotation = annotation_lookup.get(child.id)
                product_override = annotation.product_name if annotation else ""
                transcript_locus_tag = (
                    annotation.custom_locus_tag
                    if annotation is not None and annotation.custom_locus_tag
                    else locus_tag
                )
                if gene_is_pseudogenic:
                    render_generic(
                        child,
                        locus_tag=transcript_locus_tag,
                        product_override=product_override,
                        feature_key_override="misc_RNA",
                        pseudogene_override=pseudogene_type(gene_feature.attributes),
                    )
                    for descendant in parent_lookup.get(child.id, []):
                        processed.add(id(descendant))
                    continue
                emit_mrna = _mrna_adds_submission_information(
                    child,
                    parent_lookup,
                    alternative_isoforms=alternative_isoforms,
                )
                if emit_mrna:
                    render_generic(
                        child,
                        locus_tag=transcript_locus_tag,
                        product_override=product_override,
                    )
                render_descendants(
                    child.id,
                    locus_tag=transcript_locus_tag,
                    emit_transcript_structure=emit_mrna,
                )
                if any(
                    descendant.type == "CDS"
                    for descendant in parent_lookup.get(child.id, [])
                ):
                    chunks.append(
                        build_mrna_text(
                            parent_lookup=parent_lookup,
                            rna_feature=child,
                            locus_tag_prefix=locus_tag_prefix,
                            locus_tag_counter=locus_tag_counter,
                            annotation_lookup=annotation_lookup,
                            protein_lookup=protein_lookup,
                            gap_regions=gap_regions,
                            contig_sequence=contig_sequence,
                            infer_boundary=infer_boundary,
                            start_codons=start_codons,
                            stop_codons=stop_codons,
                            feature_with_gap=feature_with_gap,
                            minimum_intron_size_cutoff=minimum_intron_size_cutoff,
                            transl_table=transl_table,
                            event_counts=event_counts,
                        )
                    )
            elif child.type == "CDS":
                chunks.append(
                    build_standalone_cds_text(
                        cds_feature=child,
                        locus_tag=locus_tag,
                        annotation_lookup=annotation_lookup,
                        gap_regions=gap_regions,
                        transl_table=transl_table,
                    )
                )
                render_descendants(child.id, locus_tag=locus_tag)
            elif child.type == "rRNA":
                chunks.append(
                    build_rrna_text(
                        parent_lookup=parent_lookup,
                        rna_feature=child,
                        locus_tag_prefix=locus_tag_prefix,
                        locus_tag_counter=locus_tag_counter,
                        gap_regions=gap_regions,
                    )
                )
                event_counts["non_cds_features"] += 1
                render_descendants(child.id, locus_tag=locus_tag)
            elif child.type == "tRNA":
                chunks.append(
                    build_trna_text(
                        parent_lookup=parent_lookup,
                        rna_feature=child,
                        locus_tag_prefix=locus_tag_prefix,
                        locus_tag_counter=locus_tag_counter,
                        gap_regions=gap_regions,
                    )
                )
                event_counts["non_cds_features"] += 1
                render_descendants(child.id, locus_tag=locus_tag)
            elif child.type not in CODON_TYPES:
                render_generic(child, locus_tag=locus_tag)
                render_descendants(child.id, locus_tag=locus_tag)

    for feature in seq_lookup.get(contig_name, []):
        if id(feature) in processed:
            continue
        processed.add(id(feature))
        if feature.type in _CONTAINER_TYPES or feature.type in CODON_TYPES:
            continue
        if feature.type == "gene":
            # Gene rows were handled above. This branch is only for malformed
            # duplicate indexes and deliberately avoids duplicate MSS features.
            continue
        if feature.type == "CDS":
            locus_tag_counter += 100
            chunks.append(
                build_standalone_cds_text(
                    cds_feature=feature,
                    locus_tag=current_locus_tag(),
                    annotation_lookup=annotation_lookup,
                    gap_regions=gap_regions,
                    transl_table=transl_table,
                )
            )
        elif feature.type in CODING_TRANSCRIPT_TYPES:
            locus_tag_counter += 100
            annotation = annotation_lookup.get(feature.id)
            transcript_locus_tag = (
                annotation.custom_locus_tag
                if annotation is not None and annotation.custom_locus_tag
                else current_locus_tag()
            )
            emit_mrna = _mrna_adds_submission_information(
                feature,
                parent_lookup,
                alternative_isoforms=False,
            )
            if emit_mrna:
                render_generic(
                    feature,
                    locus_tag=transcript_locus_tag,
                    product_override=annotation.product_name if annotation else "",
                )
            render_descendants(
                feature.id,
                locus_tag=transcript_locus_tag,
                emit_transcript_structure=emit_mrna,
            )
            if any(
                descendant.type == "CDS"
                for descendant in parent_lookup.get(feature.id, [])
            ):
                chunks.append(
                    build_mrna_text(
                        parent_lookup=parent_lookup,
                        rna_feature=feature,
                        locus_tag_prefix=locus_tag_prefix,
                        locus_tag_counter=locus_tag_counter,
                        annotation_lookup=annotation_lookup,
                        protein_lookup=protein_lookup,
                        gap_regions=gap_regions,
                        contig_sequence=contig_sequence,
                        infer_boundary=infer_boundary,
                        start_codons=start_codons,
                        stop_codons=stop_codons,
                        feature_with_gap=feature_with_gap,
                        minimum_intron_size_cutoff=minimum_intron_size_cutoff,
                        transl_table=transl_table,
                        event_counts=event_counts,
                    )
                )
        else:
            needs_locus = feature.type in RNA_FEATURE_TYPES or feature.type == "pseudogene"
            if needs_locus:
                locus_tag_counter += 100
            render_generic(
                feature,
                locus_tag=current_locus_tag() if needs_locus else None,
            )
        if feature.id:
            render_descendants(
                feature.id,
                locus_tag=(
                    current_locus_tag()
                    if feature.type in RNA_FEATURE_TYPES or feature.type == "pseudogene"
                    else None
                ),
            )
    return locus_tag_counter, "".join(chunks)
