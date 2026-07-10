from __future__ import annotations

import re
from collections import Counter

from ..fasta import reverse_complement
from .gaps import GapRegions, append_position
from .models import AnnotationEntry, FeatureRecord
from .render import (
    render_cds_feature,
    render_rrna_feature,
    render_trna_feature,
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
    annotation = annotation_lookup[transcript_id]
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

    sub_features = parent_lookup.get(rna_feature.id)
    if not sub_features:
        return ""
    for sub_feature in sub_features:
        if sub_feature.type == "exon":
            count += 1
            position, joint_prefix, joint_suffix, out_gap_flag = append_position(
                sub_feature,
                count,
                position,
                gap_regions,
                rna_feature.strand,
            )
    joined_location = strand_prefix + joint_prefix + position + joint_suffix + strand_suffix
    out = render_rrna_feature(rna_feature.rna_type, joined_location)
    out += f"\t\t\tlocus_tag\t{locus_tag_prefix}{str(locus_tag_counter).zfill(9)}\n"
    out += f"\t\t\tnote\ttranscript_id:{rna_feature.name}\n"
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

    sub_features = parent_lookup.get(rna_feature.id)
    if not sub_features:
        return ""
    for sub_feature in sub_features:
        if sub_feature.type == "exon":
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
        product=rna_feature.name,
        anticodon=rna_feature.anticodon,
        note="",
    )
    if out_gap_flag:
        out += "\t\t\tartificial_location\tlow-quality sequence region\n"
    return out


def convert_contig_features(
    *,
    gene_lookup: dict[str, list[FeatureRecord]],
    parent_lookup: dict[str, list[FeatureRecord]],
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
    if not genes:
        return locus_tag_counter, ""
    chunks: list[str] = []
    for gene_feature in genes:
        locus_tag_counter += 100
        gene_children = parent_lookup.get(gene_feature.id)
        if not gene_children:
            continue
        for child in gene_children:
            if child.type in ("mRNA", "transcript"):
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
    return locus_tag_counter, "".join(chunks)
