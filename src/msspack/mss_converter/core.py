from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from Bio.Data import CodonTable

from ..fasta import iter_fasta
from ..utils import atomic_text_writer
from .features import convert_contig_features
from .gaps import detect_gap_regions
from .models import (
    FeatureRecord,
    build_gff_indexes,
    load_annotation_lookup,
    load_gff_features,
    load_protein_id_lookup,
)
from .render import render_source_feature


@dataclass(frozen=True)
class ConversionOptions:
    fasta_path: Path
    gff_path: Path
    annotation_path: Path
    output_path: Path
    locus_tag_prefix: str
    organism_name: str
    strain: str = ""
    mol_type: str = "genomic DNA"
    protein_id_path: Path | None = None
    linkage_evidence: str = "paired-ends"
    genetic_code: str = "1"
    infer_boundary: bool = False
    start_codons: tuple[str, ...] | None = None
    isolate: str = ""
    sex: str = ""
    country: str = ""
    collection_date: str = ""
    min_assembly_gap_size: int = 0
    gap_estimated_length: str = "known"
    feature_with_gap: str = "asis"
    minimum_intron_size_cutoff: int = 0


@dataclass(frozen=True)
class ContigSummary:
    contig_name: str
    event_counts: dict[str, int]


@dataclass(frozen=True)
class ConversionSummary:
    output_path: Path
    used_custom_locus_tags: bool
    contigs: list[ContigSummary]
    overall_counts: dict[str, int]


def build_gff_feature_indexes(
    rows: list[FeatureRecord] | list[object],
) -> tuple[dict[str, list[FeatureRecord]], dict[str, list[FeatureRecord]]]:
    return build_gff_indexes(rows)


def _get_stop_codons(genetic_code: str) -> set[str]:
    return {
        str(codon)
        for codon in CodonTable.unambiguous_dna_by_id[int(genetic_code)].stop_codons
    }


def _get_start_codons(genetic_code: str) -> set[str]:
    return {
        str(codon)
        for codon in CodonTable.unambiguous_dna_by_id[int(genetic_code)].start_codons
    }


def convert_gff_to_mss(options: ConversionOptions) -> ConversionSummary:
    annotation_lookup, has_custom_locus_tags = load_annotation_lookup(options.annotation_path)
    protein_lookup = load_protein_id_lookup(options.protein_id_path)
    gff_features = load_gff_features(options.gff_path)
    gene_lookup, parent_lookup = build_gff_feature_indexes(gff_features)
    overall_counts: Counter[str] = Counter()
    contig_summaries: list[ContigSummary] = []
    locus_tag_counter = 0
    start_codons = (
        {codon.upper() for codon in options.start_codons if codon}
        if options.start_codons is not None
        else _get_start_codons(options.genetic_code)
    )
    stop_codons = _get_stop_codons(options.genetic_code)

    with atomic_text_writer(options.output_path) as out_handle:
        for record in iter_fasta(options.fasta_path):
            sequence = record.sequence.upper()
            contig_counts: Counter[str] = Counter()
            chunks = [
                render_source_feature(
                    contig_name=record.id,
                    length=len(sequence),
                    organism_name=options.organism_name,
                    strain=options.strain,
                    mol_type=options.mol_type,
                    country=options.country,
                    isolate=options.isolate,
                    collection_date=options.collection_date,
                    sex=options.sex,
                )
            ]
            gap_text, gap_regions = detect_gap_regions(
                sequence,
                linkage_evidence=options.linkage_evidence,
                min_assembly_gap_size=options.min_assembly_gap_size,
                gap_estimated_length=options.gap_estimated_length,
            )
            if gap_text:
                chunks.append(gap_text)

            locus_tag_counter, cds_text = convert_contig_features(
                gene_lookup=gene_lookup,
                parent_lookup=parent_lookup,
                contig_name=record.id,
                locus_tag_prefix=options.locus_tag_prefix,
                locus_tag_counter=locus_tag_counter,
                annotation_lookup=annotation_lookup,
                protein_lookup=protein_lookup,
                gap_regions=gap_regions,
                contig_sequence=sequence,
                infer_boundary=options.infer_boundary,
                start_codons=start_codons,
                stop_codons=stop_codons,
                feature_with_gap=options.feature_with_gap,
                minimum_intron_size_cutoff=options.minimum_intron_size_cutoff,
                transl_table=options.genetic_code,
                event_counts=contig_counts,
            )
            if cds_text:
                chunks.append(cds_text)
            out_handle.write("".join(chunks))
            contig_summaries.append(
                ContigSummary(
                    contig_name=record.id,
                    event_counts=dict(contig_counts),
                )
            )
            overall_counts.update(contig_counts)

    return ConversionSummary(
        output_path=options.output_path,
        used_custom_locus_tags=has_custom_locus_tags,
        contigs=contig_summaries,
        overall_counts=dict(overall_counts),
    )
