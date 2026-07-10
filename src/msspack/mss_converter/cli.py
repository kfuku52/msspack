from __future__ import annotations

import argparse
from pathlib import Path

from .core import ConversionOptions, convert_gff_to_mss
from .render import format_event_summary


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", required=True)
    parser.add_argument("-g", "--gff", required=True)
    parser.add_argument("-a", "--ann", required=True)
    parser.add_argument("-l", "--loc", required=True)
    parser.add_argument("-n", "--nam", required=True)
    parser.add_argument("-s", "--stn", default="")
    parser.add_argument("-o", "--out", required=True)
    parser.add_argument("-m", "--mol", default="genomic DNA")
    parser.add_argument("-p", "--pid")
    parser.add_argument("-t", "--gty", default="paired-ends")
    parser.add_argument("-c", "--gct", default="1")
    parser.add_argument("--ifc", default=False, action="store_true")
    parser.add_argument("--stc")
    parser.add_argument("--iso", default="")
    parser.add_argument("--sex", default="")
    parser.add_argument("--cou", default="")
    parser.add_argument("--cod", default="")
    parser.add_argument("--mag", default=0, type=int)
    parser.add_argument("--gel", default="known", choices=["known", "unknown"])
    parser.add_argument("--fwg", default="asis", choices=["asis", "misc_feature"])
    parser.add_argument("--mis", default=0, type=int)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    options = ConversionOptions(
        fasta_path=Path(args.fasta).expanduser().resolve(),
        gff_path=Path(args.gff).expanduser().resolve(),
        annotation_path=Path(args.ann).expanduser().resolve(),
        output_path=Path(args.out).expanduser().resolve(),
        locus_tag_prefix=args.loc,
        organism_name=args.nam,
        strain=args.stn,
        mol_type=args.mol,
        protein_id_path=Path(args.pid).expanduser().resolve() if args.pid else None,
        linkage_evidence=args.gty,
        genetic_code=args.gct,
        infer_boundary=bool(args.ifc),
        start_codons=tuple(args.stc.split(",")) if args.stc else None,
        isolate=args.iso,
        sex=args.sex,
        country=args.cou,
        collection_date=args.cod,
        min_assembly_gap_size=args.mag,
        gap_estimated_length=args.gel,
        feature_with_gap=args.fwg,
        minimum_intron_size_cutoff=args.mis,
    )
    summary = convert_gff_to_mss(options)
    if summary.used_custom_locus_tags:
        print(
            'The "Locus_tag" column exists in the annotation file. '
            "User-provided custom locus_tag values will be used."
        )
    else:
        print(
            'The "Locus_tag" column does not exist in the annotation file. '
            f"locus_tag values will be generated with the provided prefix: {args.loc}."
        )
    for contig in summary.contigs:
        print(format_event_summary(contig.event_counts, f"Finished {contig.contig_name}"))
    print(format_event_summary(summary.overall_counts, "Overall summary"))
    return 0
