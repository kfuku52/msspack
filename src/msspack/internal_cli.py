from __future__ import annotations

import argparse
from typing import Callable

from .internal_cli_handlers import (
    drop_duplicate_coordinate_gene_handler,
    extract_cds_handler,
    fasta_to_mss_handler,
    fix_gff_semicolons_handler,
    gapjust_handler,
    gff2annotation_table_handler,
    mss_cds_to_misc_handler,
    pad_from_gff_handler,
    pad_handler,
    remove_trailing_ns_handler,
    select_one_mrna_handler,
    sort_gff_handler,
    trim_gff_handler,
    update_gff_to_inframe_handler,
    update_gff_with_padding_handler,
)
from .utils import MSSPackError


def _add_command(
    subparsers: argparse._SubParsersAction[argparse.ArgumentParser],
    name: str,
    *,
    help_text: str,
    configure: Callable[[argparse.ArgumentParser], object],
    handler: Callable[[argparse.Namespace], int],
) -> None:
    parser = subparsers.add_parser(name, help=help_text)
    configure(parser)
    parser.set_defaults(internal_handler=handler)


def add_internal_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser("internal", help="run internal pipeline steps directly")
    internal_subparsers = parser.add_subparsers(dest="internal_command", required=True)

    _add_command(
        internal_subparsers,
        "remove-trailing-ns",
        help_text="trim terminal N/n runs from FASTA records",
        configure=lambda p: (
            p.add_argument("--input", required=True, help="input FASTA"),
            p.add_argument("--output", required=True, help="output FASTA"),
            p.add_argument("--log", help="optional log file"),
        ),
        handler=remove_trailing_ns_handler,
    )
    _add_command(
        internal_subparsers,
        "fasta-to-mss",
        help_text="write DDBJ MSS FASTA with // separators",
        configure=lambda p: (
            p.add_argument("--input", required=True, help="input FASTA"),
            p.add_argument("--output", required=True, help="output MSS FASTA"),
            p.add_argument("--log", help="optional log file"),
        ),
        handler=fasta_to_mss_handler,
    )
    _add_command(
        internal_subparsers,
        "extract-cds",
        help_text="extract spliced CDS FASTA from genome FASTA + GFF",
        configure=lambda p: (
            p.add_argument("--fasta", required=True, help="input genome FASTA"),
            p.add_argument("--gff", required=True, help="input GFF"),
            p.add_argument("--output", required=True, help="output CDS FASTA"),
            p.add_argument("--log", help="optional log file"),
        ),
        handler=extract_cds_handler,
    )
    _add_command(
        internal_subparsers,
        "pad",
        help_text="pad CDS FASTA entries to stay in-frame",
        configure=lambda p: (
            p.add_argument("--input", required=True, help="input FASTA"),
            p.add_argument("--output", required=True, help="output FASTA"),
            p.add_argument("--genetic-code", default="1", help="NCBI translation table"),
            p.add_argument("--padchar", default="N", help="padding character"),
            p.add_argument("--log", help="optional log file"),
        ),
        handler=pad_handler,
    )
    _add_command(
        internal_subparsers,
        "pad-from-gff",
        help_text="write a padding summary log directly from genome FASTA + GFF",
        configure=lambda p: (
            p.add_argument("--fasta", required=True, help="input genome FASTA"),
            p.add_argument("--gff", required=True, help="input GFF"),
            p.add_argument("--output-log", required=True, help="output padding log"),
            p.add_argument("--extract-log", help="optional extract-cds log"),
            p.add_argument("--genetic-code", default="1", help="NCBI translation table"),
            p.add_argument("--padchar", default="N", help="padding character"),
        ),
        handler=pad_from_gff_handler,
    )
    _add_command(
        internal_subparsers,
        "gapjust",
        help_text="normalize assembly gap lengths in FASTA and optional GFF",
        configure=lambda p: (
            p.add_argument("--fasta", required=True, help="input FASTA"),
            p.add_argument("--output", required=True, help="output FASTA"),
            p.add_argument("--gap-len", required=True, type=int, help="normalized gap length"),
            p.add_argument("--gap-just-min", type=int, help="minimum gap length to normalize"),
            p.add_argument("--gap-just-max", type=int, help="maximum gap length to normalize"),
            p.add_argument("--input-gff", help="optional input GFF"),
            p.add_argument("--output-gff", help="optional output GFF"),
            p.add_argument("--log", help="optional log file"),
        ),
        handler=gapjust_handler,
    )
    _add_command(
        internal_subparsers,
        "gff3sort",
        help_text="sort GFF with the built-in precise sorter",
        configure=lambda p: (
            p.add_argument("--input", required=True, help="input GFF"),
            p.add_argument("--output", required=True, help="output GFF"),
            p.add_argument("--log", help="optional log file"),
        ),
        handler=sort_gff_handler,
    )
    _add_command(
        internal_subparsers,
        "trim-gff-coordinates",
        help_text="clip GFF coordinates to FASTA bounds",
        configure=lambda p: (
            p.add_argument("--gff", required=True, help="input GFF"),
            p.add_argument("--fasta", required=True, help="reference FASTA"),
            p.add_argument("--output", required=True, help="output GFF"),
            p.add_argument("--log", help="optional log file"),
        ),
        handler=trim_gff_handler,
    )
    _add_command(
        internal_subparsers,
        "select-one-mrna",
        help_text="keep one transcript per gene using the current tie-break rules",
        configure=lambda p: (
            p.add_argument("--input", required=True, help="input GFF"),
            p.add_argument("--output", required=True, help="output GFF"),
            p.add_argument("--log", help="optional log file"),
        ),
        handler=select_one_mrna_handler,
    )
    _add_command(
        internal_subparsers,
        "drop-duplicate-coordinate-gene",
        help_text="remove genes that duplicate an earlier gene coordinate block",
        configure=lambda p: (
            p.add_argument("--input", required=True, help="input GFF"),
            p.add_argument("--output", required=True, help="output GFF"),
            p.add_argument("--log", help="optional log file"),
        ),
        handler=drop_duplicate_coordinate_gene_handler,
    )
    _add_command(
        internal_subparsers,
        "fix-gff-semicolons",
        help_text="repair malformed semicolons in GFF attributes",
        configure=lambda p: (
            p.add_argument("--input", required=True, help="input GFF"),
            p.add_argument("--output", required=True, help="output GFF"),
            p.add_argument("--log", help="optional log file"),
        ),
        handler=fix_gff_semicolons_handler,
    )
    _add_command(
        internal_subparsers,
        "update-gff-to-inframe",
        help_text="trim CDS models so the translated CDS is in-frame",
        configure=lambda p: (
            p.add_argument("--input", required=True, help="input GFF"),
            p.add_argument("--output", required=True, help="output GFF"),
            p.add_argument("--log", help="optional log file"),
        ),
        handler=update_gff_to_inframe_handler,
    )
    _add_command(
        internal_subparsers,
        "update-gff-with-padding",
        help_text="apply padding decisions back onto GFF features",
        configure=lambda p: (
            p.add_argument("--gff", required=True, help="input GFF"),
            p.add_argument("--padding-log", required=True, help="input padding log"),
            p.add_argument("--output", required=True, help="output GFF"),
            p.add_argument("--genes-with-stops", required=True, help="output genes-with-stops list"),
            p.add_argument("--updated-genes", required=True, help="output updated-genes list"),
            p.add_argument("--log", help="optional log file"),
        ),
        handler=update_gff_with_padding_handler,
    )
    _add_command(
        internal_subparsers,
        "gff2annotation-table",
        help_text="build the annotation TSV used by the MSS converter",
        configure=lambda p: (
            p.add_argument("--gff", required=True, help="input GFF"),
            p.add_argument("--output", required=True, help="output TSV"),
            p.add_argument("--locus-tag-prefix", default="", help="optional locus_tag prefix"),
            p.add_argument("--log", help="optional log file"),
        ),
        handler=gff2annotation_table_handler,
    )
    _add_command(
        internal_subparsers,
        "mss-cds-to-misc",
        help_text="convert selected CDS features to misc_feature in MSS output",
        configure=lambda p: (
            p.add_argument("--mss-input", required=True, help="input MSS annotation file"),
            p.add_argument("--genes-input", required=True, help="genes to convert"),
            p.add_argument("--output", required=True, help="output MSS annotation file"),
            p.add_argument("--log", help="optional log file"),
        ),
        handler=mss_cds_to_misc_handler,
    )


def handle_internal(args: argparse.Namespace) -> int:
    handler = getattr(args, "internal_handler", None)
    if handler is None:
        raise MSSPackError("Unsupported internal command")
    return handler(args)
