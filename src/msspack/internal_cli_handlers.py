from __future__ import annotations

import argparse
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator

from .annotation_table import build_annotation_table
from .fasta_steps import remove_trailing_ns_fasta, write_mss_fasta
from .gap_normalization import normalize_gap_lengths
from .gff_adjustments import fix_gff_to_inframe
from .gff_cleanup import (
    drop_duplicate_coordinate_genes,
    fix_gff_semicolons_file,
    trim_gff_to_fasta_bounds,
)
from .padding_tools import pad_fasta, write_padding_log_for_gff, write_spliced_cds_fasta
from .pipeline_actions import (
    run_apply_padding_to_gff,
    run_mss_cds_to_misc,
    run_select_one_mrna,
    sort_gff,
)
from .utils import MSSPackError


def resolved_path(value: str | Path) -> Path:
    return Path(value).expanduser().resolve()


@contextmanager
def managed_log_path(path_value: str | None, default_name: str) -> Iterator[tuple[Path, bool]]:
    if path_value:
        yield resolved_path(path_value), True
        return
    with tempfile.TemporaryDirectory(prefix="msspack-internal-") as tmp_dir:
        yield Path(tmp_dir) / default_name, False


def print_paths(*paths: Path | None) -> None:
    for path in paths:
        if path is not None:
            print(path)


def remove_trailing_ns_handler(args: argparse.Namespace) -> int:
    input_path = resolved_path(args.input)
    output_path = resolved_path(args.output)
    with managed_log_path(args.log, "remove-trailing-ns.log") as (log_path, keep_log):
        remove_trailing_ns_fasta(
            input_path=input_path,
            output_path=output_path,
            log_path=log_path,
        )
    print_paths(output_path, log_path if keep_log else None)
    return 0


def fasta_to_mss_handler(args: argparse.Namespace) -> int:
    input_path = resolved_path(args.input)
    output_path = resolved_path(args.output)
    with managed_log_path(args.log, "fasta-to-mss.log") as (log_path, keep_log):
        write_mss_fasta(
            input_path=input_path,
            output_path=output_path,
            log_path=log_path,
        )
    print_paths(output_path, log_path if keep_log else None)
    return 0


def extract_cds_handler(args: argparse.Namespace) -> int:
    fasta_path = resolved_path(args.fasta)
    gff_path = resolved_path(args.gff)
    output_path = resolved_path(args.output)
    with managed_log_path(args.log, "extract-cds.log") as (log_path, keep_log):
        write_spliced_cds_fasta(
            fasta_path=fasta_path,
            gff_path=gff_path,
            output_path=output_path,
            log_path=log_path,
        )
    print_paths(output_path, log_path if keep_log else None)
    return 0


def pad_handler(args: argparse.Namespace) -> int:
    input_path = resolved_path(args.input)
    output_path = resolved_path(args.output)
    with managed_log_path(args.log, "pad.log") as (log_path, keep_log):
        pad_fasta(
            input_path=input_path,
            output_path=output_path,
            log_path=log_path,
            genetic_code=args.genetic_code,
            padchar=args.padchar,
        )
    print_paths(output_path, log_path if keep_log else None)
    return 0


def pad_from_gff_handler(args: argparse.Namespace) -> int:
    fasta_path = resolved_path(args.fasta)
    gff_path = resolved_path(args.gff)
    output_log = resolved_path(args.output_log)
    extract_log = resolved_path(args.extract_log) if args.extract_log else None
    write_padding_log_for_gff(
        fasta_path=fasta_path,
        gff_path=gff_path,
        log_path=output_log,
        genetic_code=args.genetic_code,
        padchar=args.padchar,
        extract_log_path=extract_log,
    )
    print_paths(output_log, extract_log)
    return 0


def gapjust_handler(args: argparse.Namespace) -> int:
    fasta_path = resolved_path(args.fasta)
    output_path = resolved_path(args.output)
    input_gff = resolved_path(args.input_gff) if args.input_gff else None
    output_gff = resolved_path(args.output_gff) if args.output_gff else None
    if (input_gff is None) != (output_gff is None):
        raise MSSPackError("Provide both --input-gff and --output-gff together")
    with managed_log_path(args.log, "gapjust.log") as (log_path, keep_log):
        normalize_gap_lengths(
            fasta_path=fasta_path,
            output_fasta_path=output_path,
            log_path=log_path,
            gap_len=args.gap_len,
            gap_just_min=args.gap_just_min,
            gap_just_max=args.gap_just_max,
            input_gff_path=input_gff,
            output_gff_path=output_gff,
        )
    print_paths(output_path, output_gff, log_path if keep_log else None)
    return 0


def sort_gff_handler(args: argparse.Namespace) -> int:
    input_path = resolved_path(args.input)
    output_path = resolved_path(args.output)
    with managed_log_path(args.log, "gff3sort.log") as (log_path, keep_log):
        sort_gff(
            input_path,
            output_path,
            command_override=None,
            log_path=log_path,
        )
    print_paths(output_path, log_path if keep_log else None)
    return 0


def trim_gff_handler(args: argparse.Namespace) -> int:
    gff_path = resolved_path(args.gff)
    fasta_path = resolved_path(args.fasta)
    output_path = resolved_path(args.output)
    with managed_log_path(args.log, "trim-gff.log") as (log_path, keep_log):
        trim_gff_to_fasta_bounds(
            gff_path=gff_path,
            fasta_path=fasta_path,
            output_path=output_path,
            log_path=log_path,
        )
    print_paths(output_path, log_path if keep_log else None)
    return 0


def select_one_mrna_handler(args: argparse.Namespace) -> int:
    input_path = resolved_path(args.input)
    output_path = resolved_path(args.output)
    with managed_log_path(args.log, "select-one-mrna.log") as (log_path, keep_log):
        run_select_one_mrna(
            input_path=input_path,
            output_path=output_path,
            log_path=log_path,
        )
    print_paths(output_path, log_path if keep_log else None)
    return 0


def drop_duplicate_coordinate_gene_handler(args: argparse.Namespace) -> int:
    input_path = resolved_path(args.input)
    output_path = resolved_path(args.output)
    with managed_log_path(args.log, "drop-duplicate-coordinate-gene.log") as (log_path, keep_log):
        drop_duplicate_coordinate_genes(
            input_path=input_path,
            output_path=output_path,
            log_path=log_path,
        )
    print_paths(output_path, log_path if keep_log else None)
    return 0


def fix_gff_semicolons_handler(args: argparse.Namespace) -> int:
    input_path = resolved_path(args.input)
    output_path = resolved_path(args.output)
    with managed_log_path(args.log, "fix-gff-semicolons.log") as (log_path, keep_log):
        fix_gff_semicolons_file(
            input_path=input_path,
            output_path=output_path,
            log_path=log_path,
        )
    print_paths(output_path, log_path if keep_log else None)
    return 0


def update_gff_to_inframe_handler(args: argparse.Namespace) -> int:
    input_path = resolved_path(args.input)
    output_path = resolved_path(args.output)
    with managed_log_path(args.log, "update-gff-to-inframe.log") as (log_path, keep_log):
        fix_gff_to_inframe(
            input_path=input_path,
            output_path=output_path,
            log_path=log_path,
        )
    print_paths(output_path, log_path if keep_log else None)
    return 0


def update_gff_with_padding_handler(args: argparse.Namespace) -> int:
    gff_path = resolved_path(args.gff)
    padding_log_path = resolved_path(args.padding_log)
    output_path = resolved_path(args.output)
    genes_with_stops = resolved_path(args.genes_with_stops)
    updated_genes = resolved_path(args.updated_genes)
    with managed_log_path(args.log, "update-gff-with-padding.log") as (log_path, keep_log):
        run_apply_padding_to_gff(
            gff_path=gff_path,
            padding_log_path=padding_log_path,
            output_path=output_path,
            genes_with_stops_path=genes_with_stops,
            updated_genes_path=updated_genes,
            command_log_path=log_path,
        )
    print_paths(output_path, genes_with_stops, updated_genes, log_path if keep_log else None)
    return 0


def gff2annotation_table_handler(args: argparse.Namespace) -> int:
    gff_path = resolved_path(args.gff)
    output_path = resolved_path(args.output)
    with managed_log_path(args.log, "gff2annotation-table.log") as (log_path, keep_log):
        build_annotation_table(
            gff_path=gff_path,
            output_path=output_path,
            locus_tag_prefix=args.locus_tag_prefix,
            log_path=log_path,
        )
    print_paths(output_path, log_path if keep_log else None)
    return 0


def mss_cds_to_misc_handler(args: argparse.Namespace) -> int:
    mss_input_path = resolved_path(args.mss_input)
    genes_input_path = resolved_path(args.genes_input)
    output_path = resolved_path(args.output)
    with managed_log_path(args.log, "mss-cds-to-misc.log") as (log_path, keep_log):
        run_mss_cds_to_misc(
            mss_input_path=mss_input_path,
            genes_input_path=genes_input_path,
            mss_output_path=output_path,
            log_path=log_path,
            locus_tag_prefix=args.locus_tag_prefix,
        )
    print_paths(output_path, log_path if keep_log else None)
    return 0
