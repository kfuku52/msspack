from __future__ import annotations

import gzip
import re
from datetime import datetime
from pathlib import Path

from .config import MSSPackConfig
from .fasta import iter_fasta_handle
from .fasta_steps import remove_trailing_ns_fasta, write_mss_fasta
from .gap_normalization import normalize_gap_lengths
from .gff import sort_gff_file_precise
from .gff_adjustments import apply_padding_to_gff, fix_gff_to_inframe
from .mss_converter import ConversionOptions, convert_gff_to_mss, format_event_summary
from .mss_postprocess import convert_cds_features_to_misc
from .step_logging import (
    count_fasta_records,
    count_gff_feature_records,
    count_gff_gene_records,
    count_mss_feature_blocks,
    count_nonempty_lines,
    count_reordered_feature_lines,
    count_tabular_rows,
    write_step_log,
    write_step_metrics,
)
from .submission_render import render_final_annotation
from .transcript_selection import select_one_mrna_per_gene
from .utils import MSSPackError, copy_or_decompress, link_or_copy, write_text


def run_gapjust(
    *,
    input_fasta: Path,
    output_fasta: Path,
    log_path: Path,
    gap_len: int,
    gap_just_min: int,
    gap_just_max: int,
    input_gff: Path | None = None,
    output_gff: Path | None = None,
    metrics_path: Path | None = None,
) -> None:
    normalize_gap_lengths(
        fasta_path=input_fasta,
        output_fasta_path=output_fasta,
        log_path=log_path,
        gap_len=gap_len,
        gap_just_min=gap_just_min,
        gap_just_max=gap_just_max,
        input_gff_path=input_gff,
        output_gff_path=output_gff,
        metrics_path=metrics_path,
    )


def run_remove_trailing_ns(
    *,
    input_path: Path,
    output_path: Path,
    log_path: Path,
    metrics_path: Path | None = None,
) -> None:
    remove_trailing_ns_fasta(
        input_path=input_path,
        output_path=output_path,
        log_path=log_path,
        metrics_path=metrics_path,
    )


def write_submission_fasta(
    *,
    input_path: Path,
    stage_output_path: Path,
    final_output_path: Path,
    log_path: Path,
    metrics_path: Path | None = None,
) -> None:
    write_mss_fasta(
        input_path=input_path,
        output_path=stage_output_path,
        log_path=log_path,
        metrics_path=metrics_path,
    )
    link_or_copy(stage_output_path, final_output_path)


def copy_input_fasta(
    *,
    input_path: Path,
    output_path: Path,
    log_path: Path,
    metrics_path: Path | None = None,
) -> None:
    started_at = datetime.now()
    seen_ids: set[str] = set()
    duplicate_ids: set[str] = set()
    empty_ids: list[str] = []
    opener = gzip.open if input_path.suffix == ".gz" else Path.open
    with opener(input_path, "rt", encoding="utf-8") as input_handle:
        for record in iter_fasta_handle(input_handle):
            if record.id in seen_ids:
                duplicate_ids.add(record.id)
            seen_ids.add(record.id)
            if not record.sequence:
                empty_ids.append(record.id)
    if not seen_ids:
        raise MSSPackError("FASTA input does not contain any sequence records")
    if duplicate_ids:
        raise MSSPackError(
            "Duplicate FASTA sequence IDs are not supported: "
            + ", ".join(sorted(duplicate_ids)[:10])
        )
    if empty_ids:
        raise MSSPackError(
            "FASTA records must contain sequence data: " + ", ".join(empty_ids[:10])
        )
    copy_or_decompress(input_path, output_path)
    output_total = len(seen_ids)
    write_step_log(
        log_path=log_path,
        command=f"msspack pipeline copy-input-fasta --input {input_path} --output {output_path}",
        step="copy-input-fasta",
        started_at=started_at,
        count_unit="sequences",
        input_total=output_total,
        changed_total=0,
        output_total=output_total,
        details=[f"Source file: {input_path}"],
    )
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=metrics_path,
            step="copy-input-fasta",
            count_unit="sequences",
            input_total=output_total,
            changed_total=0,
            output_total=output_total,
            details={"source_file": str(input_path)},
        )


def copy_input_gff(
    *,
    input_path: Path,
    output_path: Path,
    log_path: Path,
    metrics_path: Path | None = None,
) -> None:
    started_at = datetime.now()
    copy_or_decompress(input_path, output_path)
    output_total = count_gff_feature_records(output_path)
    write_step_log(
        log_path=log_path,
        command=f"msspack pipeline copy-input-gff --input {input_path} --output {output_path}",
        step="copy-input-gff",
        started_at=started_at,
        count_unit="GFF features",
        input_total=output_total,
        changed_total=0,
        output_total=output_total,
        details=[f"Source file: {input_path}"],
    )
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=metrics_path,
            step="copy-input-gff",
            count_unit="GFF features",
            input_total=output_total,
            changed_total=0,
            output_total=output_total,
            details={"source_file": str(input_path)},
        )


def write_gapjust_passthrough(
    *,
    input_fasta: Path,
    input_gff: Path,
    output_fasta: Path,
    output_gff: Path,
    log_path: Path,
    metrics_path: Path | None = None,
) -> None:
    started_at = datetime.now()
    link_or_copy(input_fasta, output_fasta)
    link_or_copy(input_gff, output_gff)
    input_total = count_fasta_records(output_fasta)
    write_step_log(
        log_path=log_path,
        command="msspack internal gapjust --skipped",
        step="gapjust",
        started_at=started_at,
        status="skipped",
        count_unit="sequences",
        input_total=input_total,
        changed_total=0,
        output_total=input_total,
        details=[
            "Gap normalization skipped because pipeline.run_gapjust = false.",
            f"Linked-or-copied FASTA: {input_fasta} -> {output_fasta}",
            f"Linked-or-copied GFF: {input_gff} -> {output_gff}",
            f"Input GFF feature total: {count_gff_feature_records(output_gff):,}",
        ],
    )
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=metrics_path,
            step="gapjust",
            status="skipped",
            count_unit="sequences",
            input_total=input_total,
            changed_total=0,
            output_total=input_total,
            details={
                "skipped": True,
                "linked_fasta": str(output_fasta),
                "linked_gff": str(output_gff),
                "input_gff_feature_total": count_gff_feature_records(output_gff),
            },
        )


def run_fix_gff_to_inframe(
    *,
    input_path: Path,
    output_path: Path,
    summary_log_path: Path,
    command_log_path: Path,
    updated_gene_ids_path: Path | None = None,
    metrics_path: Path | None = None,
) -> None:
    started_at = datetime.now()
    summary = fix_gff_to_inframe(
        input_path=input_path,
        output_path=output_path,
        log_path=summary_log_path,
        updated_gene_ids_path=updated_gene_ids_path,
        metrics_path=metrics_path,
    )
    input_total = summary["updated_gene_models"] + summary["unchanged_gene_models"]
    write_step_log(
        log_path=command_log_path,
        command=f"msspack internal update-gff-to-inframe --input {input_path} --output {output_path} --log {summary_log_path}",
        step="update-gff-to-inframe",
        started_at=started_at,
        count_unit="genes",
        input_total=input_total,
        changed_total=summary["updated_gene_models"],
        output_total=input_total,
        details=[
            f"Number of unchanged gene models: {summary['unchanged_gene_models']:,}",
            f"Removed features: {summary['removed_features']:,}",
            f"Summary log: {summary_log_path}",
        ],
    )


def run_apply_padding_to_gff(
    *,
    gff_path: Path,
    padding_log_path: Path,
    output_path: Path,
    genes_with_stops_path: Path,
    updated_genes_path: Path,
    command_log_path: Path,
    metrics_path: Path | None = None,
) -> None:
    started_at = datetime.now()
    summary = apply_padding_to_gff(
        gff_path=gff_path,
        padding_log_path=padding_log_path,
        output_path=output_path,
        genes_with_stops_path=genes_with_stops_path,
        updated_genes_path=updated_genes_path,
    )
    input_total = count_gff_gene_records(gff_path)
    lines = [
        f"Number of genes with stops (new_num_stop>0): {len(summary['genes_with_stops']):,}",
        f"Number of updated genes (new_num_stop=0): {len(summary['updated_genes']):,}",
    ]
    warnings = summary["forced_first_cds_warnings"]
    if warnings:
        lines.append("Forced first CDS frame to 0 for the following:")
        lines.extend("  " + warning for warning in warnings)
    write_step_log(
        log_path=command_log_path,
        command=(
            f"msspack internal update-gff-with-padding --gff {gff_path} "
            f"--padding-log {padding_log_path} --output {output_path} "
            f"--genes-with-stops {genes_with_stops_path} --updated-genes {updated_genes_path}"
        ),
        step="update-gff-with-padding",
        started_at=started_at,
        count_unit="genes",
        input_total=input_total,
        changed_total=len(summary["updated_genes"]),
        output_total=count_gff_gene_records(output_path),
        details=lines,
    )
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=metrics_path,
            step="update-gff-with-padding",
            count_unit="genes",
            input_total=input_total,
            changed_total=len(summary["updated_genes"]) + len(summary["genes_with_stops"]),
            output_total=count_gff_gene_records(output_path),
            details={
                "genes_with_stops": len(summary["genes_with_stops"]),
                "updated_genes": len(summary["updated_genes"]),
                "forced_first_cds_warnings": summary["forced_first_cds_warnings"],
                "genes_with_stops_path": str(genes_with_stops_path),
                "updated_genes_path": str(updated_genes_path),
            },
        )


def run_select_one_mrna(
    *,
    input_path: Path,
    output_path: Path,
    log_path: Path,
    changed_gene_ids_path: Path | None = None,
    metrics_path: Path | None = None,
) -> None:
    started_at = datetime.now()
    summary = select_one_mrna_per_gene(
        input_gff_path=input_path,
        output_gff_path=output_path,
        changed_gene_ids_path=changed_gene_ids_path,
    )
    input_total = count_gff_gene_records(input_path)
    output_total = count_gff_gene_records(output_path)
    write_step_log(
        log_path=log_path,
        command=f"msspack internal select-one-mrna --input {input_path} --output {output_path}",
        step="select-one-mrna",
        started_at=started_at,
        count_unit="genes",
        input_total=input_total,
        changed_total=summary["changed_genes"],
        output_total=output_total,
        details=[f"Removed mRNAs: {summary['removed_mrnas']:,}"],
    )
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=metrics_path,
            step="select-one-mrna",
            count_unit="genes",
            input_total=input_total,
            changed_total=int(summary["changed_genes"]),
            output_total=output_total,
            details={
                "removed_mrnas": int(summary["removed_mrnas"]),
                "changed_gene_ids_path": str(changed_gene_ids_path) if changed_gene_ids_path else "",
            },
        )


def run_mss_cds_to_misc(
    *,
    mss_input_path: Path,
    genes_input_path: Path,
    mss_output_path: Path,
    log_path: Path,
    converted_gene_ids_path: Path | None = None,
    metrics_path: Path | None = None,
    locus_tag_prefix: str = "",
) -> None:
    started_at = datetime.now()
    summary = convert_cds_features_to_misc(
        mss_input_path=mss_input_path,
        genes_input_path=genes_input_path,
        mss_output_path=mss_output_path,
        converted_gene_ids_path=converted_gene_ids_path,
        locus_tag_prefix=locus_tag_prefix,
    )
    input_total = count_mss_feature_blocks(mss_input_path)
    output_total = count_mss_feature_blocks(mss_output_path)
    write_step_log(
        log_path=log_path,
        command=(
            f"msspack internal mss-cds-to-misc --mss-input {mss_input_path} "
            f"--genes-input {genes_input_path} --output {mss_output_path}"
        ),
        step="mss-cds-to-misc",
        started_at=started_at,
        count_unit="feature blocks",
        input_total=input_total,
        changed_total=summary["edited_genes"],
        output_total=output_total,
        details=[
            f"Number of genes read from {genes_input_path}: {summary['genes_read']:,}",
            f"Total number of CDS in input: {summary['cds_input']:,}",
            f"Total number of CDS in output: {summary['cds_output']:,}",
            f"Total number of misc_feature in output: {summary['misc_feature_output']:,}",
        ],
    )
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=metrics_path,
            step="mss-cds-to-misc",
            count_unit="feature blocks",
            input_total=input_total,
            changed_total=int(summary["edited_genes"]),
            output_total=output_total,
            details={
                "genes_read": int(summary["genes_read"]),
                "cds_input": int(summary["cds_input"]),
                "cds_output": int(summary["cds_output"]),
                "misc_feature_output": int(summary["misc_feature_output"]),
                "converted_gene_ids_path": str(converted_gene_ids_path) if converted_gene_ids_path else "",
            },
        )


def sort_gff(
    input_path: Path,
    output_path: Path,
    *,
    command_override: str | None,
    log_path: Path,
    metrics_path: Path | None = None,
) -> None:
    started_at = datetime.now()
    sort_gff_file_precise(input_path, output_path)
    changed_total = count_reordered_feature_lines(input_path, output_path)
    notes: list[str] = []
    if command_override:
        notes.append(f"Legacy tools.gff3sort value is ignored: {command_override}")
    notes.append("Sorted with built-in precise GFF sorter.")
    write_step_log(
        log_path=log_path,
        command=f"msspack internal gff3sort --input {input_path} --output {output_path}",
        step="gff3sort",
        started_at=started_at,
        count_unit="GFF features",
        input_total=count_gff_feature_records(input_path),
        changed_total=changed_total,
        output_total=count_gff_feature_records(output_path),
        details=notes,
    )
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=metrics_path,
            step="gff3sort",
            count_unit="GFF features",
            input_total=count_gff_feature_records(input_path),
            changed_total=changed_total,
            output_total=count_gff_feature_records(output_path),
            details={
                "legacy_gff3sort_setting": command_override or "",
                "sorter": "built-in precise GFF sorter",
            },
        )


def pad_locus_tags(
    input_path: Path,
    output_path: Path,
    prefix: str,
    digits: int,
    *,
    log_path: Path,
    metrics_path: Path | None = None,
) -> None:
    started_at = datetime.now()
    pattern = re.compile(rf"\b{re.escape(prefix)}_(\d+)\b")
    text = input_path.read_text(encoding="utf-8")
    trailing_newline = text.endswith("\n")
    replacements = 0
    modified_lines = 0

    def repl(match: re.Match[str]) -> str:
        nonlocal replacements
        replacements += 1
        return f"{prefix}_{int(match.group(1)):0{digits}d}"

    output_lines: list[str] = []
    for line in text.splitlines():
        updated = pattern.sub(repl, line)
        if updated != line:
            modified_lines += 1
        output_lines.append(updated)
    output_text = "\n".join(output_lines)
    if trailing_newline:
        output_text += "\n"
    write_text(output_path, output_text)
    write_step_log(
        log_path=log_path,
        command=f"msspack internal pad-locus-tags --input {input_path} --output {output_path}",
        step="pad-locus-tags",
        started_at=started_at,
        count_unit="GFF features",
        input_total=count_gff_feature_records(input_path),
        changed_total=modified_lines,
        output_total=count_gff_feature_records(output_path),
        details=[
            f"Prefix: {prefix}",
            f"Digits: {digits}",
            f"Updated locus_tag references: {replacements:,}",
        ],
    )
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=metrics_path,
            step="pad-locus-tags",
            count_unit="GFF features",
            input_total=count_gff_feature_records(input_path),
            changed_total=modified_lines,
            output_total=count_gff_feature_records(output_path),
            details={
                "prefix": prefix,
                "digits": digits,
                "updated_locus_tag_references": replacements,
            },
        )


def replace_annotation_products(
    input_path: Path,
    output_path: Path,
    patterns: list[str],
    replacement: str,
    *,
    log_path: Path,
    metrics_path: Path | None = None,
) -> None:
    started_at = datetime.now()
    lines = input_path.read_text(encoding="utf-8").splitlines()
    input_total = count_tabular_rows(input_path)
    replaced_count = 0
    if not patterns:
        write_text(output_path, "\n".join(lines) + "\n")
        write_step_log(
            log_path=log_path,
            command=f"msspack internal normalize-annotation-products --input {input_path} --output {output_path}",
            step="normalize-annotation-products",
            started_at=started_at,
            count_unit="annotation rows",
            input_total=input_total,
            changed_total=0,
            output_total=input_total,
            details=["No replacement patterns configured."],
        )
        if metrics_path is not None:
            write_step_metrics(
                metrics_path=metrics_path,
                step="normalize-annotation-products",
                count_unit="annotation rows",
                input_total=input_total,
                changed_total=0,
                output_total=input_total,
                details={"replacement_patterns": 0, "rows_replaced": 0},
            )
        return

    matcher = re.compile(rf"^({'|'.join(pattern for pattern in patterns if pattern)})$", flags=re.IGNORECASE)
    output_lines: list[str] = []
    for index, line in enumerate(lines):
        if index == 0:
            output_lines.append(line)
            continue
        parts = line.split("\t")
        if len(parts) >= 3:
            description = parts[1].replace("\r", "").replace("\xa0", "").strip()
            if description.startswith("Belongs to the "):
                description = description[len("Belongs to the ") :]
            if matcher.match(description):
                parts[1] = replacement
                replaced_count += 1
            else:
                parts[1] = description
        output_lines.append("\t".join(parts))
    write_text(output_path, "\n".join(output_lines) + "\n")
    write_step_log(
        log_path=log_path,
        command=f"msspack internal normalize-annotation-products --input {input_path} --output {output_path}",
        step="normalize-annotation-products",
        started_at=started_at,
        count_unit="annotation rows",
        input_total=input_total,
        changed_total=replaced_count,
        output_total=count_tabular_rows(output_path),
        details=[
            f"Configured replacement patterns: {len(patterns):,}",
            f"Rows replaced: {replaced_count:,}",
        ],
    )
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=metrics_path,
            step="normalize-annotation-products",
            count_unit="annotation rows",
            input_total=input_total,
            changed_total=replaced_count,
            output_total=count_tabular_rows(output_path),
            details={
                "replacement_patterns": len(patterns),
                "rows_replaced": replaced_count,
                "replacement": replacement,
            },
        )


def run_gff2mss(
    config: MSSPackConfig,
    *,
    fasta_path: Path,
    gff_path: Path,
    ann_table_path: Path,
    output_path: Path,
    log_path: Path,
    metrics_path: Path | None = None,
) -> None:
    started_at = datetime.now()
    summary = convert_gff_to_mss(
        ConversionOptions(
            fasta_path=fasta_path,
            gff_path=gff_path,
            annotation_path=ann_table_path,
            output_path=output_path,
            locus_tag_prefix=config.sample.locus_tag,
            organism_name=config.sample.scientific_name,
            strain=config.sample.strain,
            linkage_evidence=config.sample.linkage_evidence,
            genetic_code=config.sample.genetic_code,
            infer_boundary=config.pipeline.infer_complete,
            isolate=config.sample.isolate,
            sex=config.sample.sex,
            country=config.sample.geo_loc_name,
            collection_date=config.sample.collection_date,
            min_assembly_gap_size=config.pipeline.min_assembly_gap,
            gap_estimated_length=config.pipeline.gap_assembly,
            feature_with_gap=config.pipeline.feature_with_gap,
            minimum_intron_size_cutoff=config.pipeline.min_artificial_intron_size,
        )
    )
    lines = [
        (
            'Annotation file contains "Locus_tag" column.'
            if summary.used_custom_locus_tags
            else 'Annotation file does not contain "Locus_tag" column.'
        ),
    ]
    lines.extend(
        format_event_summary(contig.event_counts, f"Finished {contig.contig_name}")
        for contig in summary.contigs
    )
    lines.append(format_event_summary(summary.overall_counts, "Overall summary"))
    write_step_log(
        log_path=log_path,
        command=(
            "msspack convert-gff-to-mss "
            f"--fasta {fasta_path} --gff {gff_path} --ann {ann_table_path} "
            f"--out {output_path} --loc {config.sample.locus_tag}"
        ),
        step="gff2mss",
        started_at=started_at,
        count_unit="contigs",
        input_total=len(summary.contigs),
        changed_total=sum(1 for contig in summary.contigs if contig.event_counts),
        output_total=len(summary.contigs),
        details=lines,
    )
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=metrics_path,
            step="gff2mss",
            count_unit="contigs",
            input_total=len(summary.contigs),
            changed_total=sum(1 for contig in summary.contigs if contig.event_counts),
            output_total=len(summary.contigs),
            details={
                "used_custom_locus_tags": summary.used_custom_locus_tags,
                "contigs_with_events": sum(1 for contig in summary.contigs if contig.event_counts),
                "summary_lines": lines,
            },
        )


def render_final_annotation_with_log(
    *,
    header: str,
    body_path: Path,
    output_path: Path,
    log_path: Path,
    metrics_path: Path | None = None,
) -> None:
    started_at = datetime.now()
    body = body_path.read_text(encoding="utf-8")
    final_text = render_final_annotation(header, body)
    write_text(output_path, final_text)
    header_rows = count_nonempty_lines(header)
    body_rows = count_nonempty_lines(body)
    changed_total = len(re.findall(r"(^\t\t\t)country(\t)", body, flags=re.MULTILINE))
    write_step_log(
        log_path=log_path,
        command=f"msspack pipeline render-final-annotation --body {body_path} --output {output_path}",
        step="render-final-annotation",
        started_at=started_at,
        count_unit="annotation rows",
        input_total=header_rows + body_rows,
        changed_total=changed_total,
        output_total=count_nonempty_lines(final_text),
        details=[
            f"Header rows: {header_rows:,}",
            f"Body rows: {body_rows:,}",
            f"Normalized geo_loc_name rows: {changed_total:,}",
        ],
    )
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=metrics_path,
            step="render-final-annotation",
            count_unit="annotation rows",
            input_total=header_rows + body_rows,
            changed_total=changed_total,
            output_total=count_nonempty_lines(final_text),
            details={
                "header_rows": header_rows,
                "body_rows": body_rows,
                "normalized_geo_loc_name_rows": changed_total,
            },
        )
