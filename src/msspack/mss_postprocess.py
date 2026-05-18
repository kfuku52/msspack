from __future__ import annotations

import re
from pathlib import Path
from typing import TypedDict

from .step_logging import write_id_list

_LOCUS_TAG_PATTERN = re.compile(r"\blocus_tag\s+(\S+)")


class MssPostprocessSummary(TypedDict):
    genes_read: int
    edited_genes: int
    cds_input: int
    cds_output: int
    misc_feature_output: int
    converted_gene_ids: list[str]


def normalize_gene_identifier(identifier: str) -> str:
    if "_" in identifier:
        return identifier.split("_")[-1]
    return identifier


def read_gene_lookup(genes_file: str | Path) -> dict[str, str]:
    genes: dict[str, str] = {}
    with Path(genes_file).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            gene = raw_line.strip()
            if gene:
                genes.setdefault(gene, gene)
                genes.setdefault(normalize_gene_identifier(gene), gene)
    return genes


def get_locus_tag_from_block(block_lines: list[str]) -> str | None:
    for line in block_lines:
        match = _LOCUS_TAG_PATTERN.search(line)
        if match:
            return match.group(1)
    return None


def process_feature_block(
    block_lines: list[str],
    gene_lookup: dict[str, str],
) -> tuple[list[str], bool, bool, str | None, str | None]:
    if not block_lines:
        return block_lines, False, False, None, None

    header_line = block_lines[0]
    header_stripped = header_line.lstrip()
    was_cds = False
    was_converted = False
    final_feature_type: str | None

    if header_stripped.startswith("CDS"):
        was_cds = True
        locus_tag = get_locus_tag_from_block(block_lines)
        matched_gene_id = None
        if locus_tag:
            matched_gene_id = gene_lookup.get(locus_tag) or gene_lookup.get(
                normalize_gene_identifier(locus_tag)
            )
        if matched_gene_id is not None:
            indent = header_line[: len(header_line) - len(header_line.lstrip())]
            block_lines[0] = indent + header_stripped.replace("CDS", "misc_feature", 1)
            was_converted = True
            new_block_lines: list[str] = []
            for line in block_lines:
                stripped = line.strip()
                if stripped.startswith("transl_table") or stripped.startswith("codon_start"):
                    continue
                indentation = line[: len(line) - len(line.lstrip())]
                content = line.lstrip()
                if content.startswith("product"):
                    content = content.replace("product", "note", 1)
                new_block_lines.append(indentation + content)
            return new_block_lines, was_cds, was_converted, "misc_feature", matched_gene_id
        final_feature_type = "CDS"
        return block_lines, was_cds, was_converted, final_feature_type, None

    if header_stripped.startswith("misc_feature"):
        final_feature_type = "misc_feature"
    else:
        final_feature_type = header_stripped.split()[0]
    return block_lines, was_cds, was_converted, final_feature_type, None


def _is_block_header(line: str) -> bool:
    if not line.strip():
        return False
    return not line.startswith("\t\t\t")


def convert_cds_features_to_misc(
    *,
    mss_input_path: str | Path,
    genes_input_path: str | Path,
    mss_output_path: str | Path,
    converted_gene_ids_path: str | Path | None = None,
) -> MssPostprocessSummary:
    gene_lookup = read_gene_lookup(genes_input_path)
    total_cds_input = 0
    total_edited_genes = 0
    total_cds_output = 0
    total_misc_feature_output = 0
    converted_gene_ids: list[str] = []

    with Path(mss_input_path).open("r", encoding="utf-8") as infile, Path(mss_output_path).open(
        "w", encoding="utf-8"
    ) as outfile:
        current_block: list[str] = []

        def write_processed_block(block: list[str]) -> None:
            nonlocal total_cds_input, total_edited_genes
            nonlocal total_cds_output, total_misc_feature_output

            processed_lines, was_cds, was_converted, final_ft_type, matched_gene_id = process_feature_block(
                block,
                gene_lookup,
            )
            if was_cds:
                total_cds_input += 1
            if was_converted:
                total_edited_genes += 1
                if matched_gene_id:
                    converted_gene_ids.append(matched_gene_id)
            if final_ft_type == "CDS":
                total_cds_output += 1
            elif final_ft_type == "misc_feature":
                total_misc_feature_output += 1
            outfile.write("".join(processed_lines))

        for line in infile:
            if not line.strip():
                current_block.append(line)
                continue
            if _is_block_header(line):
                if current_block:
                    write_processed_block(current_block)
                    current_block = []
                current_block = [line]
            else:
                current_block.append(line)

        if current_block:
            write_processed_block(current_block)

    if converted_gene_ids_path is not None:
        write_id_list(Path(converted_gene_ids_path), converted_gene_ids)
    return {
        "genes_read": len(set(gene_lookup.values())),
        "edited_genes": total_edited_genes,
        "cds_input": total_cds_input,
        "cds_output": total_cds_output,
        "misc_feature_output": total_misc_feature_output,
        "converted_gene_ids": converted_gene_ids,
    }
