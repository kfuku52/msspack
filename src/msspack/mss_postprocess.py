from __future__ import annotations

from pathlib import Path
from typing import TypedDict

from .step_logging import write_id_list
from .utils import atomic_text_writer


class MssPostprocessSummary(TypedDict):
    genes_read: int
    edited_genes: int
    cds_input: int
    cds_output: int
    misc_feature_output: int
    converted_gene_ids: list[str]


def read_gene_lookup(
    genes_file: str | Path,
    *,
    locus_tag_prefix: str = "",
) -> dict[str, str]:
    genes: dict[str, str] = {}
    with Path(genes_file).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            gene = raw_line.strip()
            if gene:
                genes.setdefault(gene, gene)
                if locus_tag_prefix:
                    prefixed = (
                        gene
                        if gene == locus_tag_prefix or gene.startswith(locus_tag_prefix + "_")
                        else f"{locus_tag_prefix}_{gene}"
                    )
                    genes.setdefault(prefixed, gene)
    return genes


def get_locus_tag_from_block(block_lines: list[str]) -> str | None:
    for line in block_lines:
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) == 5 and fields[3] == "locus_tag":
            return fields[4]
    return None


def process_feature_block(
    block_lines: list[str],
    gene_lookup: dict[str, str],
) -> tuple[list[str], bool, bool, str | None, str | None]:
    if not block_lines:
        return block_lines, False, False, None, None

    header = block_lines[0].rstrip("\r\n").split("\t")
    feature_type = header[1] if len(header) >= 3 else None
    was_cds = feature_type == "CDS"
    locus_tag = get_locus_tag_from_block(block_lines)
    matched_gene_id = gene_lookup.get(locus_tag) if locus_tag and was_cds else None
    if matched_gene_id is None:
        return block_lines, was_cds, False, feature_type, None

    new_block_lines: list[str] = []
    for index, line in enumerate(block_lines):
        fields = line.rstrip("\r\n").split("\t")
        if index == 0:
            fields[1] = "misc_feature"
        if len(fields) == 5:
            if fields[3] in {"transl_table", "codon_start"}:
                if index != 0:
                    continue
                fields[3:] = ["", ""]
            elif fields[3] == "product":
                fields[3] = "note"
        ending = line[len(line.rstrip("\r\n")):]
        new_block_lines.append("\t".join(fields) + ending)
    return new_block_lines, True, True, "misc_feature", matched_gene_id


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
    locus_tag_prefix: str = "",
) -> MssPostprocessSummary:
    gene_lookup = read_gene_lookup(genes_input_path, locus_tag_prefix=locus_tag_prefix)
    total_cds_input = 0
    total_edited_genes = 0
    total_cds_output = 0
    total_misc_feature_output = 0
    converted_gene_ids: list[str] = []

    with Path(mss_input_path).open("r", encoding="utf-8") as infile, atomic_text_writer(
        Path(mss_output_path)
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
