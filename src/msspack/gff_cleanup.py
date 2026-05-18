from __future__ import annotations

from collections import defaultdict, deque
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .fasta import iter_fasta
from .gff import parse_attributes
from .step_logging import write_id_list, write_step_log, write_step_metrics
from .utils import ensure_dir


def trim_gff_to_fasta_bounds(
    *,
    gff_path: Path,
    fasta_path: Path,
    output_path: Path,
    log_path: Path,
    metrics_path: Optional[Path] = None,
) -> None:
    started_at = datetime.now()
    seq_lengths: Dict[str, int] = {}
    for record in iter_fasta(fasta_path):
        seq_lengths[record.id] = len(record.sequence)

    input_total = 0
    updated_count = 0
    written_count = 0
    ensure_dir(output_path.parent)
    with gff_path.open("r", encoding="utf-8") as in_handle, output_path.open(
        "w", encoding="utf-8"
    ) as out_handle:
        for line in in_handle:
            if line.startswith("#") or not line.strip():
                out_handle.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            input_total += 1
            if len(fields) < 9:
                out_handle.write(line)
                continue
            seqid = fields[0]
            if seqid not in seq_lengths:
                continue
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                out_handle.write(line)
                continue
            original = (start, end)
            if start < 1:
                start = 1
            if end > seq_lengths[seqid]:
                end = seq_lengths[seqid]
            if start > end:
                continue
            if (start, end) != original:
                updated_count += 1
            fields[3] = str(start)
            fields[4] = str(end)
            out_handle.write("\t".join(fields) + "\n")
            written_count += 1

    write_step_log(
        log_path=log_path,
        command=(
            f"msspack internal trim-gff-coordinates --gff {gff_path} "
            f"--fasta {fasta_path} --output {output_path}"
        ),
        step="trim-gff-coordinates",
        started_at=started_at,
        count_unit="GFF features",
        input_total=input_total,
        changed_total=updated_count,
        output_total=written_count,
        details=[f"Skipped features: {input_total - written_count:,}"],
    )
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=metrics_path,
            step="trim-gff-coordinates",
            count_unit="GFF features",
            input_total=input_total,
            changed_total=updated_count,
            output_total=written_count,
            details={"skipped_features": input_total - written_count},
        )


def drop_duplicate_coordinate_genes(
    *,
    input_path: Path,
    output_path: Path,
    log_path: Path,
    removed_gene_ids_path: Optional[Path] = None,
    metrics_path: Optional[Path] = None,
) -> dict[str, object]:
    started_at = datetime.now()
    lines: List[Tuple[bool, str, Optional[Tuple[str, int, int, str]], Optional[str], Optional[str]]] = []
    parent_to_children: Dict[str, List[str]] = defaultdict(list)
    coords_to_gene: Dict[Tuple[str, int, int, str], str] = {}
    removed_gene_ids: List[str] = []
    input_gene_total = 0

    with input_path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if (not line) or line.startswith("#"):
                lines.append((False, line, None, None, None))
                continue
            fields = line.split("\t")
            if len(fields) < 9:
                lines.append((False, line, None, None, None))
                continue
            attr = parse_attributes(fields[8])
            rec_id = attr.get("ID")
            parent = attr.get("Parent")
            coord_key: Optional[Tuple[str, int, int, str]] = None
            if fields[2] == "gene" and rec_id is not None:
                input_gene_total += 1
                coord_key = (fields[0], int(fields[3]), int(fields[4]), fields[6])
                if coord_key in coords_to_gene:
                    removed_gene_ids.append(rec_id)
                else:
                    coords_to_gene[coord_key] = rec_id
            if rec_id and parent:
                for parent_id in parent.split(","):
                    parent_to_children[parent_id].append(rec_id)
            lines.append((True, line, coord_key, rec_id, parent))

    to_remove = set()
    queue = deque(removed_gene_ids)
    while queue:
        current = queue.popleft()
        if current in to_remove:
            continue
        to_remove.add(current)
        queue.extend(parent_to_children.get(current, []))

    kept = 0
    ensure_dir(output_path.parent)
    with output_path.open("w", encoding="utf-8") as out_handle:
        for is_data, line, _, rec_id, parent in lines:
            if not is_data:
                out_handle.write(line + "\n")
                continue
            if rec_id in to_remove:
                continue
            if parent and any(parent_id in to_remove for parent_id in parent.split(",")):
                continue
            out_handle.write(line + "\n")
            kept += 1

    write_step_log(
        log_path=log_path,
        command=f"msspack internal drop-duplicate-coordinate-gene --input {input_path} --output {output_path}",
        step="drop-duplicate-coordinate-gene",
        started_at=started_at,
        count_unit="genes",
        input_total=input_gene_total,
        changed_total=len(removed_gene_ids),
        output_total=input_gene_total - len(removed_gene_ids),
        details=[f"Output feature count: {kept:,}"],
    )
    if removed_gene_ids_path is not None:
        write_id_list(removed_gene_ids_path, removed_gene_ids)
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=metrics_path,
            step="drop-duplicate-coordinate-gene",
            count_unit="genes",
            input_total=input_gene_total,
            changed_total=len(removed_gene_ids),
            output_total=input_gene_total - len(removed_gene_ids),
            details={
                "output_feature_count": kept,
                "removed_gene_ids_path": str(removed_gene_ids_path) if removed_gene_ids_path else "",
            },
        )
    return {
        "removed_gene_ids": removed_gene_ids,
        "input_gene_total": input_gene_total,
        "output_gene_total": input_gene_total - len(removed_gene_ids),
        "output_feature_count": kept,
    }


def _fix_attributes(attribute_string: str) -> Tuple[str, bool, bool]:
    trimmed_attr = attribute_string.rstrip(";")
    trailing_semicolons_removed = trimmed_attr != attribute_string
    parts = trimmed_attr.split(";")
    new_attributes: List[str] = []
    current_key: Optional[str] = None
    current_value: List[str] = []
    semicolon_value_fixed = False
    for chunk in parts:
        if "=" in chunk:
            if current_key is not None:
                new_attributes.append(f"{current_key}={'.'.join(current_value)}")
            key, value = chunk.split("=", 1)
            current_key = key
            current_value = [value]
        else:
            if current_value:
                current_value[-1] = current_value[-1] + "." + chunk
            else:
                current_value = [chunk]
            semicolon_value_fixed = True
    if current_key is not None:
        new_attributes.append(f"{current_key}={'.'.join(current_value)}")
    return ";".join(new_attributes), semicolon_value_fixed, trailing_semicolons_removed


def fix_gff_semicolons_file(
    *,
    input_path: Path,
    output_path: Path,
    log_path: Path,
    metrics_path: Optional[Path] = None,
) -> None:
    started_at = datetime.now()
    input_total = 0
    total_modified = 0
    internal_fixed = 0
    trailing_fixed = 0
    ensure_dir(output_path.parent)
    with input_path.open("r", encoding="utf-8") as in_handle, output_path.open(
        "w", encoding="utf-8"
    ) as out_handle:
        for raw in in_handle:
            line = raw.rstrip("\n")
            if (not line) or line.startswith("#"):
                out_handle.write(line + "\n")
                continue
            fields = line.split("\t")
            input_total += 1
            if len(fields) < 9:
                out_handle.write(line + "\n")
                continue
            fixed_attr, semicolon_value_fixed, trailing_semicolons_removed = _fix_attributes(fields[8])
            if fixed_attr != fields[8]:
                total_modified += 1
                if semicolon_value_fixed:
                    internal_fixed += 1
                if trailing_semicolons_removed:
                    trailing_fixed += 1
            fields[8] = fixed_attr
            out_handle.write("\t".join(fields) + "\n")

    write_step_log(
        log_path=log_path,
        command=f"msspack internal fix-gff-semicolons --input {input_path} --output {output_path}",
        step="fix-gff-semicolons",
        started_at=started_at,
        count_unit="GFF features",
        input_total=input_total,
        changed_total=total_modified,
        output_total=input_total,
        details=[
            f"Lines with internal fix: {internal_fixed:,}",
            f"Lines with trailing fix: {trailing_fixed:,}",
        ],
    )
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=metrics_path,
            step="fix-gff-semicolons",
            count_unit="GFF features",
            input_total=input_total,
            changed_total=total_modified,
            output_total=input_total,
            details={
                "internal_fixed": internal_fixed,
                "trailing_fixed": trailing_fixed,
            },
        )
