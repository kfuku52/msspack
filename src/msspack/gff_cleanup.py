from __future__ import annotations

from collections import defaultdict, deque
from datetime import datetime
from pathlib import Path

from .fasta import iter_fasta
from .gff import parse_attributes
from .step_logging import write_id_list, write_step_log, write_step_metrics
from .utils import MSSPackError, atomic_text_writer, ensure_dir


def trim_gff_to_fasta_bounds(
    *,
    gff_path: Path,
    fasta_path: Path,
    output_path: Path,
    log_path: Path,
    metrics_path: Path | None = None,
) -> None:
    started_at = datetime.now()
    seq_lengths: dict[str, int] = {}
    for record in iter_fasta(fasta_path):
        seq_lengths[record.id] = len(record.sequence)

    unknown_seqids: dict[str, int] = defaultdict(int)
    with gff_path.open("r", encoding="utf-8") as validation_handle:
        for line_number, line in enumerate(validation_handle, start=1):
            if line.rstrip("\n") == "##FASTA":
                break
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise MSSPackError(
                    f"Invalid GFF record at {gff_path}:{line_number}: expected 9 columns, "
                    f"found {len(fields)}"
                )
            if fields[0] not in seq_lengths:
                unknown_seqids[fields[0]] += 1
    if unknown_seqids:
        shown = ", ".join(
            f"{seqid} ({count} features)"
            for seqid, count in sorted(unknown_seqids.items())[:10]
        )
        if len(unknown_seqids) > 10:
            shown += ", ..."
        raise MSSPackError(
            "GFF seqids are missing from the FASTA; refusing to discard annotations: " + shown
        )

    input_total = 0
    updated_count = 0
    written_count = 0
    ensure_dir(output_path.parent)
    with gff_path.open("r", encoding="utf-8") as in_handle, atomic_text_writer(
        output_path
    ) as out_handle:
        in_fasta = False
        for line in in_handle:
            if in_fasta:
                out_handle.write(line)
                continue
            if line.rstrip("\n") == "##FASTA":
                in_fasta = True
                out_handle.write(line)
                continue
            if line.startswith("#") or not line.strip():
                out_handle.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            input_total += 1
            if len(fields) != 9:  # validated above; retained as a defensive guard
                raise MSSPackError(f"Invalid GFF record encountered while writing: {line!r}")
            seqid = fields[0]
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
    removed_gene_ids_path: Path | None = None,
    metrics_path: Path | None = None,
) -> dict[str, object]:
    started_at = datetime.now()
    lines: list[tuple[bool, str, tuple[str, int, int, str] | None, str | None, str | None]] = []
    parent_to_children: dict[str, list[str]] = defaultdict(list)
    coords_to_gene: dict[tuple[str, int, int, str], str] = {}
    removed_gene_ids: list[str] = []
    input_gene_total = 0

    with input_path.open("r", encoding="utf-8") as handle:
        in_fasta = False
        for raw in handle:
            line = raw.rstrip("\n")
            if in_fasta:
                lines.append((False, line, None, None, None))
                continue
            if line == "##FASTA":
                in_fasta = True
                lines.append((False, line, None, None, None))
                continue
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
            coord_key: tuple[str, int, int, str] | None = None
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
    with atomic_text_writer(output_path) as out_handle:
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


def _fix_attributes(attribute_string: str) -> tuple[str, bool, bool]:
    if attribute_string.strip() in ("", "."):
        return attribute_string, False, False
    trimmed_attr = attribute_string.rstrip(";")
    trailing_semicolons_removed = trimmed_attr != attribute_string
    parts = trimmed_attr.split(";")
    new_attributes: list[str] = []
    current_key: str | None = None
    current_value: list[str] = []
    semicolon_value_fixed = False
    for chunk in parts:
        if "=" in chunk:
            if current_key is not None:
                new_attributes.append(f"{current_key}={'.'.join(current_value)}")
            key, value = chunk.split("=", 1)
            current_key = key
            current_value = [value]
        else:
            if current_key is None or not current_value:
                raise MSSPackError(
                    f"Cannot repair GFF3 attribute fragment without a preceding key: {chunk!r}"
                )
            current_value[-1] = current_value[-1] + "." + chunk
            semicolon_value_fixed = True
    if current_key is not None:
        new_attributes.append(f"{current_key}={'.'.join(current_value)}")
    return ";".join(new_attributes), semicolon_value_fixed, trailing_semicolons_removed


def fix_gff_semicolons_file(
    *,
    input_path: Path,
    output_path: Path,
    log_path: Path,
    metrics_path: Path | None = None,
) -> None:
    started_at = datetime.now()
    input_total = 0
    total_modified = 0
    internal_fixed = 0
    trailing_fixed = 0
    ensure_dir(output_path.parent)
    with input_path.open("r", encoding="utf-8") as in_handle, atomic_text_writer(
        output_path
    ) as out_handle:
        in_fasta = False
        for raw in in_handle:
            line = raw.rstrip("\n")
            if in_fasta:
                out_handle.write(line + "\n")
                continue
            if line == "##FASTA":
                in_fasta = True
                out_handle.write(line + "\n")
                continue
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
