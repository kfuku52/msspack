from __future__ import annotations

import re
from collections import Counter
from collections.abc import Iterable
from datetime import datetime
from pathlib import Path

from .fasta import iter_fasta, write_fasta_record
from .gff import GFFDocument, read_gff_document, write_gff_document
from .step_logging import write_step_log, write_step_metrics
from .utils import atomic_text_writer, ensure_dir, link_or_copy

# Some of the gap-normalization logic below is carried forward from earlier
# internal tooling used before msspack unified the MSS packaging workflow.

GapJustification = tuple[int, int]


def _justify_gap(
    gap_length: int,
    target_length: int,
    gap_just_min: int | None,
    gap_just_max: int | None,
) -> bool:
    if gap_length == target_length:
        return False
    if gap_length < target_length:
        if gap_just_min is not None and gap_length < gap_just_min:
            return False
        return True
    if gap_just_max is not None and gap_length > gap_just_max:
        return False
    return True


def _apply_gapjust_to_sequence(
    seq: str,
    *,
    target_gap_length: int,
    gap_just_min: int | None,
    gap_just_max: int | None,
) -> tuple[str, list[GapJustification], int, int | None, int]:
    justifications: list[GapJustification] = []
    rebuilt: list[str] = []
    cursor = 0
    num_justifications = 0
    min_original_gap_length: int | None = None
    max_original_gap_length = 0

    for match in re.finditer("N+", seq.replace("n", "N")):
        start = match.start()
        end = match.end()
        gap_length = end - start
        rebuilt.append(seq[cursor:start])
        if _justify_gap(gap_length, target_gap_length, gap_just_min, gap_just_max):
            rebuilt.append("N" * target_gap_length)
            edit_len = target_gap_length - gap_length
            justifications.append((start, edit_len))
            num_justifications += 1
            if (min_original_gap_length is None) or (gap_length < min_original_gap_length):
                min_original_gap_length = gap_length
            if gap_length > max_original_gap_length:
                max_original_gap_length = gap_length
        else:
            rebuilt.append(seq[start:end])
        cursor = end
    rebuilt.append(seq[cursor:])
    return "".join(rebuilt), justifications, num_justifications, min_original_gap_length, max_original_gap_length


def _shift_coordinate(value: int, justifications: Iterable[GapJustification]) -> int:
    updated = value
    cumulative_offset = 0
    for original_start_zero_based, edit_len in sorted(justifications, key=lambda item: item[0]):
        if edit_len == 0:
            continue
        actual_edit_start = original_start_zero_based + 1 + cumulative_offset
        if updated > actual_edit_start:
            updated += edit_len
        cumulative_offset += edit_len
    return updated


def _apply_gapjust_to_gff(
    document: GFFDocument,
    justifications_by_seq: dict[str, list[GapJustification]],
) -> tuple[int, int, int]:
    changed_starts = 0
    changed_ends = 0
    changed_genes = 0
    for record in document.records:
        justifications = justifications_by_seq.get(record.seqid)
        if not justifications:
            continue
        old_start = record.start
        old_end = record.end
        record.start = _shift_coordinate(record.start, justifications)
        record.end = _shift_coordinate(record.end, justifications)
        start_changed = old_start != record.start
        end_changed = old_end != record.end
        if start_changed:
            changed_starts += 1
        if end_changed:
            changed_ends += 1
        if record.type == "gene" and (start_changed or end_changed):
            changed_genes += 1
    return changed_starts, changed_ends, changed_genes


def normalize_gap_lengths(
    *,
    fasta_path: Path,
    output_fasta_path: Path,
    log_path: Path,
    gap_len: int,
    gap_just_min: int | None = None,
    gap_just_max: int | None = None,
    input_gff_path: Path | None = None,
    output_gff_path: Path | None = None,
    metrics_path: Path | None = None,
) -> None:
    if gap_len < 0:
        raise ValueError("--gap_len must be >= 0")
    if gap_just_min is not None and gap_just_min < 0:
        raise ValueError("--gap_just_min must be >= 0")
    if gap_just_max is not None and gap_just_max < 0:
        raise ValueError("--gap_just_max must be >= 0")

    command = [
        "msspack",
        "internal",
        "gapjust",
        "--gap_len",
        str(gap_len),
        "--seqfile",
        str(fasta_path),
        "--outfile",
        str(output_fasta_path),
    ]
    if gap_just_min is not None:
        command.extend(["--gap_just_min", str(gap_just_min)])
    if gap_just_max is not None:
        command.extend(["--gap_just_max", str(gap_just_max)])
    if input_gff_path is not None and output_gff_path is not None:
        command.extend(["--ingff", str(input_gff_path), "--outgff", str(output_gff_path)])

    started_at = datetime.now()
    justifications_by_seq: dict[str, list[GapJustification]] = {}
    num_justifications = 0
    min_original_gap_length: int | None = None
    max_original_gap_length = 0
    seen_ids: Counter[str] = Counter()
    records: list[tuple[str, str, str]] = []

    for record in iter_fasta(fasta_path):
        seen_ids[record.id] += 1
        if seen_ids[record.id] > 1 and input_gff_path is not None:
            shown = ",".join(sorted(seq_id for seq_id, count in seen_ids.items() if count > 1)[:10])
            if sum(count > 1 for count in seen_ids.values()) > 10:
                shown += ",..."
            raise ValueError(
                "Duplicate sequence IDs are not supported with GFF gapjust: " + shown
            )
        (
            new_seq,
            justifications,
            record_num_justifications,
            record_min_original_gap_length,
            record_max_original_gap_length,
        ) = _apply_gapjust_to_sequence(
            record.sequence,
            target_gap_length=gap_len,
            gap_just_min=gap_just_min,
            gap_just_max=gap_just_max,
        )
        records.append((record.id, record.description, new_seq))
        if justifications:
            justifications_by_seq[record.id] = justifications
        num_justifications += record_num_justifications
        if record_min_original_gap_length is not None:
            if (min_original_gap_length is None) or (
                record_min_original_gap_length < min_original_gap_length
            ):
                min_original_gap_length = record_min_original_gap_length
            max_original_gap_length = max(
                max_original_gap_length,
                record_max_original_gap_length,
            )

    changed_sequences = len(justifications_by_seq)
    detail_lines = [f"Number of gap justifications: {num_justifications:,}"]
    metrics_details: dict[str, object] = {
        "gap_len": gap_len,
        "gap_just_min": gap_just_min,
        "gap_just_max": gap_just_max,
        "gap_justifications": num_justifications,
        "min_original_gap_length": min_original_gap_length,
        "max_original_gap_length": max_original_gap_length,
    }
    if num_justifications:
        detail_lines.append(
            f"Minimum and maximum original gap lengths: {min_original_gap_length} and {max_original_gap_length}"
        )
    else:
        detail_lines.append("No gap edits were made.")
        link_or_copy(fasta_path, output_fasta_path)
    if num_justifications:
        ensure_dir(output_fasta_path.parent)
        with atomic_text_writer(output_fasta_path) as handle:
            for record_id, description, sequence in records:
                write_fasta_record(
                    handle,
                    identifier=record_id,
                    description=description,
                    sequence=sequence,
                )

    if input_gff_path is not None and output_gff_path is not None:
        document = read_gff_document(input_gff_path)
        detail_lines.append(f"Number of input GFF header lines: {len(document.header_lines):,}")
        detail_lines.append(f"Number of input GFF records: {len(document.records):,}")
        detail_lines.append(
            f"Number of input GFF unique seqids: {len({record.seqid for record in document.records}):,}"
        )
        if num_justifications:
            changed_starts, changed_ends, changed_genes = _apply_gapjust_to_gff(
                document,
                justifications_by_seq,
            )
        else:
            changed_starts = changed_ends = changed_genes = 0
        detail_lines.append(f"Number of justified GFF start coordinates: {changed_starts:,}")
        detail_lines.append(f"Number of justified GFF end coordinates: {changed_ends:,}")
        detail_lines.append(f"Number of justified GFF gene features: {changed_genes:,}")
        metrics_details.update(
            {
                "input_gff_records": len(document.records),
                "output_gff_records": len(document.records),
                "input_gff_header_lines": len(document.header_lines),
                "output_gff_header_lines": len(document.header_lines),
                "justified_gff_start_coordinates": changed_starts,
                "justified_gff_end_coordinates": changed_ends,
                "justified_gff_gene_features": changed_genes,
            }
        )
        if num_justifications:
            write_gff_document(output_gff_path, document)
        else:
            link_or_copy(input_gff_path, output_gff_path)
        detail_lines.append(f"Number of output GFF header lines: {len(document.header_lines):,}")
        detail_lines.append(f"Number of output GFF records: {len(document.records):,}")
        detail_lines.append(
            f"Number of output GFF unique seqids: {len({record.seqid for record in document.records}):,}"
        )

    write_step_log(
        log_path=log_path,
        command=command,
        step="gapjust",
        started_at=started_at,
        count_unit="sequences",
        input_total=len(records),
        changed_total=changed_sequences,
        output_total=len(records),
        details=detail_lines,
    )
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=metrics_path,
            step="gapjust",
            count_unit="sequences",
            input_total=len(records),
            changed_total=changed_sequences,
            output_total=len(records),
            details=metrics_details,
        )
