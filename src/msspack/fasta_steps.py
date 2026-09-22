from __future__ import annotations

from datetime import datetime
from pathlib import Path

from .fasta import iter_fasta, write_fasta_record
from .step_logging import write_step_log, write_step_metrics
from .utils import MSSPackError, atomic_text_writer, ensure_dir


def remove_trailing_ns_fasta(
    *,
    input_path: Path,
    output_path: Path,
    log_path: Path,
    metrics_path: Path | None = None,
) -> None:
    started_at = datetime.now()
    removed_total = 0
    input_total = 0
    changed_total = 0
    ensure_dir(output_path.parent)
    with atomic_text_writer(output_path) as out_handle:
        for record in iter_fasta(input_path):
            sequence = record.sequence.rstrip("Nn")
            if not sequence:
                raise MSSPackError(f"No sequence remains after trimming terminal Ns: {record.id}")
            input_total += 1
            removed = len(record.sequence) - len(sequence)
            removed_total += removed
            changed_total += int(removed > 0)
            write_fasta_record(out_handle, identifier=record.id,
                               description=record.description, sequence=sequence)

    write_step_log(
        log_path=log_path,
        command=f"msspack internal remove-trailing-ns --input {input_path} --output {output_path}",
        step="remove-trailing-ns",
        started_at=started_at,
        count_unit="sequences",
        input_total=input_total,
        changed_total=changed_total,
        output_total=input_total,
        details=[f"Total Ns removed: {removed_total:,}"],
    )
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=metrics_path,
            step="remove-trailing-ns",
            count_unit="sequences",
            input_total=input_total,
            changed_total=changed_total,
            output_total=input_total,
            details={"removed_bases": removed_total},
        )


def write_mss_fasta(
    *,
    input_path: Path,
    output_path: Path,
    log_path: Path,
    metrics_path: Path | None = None,
) -> None:
    started_at = datetime.now()
    entries = 0
    ensure_dir(output_path.parent)
    with atomic_text_writer(output_path) as out_handle:
        for record in iter_fasta(input_path):
            if not record.sequence:
                raise MSSPackError(f"Cannot write an empty submission sequence: {record.id}")
            out_handle.write(f">{record.id}\n")
            seq = record.sequence
            for start in range(0, len(seq), 60):
                out_handle.write(seq[start : start + 60] + "\n")
            out_handle.write("//\n")
            entries += 1

    write_step_log(
        log_path=log_path,
        command=f"msspack internal fasta-to-mss --input {input_path} --output {output_path}",
        step="fasta-to-mss",
        started_at=started_at,
        count_unit="sequences",
        input_total=entries,
        changed_total=entries,
        output_total=entries,
        details=[f"Entries written: {entries:,}"],
    )
    if metrics_path is not None:
        write_step_metrics(
            metrics_path=metrics_path,
            step="fasta-to-mss",
            count_unit="sequences",
            input_total=entries,
            changed_total=entries,
            output_total=entries,
            details={"entries_written": entries},
        )
