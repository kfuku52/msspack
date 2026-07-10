from __future__ import annotations

from datetime import datetime
from pathlib import Path

from .fasta import iter_fasta
from .step_logging import write_step_log, write_step_metrics
from .utils import atomic_text_writer, ensure_dir


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
    with input_path.open("r", encoding="utf-8") as in_handle, atomic_text_writer(
        output_path
    ) as out_handle:
        header: str | None = None
        seq_lines: list[str] = []
        for raw_line in in_handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    input_total += 1
                    record_changed = False
                    idx = len(seq_lines) - 1
                    while idx >= 0:
                        original = seq_lines[idx]
                        trimmed = original.rstrip("Nn")
                        removed_total += len(original) - len(trimmed)
                        if original != trimmed:
                            record_changed = True
                        if trimmed:
                            seq_lines[idx] = trimmed
                            break
                        seq_lines.pop()
                        idx -= 1
                    if record_changed:
                        changed_total += 1
                    out_handle.write(header + "\n")
                    for seq_line in seq_lines:
                        out_handle.write(seq_line + "\n")
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            input_total += 1
            record_changed = False
            idx = len(seq_lines) - 1
            while idx >= 0:
                original = seq_lines[idx]
                trimmed = original.rstrip("Nn")
                removed_total += len(original) - len(trimmed)
                if original != trimmed:
                    record_changed = True
                if trimmed:
                    seq_lines[idx] = trimmed
                    break
                seq_lines.pop()
                idx -= 1
            if record_changed:
                changed_total += 1
            out_handle.write(header + "\n")
            for seq_line in seq_lines:
                out_handle.write(seq_line + "\n")

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
