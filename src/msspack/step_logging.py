from __future__ import annotations

import json
from collections.abc import Iterable, Sequence
from datetime import datetime
from pathlib import Path

from .fasta import iter_fasta
from .utils import shell_join, write_text


def _format_count(value: int | None) -> str:
    if value is None:
        return "n/a"
    return f"{value:,}"


def _format_command(command: str | Sequence[str]) -> str:
    if isinstance(command, str):
        return command if command.startswith("$ ") else f"$ {command}"
    return f"$ {shell_join(command)}"


def _format_table_value(value: str) -> str:
    return value.replace("|", r"\|").replace("\n", "<br>")


def _summary_table_lines(
    *,
    step: str,
    status: str,
    started_at: datetime,
    finished_at: datetime,
    count_unit: str | None,
    input_total: int | None,
    changed_total: int | None,
    output_total: int | None,
) -> list[str]:
    rows = [
        ("Step", step),
        ("Status", status),
        ("Started at", started_at.isoformat(sep=" ")),
        ("Ended at", finished_at.isoformat(sep=" ")),
        ("Duration seconds", f"{(finished_at - started_at).total_seconds():.3f}"),
        ("Count unit", count_unit or "n/a"),
        ("Input total", _format_count(input_total)),
        ("Changed total", _format_count(changed_total)),
        ("Output total", _format_count(output_total)),
    ]
    lines = [
        "Summary table:",
        "| Metric | Value |",
        "| --- | --- |",
    ]
    lines.extend(
        f"| {_format_table_value(metric)} | {_format_table_value(value)} |"
        for metric, value in rows
    )
    return lines


def write_step_log(
    *,
    log_path: Path,
    command: str | Sequence[str],
    step: str,
    started_at: datetime,
    ended_at: datetime | None = None,
    status: str = "completed",
    count_unit: str | None = None,
    input_total: int | None = None,
    changed_total: int | None = None,
    output_total: int | None = None,
    details: Iterable[str] = (),
) -> Path:
    finished_at = ended_at or datetime.now()
    lines = [
        _format_command(command),
        "",
        *_summary_table_lines(
            step=step,
            status=status,
            started_at=started_at,
            finished_at=finished_at,
            count_unit=count_unit,
            input_total=input_total,
            changed_total=changed_total,
            output_total=output_total,
        ),
        "",
        f"Step: {step}",
        f"Status: {status}",
        f"Started at: {started_at.isoformat(sep=' ')}",
        f"Ended at: {finished_at.isoformat(sep=' ')}",
        f"Duration seconds: {(finished_at - started_at).total_seconds():.3f}",
    ]
    if count_unit:
        lines.append(f"Count unit: {count_unit}")
    lines.extend(
        [
            f"Input total: {_format_count(input_total)}",
            f"Changed total: {_format_count(changed_total)}",
            f"Output total: {_format_count(output_total)}",
        ]
    )
    detail_lines = [line for line in details if line]
    if detail_lines:
        lines.append("")
        lines.extend(detail_lines)
    lines.append("")
    write_text(log_path, "\n".join(lines))
    return log_path


def write_step_metrics(
    *,
    metrics_path: Path,
    step: str,
    status: str = "completed",
    count_unit: str | None = None,
    input_total: int | None = None,
    changed_total: int | None = None,
    output_total: int | None = None,
    details: dict[str, object] | None = None,
) -> Path:
    payload: dict[str, object] = {
        "schema_version": 1,
        "step": step,
        "status": status,
        "count_unit": count_unit or "n/a",
        "input_total": input_total,
        "changed_total": changed_total,
        "output_total": output_total,
        "details": details or {},
    }
    write_text(metrics_path, json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return metrics_path


def write_id_list(path: Path, identifiers: Iterable[str]) -> Path:
    ordered = [identifier for identifier in identifiers if identifier]
    write_text(path, "\n".join(ordered) + ("\n" if ordered else ""))
    return path


def count_fasta_records(path: str | Path) -> int:
    return sum(1 for _ in iter_fasta(path))


def count_gff_feature_records(path: str | Path) -> int:
    count = 0
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line == "##FASTA":
                break
            if line and not line.startswith("#"):
                count += 1
    return count


def count_gff_gene_records(path: str | Path) -> int:
    count = 0
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line == "##FASTA":
                break
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) >= 3 and fields[2] == "gene":
                count += 1
    return count


def count_tabular_rows(path: str | Path, *, has_header: bool = True) -> int:
    with Path(path).open("r", encoding="utf-8") as handle:
        rows = [line for line in handle.read().splitlines() if line.strip()]
    if has_header and rows:
        return max(0, len(rows) - 1)
    return len(rows)


def count_nonempty_lines(text: str) -> int:
    return sum(1 for line in text.splitlines() if line.strip())


def count_mss_feature_blocks(path: str | Path) -> int:
    count = 0
    with Path(path).open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.strip() and not line.startswith("\t\t\t"):
                count += 1
    return count


def count_reordered_feature_lines(
    before_path: str | Path,
    after_path: str | Path,
) -> int:
    def feature_lines(path: str | Path) -> list[str]:
        lines: list[str] = []
        with Path(path).open("r", encoding="utf-8") as handle:
            for raw_line in handle:
                line = raw_line.strip()
                if line == "##FASTA":
                    break
                if line and not line.startswith("#"):
                    lines.append(line)
        return lines

    before_lines = feature_lines(before_path)
    after_lines = feature_lines(after_path)
    changed = sum(
        1
        for before, after in zip(before_lines, after_lines, strict=False)
        if before != after
    )
    changed += abs(len(before_lines) - len(after_lines))
    return changed
