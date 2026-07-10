from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Optional, Sequence, TypedDict
from xml.sax.saxutils import escape

from .chart_primitives import (
    CHART_FONT_SIZE_PT,
    GRID_RGB,
    MUTED_RGB,
    SVG_FONT_SIZE,
    TEXT_RGB,
    write_single_page_pdf,
)
from .chart_primitives import (
    pdf_text_command as _pdf_text_command,
)
from .chart_primitives import (
    pdf_top_to_bottom as _pdf_top_to_bottom,
)
from .config import BuscoConfig, MSSPackConfig, load_config
from .config_validation import validate_busco_config
from .execution import module_origin, run_if_needed
from .padding_tools import write_spliced_cds_fasta
from .pipeline import PipelineOutputs, run_pipeline
from .utils import (
    MSSPackError,
    default_cache_dir,
    ensure_dir,
    link_or_copy,
    run_command,
    which,
    write_text,
)

SUMMARY_RE = re.compile(
    r"C:(?P<complete>\d+(?:\.\d+)?)%\[S:(?P<single>\d+(?:\.\d+)?)%,D:(?P<duplicated>\d+(?:\.\d+)?)%\],"
    r"F:(?P<fragmented>\d+(?:\.\d+)?)%,M:(?P<missing>\d+(?:\.\d+)?)%,n:(?P<total>\d+)"
)
COUNT_PATTERNS = {
    "complete": re.compile(r"^(?P<count>\d+)\s+Complete BUSCOs \(C\)$"),
    "single_copy": re.compile(r"^(?P<count>\d+)\s+Complete and single-copy BUSCOs \(S\)$"),
    "duplicated": re.compile(r"^(?P<count>\d+)\s+Complete and duplicated BUSCOs \(D\)$"),
    "fragmented": re.compile(r"^(?P<count>\d+)\s+Fragmented BUSCOs \(F\)$"),
    "missing": re.compile(r"^(?P<count>\d+)\s+Missing BUSCOs \(M\)$"),
}
SEGMENT_COLORS_HEX = {
    "single_copy": "#2ca25f",
    "duplicated": "#3b82f6",
    "fragmented": "#f59e0b",
    "missing": "#ef4444",
}
SEGMENT_COLORS_RGB = {
    key: tuple(int(value[index : index + 2], 16) / 255.0 for index in (1, 3, 5))
    for key, value in SEGMENT_COLORS_HEX.items()
}
PDF_POINTS_PER_INCH = 72.0
BUSCO_COMPARISON_WIDTH_IN = 3.6
BUSCO_COMPARISON_WIDTH_PT = BUSCO_COMPARISON_WIDTH_IN * PDF_POINTS_PER_INCH


@dataclass(frozen=True)
class BuscoSummary:
    label: str
    input_fasta: Path
    raw_output_dir: Path
    short_summary_path: Path
    mode: str
    lineage_dataset: str
    busco_version: str
    complete_pct: float
    single_copy_pct: float
    duplicated_pct: float
    fragmented_pct: float
    missing_pct: float
    total_buscos: int
    complete_count: Optional[int] = None
    single_copy_count: Optional[int] = None
    duplicated_count: Optional[int] = None
    fragmented_count: Optional[int] = None
    missing_count: Optional[int] = None
    selection_strategy: str = ""

    def to_dict(self) -> dict[str, object]:
        return {
            "label": self.label,
            "input_fasta": str(self.input_fasta),
            "raw_output_dir": str(self.raw_output_dir),
            "short_summary_path": str(self.short_summary_path),
            "mode": self.mode,
            "lineage_dataset": self.lineage_dataset,
            "busco_version": self.busco_version,
            "selection_strategy": self.selection_strategy,
            "percentages": {
                "complete": self.complete_pct,
                "single_copy": self.single_copy_pct,
                "duplicated": self.duplicated_pct,
                "fragmented": self.fragmented_pct,
                "missing": self.missing_pct,
            },
            "counts": {
                "complete": self.complete_count,
                "single_copy": self.single_copy_count,
                "duplicated": self.duplicated_count,
                "fragmented": self.fragmented_count,
                "missing": self.missing_count,
                "total_buscos": self.total_buscos,
            },
        }


class TickLayout(TypedDict):
    value: int
    x: float


class LegendItemLayout(TypedDict):
    label: str
    key: str


class SegmentLayout(TypedDict):
    key: str
    value: float
    x: float
    y: float
    width: float
    height: float


class RowLayout(TypedDict):
    summary: BuscoSummary
    y: float
    label_y: float
    summary_y: float
    segments: list[SegmentLayout]
    summary_lines: list[str]


class ComparisonChartLayout(TypedDict):
    comparison_name: str
    width: float
    height: float
    left: float
    label_x: float
    top: float
    bar_height: float
    bar_gap: float
    bar_width: float
    axis_bottom: float
    tick_label_y: float
    title_y: float
    subtitle_y: float
    legend_y: float
    legend_gap: float
    legend_marker_size: float
    ticks: list[TickLayout]
    legend: list[LegendItemLayout]
    rows: list[RowLayout]
    lineage_note: str


@dataclass(frozen=True)
class BuscoComparisonArtifacts:
    root: Path
    logs_dir: Path
    raw_root: Path
    input_fasta: Path
    processed_fasta: Path
    input_summary_json: Path
    processed_summary_json: Path
    comparison_json: Path
    comparison_tsv: Path
    comparison_svg: Path
    comparison_pdf: Path


@dataclass(frozen=True)
class BuscoArtifacts:
    root: Path
    cds: Optional[BuscoComparisonArtifacts]
    genome: Optional[BuscoComparisonArtifacts]


def _build_comparison_artifacts(root: Path, name: str) -> BuscoComparisonArtifacts:
    group_root = ensure_dir(root / name)
    if name == "genome":
        input_fasta = group_root / "input.genome.fasta"
        processed_fasta = group_root / "processed.genome.fasta"
    else:
        input_fasta = group_root / "input.cds.fasta"
        processed_fasta = group_root / "processed.cds.fasta"
    return BuscoComparisonArtifacts(
        root=group_root,
        logs_dir=ensure_dir(group_root / "logs"),
        raw_root=ensure_dir(group_root / "raw"),
        input_fasta=input_fasta,
        processed_fasta=processed_fasta,
        input_summary_json=group_root / "input.summary.json",
        processed_summary_json=group_root / "processed.summary.json",
        comparison_json=group_root / "comparison.json",
        comparison_tsv=group_root / "comparison.tsv",
        comparison_svg=group_root / "comparison.svg",
        comparison_pdf=group_root / "comparison.pdf",
    )


def _build_busco_artifacts(outputs: PipelineOutputs, busco: BuscoConfig) -> BuscoArtifacts:
    root = ensure_dir(outputs.root / "busco")
    return BuscoArtifacts(
        root=root,
        cds=_build_comparison_artifacts(root, "cds") if busco.run_cds else None,
        genome=_build_comparison_artifacts(root, "genome") if busco.run_genome else None,
    )


def _busco_scope_flag(scope: str) -> str:
    if scope == "all":
        return "--auto-lineage"
    if scope == "euk":
        return "--auto-lineage-euk"
    if scope == "prok":
        return "--auto-lineage-prok"
    raise MSSPackError(f"Unsupported BUSCO auto-lineage scope: {scope}")


def _build_busco_command(
    *,
    busco: BuscoConfig,
    input_fasta: Path,
    out_name: str,
    out_path: Path,
    lineage_dataset: str = "",
    use_auto_lineage: bool = False,
) -> list[str]:
    command = [
        busco.command,
        "-i",
        str(input_fasta),
        "-m",
        busco.mode,
        "-o",
        out_name,
        "--out_path",
        str(out_path),
        "-c",
        str(busco.threads),
    ]
    if lineage_dataset:
        command.extend(["-l", lineage_dataset])
    elif use_auto_lineage:
        command.append(_busco_scope_flag(busco.auto_lineage_scope))
    if busco.download_path:
        command.extend(["--download_path", busco.download_path])
    if busco.offline:
        command.append("--offline")
    if busco.force:
        command.append("--force")
    if busco.opt_out_run_stats:
        command.append("--opt-out-run-stats")
    return command


def _busco_env(command0: str) -> Optional[dict[str, str]]:
    command_path = Path(command0)
    if not command_path.is_absolute():
        return None
    env = dict(os.environ)
    current_path = env.get("PATH", "")
    env["PATH"] = str(command_path.parent) + (os.pathsep + current_path if current_path else "")
    return env


def _remove_busco_output_dir(path: Path) -> None:
    _remove_busco_path(path)


def _remove_busco_path(path: Path) -> None:
    if path.is_symlink():
        path.unlink()
    elif path.exists():
        shutil.rmtree(path)


def _safe_workspace_name(path: Path) -> str:
    digest = hashlib.sha1(str(path).encode("utf-8")).hexdigest()[:12]
    return f"{path.name}-{digest}"


def _prepare_busco_paths(input_fasta: Path, raw_dir: Path) -> tuple[Path, Path]:
    if " " not in str(input_fasta) and " " not in str(raw_dir):
        return input_fasta, raw_dir
    workspace_root = ensure_dir(busco_workspace_root())
    safe_raw_dir = workspace_root / _safe_workspace_name(raw_dir)
    _remove_busco_path(safe_raw_dir)
    ensure_dir(safe_raw_dir)
    safe_input_fasta = safe_raw_dir / input_fasta.name.replace(" ", "_")
    link_or_copy(input_fasta, safe_input_fasta)
    return safe_input_fasta, safe_raw_dir


def _publish_busco_workspace(work_dir: Path, raw_dir: Path) -> None:
    if work_dir == raw_dir:
        return
    ensure_dir(raw_dir.parent)
    _remove_busco_path(raw_dir)
    shutil.copytree(work_dir, raw_dir)
    _remove_busco_path(work_dir)


def busco_workspace_root() -> Path:
    return default_cache_dir() / "busco-work"


def cleanup_busco_cache() -> Optional[Path]:
    workspace_root = busco_workspace_root()
    if not workspace_root.exists():
        return None
    _remove_busco_path(workspace_root)
    return workspace_root


def _discover_short_summary(raw_dir: Path) -> Path:
    candidates = sorted(
        raw_dir.rglob("short_summary*.txt"),
        key=lambda path: (path.stat().st_mtime, str(path)),
        reverse=True,
    )
    if not candidates:
        raise MSSPackError(f"BUSCO short summary was not found under {raw_dir}")
    return candidates[0]


def _parse_prefixed_value(line: str, prefix: str) -> str:
    return line.split(prefix, 1)[1].strip()


def parse_short_summary(
    summary_path: Path,
    *,
    label: str,
    input_fasta: Path,
    raw_output_dir: Path,
    selection_strategy: str = "",
) -> BuscoSummary:
    busco_version = ""
    lineage_dataset = ""
    mode = ""
    counts: dict[str, int] = {}
    summary_match: Optional[re.Match[str]] = None
    for raw_line in summary_path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if "BUSCO version is:" in line:
            busco_version = _parse_prefixed_value(line, "BUSCO version is:")
            continue
        if "The lineage dataset is:" in line:
            lineage_text = _parse_prefixed_value(line, "The lineage dataset is:")
            lineage_dataset = lineage_text.split("(", 1)[0].strip()
            continue
        if "BUSCO was run in mode:" in line:
            mode = _parse_prefixed_value(line, "BUSCO was run in mode:")
            continue
        candidate = SUMMARY_RE.search(line)
        if candidate:
            summary_match = candidate
            continue
        for count_name, pattern in COUNT_PATTERNS.items():
            count_match = pattern.match(line)
            if count_match:
                counts[count_name] = int(count_match.group("count"))
                break
    if summary_match is None:
        raise MSSPackError(f"BUSCO summary metrics were not found in {summary_path}")
    if not lineage_dataset:
        raise MSSPackError(f"BUSCO lineage dataset was not found in {summary_path}")
    if not mode:
        raise MSSPackError(f"BUSCO mode was not found in {summary_path}")
    if not busco_version:
        busco_version = "unknown"
    return BuscoSummary(
        label=label,
        input_fasta=input_fasta,
        raw_output_dir=raw_output_dir,
        short_summary_path=summary_path,
        mode=mode,
        lineage_dataset=lineage_dataset,
        busco_version=busco_version,
        complete_pct=float(summary_match.group("complete")),
        single_copy_pct=float(summary_match.group("single")),
        duplicated_pct=float(summary_match.group("duplicated")),
        fragmented_pct=float(summary_match.group("fragmented")),
        missing_pct=float(summary_match.group("missing")),
        total_buscos=int(summary_match.group("total")),
        complete_count=counts.get("complete"),
        single_copy_count=counts.get("single_copy"),
        duplicated_count=counts.get("duplicated"),
        fragmented_count=counts.get("fragmented"),
        missing_count=counts.get("missing"),
        selection_strategy=selection_strategy,
    )


def _write_summary_json(summary: BuscoSummary, output_path: Path) -> Path:
    return write_text(output_path, json.dumps(summary.to_dict(), indent=2, sort_keys=True) + "\n")


def _read_summary_json(path: Path) -> BuscoSummary:
    data = json.loads(path.read_text(encoding="utf-8"))
    percentages = data["percentages"]
    counts = data["counts"]
    return BuscoSummary(
        label=str(data["label"]),
        input_fasta=Path(str(data["input_fasta"])),
        raw_output_dir=Path(str(data["raw_output_dir"])),
        short_summary_path=Path(str(data["short_summary_path"])),
        mode=str(data["mode"]),
        lineage_dataset=str(data["lineage_dataset"]),
        busco_version=str(data["busco_version"]),
        selection_strategy=str(data.get("selection_strategy", "")),
        complete_pct=float(percentages["complete"]),
        single_copy_pct=float(percentages["single_copy"]),
        duplicated_pct=float(percentages["duplicated"]),
        fragmented_pct=float(percentages["fragmented"]),
        missing_pct=float(percentages["missing"]),
        total_buscos=int(counts["total_buscos"]),
        complete_count=_optional_int(counts.get("complete")),
        single_copy_count=_optional_int(counts.get("single_copy")),
        duplicated_count=_optional_int(counts.get("duplicated")),
        fragmented_count=_optional_int(counts.get("fragmented")),
        missing_count=_optional_int(counts.get("missing")),
    )


def _optional_int(value: object) -> Optional[int]:
    if value is None:
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        stripped = value.strip()
        if not stripped:
            return None
        return int(stripped)
    raise MSSPackError(f"Expected BUSCO count to be an int-compatible value, got {type(value)!r}")


def _run_busco_once(
    *,
    busco: BuscoConfig,
    label: str,
    input_fasta: Path,
    summary_json_path: Path,
    raw_dir: Path,
    log_path: Path,
    dependencies: Sequence[Path],
    lineage_dataset: str = "",
    use_auto_lineage: bool = False,
    selection_strategy: str = "",
) -> BuscoSummary:
    if busco.force:
        if summary_json_path.exists():
            summary_json_path.unlink()
        _remove_busco_output_dir(raw_dir)
    ensure_dir(raw_dir)
    ran = run_if_needed(
        outputs=[summary_json_path],
        dependencies=list(dependencies),
        cache_key={
            "command": busco.command,
            "mode": busco.mode,
            "threads": busco.threads,
            "download_path": busco.download_path,
            "offline": busco.offline,
            "opt_out_run_stats": busco.opt_out_run_stats,
            "lineage_dataset": lineage_dataset,
            "use_auto_lineage": use_auto_lineage,
            "auto_lineage_scope": busco.auto_lineage_scope,
            "selection_strategy": selection_strategy,
        },
        action=lambda: _execute_and_capture_summary(
            busco=busco,
            label=label,
            input_fasta=input_fasta,
            summary_json_path=summary_json_path,
            raw_dir=raw_dir,
            log_path=log_path,
            lineage_dataset=lineage_dataset,
            use_auto_lineage=use_auto_lineage,
            selection_strategy=selection_strategy,
        ),
    )
    if not ran and not summary_json_path.exists():
        raise MSSPackError(f"BUSCO summary JSON was not created: {summary_json_path}")
    return _read_summary_json(summary_json_path)


def _execute_and_capture_summary(
    *,
    busco: BuscoConfig,
    label: str,
    input_fasta: Path,
    summary_json_path: Path,
    raw_dir: Path,
    log_path: Path,
    lineage_dataset: str,
    use_auto_lineage: bool,
    selection_strategy: str,
) -> None:
    _remove_busco_output_dir(raw_dir)
    safe_input_fasta, work_raw_dir = _prepare_busco_paths(input_fasta, raw_dir)
    out_name = "run"
    command = _build_busco_command(
        busco=busco,
        input_fasta=safe_input_fasta,
        out_name=out_name,
        out_path=work_raw_dir,
        lineage_dataset=lineage_dataset,
        use_auto_lineage=use_auto_lineage,
    )
    run_command(command, log_path=log_path, env=_busco_env(command[0]))
    _publish_busco_workspace(work_raw_dir, raw_dir)
    short_summary_path = _discover_short_summary(raw_dir)
    summary = parse_short_summary(
        short_summary_path,
        label=label,
        input_fasta=input_fasta,
        raw_output_dir=raw_dir,
        selection_strategy=selection_strategy,
    )
    _write_summary_json(summary, summary_json_path)


def _write_comparison_json(
    summaries: Sequence[BuscoSummary],
    output_path: Path,
    *,
    comparison_name: str,
) -> Path:
    payload = {
        "comparison_name": comparison_name,
        "labels": [summary.label for summary in summaries],
        "lineage_datasets": {summary.label: summary.lineage_dataset for summary in summaries},
        "mode": summaries[0].mode if summaries else "",
        "summaries": [summary.to_dict() for summary in summaries],
    }
    return write_text(output_path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _stringify_tsv_value(value: object) -> str:
    if value is None:
        return ""
    return str(value)


def _write_comparison_tsv(
    summaries: Sequence[BuscoSummary],
    output_path: Path,
    *,
    comparison_name: str,
) -> Path:
    lines = [
        f"# comparison_name\t{comparison_name}",
        "\t".join(
            [
                "label",
                "input_fasta",
                "lineage_dataset",
                "mode",
                "busco_version",
                "complete_pct",
                "single_copy_pct",
                "duplicated_pct",
                "fragmented_pct",
                "missing_pct",
                "total_buscos",
                "complete_count",
                "single_copy_count",
                "duplicated_count",
                "fragmented_count",
                "missing_count",
                "selection_strategy",
            ]
        )
    ]
    for summary in summaries:
        lines.append(
            "\t".join(
                _stringify_tsv_value(value)
                for value in (
                    summary.label,
                    summary.input_fasta,
                    summary.lineage_dataset,
                    summary.mode,
                    summary.busco_version,
                    summary.complete_pct,
                    summary.single_copy_pct,
                    summary.duplicated_pct,
                    summary.fragmented_pct,
                    summary.missing_pct,
                    summary.total_buscos,
                    summary.complete_count,
                    summary.single_copy_count,
                    summary.duplicated_count,
                    summary.fragmented_count,
                    summary.missing_count,
                    summary.selection_strategy,
                )
            )
        )
    return write_text(output_path, "\n".join(lines) + "\n")


def _comparison_chart_layout(
    summaries: Sequence[BuscoSummary],
    *,
    comparison_name: str,
) -> ComparisonChartLayout:
    if len(summaries) != 2:
        raise MSSPackError("BUSCO comparison plot expects exactly two summaries")
    width = BUSCO_COMPARISON_WIDTH_PT
    height = 187.2
    left = 70.0
    label_x = left - 4.0
    top = 64.0
    bar_height = 12.0
    bar_gap = 54.0
    right_margin = 12.0
    bar_width = width - left - right_margin
    axis_bottom = top + bar_gap + bar_height
    tick_label_y = axis_bottom + 13.0
    legend: list[LegendItemLayout] = [
        {"label": "Single", "key": "single_copy"},
        {"label": "Dup", "key": "duplicated"},
        {"label": "Frag", "key": "fragmented"},
        {"label": "Miss", "key": "missing"},
    ]
    lineage_datasets = {summary.lineage_dataset for summary in summaries}
    if len(lineage_datasets) == 1:
        lineage_note = summaries[0].lineage_dataset
    else:
        lineage_note = " | ".join(
            f"{summary.label}:{summary.lineage_dataset}" for summary in summaries
        )
    rows: list[RowLayout] = []
    for index, summary in enumerate(summaries):
        y = top + index * bar_gap
        cursor = float(left)
        segments: list[SegmentLayout] = []
        for key, value in (
            ("single_copy", summary.single_copy_pct),
            ("duplicated", summary.duplicated_pct),
            ("fragmented", summary.fragmented_pct),
            ("missing", summary.missing_pct),
        ):
            segment_width = bar_width * value / 100.0
            if segment_width <= 0:
                continue
            segments.append(
                {
                    "key": key,
                    "value": value,
                    "x": cursor,
                    "y": y,
                    "width": segment_width,
                    "height": bar_height,
                }
            )
            cursor += segment_width
        summary_y = y + bar_height + 10.0
        if index == len(summaries) - 1:
            summary_y = tick_label_y + 19.0
        summary_lines = [
            f"C{summary.complete_pct:.1f}% S{summary.single_copy_pct:.1f}% "
            f"D{summary.duplicated_pct:.1f}%",
            f"F{summary.fragmented_pct:.1f}% M{summary.missing_pct:.1f}% "
            f"n{summary.total_buscos}",
        ]
        rows.append(
            {
                "summary": summary,
                "y": y,
                "label_y": y + 9,
                "summary_y": summary_y,
                "segments": segments,
                "summary_lines": summary_lines,
            }
        )
    return {
        "comparison_name": comparison_name,
        "width": width,
        "height": height,
        "left": left,
        "label_x": label_x,
        "top": top,
        "bar_height": bar_height,
        "bar_gap": bar_gap,
        "bar_width": bar_width,
        "axis_bottom": axis_bottom,
        "tick_label_y": tick_label_y,
        "title_y": 14.0,
        "subtitle_y": 27.0,
        "legend_y": 42.0,
        "legend_gap": 46.0,
        "legend_marker_size": 5.5,
        "ticks": [
            {"value": tick, "x": left + (bar_width * tick / 100.0)}
            for tick in (0, 25, 50, 75, 100)
        ],
        "legend": legend,
        "rows": rows,
        "lineage_note": lineage_note,
    }


def _write_comparison_svg(
    summaries: Sequence[BuscoSummary],
    output_path: Path,
    *,
    comparison_name: str,
) -> Path:
    layout = _comparison_chart_layout(summaries, comparison_name=comparison_name)
    width = float(layout["width"])
    height = float(layout["height"])
    left = float(layout["left"])
    label_x = float(layout["label_x"])
    top = float(layout["top"])
    bar_width = float(layout["bar_width"])
    bar_height = float(layout["bar_height"])
    axis_bottom = float(layout["axis_bottom"])
    tick_label_y = float(layout["tick_label_y"])
    svg_width = f"{width / PDF_POINTS_PER_INCH:g}in"
    svg_height = f"{height / PDF_POINTS_PER_INCH:g}in"
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{svg_width}" height="{svg_height}" viewBox="0 0 {width:.2f} {height:.2f}">',
        f'<style>text{{font-family:Helvetica,Arial,sans-serif;fill:#111827}} .muted{{fill:#4b5563;font-size:{SVG_FONT_SIZE}}} .label{{font-size:{SVG_FONT_SIZE};font-weight:700}} .title{{font-size:{SVG_FONT_SIZE};font-weight:700}} .tick{{font-size:{SVG_FONT_SIZE};fill:#6b7280}}</style>',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="14" y="{float(layout["title_y"]):.1f}" class="title">BUSCO comparison: {escape(comparison_name)}</text>',
    ]
    parts.append(
        f'<text x="14" y="{float(layout["subtitle_y"]):.1f}" class="muted">{escape(summaries[0].mode)} | {escape(str(layout["lineage_note"]))}</text>'
    )
    for tick in layout["ticks"]:
        x = float(tick["x"])
        parts.append(
            f'<line x1="{x:.1f}" y1="{top:.1f}" x2="{x:.1f}" y2="{axis_bottom:.1f}" stroke="#e5e7eb" stroke-width="0.6"/>'
        )
        parts.append(
            f'<text x="{x:.1f}" y="{tick_label_y:.1f}" text-anchor="middle" class="tick">{tick["value"]}%</text>'
        )
    legend_x = left
    legend_y = float(layout["legend_y"])
    legend_gap = float(layout["legend_gap"])
    legend_marker_size = float(layout["legend_marker_size"])
    for index, item in enumerate(layout["legend"]):
        x = legend_x + index * legend_gap
        color = SEGMENT_COLORS_HEX[str(item["key"])]
        parts.append(
            f'<rect x="{x:.1f}" y="{legend_y:.1f}" width="{legend_marker_size:.1f}" height="{legend_marker_size:.1f}" rx="1" fill="{color}"/>'
        )
        parts.append(
            f'<text x="{x + legend_marker_size + 3:.1f}" y="{legend_y + legend_marker_size + 1.6:.1f}" class="muted">{escape(str(item["label"]))}</text>'
        )
    for row in layout["rows"]:
        summary = row["summary"]
        y = float(row["y"])
        parts.append(
            f'<text x="{label_x:.1f}" y="{float(row["label_y"]):.1f}" text-anchor="end" class="label">{escape(summary.label)}</text>'
        )
        for segment in row["segments"]:
            parts.append(
                f'<rect x="{float(segment["x"]):.2f}" y="{float(segment["y"]):.2f}" width="{float(segment["width"]):.2f}" height="{float(segment["height"]):.2f}" fill="{SEGMENT_COLORS_HEX[str(segment["key"])]}"/>'
            )
        parts.append(
            f'<rect x="{left:.1f}" y="{y:.1f}" width="{bar_width:.1f}" height="{bar_height:.1f}" fill="none" stroke="#111827" stroke-width="0.6"/>'
        )
        for line_index, line in enumerate(row["summary_lines"]):
            summary_y = float(row["summary_y"]) + line_index * 10.0
            parts.append(
                f'<text x="{left}" y="{summary_y:.1f}" class="muted">{escape(line)}</text>'
            )
    parts.append("</svg>")
    return write_text(output_path, "\n".join(parts) + "\n")


def _write_comparison_pdf(
    summaries: Sequence[BuscoSummary],
    output_path: Path,
    *,
    comparison_name: str,
) -> Path:
    layout = _comparison_chart_layout(summaries, comparison_name=comparison_name)
    width = float(layout["width"])
    height = float(layout["height"])
    left = float(layout["left"])
    label_x = float(layout["label_x"])
    top = float(layout["top"])
    bar_width = float(layout["bar_width"])
    bar_height = float(layout["bar_height"])
    axis_bottom = float(layout["axis_bottom"])
    tick_label_y = float(layout["tick_label_y"])
    commands = [
        f"1 1 1 rg 0 0 {width:.2f} {height:.2f} re f",
        _pdf_text_command(
            page_height=height,
            x=14,
            y_top=float(layout["title_y"]),
            text=f"BUSCO comparison: {comparison_name}",
            font="F2",
            size=CHART_FONT_SIZE_PT,
            color=TEXT_RGB,
        ),
        _pdf_text_command(
            page_height=height,
            x=14,
            y_top=float(layout["subtitle_y"]),
            text=f"{summaries[0].mode} | {str(layout['lineage_note'])}",
            font="F1",
            size=CHART_FONT_SIZE_PT,
            color=MUTED_RGB,
        ),
    ]
    for tick in layout["ticks"]:
        x = float(tick["x"])
        y1 = _pdf_top_to_bottom(height, top)
        y2 = _pdf_top_to_bottom(height, axis_bottom)
        commands.append(f"{GRID_RGB[0]:.3f} {GRID_RGB[1]:.3f} {GRID_RGB[2]:.3f} RG 0.6 w {x:.2f} {y1:.2f} m {x:.2f} {y2:.2f} l S")
        tick_text = f'{tick["value"]}%'
        tick_offset = 2.1 * len(tick_text)
        commands.append(
            _pdf_text_command(
                page_height=height,
                x=x - tick_offset,
                y_top=tick_label_y,
                text=tick_text,
                font="F1",
                size=CHART_FONT_SIZE_PT,
                color=MUTED_RGB,
            )
        )
    legend_y = float(layout["legend_y"])
    legend_gap = float(layout["legend_gap"])
    legend_marker_size = float(layout["legend_marker_size"])
    for index, item in enumerate(layout["legend"]):
        x = left + index * legend_gap
        r, g, b = SEGMENT_COLORS_RGB[str(item["key"])]
        rect_y = _pdf_top_to_bottom(height, legend_y, legend_marker_size)
        commands.append(
            f"{r:.3f} {g:.3f} {b:.3f} rg {x:.2f} {rect_y:.2f} {legend_marker_size:.2f} {legend_marker_size:.2f} re f"
        )
        commands.append(
            _pdf_text_command(
                page_height=height,
                x=x + legend_marker_size + 3,
                y_top=legend_y + legend_marker_size + 1.6,
                text=str(item["label"]),
                font="F1",
                size=CHART_FONT_SIZE_PT,
                color=MUTED_RGB,
            )
        )
    for row in layout["rows"]:
        summary = row["summary"]
        label_offset = CHART_FONT_SIZE_PT * 0.52 * len(summary.label)
        commands.append(
            _pdf_text_command(
                page_height=height,
                x=max(2.0, label_x - label_offset),
                y_top=float(row["label_y"]),
                text=summary.label,
                font="F2",
                size=CHART_FONT_SIZE_PT,
                color=TEXT_RGB,
            )
        )
        for segment in row["segments"]:
            r, g, b = SEGMENT_COLORS_RGB[str(segment["key"])]
            rect_y = _pdf_top_to_bottom(
                height,
                float(segment["y"]),
                float(segment["height"]),
            )
            commands.append(
                f"{r:.3f} {g:.3f} {b:.3f} rg {float(segment['x']):.2f} {rect_y:.2f} {float(segment['width']):.2f} {float(segment['height']):.2f} re f"
            )
        border_y = _pdf_top_to_bottom(height, float(row["y"]), bar_height)
        commands.append(
            f"{TEXT_RGB[0]:.3f} {TEXT_RGB[1]:.3f} {TEXT_RGB[2]:.3f} RG 0.6 w {left:.2f} {border_y:.2f} {bar_width:.2f} {bar_height:.2f} re S"
        )
        for line_index, line in enumerate(row["summary_lines"]):
            commands.append(
                _pdf_text_command(
                    page_height=height,
                    x=left,
                    y_top=float(row["summary_y"]) + line_index * 10.0,
                    text=line,
                    font="F1",
                    size=CHART_FONT_SIZE_PT,
                    color=MUTED_RGB,
                )
            )

    return write_single_page_pdf(
        width=width,
        height=height,
        commands=commands,
        output_path=output_path,
    )


def _stage_fastas(outputs: PipelineOutputs) -> tuple[Path, Path]:
    input_fasta = outputs.intermediate / "00.input.genome.fasta"
    processed_fasta = outputs.intermediate / "02.gap-normalized.genome.fasta"
    if not input_fasta.exists():
        raise MSSPackError(f"BUSCO input FASTA is missing: {input_fasta}")
    if not processed_fasta.exists():
        raise MSSPackError(f"BUSCO processed FASTA is missing: {processed_fasta}")
    return input_fasta, processed_fasta


def _stage_gffs(outputs: PipelineOutputs) -> tuple[Path, Path]:
    input_gff = outputs.intermediate / "00.input.annotations.gff"
    processed_gff = outputs.intermediate / "12.gff.final-sorted.gff"
    if not input_gff.exists():
        raise MSSPackError(f"BUSCO input GFF is missing: {input_gff}")
    if not processed_gff.exists():
        raise MSSPackError(f"BUSCO processed GFF is missing: {processed_gff}")
    return input_gff, processed_gff


def _resolved_busco_config(
    config: MSSPackConfig,
    *,
    lineage_dataset: str = "",
    threads: Optional[int] = None,
    force: bool = False,
    auto_lineage: Optional[bool] = None,
    auto_lineage_scope: str = "",
    run_genome: Optional[bool] = None,
    run_cds: Optional[bool] = None,
) -> BuscoConfig:
    busco = config.busco
    return BuscoConfig(
        command=busco.command,
        run_cds=busco.run_cds if run_cds is None else run_cds,
        run_genome=busco.run_genome if run_genome is None else run_genome,
        mode=busco.mode,
        cds_mode=busco.cds_mode,
        lineage_dataset=lineage_dataset or busco.lineage_dataset,
        auto_lineage=busco.auto_lineage if auto_lineage is None else auto_lineage,
        auto_lineage_scope=auto_lineage_scope or busco.auto_lineage_scope,
        threads=threads if threads is not None else busco.threads,
        download_path=busco.download_path,
        offline=busco.offline,
        force=force or busco.force,
        opt_out_run_stats=busco.opt_out_run_stats,
    )


def _extract_cds_fastas(
    *,
    artifacts: BuscoComparisonArtifacts,
    input_genome_fasta: Path,
    input_gff: Path,
    processed_genome_fasta: Path,
    processed_gff: Path,
    dependencies: Sequence[Path],
) -> tuple[Path, Path]:
    run_if_needed(
        outputs=[artifacts.input_fasta, artifacts.logs_dir / "extract-input-cds.log"],
        dependencies=[input_genome_fasta, input_gff, *dependencies],
        action=lambda: write_spliced_cds_fasta(
            fasta_path=input_genome_fasta,
            gff_path=input_gff,
            output_path=artifacts.input_fasta,
            log_path=artifacts.logs_dir / "extract-input-cds.log",
        ),
    )
    run_if_needed(
        outputs=[artifacts.processed_fasta, artifacts.logs_dir / "extract-processed-cds.log"],
        dependencies=[processed_genome_fasta, processed_gff, *dependencies],
        action=lambda: write_spliced_cds_fasta(
            fasta_path=processed_genome_fasta,
            gff_path=processed_gff,
            output_path=artifacts.processed_fasta,
            log_path=artifacts.logs_dir / "extract-processed-cds.log",
        ),
    )
    return artifacts.input_fasta, artifacts.processed_fasta


def _materialize_genome_fastas(
    *,
    artifacts: BuscoComparisonArtifacts,
    input_genome_fasta: Path,
    processed_genome_fasta: Path,
    dependencies: Sequence[Path],
) -> tuple[Path, Path]:
    run_if_needed(
        outputs=[artifacts.input_fasta],
        dependencies=[input_genome_fasta, *dependencies],
        action=lambda: link_or_copy(input_genome_fasta, artifacts.input_fasta),
    )
    run_if_needed(
        outputs=[artifacts.processed_fasta],
        dependencies=[processed_genome_fasta, *dependencies],
        action=lambda: link_or_copy(processed_genome_fasta, artifacts.processed_fasta),
    )
    return artifacts.input_fasta, artifacts.processed_fasta


def _write_group_outputs(
    artifacts: BuscoComparisonArtifacts,
    summaries: Sequence[BuscoSummary],
    *,
    comparison_name: str,
    dependency_paths: Sequence[Path],
) -> None:
    run_if_needed(
        outputs=[artifacts.comparison_json],
        dependencies=list(dependency_paths),
        action=lambda: _write_comparison_json(
            summaries,
            artifacts.comparison_json,
            comparison_name=comparison_name,
        ),
    )
    run_if_needed(
        outputs=[artifacts.comparison_tsv],
        dependencies=list(dependency_paths),
        action=lambda: _write_comparison_tsv(
            summaries,
            artifacts.comparison_tsv,
            comparison_name=comparison_name,
        ),
    )
    run_if_needed(
        outputs=[artifacts.comparison_svg],
        dependencies=list(dependency_paths),
        action=lambda: _write_comparison_svg(
            summaries,
            artifacts.comparison_svg,
            comparison_name=comparison_name,
        ),
    )
    run_if_needed(
        outputs=[artifacts.comparison_pdf],
        dependencies=list(dependency_paths),
        action=lambda: _write_comparison_pdf(
            summaries,
            artifacts.comparison_pdf,
            comparison_name=comparison_name,
        ),
    )


def _load_comparison_pair(
    artifacts: BuscoComparisonArtifacts,
) -> tuple[BuscoSummary, BuscoSummary]:
    return (
        _read_summary_json(artifacts.input_summary_json),
        _read_summary_json(artifacts.processed_summary_json),
    )


def _comparison_manifest_entry(artifacts: BuscoComparisonArtifacts) -> dict[str, object]:
    input_summary, processed_summary = _load_comparison_pair(artifacts)
    return {
        "input_fasta": str(artifacts.input_fasta),
        "processed_fasta": str(artifacts.processed_fasta),
        "input_summary_json": str(artifacts.input_summary_json),
        "processed_summary_json": str(artifacts.processed_summary_json),
        "comparison_json": str(artifacts.comparison_json),
        "comparison_tsv": str(artifacts.comparison_tsv),
        "comparison_svg": str(artifacts.comparison_svg),
        "comparison_pdf": str(artifacts.comparison_pdf),
        "logs": {
            "input": str(artifacts.logs_dir / "input.log"),
            "processed": str(artifacts.logs_dir / "processed.log"),
        },
        "lineage_dataset": processed_summary.lineage_dataset,
        "mode": processed_summary.mode,
        "busco_version": processed_summary.busco_version,
        "selection_strategy": {
            "input": input_summary.selection_strategy,
            "processed": processed_summary.selection_strategy,
        },
        "metrics": {
            "input": input_summary.to_dict()["percentages"],
            "processed": processed_summary.to_dict()["percentages"],
            "delta": {
                "complete": round(processed_summary.complete_pct - input_summary.complete_pct, 4),
                "single_copy": round(
                    processed_summary.single_copy_pct - input_summary.single_copy_pct,
                    4,
                ),
                "duplicated": round(
                    processed_summary.duplicated_pct - input_summary.duplicated_pct,
                    4,
                ),
                "fragmented": round(
                    processed_summary.fragmented_pct - input_summary.fragmented_pct,
                    4,
                ),
                "missing": round(processed_summary.missing_pct - input_summary.missing_pct, 4),
            },
        },
        "counts": {
            "input": input_summary.to_dict()["counts"],
            "processed": processed_summary.to_dict()["counts"],
        },
    }


def summarize_busco_artifacts(artifacts: BuscoArtifacts) -> list[str]:
    lines: list[str] = []
    for name, comparison in (("cds", artifacts.cds), ("genome", artifacts.genome)):
        if comparison is None:
            continue
        input_summary, processed_summary = _load_comparison_pair(comparison)
        lines.append(
            (
                f"{name}\tlineage={processed_summary.lineage_dataset}"
                f"\tcomplete={input_summary.complete_pct:.1f}->{processed_summary.complete_pct:.1f}"
                f"\tsingle_copy={input_summary.single_copy_pct:.1f}->{processed_summary.single_copy_pct:.1f}"
                f"\tduplicated={input_summary.duplicated_pct:.1f}->{processed_summary.duplicated_pct:.1f}"
                f"\tfragmented={input_summary.fragmented_pct:.1f}->{processed_summary.fragmented_pct:.1f}"
                f"\tmissing={input_summary.missing_pct:.1f}->{processed_summary.missing_pct:.1f}"
            )
        )
    return lines


def _update_busco_manifest(
    manifest_path: Path,
    *,
    busco: BuscoConfig,
    artifacts: BuscoArtifacts,
) -> None:
    payload: dict[str, object]
    if manifest_path.exists():
        payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    else:
        payload = {}
    comparisons: dict[str, object] = {}
    if artifacts.cds is not None:
        comparisons["cds"] = _comparison_manifest_entry(artifacts.cds)
    if artifacts.genome is not None:
        comparisons["genome"] = _comparison_manifest_entry(artifacts.genome)
    payload["busco"] = {
        "enabled": True,
        "command": busco.command,
        "run_cds": busco.run_cds,
        "run_genome": busco.run_genome,
        "mode": busco.mode,
        "cds_mode": busco.cds_mode,
        "lineage_dataset": busco.lineage_dataset,
        "auto_lineage": busco.auto_lineage,
        "auto_lineage_scope": busco.auto_lineage_scope,
        "threads": busco.threads,
        "comparisons": comparisons,
    }
    write_text(manifest_path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def run_busco_comparison(
    config_file: str | Path,
    *,
    lineage_dataset: str = "",
    threads: Optional[int] = None,
    force: bool = False,
    clean_cache: bool = False,
    auto_lineage: Optional[bool] = None,
    auto_lineage_scope: str = "",
    run_genome: Optional[bool] = None,
    run_cds: Optional[bool] = None,
) -> BuscoArtifacts:
    config = load_config(config_file)
    busco = _resolved_busco_config(
        config,
        lineage_dataset=lineage_dataset,
        threads=threads,
        force=force,
        auto_lineage=auto_lineage,
        auto_lineage_scope=auto_lineage_scope,
        run_genome=run_genome,
        run_cds=run_cds,
    )
    validate_busco_config(busco)
    if clean_cache:
        cleanup_busco_cache()
    outputs = run_pipeline(config_file, validate=False)
    artifacts = _build_busco_artifacts(outputs, busco)
    input_fasta, processed_fasta = _stage_fastas(outputs)
    input_gff, processed_gff = _stage_gffs(outputs)
    busco_module = module_origin("msspack.busco")
    padding_tools_module = module_origin("msspack.padding_tools")
    config_path = Path(config_file).expanduser().resolve()
    busco_binary = which(busco.command)
    busco_dependencies = [config_path, busco_module]
    if busco_binary:
        busco_dependencies.append(Path(busco_binary).resolve())

    shared_lineage = busco.lineage_dataset
    if artifacts.genome is not None:
        genome_input_materialized, genome_processed_materialized = _materialize_genome_fastas(
            artifacts=artifacts.genome,
            input_genome_fasta=input_fasta,
            processed_genome_fasta=processed_fasta,
            dependencies=[config_path, busco_module],
        )
        genome_input_summary = _run_busco_once(
            busco=busco,
            label="input",
            input_fasta=genome_input_materialized,
            summary_json_path=artifacts.genome.input_summary_json,
            raw_dir=artifacts.genome.raw_root / "input",
            log_path=artifacts.genome.logs_dir / "input.log",
            dependencies=[genome_input_materialized, *busco_dependencies],
            lineage_dataset=busco.lineage_dataset,
            use_auto_lineage=busco.auto_lineage and not busco.lineage_dataset,
            selection_strategy="configured-lineage"
            if busco.lineage_dataset
            else ("auto-lineage" if busco.auto_lineage else "configured-lineage"),
        )
        shared_lineage = busco.lineage_dataset or genome_input_summary.lineage_dataset
        genome_processed_summary = _run_busco_once(
            busco=busco,
            label="processed",
            input_fasta=genome_processed_materialized,
            summary_json_path=artifacts.genome.processed_summary_json,
            raw_dir=artifacts.genome.raw_root / "processed",
            log_path=artifacts.genome.logs_dir / "processed.log",
            dependencies=[
                genome_processed_materialized,
                artifacts.genome.input_summary_json,
                *busco_dependencies,
            ],
            lineage_dataset=shared_lineage,
            use_auto_lineage=False,
            selection_strategy="configured-lineage"
            if busco.lineage_dataset
            else "auto-lineage-from-genome-input",
        )
        _write_group_outputs(
            artifacts.genome,
            [genome_input_summary, genome_processed_summary],
            comparison_name="genome",
            dependency_paths=[
                artifacts.genome.input_summary_json,
                artifacts.genome.processed_summary_json,
                busco_module,
            ],
        )

    if artifacts.cds is not None:
        input_cds_fasta, processed_cds_fasta = _extract_cds_fastas(
            artifacts=artifacts.cds,
            input_genome_fasta=input_fasta,
            input_gff=input_gff,
            processed_genome_fasta=processed_fasta,
            processed_gff=processed_gff,
            dependencies=[config_path, padding_tools_module],
        )
        cds_busco = replace(busco, mode=busco.cds_mode)
        cds_input_summary = _run_busco_once(
            busco=cds_busco,
            label="input",
            input_fasta=input_cds_fasta,
            summary_json_path=artifacts.cds.input_summary_json,
            raw_dir=artifacts.cds.raw_root / "input",
            log_path=artifacts.cds.logs_dir / "input.log",
            dependencies=[input_cds_fasta, *busco_dependencies],
            lineage_dataset=shared_lineage,
            use_auto_lineage=busco.auto_lineage and not shared_lineage,
            selection_strategy="configured-lineage"
            if shared_lineage
            else ("auto-lineage" if busco.auto_lineage else "configured-lineage"),
        )
        shared_lineage = shared_lineage or cds_input_summary.lineage_dataset
        cds_processed_summary = _run_busco_once(
            busco=cds_busco,
            label="processed",
            input_fasta=processed_cds_fasta,
            summary_json_path=artifacts.cds.processed_summary_json,
            raw_dir=artifacts.cds.raw_root / "processed",
            log_path=artifacts.cds.logs_dir / "processed.log",
            dependencies=[
                processed_cds_fasta,
                artifacts.cds.input_summary_json,
                *busco_dependencies,
            ],
            lineage_dataset=shared_lineage,
            use_auto_lineage=False,
            selection_strategy="configured-lineage"
            if busco.lineage_dataset
            else (
                "auto-lineage-from-genome-input"
                if artifacts.genome is not None
                else "auto-lineage-from-cds-input"
            ),
        )
        _write_group_outputs(
            artifacts.cds,
            [cds_input_summary, cds_processed_summary],
            comparison_name="cds",
            dependency_paths=[
                artifacts.cds.input_summary_json,
                artifacts.cds.processed_summary_json,
                busco_module,
            ],
        )
    _update_busco_manifest(outputs.manifest_path, busco=busco, artifacts=artifacts)
    return artifacts
