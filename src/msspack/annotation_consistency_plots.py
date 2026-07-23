from __future__ import annotations

import csv
from pathlib import Path
from xml.sax.saxutils import escape

from .chart_primitives import (
    CHART_FONT_SIZE_PT,
    GRID_RGB,
    MUTED_RGB,
    SVG_FONT_SIZE,
    TEXT_RGB,
    hex_to_rgb,
    pdf_helvetica_text_width,
    pdf_text_command,
    pdf_top_to_bottom,
    write_single_page_pdf,
)
from .utils import MSSPackError, atomic_text_writer, write_text

CONSISTENCY_COLORS = {
    "exact": "#1d4ed8",
    "canonical_equivalent": "#60a5fa",
    "compatible_granularity": "#dbeafe",
    "review": "#d97706",
}

CONSISTENCY_LABELS = {
    "exact": "Exact",
    "canonical_equivalent": "Canonical equivalent",
    "compatible_granularity": "Compatible granularity",
    "review": "Name conflict",
}

TIER_LABELS = {
    "near_identical": "Near-identical",
    "family": "Close family peer",
    "broad": "Broad homolog",
}

TIER_ORDER = {key: index for index, key in enumerate(TIER_LABELS)}
SOURCE_ORDER = {
    "existing": 0,
    "reference": 1,
    "swissprot": 2,
    "uniref90": 3,
    "pfam": 4,
    "cdd": 5,
    "none": 9,
}
SOURCE_DISPLAY_LABELS = {
    "swissprot": "Swiss",
    "uniref90": "UniRef90",
    "pfam": "Pfam",
    "cdd": "CDD",
}

PDF_POINTS_PER_INCH = 72.0
COMPACT_FIGURE_WIDTH_IN = 3.6
COMPACT_FIGURE_WIDTH = COMPACT_FIGURE_WIDTH_IN * PDF_POINTS_PER_INCH


def _object_int(value: object) -> int:
    if isinstance(value, bool):
        raise TypeError("Boolean is not an integer count")
    if isinstance(value, int | float | str):
        return int(value)
    raise TypeError(f"Expected an integer-compatible value, got {type(value).__name__}")


def _object_float(value: object) -> float:
    if isinstance(value, bool):
        raise TypeError("Boolean is not a floating-point measurement")
    if isinstance(value, int | float | str):
        return float(value)
    raise TypeError(f"Expected a float-compatible value, got {type(value).__name__}")


def _read_tsv(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise MSSPackError(
                f"Plot input is missing required columns {sorted(required)}: {path}"
            )
        return list(reader)


def _write_rows(path: Path, fieldnames: tuple[str, ...], rows: list[dict[str, object]]) -> Path:
    with atomic_text_writer(path) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    return path


def prepare_name_consistency_plot_data(
    summary_path: Path,
    output_path: Path,
) -> list[dict[str, object]]:
    raw_rows = _read_tsv(
        summary_path,
        {
            "tier",
            "threshold",
            "evaluable_clusters",
            "exact_clusters",
            "canonical_equivalent_clusters",
            "compatible_granularity_clusters",
            "review_clusters",
            "strict_equivalence_rate",
            "biological_compatibility_rate",
        },
    )
    rows: list[dict[str, object]] = []
    for row in sorted(raw_rows, key=lambda item: TIER_ORDER.get(item["tier"], 99)):
        evaluable = int(row["evaluable_clusters"])
        values = {
            "exact": int(row["exact_clusters"]),
            "canonical_equivalent": int(row["canonical_equivalent_clusters"]),
            "compatible_granularity": int(row["compatible_granularity_clusters"]),
            "review": int(row["review_clusters"]),
        }
        rows.append(
            {
                "tier": row["tier"],
                "label": TIER_LABELS.get(row["tier"], row["tier"]),
                "threshold": row["threshold"],
                "evaluable_clusters": evaluable,
                **{f"{key}_clusters": value for key, value in values.items()},
                **{
                    f"{key}_fraction": value / evaluable if evaluable else 0.0
                    for key, value in values.items()
                },
                "strict_equivalence_rate": float(row["strict_equivalence_rate"]),
                "biological_compatibility_rate": float(
                    row["biological_compatibility_rate"]
                ),
            }
        )
    fields = (
        "tier",
        "label",
        "threshold",
        "evaluable_clusters",
        "exact_clusters",
        "canonical_equivalent_clusters",
        "compatible_granularity_clusters",
        "review_clusters",
        "exact_fraction",
        "canonical_equivalent_fraction",
        "compatible_granularity_fraction",
        "review_fraction",
        "strict_equivalence_rate",
        "biological_compatibility_rate",
    )
    _write_rows(output_path, fields, rows)
    return rows


def load_name_consistency_plot_data(path: Path) -> list[dict[str, object]]:
    raw_rows = _read_tsv(
        path,
        {
            "tier",
            "label",
            "evaluable_clusters",
            "exact_clusters",
            "canonical_equivalent_clusters",
            "compatible_granularity_clusters",
            "review_clusters",
            "exact_fraction",
            "canonical_equivalent_fraction",
            "compatible_granularity_fraction",
            "review_fraction",
        },
    )
    rows: list[dict[str, object]] = []
    for row in raw_rows:
        converted: dict[str, object] = dict(row)
        converted["evaluable_clusters"] = int(row["evaluable_clusters"])
        for key in CONSISTENCY_LABELS:
            converted[f"{key}_clusters"] = int(row[f"{key}_clusters"])
            converted[f"{key}_fraction"] = float(row[f"{key}_fraction"])
        rows.append(converted)
    return rows


def prepare_source_consistency_plot_data(
    source_pair_path: Path,
    output_path: Path,
) -> list[dict[str, object]]:
    raw_rows = _read_tsv(
        source_pair_path,
        {
            "tier",
            "source_a",
            "source_b",
            "annotated_direct_pairs",
            "review_pairs",
            "review_rate",
            "meets_minimum_pairs",
        },
    )
    rows: list[dict[str, object]] = []
    for row in raw_rows:
        rows.append(
            {
                "tier": row["tier"],
                "source_a": row["source_a"],
                "source_b": row["source_b"],
                "annotated_direct_pairs": int(row["annotated_direct_pairs"]),
                "review_pairs": int(row["review_pairs"]),
                "review_rate": float(row["review_rate"]),
                "meets_minimum_pairs": row["meets_minimum_pairs"],
            }
        )
    rows.sort(
        key=lambda row: (
            SOURCE_ORDER.get(str(row["source_a"]).casefold(), 8),
            SOURCE_ORDER.get(str(row["source_b"]).casefold(), 8),
            str(row["source_a"]).casefold(),
            str(row["source_b"]).casefold(),
        )
    )
    _write_rows(
        output_path,
        (
            "tier",
            "source_a",
            "source_b",
            "annotated_direct_pairs",
            "review_pairs",
            "review_rate",
            "meets_minimum_pairs",
        ),
        rows,
    )
    return rows


def load_source_consistency_plot_data(path: Path) -> list[dict[str, object]]:
    raw_rows = _read_tsv(
        path,
        {
            "tier",
            "source_a",
            "source_b",
            "annotated_direct_pairs",
            "review_pairs",
            "review_rate",
            "meets_minimum_pairs",
        },
    )
    rows: list[dict[str, object]] = []
    for row in raw_rows:
        rows.append(
            {
                **row,
                "annotated_direct_pairs": int(row["annotated_direct_pairs"]),
                "review_pairs": int(row["review_pairs"]),
                "review_rate": float(row["review_rate"]),
            }
        )
    return rows


def _svg_document(width: float, height: float, body: str) -> str:
    width_in = width / PDF_POINTS_PER_INCH
    height_in = height / PDF_POINTS_PER_INCH
    return "\n".join(
        (
            '<?xml version="1.0" encoding="UTF-8"?>',
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width_in:g}in" '
            f'height="{height_in:g}in" viewBox="0 0 {width:.2f} {height:.2f}">',
            "<style>"
            f"text{{font-family:Helvetica,Arial,sans-serif;font-size:{SVG_FONT_SIZE};fill:#111827}}"
            ".title{font-weight:700}.subtitle{fill:#4b5563}"
            ".axis{fill:#4b5563}.label{font-weight:700}"
            "</style>",
            '<rect width="100%" height="100%" fill="#ffffff"/>',
            body,
            "</svg>",
            "",
        )
    )


def _threshold_lines(value: object) -> tuple[str, str]:
    text = str(value)
    identity = "?"
    coverage = "?"
    for component in text.split(";"):
        key, _, threshold = component.partition(">=")
        if not threshold:
            continue
        if key.strip() == "identity":
            identity = threshold.strip()
        elif key.strip() == "mutual coverage":
            coverage = threshold.strip()
    return f"id >={identity}", f"mutual cov >={coverage}"


def _source_display_label(source: str) -> str:
    return SOURCE_DISPLAY_LABELS.get(source.casefold(), source)


def _legend_positions() -> tuple[tuple[str, float, float], ...]:
    return (
        ("exact", 10.0, 48.0),
        ("canonical_equivalent", 112.0, 48.0),
        ("compatible_granularity", 10.0, 62.0),
        ("review", 153.0, 62.0),
    )


def write_name_consistency_svg(rows: list[dict[str, object]], output_path: Path) -> Path:
    width = COMPACT_FIGURE_WIDTH
    height = 260.0
    left = 96.0
    right = 48.0
    plot_width = width - left - right
    bar_height = 30.0
    bar_tops = [108.0 + index * 54.0 for index in range(len(rows))]
    parts = [
        '<text x="10" y="15" class="title">Functional annotation name consistency</text>',
        '<text x="10" y="29" class="subtitle">Evaluable clusters; each bar totals 100%</text>',
        '<text x="10" y="41" class="subtitle">Labels: identity / mutual-coverage thresholds</text>',
    ]
    for key, legend_x, legend_y in _legend_positions():
        parts.append(
            f'<rect x="{legend_x:.1f}" y="{legend_y:.1f}" width="7" height="7" '
            f'fill="{CONSISTENCY_COLORS[key]}" stroke="#334155" stroke-width="0.5"/>'
        )
        parts.append(
            f'<text x="{legend_x + 10:.1f}" y="{legend_y + 7:.1f}">'
            f'{escape(CONSISTENCY_LABELS[key])}</text>'
        )
    for tick in (0, 25, 50, 75, 100):
        x = left + plot_width * tick / 100.0
        parts.append(
            f'<line x1="{x:.2f}" y1="94" x2="{x:.2f}" y2="{height - 8:.2f}" '
            'stroke="#e5e7eb" stroke-width="0.8"/>'
        )
        parts.append(
            f'<text x="{x:.2f}" y="90" text-anchor="middle" class="axis">{tick}%</text>'
        )
    for row, top in zip(rows, bar_tops, strict=False):
        label = str(row["label"])
        denominator = _object_int(row["evaluable_clusters"])
        parts.append(
            f'<text x="{left - 5:.2f}" y="{top + 7:.2f}" text-anchor="end" '
            f'class="label">{escape(label)}</text>'
        )
        threshold_line1, threshold_line2 = _threshold_lines(row.get("threshold", ""))
        parts.append(
            f'<text x="{left - 5:.2f}" y="{top + 18:.2f}" text-anchor="end" '
            f'class="axis">{escape(threshold_line1)}</text>'
        )
        parts.append(
            f'<text x="{left - 5:.2f}" y="{top + 29:.2f}" text-anchor="end" '
            f'class="axis">{escape(threshold_line2)}</text>'
        )
        cursor = left
        for key in CONSISTENCY_LABELS:
            fraction = _object_float(row[f"{key}_fraction"])
            count = _object_int(row[f"{key}_clusters"])
            segment_width = plot_width * fraction
            parts.append(
                f'<rect x="{cursor:.2f}" y="{top:.2f}" width="{segment_width:.2f}" '
                f'height="{bar_height:.2f}" fill="{CONSISTENCY_COLORS[key]}" '
                'stroke="#ffffff" stroke-width="1"/>'
            )
            if segment_width >= 23.0 and count > 0:
                text_color = "#ffffff" if key in {"exact", "review"} else "#111827"
                parts.append(
                    f'<text x="{cursor + segment_width / 2:.2f}" y="{top + 19:.2f}" '
                    f'text-anchor="middle" class="value" style="fill:{text_color}">{count:,}</text>'
                )
            cursor += segment_width
        parts.append(
            f'<rect x="{left:.2f}" y="{top:.2f}" width="{plot_width:.2f}" '
            f'height="{bar_height:.2f}" fill="none" stroke="#334155" stroke-width="0.7"/>'
        )
        parts.append(
            f'<text x="{width - 2.0:.2f}" y="{top + 19:.2f}" '
            f'text-anchor="end" class="axis">n={denominator:,}</text>'
        )
    return write_text(output_path, _svg_document(width, height, "\n".join(parts)))


def _pdf_rect(
    *,
    page_height: float,
    x: float,
    y_top: float,
    width: float,
    height: float,
    color: tuple[float, float, float],
    stroke: tuple[float, float, float] | None = None,
) -> str:
    y = pdf_top_to_bottom(page_height, y_top, height)
    command = (
        f"{color[0]:.3f} {color[1]:.3f} {color[2]:.3f} rg "
        f"{x:.2f} {y:.2f} {width:.2f} {height:.2f} re f"
    )
    if stroke is not None:
        command += (
            f" {stroke[0]:.3f} {stroke[1]:.3f} {stroke[2]:.3f} RG 0.5 w "
            f"{x:.2f} {y:.2f} {width:.2f} {height:.2f} re S"
        )
    return command


def write_name_consistency_pdf(rows: list[dict[str, object]], output_path: Path) -> Path:
    width = COMPACT_FIGURE_WIDTH
    height = 260.0
    left = 96.0
    right = 48.0
    plot_width = width - left - right
    commands = [
        pdf_text_command(
            page_height=height,
            x=10,
            y_top=15,
            text="Functional annotation name consistency",
            font="F2",
            size=CHART_FONT_SIZE_PT,
            color=TEXT_RGB,
        ),
        pdf_text_command(
            page_height=height,
            x=10,
            y_top=29,
            text="Evaluable clusters; each bar totals 100%",
            font="F1",
            size=CHART_FONT_SIZE_PT,
            color=MUTED_RGB,
        ),
        pdf_text_command(
            page_height=height,
            x=10,
            y_top=41,
            text="Labels: identity / mutual-coverage thresholds",
            font="F1",
            size=CHART_FONT_SIZE_PT,
            color=MUTED_RGB,
        ),
    ]
    for key, legend_x, legend_y in _legend_positions():
        label = CONSISTENCY_LABELS[key]
        commands.append(
            _pdf_rect(
                page_height=height,
                x=legend_x,
                y_top=legend_y,
                width=7,
                height=7,
                color=hex_to_rgb(CONSISTENCY_COLORS[key]),
                stroke=TEXT_RGB,
            )
        )
        commands.append(
            pdf_text_command(
                page_height=height,
                x=legend_x + 10,
                y_top=legend_y + 7,
                text=label,
                font="F1",
                size=CHART_FONT_SIZE_PT,
                color=TEXT_RGB,
            )
        )
    for tick in (0, 25, 50, 75, 100):
        x = left + plot_width * tick / 100.0
        y1 = pdf_top_to_bottom(height, 94)
        y2 = pdf_top_to_bottom(height, height - 8)
        commands.append(
            f"{GRID_RGB[0]:.3f} {GRID_RGB[1]:.3f} {GRID_RGB[2]:.3f} RG 0.5 w "
            f"{x:.2f} {y1:.2f} m {x:.2f} {y2:.2f} l S"
        )
        tick_label = f"{tick}%"
        commands.append(
            pdf_text_command(
                page_height=height,
                x=x
                - pdf_helvetica_text_width(
                    tick_label,
                    size=CHART_FONT_SIZE_PT,
                )
                / 2.0,
                y_top=90,
                text=tick_label,
                font="F1",
                size=CHART_FONT_SIZE_PT,
                color=MUTED_RGB,
            )
        )
    for index, row in enumerate(rows):
        top = 108.0 + index * 54.0
        label = str(row["label"])
        threshold_line1, threshold_line2 = _threshold_lines(row.get("threshold", ""))
        commands.append(
            pdf_text_command(
                page_height=height,
                x=left
                - 5
                - pdf_helvetica_text_width(
                    label,
                    size=CHART_FONT_SIZE_PT,
                    bold=True,
                ),
                y_top=top + 7,
                text=label,
                font="F2",
                size=CHART_FONT_SIZE_PT,
                color=TEXT_RGB,
            )
        )
        commands.append(
            pdf_text_command(
                page_height=height,
                x=left
                - 5
                - pdf_helvetica_text_width(
                    threshold_line1,
                    size=CHART_FONT_SIZE_PT,
                ),
                y_top=top + 18,
                text=threshold_line1,
                font="F1",
                size=CHART_FONT_SIZE_PT,
                color=MUTED_RGB,
            )
        )
        commands.append(
            pdf_text_command(
                page_height=height,
                x=left
                - 5
                - pdf_helvetica_text_width(
                    threshold_line2,
                    size=CHART_FONT_SIZE_PT,
                ),
                y_top=top + 29,
                text=threshold_line2,
                font="F1",
                size=CHART_FONT_SIZE_PT,
                color=MUTED_RGB,
            )
        )
        cursor = left
        for key in CONSISTENCY_LABELS:
            fraction = _object_float(row[f"{key}_fraction"])
            count = _object_int(row[f"{key}_clusters"])
            segment_width = plot_width * fraction
            commands.append(
                _pdf_rect(
                    page_height=height,
                    x=cursor,
                    y_top=top,
                    width=segment_width,
                    height=30,
                    color=hex_to_rgb(CONSISTENCY_COLORS[key]),
                )
            )
            if segment_width >= 23.0 and count > 0:
                count_label = f"{count:,}"
                text_color = (
                    (1.0, 1.0, 1.0)
                    if key in {"exact", "review"}
                    else TEXT_RGB
                )
                commands.append(
                    pdf_text_command(
                        page_height=height,
                        x=cursor
                        + segment_width / 2.0
                        - pdf_helvetica_text_width(
                            count_label,
                            size=CHART_FONT_SIZE_PT,
                        )
                        / 2.0,
                        y_top=top + 19,
                        text=count_label,
                        font="F1",
                        size=CHART_FONT_SIZE_PT,
                        color=text_color,
                    )
                )
            cursor += segment_width
        pdf_label = f"n={_object_int(row['evaluable_clusters']):,}"
        commands.append(
            pdf_text_command(
                page_height=height,
                x=width
                - 2.0
                - pdf_helvetica_text_width(
                    pdf_label,
                    size=CHART_FONT_SIZE_PT,
                ),
                y_top=top + 19,
                text=pdf_label,
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


def _heat_color(rate: float) -> str:
    rate = max(0.0, min(1.0, rate))
    start = (255, 247, 237)
    end = (194, 65, 12)
    rgb = tuple(round(start[index] + (end[index] - start[index]) * rate) for index in range(3))
    return "#" + "".join(f"{channel:02x}" for channel in rgb)


def _source_matrix(
    rows: list[dict[str, object]],
) -> tuple[list[str], dict[tuple[str, str], dict[str, object]]]:
    sources = sorted(
        {
            str(row[key])
            for row in rows
            for key in ("source_a", "source_b")
        },
        key=lambda source: (SOURCE_ORDER.get(source.casefold(), 8), source.casefold()),
    )
    matrix: dict[tuple[str, str], dict[str, object]] = {}
    for row in rows:
        left = str(row["source_a"])
        right = str(row["source_b"])
        matrix[(left, right)] = row
        matrix[(right, left)] = row
    return sources, matrix


def write_source_consistency_svg(rows: list[dict[str, object]], output_path: Path) -> Path:
    sources, matrix = _source_matrix(rows)
    width = COMPACT_FIGURE_WIDTH
    left = 72.0
    top = 92.0
    cell = min(40.0, (width - left - 12.0) / max(1, len(sources)))
    height = top + cell * len(sources) + 14.0
    parts = [
        '<text x="10" y="15" class="title">Near-identical name conflict rate</text>',
        '<text x="10" y="27" class="title">by evidence source</text>',
        '<text x="10" y="43" class="subtitle">Cells: conflicts / annotated pairs and conflict rate.</text>',
        '<text x="10" y="55" class="subtitle">grey cells: below the configured minimum n.</text>',
    ]
    for index, source in enumerate(sources):
        display_source = _source_display_label(source)
        x = left + index * cell + cell / 2
        y = top + index * cell + cell / 2
        parts.append(
            f'<text x="{x:.2f}" y="{top - 10:.2f}" text-anchor="middle" '
            f'class="label">{escape(display_source)}</text>'
        )
        parts.append(
            f'<text x="{left - 10:.2f}" y="{y + 3:.2f}" text-anchor="end" '
            f'class="label">{escape(display_source)}</text>'
        )
    for row_index, source_a in enumerate(sources):
        for column_index, source_b in enumerate(sources):
            x = left + column_index * cell
            y = top + row_index * cell
            data = matrix.get((source_a, source_b))
            if data is None:
                color = "#ffffff"
                line1, line2 = "-", ""
            else:
                total = _object_int(data["annotated_direct_pairs"])
                review = _object_int(data["review_pairs"])
                rate = _object_float(data["review_rate"])
                sufficient = str(data["meets_minimum_pairs"]).casefold() == "true"
                color = _heat_color(rate) if sufficient else "#e5e7eb"
                line1 = f"{review}/{total}"
                line2 = f"{100 * rate:.1f}%" if sufficient else "low n"
            parts.append(
                f'<rect x="{x:.2f}" y="{y:.2f}" width="{cell:.2f}" height="{cell:.2f}" '
                f'fill="{color}" stroke="#cbd5e1" stroke-width="0.8"/>'
            )
            parts.append(
                f'<text x="{x + cell / 2:.2f}" y="{y + cell / 2 - 2:.2f}" '
                f'text-anchor="middle" class="value">{line1}</text>'
            )
            if line2:
                parts.append(
                    f'<text x="{x + cell / 2:.2f}" y="{y + cell / 2 + 9:.2f}" '
                    f'text-anchor="middle" class="axis">{line2}</text>'
                )
    return write_text(output_path, _svg_document(width, height, "\n".join(parts)))


def write_source_consistency_pdf(rows: list[dict[str, object]], output_path: Path) -> Path:
    sources, matrix = _source_matrix(rows)
    width = COMPACT_FIGURE_WIDTH
    left = 72.0
    top = 92.0
    cell = min(40.0, (width - left - 12.0) / max(1, len(sources)))
    height = top + cell * len(sources) + 14.0
    commands = [
        pdf_text_command(
            page_height=height,
            x=10,
            y_top=15,
            text="Near-identical name conflict rate",
            font="F2",
            size=CHART_FONT_SIZE_PT,
            color=TEXT_RGB,
        ),
        pdf_text_command(
            page_height=height,
            x=10,
            y_top=27,
            text="by evidence source",
            font="F2",
            size=CHART_FONT_SIZE_PT,
            color=TEXT_RGB,
        ),
        pdf_text_command(
            page_height=height,
            x=10,
            y_top=43,
            text="Cells: conflicts / annotated pairs and conflict rate.",
            font="F1",
            size=CHART_FONT_SIZE_PT,
            color=MUTED_RGB,
        ),
        pdf_text_command(
            page_height=height,
            x=10,
            y_top=55,
            text="grey cells: below the configured minimum n.",
            font="F1",
            size=CHART_FONT_SIZE_PT,
            color=MUTED_RGB,
        ),
    ]
    for index, source in enumerate(sources):
        display_source = _source_display_label(source)
        column_center = left + index * cell + cell / 2.0
        source_width = pdf_helvetica_text_width(
            display_source,
            size=CHART_FONT_SIZE_PT,
            bold=True,
        )
        commands.append(
            pdf_text_command(
                page_height=height,
                x=column_center - source_width / 2.0,
                y_top=top - 10,
                text=display_source,
                font="F2",
                size=CHART_FONT_SIZE_PT,
                color=TEXT_RGB,
            )
        )
        commands.append(
            pdf_text_command(
                page_height=height,
                x=left
                - 10
                - pdf_helvetica_text_width(
                    display_source,
                    size=CHART_FONT_SIZE_PT,
                    bold=True,
                ),
                y_top=top + index * cell + cell / 2 + 3,
                text=display_source,
                font="F2",
                size=CHART_FONT_SIZE_PT,
                color=TEXT_RGB,
            )
        )
    for row_index, source_a in enumerate(sources):
        for column_index, source_b in enumerate(sources):
            x = left + column_index * cell
            y = top + row_index * cell
            data = matrix.get((source_a, source_b))
            if data is None:
                color = (1.0, 1.0, 1.0)
                line1, line2 = "-", ""
            else:
                total = _object_int(data["annotated_direct_pairs"])
                review = _object_int(data["review_pairs"])
                rate = _object_float(data["review_rate"])
                sufficient = str(data["meets_minimum_pairs"]).casefold() == "true"
                color = (
                    hex_to_rgb(_heat_color(rate))
                    if sufficient
                    else hex_to_rgb("#e5e7eb")
                )
                line1 = f"{review}/{total}"
                line2 = f"{100 * rate:.1f}%" if sufficient else "low n"
            commands.append(
                _pdf_rect(
                    page_height=height,
                    x=x,
                    y_top=y,
                    width=cell,
                    height=cell,
                    color=color,
                    stroke=hex_to_rgb("#cbd5e1"),
                )
            )
            line1_width = pdf_helvetica_text_width(
                line1,
                size=CHART_FONT_SIZE_PT,
            )
            commands.append(
                pdf_text_command(
                    page_height=height,
                    x=x + cell / 2 - line1_width / 2.0,
                    y_top=y + cell / 2 - 2,
                    text=line1,
                    font="F1",
                    size=CHART_FONT_SIZE_PT,
                    color=TEXT_RGB,
                )
            )
            if line2:
                commands.append(
                    pdf_text_command(
                        page_height=height,
                        x=x
                        + cell / 2
                        - pdf_helvetica_text_width(
                            line2,
                            size=CHART_FONT_SIZE_PT,
                        )
                        / 2.0,
                        y_top=y + cell / 2 + 9,
                        text=line2,
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


__all__ = [
    "load_name_consistency_plot_data",
    "load_source_consistency_plot_data",
    "prepare_name_consistency_plot_data",
    "prepare_source_consistency_plot_data",
    "write_name_consistency_pdf",
    "write_name_consistency_svg",
    "write_source_consistency_pdf",
    "write_source_consistency_svg",
]
