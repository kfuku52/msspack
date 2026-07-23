from __future__ import annotations

import json
import os
import re
from collections.abc import Mapping
from dataclasses import dataclass
from html import escape
from pathlib import Path
from typing import Any

from .config import load_config
from .pipeline_plots import run_pipeline_plots
from .utils import MSSPackError, ensure_dir, write_text


@dataclass(frozen=True)
class ReportArtifacts:
    root: Path
    index_html: Path


def _read_manifest(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise MSSPackError(f"Build manifest not found: {path}. Run 'msspack pack' first.")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise MSSPackError(f"Invalid build manifest: {path}")
    return payload


def _relative_href(report_root: Path, target: str | Path) -> str:
    target_path = Path(target)
    if not target_path.is_absolute():
        target_path = (report_root.parent / target_path).resolve()
    return os.path.relpath(target_path, report_root).replace(os.sep, "/")


def _format_int(value: object) -> str:
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, int):
        return f"{value:,}"
    if isinstance(value, float) and value.is_integer():
        return f"{int(value):,}"
    if value in (None, "", "n/a"):
        return "n/a"
    return escape(str(value))


def _format_seconds(value: object) -> str:
    if isinstance(value, int | float):
        return f"{value:.2f}s"
    return "n/a"


def _render_link_list(report_root: Path, entries: Mapping[str, object]) -> str:
    items: list[str] = []
    for label, raw_path in entries.items():
        if not raw_path:
            continue
        items.append(
            f'<li><a href="{escape(_relative_href(report_root, str(raw_path)))}">{escape(label)}</a></li>'
        )
    if not items:
        return "<p class='muted'>No files recorded.</p>"
    return "<ul>" + "".join(items) + "</ul>"


def _render_plot_block(report_root: Path, title: str, svg_path: str | None, pdf_path: str | None, tsv_path: str | None) -> str:
    if not svg_path and not pdf_path and not tsv_path:
        return ""
    links: list[str] = []
    if svg_path:
        links.append(f'<a href="{escape(_relative_href(report_root, svg_path))}">SVG</a>')
    if pdf_path:
        links.append(f'<a href="{escape(_relative_href(report_root, pdf_path))}">PDF</a>')
    if tsv_path:
        links.append(f'<a href="{escape(_relative_href(report_root, tsv_path))}">TSV</a>')
    image = ""
    if svg_path:
        image = f'<div class="figure"><img src="{escape(_relative_href(report_root, svg_path))}" alt="{escape(title)}"></div>'
    return (
        '<section class="card">'
        f"<h3>{escape(title)}</h3>"
        f"<p class='links'>{' | '.join(links)}</p>"
        f"{image}"
        "</section>"
    )


def _stage_sort_key(path: Path) -> tuple[int, str]:
    match = re.match(r"(?P<num>\d+)", path.name)
    if match is None:
        return (9999, path.name)
    return (int(match.group("num")), path.name)


def _load_stage_metrics(logs_dir: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for metrics_path in sorted(logs_dir.glob("*.metrics.json"), key=_stage_sort_key):
        payload = json.loads(metrics_path.read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            continue
        log_name = metrics_path.name.replace(".metrics.json", ".log")
        log_path = logs_dir / log_name
        details = payload.get("details")
        detail_preview = ""
        if isinstance(details, dict) and details:
            pairs = []
            for key in sorted(details)[:4]:
                value = details[key]
                if isinstance(value, list | dict):
                    continue
                pairs.append(f"{key}={value}")
            detail_preview = ", ".join(pairs)
        rows.append(
            {
                "stage_file": metrics_path.name.replace(".metrics.json", ""),
                "step": payload.get("step", ""),
                "status": payload.get("status", ""),
                "count_unit": payload.get("count_unit", ""),
                "input_total": payload.get("input_total"),
                "changed_total": payload.get("changed_total"),
                "output_total": payload.get("output_total"),
                "detail_preview": detail_preview,
                "metrics_path": metrics_path,
                "log_path": log_path if log_path.exists() else None,
            }
        )
    return rows


def _render_stage_table(report_root: Path, logs_dir: Path) -> str:
    rows = _load_stage_metrics(logs_dir)
    if not rows:
        return "<section class='card'><h2>Stage metrics</h2><p class='muted'>No stage metrics were found.</p></section>"
    body_rows: list[str] = []
    for row in rows:
        links = [f'<a href="{escape(_relative_href(report_root, str(row["metrics_path"])))}">metrics</a>']
        if row["log_path"] is not None:
            links.append(
                f'<a href="{escape(_relative_href(report_root, str(row["log_path"])))}">log</a>'
            )
        body_rows.append(
            "<tr>"
            f"<td>{escape(str(row['stage_file']))}</td>"
            f"<td>{escape(str(row['step']))}</td>"
            f"<td>{escape(str(row['count_unit']))}</td>"
            f"<td>{_format_int(row['input_total'])}</td>"
            f"<td>{_format_int(row['changed_total'])}</td>"
            f"<td>{_format_int(row['output_total'])}</td>"
            f"<td>{escape(str(row['status']))}</td>"
            f"<td>{escape(str(row['detail_preview']))}</td>"
            f"<td>{' | '.join(links)}</td>"
            "</tr>"
        )
    return (
        "<section class='card'>"
        "<h2>Stage metrics</h2>"
        "<table>"
        "<thead><tr><th>Stage</th><th>Step</th><th>Unit</th><th>Input</th><th>Changed</th><th>Output</th><th>Status</th><th>Details</th><th>Files</th></tr></thead>"
        f"<tbody>{''.join(body_rows)}</tbody>"
        "</table>"
        "</section>"
    )


def _render_busco_section(report_root: Path, payload: dict[str, Any]) -> str:
    busco = payload.get("busco")
    if not isinstance(busco, dict) or not busco.get("enabled"):
        return ""
    comparisons = busco.get("comparisons")
    if not isinstance(comparisons, dict):
        return ""
    blocks: list[str] = []
    taxonomy_crosscheck = busco.get("taxonomy_crosscheck")
    if taxonomy_crosscheck:
        blocks.append(
            "<section class='card'>"
            "<h3>BUSCO taxonomy cross-check</h3>"
            f"<p class='links'><a href='{escape(_relative_href(report_root, str(taxonomy_crosscheck)))}'>JSON</a></p>"
            "</section>"
        )
    for name in ("cds", "genome"):
        comparison = comparisons.get(name)
        if not isinstance(comparison, dict):
            continue
        metrics = comparison.get("metrics")
        counts = comparison.get("counts")
        table_html = ""
        if isinstance(metrics, dict):
            input_metrics = metrics.get("input", {})
            processed_metrics = metrics.get("processed", {})
            delta_metrics = metrics.get("delta", {})
            rows = []
            for metric_name in ("complete", "single_copy", "duplicated", "fragmented", "missing"):
                rows.append(
                    "<tr>"
                    f"<td>{escape(metric_name)}</td>"
                    f"<td>{_format_int(input_metrics.get(metric_name))}</td>"
                    f"<td>{_format_int(processed_metrics.get(metric_name))}</td>"
                    f"<td>{_format_int(delta_metrics.get(metric_name))}</td>"
                    "</tr>"
                )
            table_html = (
                "<table>"
                "<thead><tr><th>Metric</th><th>Input</th><th>Processed</th><th>Delta</th></tr></thead>"
                f"<tbody>{''.join(rows)}</tbody>"
                "</table>"
            )
        count_html = ""
        if isinstance(counts, dict):
            count_html = (
                "<p class='muted'>"
                f"Counts: input={escape(json.dumps(counts.get('input', {}), sort_keys=True))}, "
                f"processed={escape(json.dumps(counts.get('processed', {}), sort_keys=True))}"
                "</p>"
            )
        blocks.append(
            "<section class='card'>"
            f"<h3>BUSCO {escape(name)}</h3>"
            f"<p>Lineage: <strong>{escape(str(comparison.get('lineage_dataset', 'n/a')))}</strong> | Mode: {escape(str(comparison.get('mode', 'n/a')))}</p>"
            f"<p class='links'><a href='{escape(_relative_href(report_root, str(comparison.get('comparison_tsv', ''))))}'>TSV</a> | "
            f"<a href='{escape(_relative_href(report_root, str(comparison.get('comparison_svg', ''))))}'>SVG</a> | "
            f"<a href='{escape(_relative_href(report_root, str(comparison.get('comparison_pdf', ''))))}'>PDF</a></p>"
            f"<div class='figure'><img src='{escape(_relative_href(report_root, str(comparison.get('comparison_svg', ''))))}' alt='BUSCO {escape(name)}'></div>"
            f"{table_html}"
            f"{count_html}"
            "</section>"
        )
    if not blocks:
        return ""
    return "<section><h2>BUSCO</h2>" + "".join(blocks) + "</section>"


def _render_validation_section(report_root: Path, payload: dict[str, Any]) -> str:
    validation = payload.get("validation")
    if not isinstance(validation, dict):
        return ""
    options = validation.get("options")
    outputs = validation.get("outputs")
    parts = ["<section class='card'><h2>Validation</h2>"]
    parts.append(f"<p>Status: <strong>{'enabled' if validation.get('enabled') else 'disabled'}</strong></p>")
    if isinstance(options, dict):
        parts.append(
            "<p class='muted'>"
            f"parser={escape(str(options.get('run_parser', 'n/a')))}, "
            f"transChecker={escape(str(options.get('run_transchecker', 'n/a')))}, "
            f"parallel={escape(str(options.get('parallel', 'n/a')))}, "
            f"heap={escape(str(options.get('heap', 'n/a')))}"
            "</p>"
        )
    if isinstance(outputs, dict):
        parts.append(_render_link_list(report_root, outputs))
    parts.append("</section>")
    return "".join(parts)


def _render_summary_cards(report_root: Path, payload: dict[str, Any]) -> str:
    outputs = payload.get("outputs", {})
    stage_summary = payload.get("stage_summary", {})
    runtime = payload.get("runtime", {})
    cards = [
        ("Project", payload.get("project", "n/a")),
        ("Status", payload.get("status", "n/a")),
        ("Started", payload.get("started_at", "n/a")),
        ("Finished", payload.get("finished_at", "n/a")),
        ("Duration", _format_seconds(payload.get("duration_seconds"))),
        ("Stages", f"{_format_int(stage_summary.get('count'))} total / {_format_int(stage_summary.get('ran'))} ran / {_format_int(stage_summary.get('reused'))} reused"),
        ("Python", runtime.get("python_version", "n/a")),
    ]
    output_links = ""
    if isinstance(outputs, dict):
        link_targets = {}
        annotation = outputs.get("annotation")
        fasta = outputs.get("fasta")
        functional_taxonomy = outputs.get("functional_annotation_taxonomy")
        functional_evidence = outputs.get("functional_annotation_evidence")
        name_standardization = outputs.get(
            "functional_annotation_name_standardization"
        )
        domain_comparison = outputs.get("functional_domain_search_comparison")
        consistency_audit = outputs.get("functional_annotation_consistency")
        consistency_families = outputs.get("functional_annotation_families")
        consistency_review = outputs.get("functional_annotation_conflicts") or outputs.get(
            "functional_annotation_review"
        )
        if isinstance(annotation, dict) and annotation.get("path"):
            link_targets["annotation"] = str(annotation["path"])
        if isinstance(fasta, dict) and fasta.get("path"):
            link_targets["fasta"] = str(fasta["path"])
        if isinstance(functional_taxonomy, dict) and functional_taxonomy.get("path"):
            link_targets["functional annotation taxonomy"] = str(
                functional_taxonomy["path"]
            )
        if isinstance(functional_evidence, dict) and functional_evidence.get("path"):
            link_targets["functional annotation evidence"] = str(functional_evidence["path"])
        if isinstance(name_standardization, dict) and name_standardization.get("path"):
            link_targets["product-name standardization"] = str(
                name_standardization["path"]
            )
        if isinstance(domain_comparison, dict) and domain_comparison.get("path"):
            link_targets["Pfam/CDD search comparison"] = str(domain_comparison["path"])
        if isinstance(consistency_audit, dict) and consistency_audit.get("path"):
            link_targets["annotation name consistency"] = str(consistency_audit["path"])
        if isinstance(consistency_families, dict) and consistency_families.get("path"):
            link_targets["annotation family summary"] = str(consistency_families["path"])
        if isinstance(consistency_review, dict) and consistency_review.get("path"):
            link_targets["annotation conflict diagnostics"] = str(consistency_review["path"])
        output_links = _render_link_list(report_root, link_targets)
    body = "".join(
        f"<div class='summary-card'><h3>{escape(str(label))}</h3><p>{escape(str(value))}</p></div>"
        for label, value in cards
    )
    return "<section><h2>Run summary</h2><div class='summary-grid'>" + body + "</div>" + output_links + "</section>"


def _render_pipeline_plots(report_root: Path, payload: dict[str, Any]) -> str:
    plots = payload.get("plots")
    if not isinstance(plots, dict):
        return ""
    pipeline = plots.get("pipeline")
    if not isinstance(pipeline, dict):
        return ""
    blocks = [
        _render_plot_block(
            report_root,
            "Pipeline gene flow",
            pipeline.get("gene_flow_svg"),
            pipeline.get("gene_flow_pdf"),
            pipeline.get("gene_flow_tsv"),
        ),
        _render_plot_block(
            report_root,
            "Pipeline event counts",
            pipeline.get("event_counts_svg"),
            pipeline.get("event_counts_pdf"),
            pipeline.get("event_counts_tsv"),
        ),
    ]
    annotation_consistency = pipeline.get("annotation_consistency")
    if isinstance(annotation_consistency, dict):
        blocks.extend(
            [
                _render_plot_block(
                    report_root,
                    "Functional annotation name consistency",
                    annotation_consistency.get("name_consistency_svg"),
                    annotation_consistency.get("name_consistency_pdf"),
                    annotation_consistency.get("name_consistency_tsv"),
                ),
                _render_plot_block(
                    report_root,
                    "Name review rate by evidence source",
                    annotation_consistency.get("source_consistency_svg"),
                    annotation_consistency.get("source_consistency_pdf"),
                    annotation_consistency.get("source_consistency_tsv"),
                ),
            ]
        )
    return "<section><h2>Plots</h2>" + "".join(block for block in blocks if block) + "</section>"


def _render_html(report_root: Path, payload: dict[str, Any], logs_dir: Path) -> str:
    return "\n".join(
        [
            "<!DOCTYPE html>",
            "<html lang='en'>",
            "<head>",
            "<meta charset='utf-8'>",
            "<meta name='viewport' content='width=device-width, initial-scale=1'>",
            f"<title>{escape(str(payload.get('project', 'msspack')))} report</title>",
            "<style>",
            "body{font-family:Helvetica,Arial,sans-serif;margin:0;background:#f8fafc;color:#0f172a;line-height:1.5}",
            "main{max-width:1320px;margin:0 auto;padding:32px 24px 64px}",
            "h1,h2,h3{margin:0 0 12px}",
            "h1{font-size:32px} h2{font-size:22px;margin-top:28px} h3{font-size:18px}",
            "p{margin:0 0 12px} .muted{color:#475569}",
            ".summary-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:16px;margin-bottom:18px}",
            ".summary-card,.card{background:white;border:1px solid #dbe3ee;border-radius:14px;padding:18px;box-shadow:0 1px 2px rgba(15,23,42,.05)}",
            ".card{margin-bottom:18px}",
            "table{width:100%;border-collapse:collapse;font-size:14px;background:white}",
            "th,td{padding:10px 12px;border-bottom:1px solid #e2e8f0;text-align:left;vertical-align:top}",
            "th{background:#f1f5f9;font-size:13px;color:#334155}",
            ".links{font-size:14px;color:#1d4ed8}",
            ".figure{margin-top:12px;border:1px solid #e2e8f0;border-radius:10px;overflow:hidden;background:#fff}",
            ".figure img{display:block;max-width:100%;height:auto}",
            "ul{margin:8px 0 0 18px;padding:0}",
            "code{background:#e2e8f0;padding:2px 6px;border-radius:6px}",
            "section{margin-bottom:20px}",
            "</style>",
            "</head>",
            "<body>",
            "<main>",
            f"<h1>{escape(str(payload.get('project', 'msspack')))} packaging report</h1>",
            "<p class='muted'>This report combines build metadata, structured step metrics, validation outputs, pipeline plots, and BUSCO comparisons for one msspack run.</p>",
            _render_summary_cards(report_root, payload),
            _render_validation_section(report_root, payload),
            _render_pipeline_plots(report_root, payload),
            _render_busco_section(report_root, payload),
            _render_stage_table(report_root, logs_dir),
            "</main>",
            "</body>",
            "</html>",
        ]
    ) + "\n"


def _update_report_manifest(manifest_path: Path, artifacts: ReportArtifacts) -> None:
    payload = _read_manifest(manifest_path)
    payload["report"] = {
        "enabled": True,
        "root": str(artifacts.root),
        "index_html": str(artifacts.index_html),
    }
    write_text(manifest_path, json.dumps(payload, indent=2, sort_keys=True) + "\n")


def run_html_report(
    config_file: str | Path,
    *,
    force: bool = False,
) -> ReportArtifacts:
    config = load_config(config_file)
    output_root = config.output_dir
    logs_dir = output_root / "logs"
    manifest_path = output_root / "build-manifest.json"
    if not logs_dir.exists():
        raise MSSPackError(
            f"Build logs are missing under {logs_dir}. Run 'msspack pack --config {config_file}' first."
        )
    run_pipeline_plots(config_file, force=force)
    payload = _read_manifest(manifest_path)
    artifacts = ReportArtifacts(
        root=ensure_dir(output_root / "report"),
        index_html=output_root / "report" / "index.html",
    )
    html = _render_html(artifacts.root, payload, logs_dir)
    write_text(artifacts.index_html, html)
    _update_report_manifest(manifest_path, artifacts)
    return artifacts
