"""Import external organelle classification/QC without predicting annotation."""
from __future__ import annotations

import csv
import hashlib
import html
import json
import os
from collections import Counter
from pathlib import Path
from typing import Any

from .chart_primitives import MUTED_RGB, TEXT_RGB, pdf_text_command, write_single_page_pdf
from .config import MSSPackConfig
from .config_models import OrganelleConfig
from .pipeline_plot_models import EventCount
from .pipeline_plot_render import write_event_counts_pdf, write_event_counts_svg
from .utils import MSSPackError, expand_path, write_text

STATES = {
    "complete": ("Complete models", "#059669"),
    "fragmented": ("Fragmented / review", "#d97706"),
    "not_detected": ("Not detected", "#dc2626"),
    "not_evaluated": ("Not assessed", "#94a3b8"),
}
COLORS = ("#059669", "#2563eb", "#7c3aed", "#d97706", "#0d9488", "#c026d3")


def _require(value: object, message: str) -> None:
    if not value:
        raise MSSPackError("Organelle assessment: " + message)


def _text(value: Any, label: str) -> str:
    _require(isinstance(value, str) and value.strip() and not any(c in value for c in "\x00\t\r\n"),
             "invalid " + label)
    return str(value)


def _validate_swissprot(data: Any, identities: set[tuple[str, str]]) -> None:
    if data is None:
        return
    _require(isinstance(data, dict), "invalid Swiss-Prot assessment")
    _text(data.get("method"), "Swiss-Prot method")
    database = data.get("database")
    _require(isinstance(database, dict), "Swiss-Prot database provenance is required")
    _require(database.get("name") == "Swiss-Prot", "Swiss-Prot database name must be Swiss-Prot")
    _text(database.get("version"), "Swiss-Prot database version")
    digest = _text(database.get("fasta_sha256"), "Swiss-Prot FASTA SHA-256")
    _require(len(digest) == 64 and all(c in "0123456789abcdef" for c in digest), "invalid Swiss-Prot FASTA SHA-256")
    _require(isinstance(data.get("thresholds"), dict) and bool(data["thresholds"]), "Swiss-Prot thresholds are required")
    results = data.get("results")
    _require(isinstance(results, list), "Swiss-Prot results must be a list")
    seen: set[tuple[str, str]] = set()
    for row in results:
        _require(isinstance(row, dict), "invalid Swiss-Prot result")
        key = (_text(row.get("entry"), "Swiss-Prot entry"), _text(row.get("locus_tag"), "Swiss-Prot locus tag"))
        _require(key in identities and key not in seen, "unknown or duplicate Swiss-Prot CDS: " + repr(key))
        seen.add(key)
        _require(row.get("status") in ("hit", "no_hit", "not_evaluated"), "invalid Swiss-Prot result status")
        if row["status"] == "hit":
            _text(row.get("accession"), "Swiss-Prot accession")
        else:
            _require(not row.get("accession"), "Swiss-Prot result without a hit has an accession")
        _text(row.get("evidence"), "Swiss-Prot evidence")
    _require(seen == identities, "Swiss-Prot results must cover every supplied CDS; use not_evaluated explicitly")


def _load_external(
    config: MSSPackConfig, item: OrganelleConfig, records: list[dict[str, Any]],
    identities: set[tuple[str, str]],
    partial: set[tuple[str, str]],
) -> dict[str, Any] | None:
    if not item.assessment_json:
        return None
    path = expand_path(item.assessment_json, config.base_dir)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise MSSPackError(f"Could not read organelle assessment: {path}") from exc
    _require(isinstance(data, dict) and type(data.get("schema_version")) is int and data["schema_version"] == 1, "unsupported schema")
    _require(data.get("name") == item.name and data.get("organelle") == item.organelle, "input identity differs")
    expected = sorted(hashlib.sha256(expand_path(source, config.base_dir).read_bytes()).hexdigest()
                      for source in (item.genbank, item.fasta, item.gff) if source)
    _require(data.get("input_sha256") == expected, "source annotation changed; regenerate external assessment")
    _require(data.get("sequence_sha256") == {r["input_id"]: r["sequence_sha256"] for r in records},
             "source sequences differ")
    _text(data.get("method"), "method")
    functional = data.get("functional_classification", {})
    _require(isinstance(functional, dict), "invalid functional classification")
    assignments = functional.get("assignments", [])
    _require(isinstance(assignments, list), "functional assignments must be a list")
    seen: set[tuple[str, str]] = set()
    for row in assignments:
        _require(isinstance(row, dict), "invalid functional assignment")
        key = (_text(row.get("entry"), "entry"), _text(row.get("locus_tag"), "locus tag"))
        _require(key in identities and key not in seen, "unknown or duplicate functional feature: " + repr(key))
        seen.add(key)
        _text(row.get("category"), "functional category")
        _text(row.get("evidence"), "functional evidence")
    content = data.get("gene_content")
    if content is not None:
        _require(isinstance(content, dict), "invalid gene-content section")
        panel = content.get("marker_set")
        _require(isinstance(panel, dict), "marker-set provenance is required")
        for field in ("id", "version", "definition"):
            _text(panel.get(field), "marker-set " + field)
        _require(isinstance(panel.get("references"), list) and bool(panel["references"]), "reference list is required")
        for reference in panel["references"]:
            _text(reference, "reference identifier")
        _require(content.get("scope") in {"supplied_annotations", "assembly_and_annotation"}, "invalid assessment scope")
        markers = content.get("markers")
        _require(isinstance(markers, list) and bool(markers), "empty marker set cannot yield a score")
        names: set[str] = set()
        for marker in markers:
            _require(isinstance(marker, dict), "invalid marker")
            name = _text(marker.get("gene"), "marker gene")
            _require(name.casefold() not in names, "duplicate marker: " + name)
            names.add(name.casefold())
            status = marker.get("status")
            _require(isinstance(status, str) and status in STATES, "invalid marker status")
            matches = marker.get("matches", [])
            _require(isinstance(matches, list), "invalid marker matches")
            for match in matches:
                _require(isinstance(match, dict), "invalid marker match")
                key = (_text(match.get("entry"), "match entry"), _text(match.get("locus_tag"), "match locus tag"))
                _require(key in identities, "marker refers to an unknown CDS")
            if status == "complete":
                _require(bool(matches), "complete marker requires a supplied CDS")
                _require(any((match["entry"], match["locus_tag"]) not in partial for match in matches),
                         "complete marker has only partial/pseudo CDS models")
            if status == "not_detected":
                _require(not matches, "not-detected marker has a matched CDS")
            _text(marker.get("evidence"), "marker evidence")
    _require(isinstance(data.get("limitations", []), list), "limitations must be a list")
    for note in data.get("limitations", []):
        _text(note, "limitation")
    _validate_swissprot(data.get("swissprot"), identities)
    return {**data, "assessment_file": str(path), "assessment_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def _status_chart(svg: Path, pdf: Path, title: str, message: str) -> None:
    write_text(svg, '<svg xmlns="http://www.w3.org/2000/svg" width="7.2in" height="1.2in" '
               'viewBox="0 0 518.4 86.4"><rect width="100%" height="100%" fill="white"/>'
               '<g font-family="Helvetica,Arial,sans-serif" font-size="10" fill="#475569">'
               '<text x="16" y="22">' + html.escape(title) + '</text><text x="16" y="48">'
               + html.escape(message) + '</text></g></svg>')
    commands = [pdf_text_command(page_height=86.4, x=16, y_top=y, text=value, font="F1", size=8, color=color)
                for y, value, color in [(22, title, TEXT_RGB), (48, message, MUTED_RGB)]]
    write_single_page_pdf(width=518.4, height=86.4, commands=commands, output_path=pdf)


def write_organelle_assessment(
    config: MSSPackConfig, item: OrganelleConfig, records: list[dict[str, Any]],
    audit: list[dict[str, Any]] | None, destination: Path,
) -> dict[str, Any]:
    _require(audit is not None or not item.assessment_json, "feature audit required for external assessment")
    cds = [row for row in audit or [] if row["organelle"] == item.name and row["type"] == "CDS"]
    identities = {(item.name + "_" + row["source_record"], tag)
                  for row in cds for tag in row["output_qualifiers"].get("locus_tag", [])}
    partial = {(item.name + "_" + row["source_record"], tag) for row in cds
               if any(char in row["location"] for char in "<>") or {"pseudo", "pseudogene"} & set(row["output_qualifiers"])
               for tag in row["output_qualifiers"].get("locus_tag", [])}
    external = _load_external(config, item, records, identities, partial)
    assignments = {(row["entry"], row["locus_tag"]): row["category"]
                   for row in (external or {}).get("functional_classification", {}).get("assignments", [])}
    categories: Counter[str] = Counter()
    for row in cds:
        keys = [(item.name + "_" + row["source_record"], tag)
                for tag in row["output_qualifiers"].get("locus_tag", [])]
        values = {assignments[key] for key in keys if key in assignments}
        _require(len(values) <= 1, "conflicting functional classes for a CDS")
        categories[next(iter(values)) if values else "Not classified"] += 1
    supplied_count = sum(row["feature_counts"].get("CDS", 0) for row in records)
    if audit is None and supplied_count:
        categories["Not classified"] = supplied_count
    if audit is not None:
        _require(len(cds) == supplied_count, "CDS feature audit does not match the conversion summary")
    swissprot = (external or {}).get("swissprot")
    swiss_counts = Counter(row["status"] for row in swissprot["results"]) if swissprot else Counter(not_evaluated=supplied_count)
    _require(sum(swiss_counts.values()) == supplied_count, "Swiss-Prot CDS counts do not match supplied CDS features")
    swissprot = {**(swissprot or {}), "state_counts": {state: swiss_counts[state] for state in ("hit", "no_hit", "not_evaluated")},
                 "assessed_CDS": swiss_counts["hit"] + swiss_counts["no_hit"]}
    stem = destination / f"organelle-{item.name}"
    functional_svg, functional_pdf = Path(str(stem) + "-functions.svg"), Path(str(stem) + "-functions.pdf")
    title = f"{item.name}: supplied functional annotation"
    if categories:
        events = [EventCount(category, category if len(category) <= 36 else category[:33] + "...", count, "CDS features", COLORS[i % len(COLORS)])
                  for i, (category, count) in enumerate(sorted(categories.items(), key=lambda p: (-p[1], p[0])))]
        for write_chart, path in [(write_event_counts_svg, functional_svg), (write_event_counts_pdf, functional_pdf)]:
            write_chart(events, path, title=title, subtitle="Classes from external annotation evidence; no genes predicted here.")
    else:
        _status_chart(functional_svg, functional_pdf, title, "No annotated CDS supplied; functional classes are not assessed.")
    content = (external or {}).get("gene_content")
    markers = content["markers"] if content else []
    counts = Counter(marker["status"] for marker in markers)
    score = 100 * counts["complete"] / len(markers) if markers and not counts["not_evaluated"] else None
    detected = 100 * (counts["complete"] + counts["fragmented"]) / len(markers) if score is not None else None
    completeness_svg, completeness_pdf = Path(str(stem) + "-gene-content.svg"), Path(str(stem) + "-gene-content.pdf")
    if content:
        events = [EventCount(state, label, counts[state], "reference marker genes", color) for state, (label, color) in STATES.items()]
        title = f"{item.name}: complete reference-gene models " + (f"{score:.1f}%" if score is not None else "NOT ASSESSED")
        for write_chart, path in [(write_event_counts_svg, completeness_svg), (write_event_counts_pdf, completeness_pdf)]:
            write_chart(events, path, title=title, subtitle=f"Reference set: {len(markers)} genes. This is not assembly completeness.")
    else:
        _status_chart(completeness_svg, completeness_pdf, f"{item.name}: reference gene set NOT ASSESSED",
                      "No external marker assessment supplied; score is unavailable.")
    tsv = Path(str(stem) + "-gene-content.tsv")
    with tsv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["gene", "status", "matched_locus_tags", "evidence"])
        for marker in markers:
            writer.writerow([marker["gene"], marker["status"], ",".join(m["locus_tag"] for m in marker.get("matches", [])), marker["evidence"]])
    swiss_tsv = Path(str(stem) + "-swissprot.tsv")
    with swiss_tsv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["entry", "locus_tag", "status", "Swiss-Prot_accession", "evidence"])
        for result in swissprot.get("results", []):
            writer.writerow([result["entry"], result["locus_tag"], result["status"], result.get("accession", ""), result["evidence"]])
    data = dict(status="assessed" if external else "not_assessed", functional_counts=dict(categories),
                supplied_CDS=supplied_count, supplied_product_count=(sum(bool(row["output_qualifiers"].get("product")) for row in cds) if audit is not None else None),
                functional_svg=functional_svg.name, functional_pdf=functional_pdf.name,
                gene_content_svg=completeness_svg.name, gene_content_pdf=completeness_pdf.name,
                gene_content_tsv=tsv.name, complete_percent=score, detected_percent=detected,
                marker_count=len(markers) if content else None, state_counts=dict(counts),
                marker_set=content["marker_set"] if content else None, scope=content["scope"] if content else None,
                thresholds=content.get("thresholds", {}) if content else {},
                functional_classification=(external or {}).get("functional_classification"),
                markers=markers, method=(external or {}).get("method"),
                assessment_sha256=(external or {}).get("assessment_sha256"),
                swissprot=swissprot, swissprot_tsv=swiss_tsv.name,
                limitations=(external or {}).get("limitations", []), annotation_models_changed=False)
    write_text(Path(str(stem) + "-assessment.json"), json.dumps(data, indent=2) + "\n")
    return data


def render_assessment_html(data: dict[str, Any], root: Path) -> str:
    def path(key: str) -> Path:
        value = Path(data[key])
        return value if value.is_absolute() else root / value
    charts = ('<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(min(340px,100%),1fr));gap:16px">'
              + "".join('<div class="figure">' + path(key).read_text(encoding="utf-8").replace('<svg ', '<svg style="max-width:100%;height:auto" ', 1) + '</div>'
                        for key in ("functional_svg", "gene_content_svg")) + '</div>')
    score = data["complete_percent"]
    label = f"{score:.1f}% ({data['state_counts'].get('complete', 0)}/{data['marker_count']})" if score is not None else "Not assessed"
    detected = data["detected_percent"]
    detection = (f"Any qualifying marker evidence (including partial hits): {detected:.1f}% "
                 f"({data['state_counts'].get('complete', 0) + data['state_counts'].get('fragmented', 0)}/{data['marker_count']})."
                 if detected is not None else "Marker detection percentage is not assessed.")
    panel = data.get("marker_set")
    description: str = html.escape(panel["definition"] + " [" + panel["id"] + ", " + panel["version"] + "]") if panel else "No reference set supplied."
    rows = "".join('<tr>' + ''.join('<td>' + html.escape(str(value)) + '</td>' for value in
                   (marker["gene"], marker["status"], ', '.join(m["locus_tag"] for m in marker.get("matches", [])), marker["evidence"])) + '</tr>'
                   for marker in data["markers"])
    links = ' | '.join('<a href="' + html.escape(os.path.relpath(path(key), root)) + '">' + label + '</a>'
                       for key, label in [("functional_pdf", "Functions PDF"), ("gene_content_pdf", "Gene-content PDF"), ("gene_content_tsv", "Per-gene TSV")])
    swiss = data.get("swissprot", {})
    swiss_html: str
    if swiss.get("database"):
        counts = swiss["state_counts"]
        swiss_html = ('<p>Swiss-Prot: ' + str(counts["hit"]) + ' hit, ' + str(counts["no_hit"]) + ' no qualifying hit, '
                      + str(counts["not_evaluated"]) + ' not assessed (CDS features; RNA is not applicable). Database: '
                      + html.escape(swiss["database"]["version"]) + '. '
                      + html.escape(swiss["method"]) + ' Product names were not reassigned.</p>')
    else:
        swiss_html = '<p>Swiss-Prot: not assessed. Supplied gene/product labels do not establish Swiss-Prot matches.</p>'
    if data.get("swissprot_tsv"):
        links += ' | <a href="' + html.escape(os.path.relpath(path("swissprot_tsv"), root)) + '">Swiss-Prot CDS evidence TSV</a>'
    method: str = data.get("method") or "External functional classification and gene-content evaluation were not supplied."
    return ('<section class="organelle-assessment"><h4>Functional annotation and reference gene-set recovery</h4>'
            + swiss_html + charts + '<p>' + links + '</p><p>Complete reference-gene models: <strong>' + label + '</strong>. '
            'This does not establish complete assembly structure or validated mature proteins.</p><p>' + description + '</p>'
            + '<p>' + html.escape(detection) + '</p>'
            + '<p>' + html.escape(method) + '</p>'
            + ''.join('<p>' + html.escape(str(note)) + '</p>' for note in data.get("limitations", []))
            + '<details><summary>Reference genes and evidence</summary><table><tr><th>Gene</th><th>Status</th><th>Matched loci</th><th>Evidence</th></tr>'
            + rows + '</table></details></section>')
