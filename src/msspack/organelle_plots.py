"""CDS/RNA feature flow from supplied organelle annotation to published MSS.

Uses conversion summaries, not nuclear cleanup metrics or inferred annotation.
Gene container rows and exon/intron/repeat rows do not count as extra loci.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from functools import partial
from pathlib import Path
from typing import Any, cast

from .config import MSSPackConfig
from .execution import run_if_needed
from .organelle_assessment import write_organelle_assessment
from .pipeline_plot_models import SankeyLink, SankeyNode
from .pipeline_plot_render import write_gene_flow_tsv, write_sankey_pdf, write_sankey_svg
from .utils import MSSPackError, ensure_dir, expand_path, write_text

FEATURE_COLORS = {
    "CDS": "#059669", "tRNA": "#2563eb", "rRNA": "#7c3aed",
    "ncRNA": "#d97706", "misc_RNA": "#0d9488", "tmRNA": "#c026d3",
}


def _write_flow_charts(
    stages: list[str], nodes: list[SankeyNode], links: list[SankeyLink],
    svg: Path, pdf: Path, tsv: Path, title: str, subtitle: str,
) -> None:
    write_gene_flow_tsv(stages, nodes, links, tsv)
    write_sankey_svg(stages, nodes, links, svg, title=title, subtitle=subtitle)
    write_sankey_pdf(stages, nodes, links, pdf, title=title, subtitle=subtitle)


def _counts(records: list[dict[str, Any]], key: str) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for record in records:
        values = record.get(key)
        if not isinstance(values, dict):
            raise MSSPackError(f"Organelle flow needs {key}; regenerate the organelle submission")
        for feature in FEATURE_COLORS:
            count = values.get(feature, 0)
            if type(count) is not int or count < 0:
                raise MSSPackError(f"Invalid organelle {feature} count in {key}")
            counts[feature] += count
    return {feature: count for feature, count in counts.items() if count}


def write_organelle_gene_flows(
    config: MSSPackConfig, summary: dict[str, Any], destination: Path,
    *, dependencies: list[Path] | None = None, force: bool = False,
    audit: list[dict[str, Any]] | None = None,
) -> list[dict[str, object]]:
    """Render one independent panel per configured organelle, in config order."""
    records = summary.get("records")
    if not isinstance(records, list) or not all(isinstance(row, dict) for row in records):
        raise MSSPackError("Invalid organelle conversion summary: records are missing")
    ensure_dir(destination)
    flows: list[dict[str, object]] = []
    for item in config.organelles:
        group = [row for row in records
                 if row.get("entry") == f"{item.name}_{row.get('input_id')}"
                 and row.get("organelle") == item.organelle]
        if not group:
            raise MSSPackError(f"No prepared organelle records for {item.name}; run 'msspack organelles' first")
        before, after = _counts(group, "input_feature_counts"), _counts(group, "feature_counts")
        if before != after:
            raise MSSPackError(f"Organelle CDS/RNA input/output counts differ for {item.name}; review the conversion")
        assessment = write_organelle_assessment(config, item, group, audit, destination)
        swiss_counts = assessment["swissprot"]["state_counts"]
        total = sum(before.values())
        kind = "Mitochondrion" if item.organelle == "mitochondrion" else "Plastid"
        title = f"{kind} ({item.name}): annotation to MSS"
        nodes, links = [], []
        if total:
            nodes = [SankeyNode("org_input", "Annotated\nfeatures", 0, total, "#334155")]
            for feature, count in sorted(before.items(), key=lambda pair: pair[0].casefold()):
                color = FEATURE_COLORS[feature]
                nodes.append(SankeyNode(feature, feature, 1, count, color))
                links.append(SankeyLink("org_input", feature, count, color))
                if feature != "CDS":
                    links.append(SankeyLink(feature, "org_annotation_RNA", count, color))
            for state, label, color in [("hit", "Swiss-Prot", "#2563eb"),
                                        ("no_hit", "No Swiss-Prot hit", "#94a3b8"),
                                        ("not_evaluated", "Swiss-Prot\nnot assessed", "#64748b")]:
                count = swiss_counts[state]
                if count:
                    node_id = "org_swissprot_" + state
                    nodes.append(SankeyNode(node_id, label, 2, count, color))
                    links.extend([SankeyLink("CDS", node_id, count, color), SankeyLink(node_id, "org_output", count, color)])
            rna_count = total - before.get("CDS", 0)
            if rna_count:
                nodes.append(SankeyNode("org_annotation_RNA", "RNA: not\napplicable", 2, rna_count, "#7c3aed"))
                links.append(SankeyLink("org_annotation_RNA", "org_output", rna_count, "#7c3aed"))
            nodes.append(SankeyNode("org_output", "Retained\nMSS features", 3, total, "#334155"))
        subtitle = ("CDS/RNA feature counts; input gene rows are not counted twice." if total else
                    "No CDS/RNA annotation supplied; no gene flow can be shown.")
        stages = ["Supplied\nannotation", "Feature\ntype", "Swiss-Prot\nannotation", "Output\nann.txt"]
        stem = destination / f"organelle-{item.name}-gene-flow"
        svg, pdf, tsv = (Path(str(stem) + suffix) for suffix in [".sankey.svg", ".sankey.pdf", ".tsv"])
        data: dict[str, object] = dict(
            name=item.name, title=title, organelle=item.organelle,
            entries=[row["entry"] for row in group],
            input_counts=before, output_counts=after, total=total,
            count_unit="CDS/RNA features", source_only=not total,
            svg=svg.name, pdf=pdf.name, tsv=tsv.name,
            swissprot_counts=swiss_counts, assessment_sha256=assessment["assessment_sha256"],
        )
        if force:
            for path in (svg, pdf, tsv):
                path.unlink(missing_ok=True)
        run_if_needed(outputs=[svg, pdf, tsv], dependencies=dependencies or [],
                      cache_key=data, action=partial(_write_flow_charts, stages, nodes, links,
                                                    svg, pdf, tsv, title, subtitle))
        data["assessment"] = assessment
        flows.append(data)
    write_text(destination / "organelle-gene-flows.json", json.dumps(flows, indent=2) + "\n")
    return flows


def prepare_organelle_gene_flows(config: MSSPackConfig, *, force: bool = False) -> list[dict[str, object]]:
    if not config.organelles:
        return []
    final = config.output_dir / "organelles" / "final"
    if not (final / "summary.json").is_file():
        raise MSSPackError("Organelle results are missing; run 'msspack organelles --config ...' first")
    # Pin one published generation while reading its summary and inputs.
    final = final.resolve(strict=True)
    summary = json.loads((final / "summary.json").read_text(encoding="utf-8"))
    expected = {(item.name, str(path), hashlib.sha256(path.read_bytes()).hexdigest())
                for item in config.organelles for source in (item.genbank, item.fasta, item.gff)
                if source for path in [expand_path(source, config.base_dir)]}
    observed = {(row["organelle"], row["path"], row["sha256"]) for row in summary.get("inputs", [])}
    if expected != observed:
        raise MSSPackError("Organelle inputs changed since preparation; regenerate organelles before plotting")
    destination = config.output_dir / "plots"
    audit_path = final / "feature-audit.json"
    audit = json.loads(audit_path.read_text()) if audit_path.exists() else None
    flows = write_organelle_gene_flows(
        config, summary, destination,
        dependencies=[final / "summary.json"], force=force, audit=audit,
    )
    result = []
    for flow in flows:
        assessment = dict(cast(dict[str, Any], flow["assessment"]))
        for key in ("functional_svg", "functional_pdf", "gene_content_svg", "gene_content_pdf", "gene_content_tsv", "swissprot_tsv"):
            assessment[key] = str(destination / assessment[key])
        result.append({**flow, **{key: str(destination / str(flow[key])) for key in ("svg", "pdf", "tsv")},
                       "assessment": assessment})
    return result
