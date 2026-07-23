from __future__ import annotations

from pathlib import Path

from .annotation_consistency_plots import (
    load_name_consistency_plot_data,
    load_source_consistency_plot_data,
    prepare_name_consistency_plot_data,
    prepare_source_consistency_plot_data,
    write_name_consistency_pdf,
    write_name_consistency_svg,
    write_source_consistency_pdf,
    write_source_consistency_svg,
)
from .config import load_config
from .execution import module_origin, run_if_needed
from .pipeline_plot_data import (
    collect_pipeline_plot_data,
    parse_pipeline_plot_metrics,
    parse_step_log,
)
from .pipeline_plot_models import PipelinePlotArtifacts
from .pipeline_plot_render import (
    build_event_counts,
    build_plot_artifacts,
    build_sankey,
    load_sankey_busco_summaries,
    summarize_pipeline_plots,
    update_plot_manifest,
    write_event_counts_pdf,
    write_event_counts_svg,
    write_event_counts_tsv,
    write_gene_flow_tsv,
    write_overlap_pdf,
    write_overlap_svg,
    write_overlap_tsv,
    write_sankey_pdf,
    write_sankey_svg,
    write_summary_json,
    write_summary_tsv,
)
from .utils import MSSPackError


def run_pipeline_plots(
    config_file: str | Path,
    *,
    force: bool = False,
) -> PipelinePlotArtifacts:
    config = load_config(config_file)
    output_root = config.output_dir
    log_dir = output_root / "logs"
    manifest_path = output_root / "build-manifest.json"
    if not log_dir.exists():
        raise MSSPackError(
            f"Pipeline logs are missing under {log_dir}. Run 'msspack pack --config {config_file}' first."
        )

    artifacts = build_plot_artifacts(output_root)
    module_paths = [
        module_origin("msspack.pipeline_plots"),
        module_origin("msspack.pipeline_plot_data"),
        module_origin("msspack.pipeline_plot_models"),
        module_origin("msspack.pipeline_plot_render"),
        module_origin("msspack.annotation_consistency_plots"),
    ]
    required_logs = [
        log_dir / "06.drop-duplicate-coordinate-gene.log",
        log_dir / "07.select-one-mrna.log",
        log_dir / "09.update-gff-to-inframe.log",
        log_dir / "11.update-gff-with-padding.log",
        log_dir / "16.mss-cds-to-misc.log",
    ]
    if force:
        for path in (
            artifacts.summary_json,
            artifacts.summary_tsv,
            artifacts.gene_flow_tsv,
            artifacts.gene_flow_svg,
            artifacts.gene_flow_pdf,
            artifacts.event_counts_tsv,
            artifacts.event_counts_svg,
            artifacts.event_counts_pdf,
            artifacts.overlap_tsv,
            artifacts.overlap_svg,
            artifacts.overlap_pdf,
            artifacts.name_consistency_tsv,
            artifacts.name_consistency_svg,
            artifacts.name_consistency_pdf,
            artifacts.source_consistency_tsv,
            artifacts.source_consistency_svg,
            artifacts.source_consistency_pdf,
        ):
            if path.exists():
                path.unlink()

    bundle = collect_pipeline_plot_data(output_root, log_dir)
    stage_labels, nodes, links = build_sankey(
        bundle.metrics,
        bundle.gene_sets,
        bundle.functional_annotation,
        bundle.annotation_consistency,
    )
    busco_summaries = load_sankey_busco_summaries(output_root)
    events = build_event_counts(bundle.metrics)
    busco_comparison_path = output_root / "busco" / "cds" / "comparison.json"
    busco_dependencies = [busco_comparison_path] if busco_comparison_path.exists() else []
    dependency_paths = [
        *required_logs,
        *(record.path for record in bundle.records.values()),
        *(gene_set.path for gene_set in bundle.gene_sets),
        *([bundle.functional_annotation.path] if bundle.functional_annotation is not None else []),
        *(
            [
                bundle.annotation_consistency.path,
                bundle.annotation_consistency.summary_path,
                bundle.annotation_consistency.source_pair_path,
            ]
            if bundle.annotation_consistency is not None
            else []
        ),
        *module_paths,
    ]

    run_if_needed(
        outputs=[
            artifacts.summary_json,
            artifacts.summary_tsv,
            artifacts.gene_flow_tsv,
            artifacts.event_counts_tsv,
            artifacts.overlap_tsv,
        ],
        dependencies=dependency_paths,
        action=lambda: (
            write_summary_json(bundle.summary_payload, artifacts.summary_json),
            write_summary_tsv(
                bundle.metrics,
                artifacts.summary_tsv,
                bundle.functional_annotation,
                bundle.annotation_consistency,
            ),
            write_gene_flow_tsv(stage_labels, nodes, links, artifacts.gene_flow_tsv),
            write_event_counts_tsv(events, artifacts.event_counts_tsv),
            write_overlap_tsv(bundle.gene_sets, bundle.overlap_rows, artifacts.overlap_tsv),
        ),
    )
    run_if_needed(
        outputs=[artifacts.gene_flow_svg],
        dependencies=[
            artifacts.summary_json,
            artifacts.gene_flow_tsv,
            *busco_dependencies,
            *module_paths,
        ],
        action=lambda: write_sankey_svg(
            stage_labels,
            nodes,
            links,
            artifacts.gene_flow_svg,
            busco_summaries=busco_summaries,
            annotation_consistency=bundle.annotation_consistency,
        ),
    )
    run_if_needed(
        outputs=[artifacts.gene_flow_pdf],
        dependencies=[
            artifacts.summary_json,
            artifacts.gene_flow_tsv,
            *busco_dependencies,
            *module_paths,
        ],
        action=lambda: write_sankey_pdf(
            stage_labels,
            nodes,
            links,
            artifacts.gene_flow_pdf,
            busco_summaries=busco_summaries,
            annotation_consistency=bundle.annotation_consistency,
        ),
    )
    run_if_needed(
        outputs=[artifacts.event_counts_svg],
        dependencies=[artifacts.summary_json, artifacts.event_counts_tsv, *module_paths],
        action=lambda: write_event_counts_svg(events, artifacts.event_counts_svg),
    )
    run_if_needed(
        outputs=[artifacts.event_counts_pdf],
        dependencies=[artifacts.summary_json, artifacts.event_counts_tsv, *module_paths],
        action=lambda: write_event_counts_pdf(events, artifacts.event_counts_pdf),
    )
    run_if_needed(
        outputs=[artifacts.overlap_svg],
        dependencies=[artifacts.summary_json, artifacts.overlap_tsv, *module_paths],
        action=lambda: write_overlap_svg(
            bundle.gene_sets, bundle.overlap_rows, artifacts.overlap_svg
        ),
    )
    run_if_needed(
        outputs=[artifacts.overlap_pdf],
        dependencies=[artifacts.summary_json, artifacts.overlap_tsv, *module_paths],
        action=lambda: write_overlap_pdf(
            bundle.gene_sets, bundle.overlap_rows, artifacts.overlap_pdf
        ),
    )
    if bundle.annotation_consistency is not None:
        consistency = bundle.annotation_consistency
        run_if_needed(
            outputs=[artifacts.name_consistency_tsv],
            dependencies=[consistency.summary_path, *module_paths],
            action=lambda: prepare_name_consistency_plot_data(
                consistency.summary_path,
                artifacts.name_consistency_tsv,
            ),
        )
        run_if_needed(
            outputs=[artifacts.source_consistency_tsv],
            dependencies=[consistency.source_pair_path, *module_paths],
            action=lambda: prepare_source_consistency_plot_data(
                consistency.source_pair_path,
                artifacts.source_consistency_tsv,
            ),
        )
        name_rows = load_name_consistency_plot_data(artifacts.name_consistency_tsv)
        source_rows = load_source_consistency_plot_data(artifacts.source_consistency_tsv)
        run_if_needed(
            outputs=[artifacts.name_consistency_svg],
            dependencies=[artifacts.name_consistency_tsv, *module_paths],
            action=lambda: write_name_consistency_svg(
                name_rows,
                artifacts.name_consistency_svg,
            ),
        )
        run_if_needed(
            outputs=[artifacts.name_consistency_pdf],
            dependencies=[artifacts.name_consistency_tsv, *module_paths],
            action=lambda: write_name_consistency_pdf(
                name_rows,
                artifacts.name_consistency_pdf,
            ),
        )
        run_if_needed(
            outputs=[artifacts.source_consistency_svg],
            dependencies=[artifacts.source_consistency_tsv, *module_paths],
            action=lambda: write_source_consistency_svg(
                source_rows,
                artifacts.source_consistency_svg,
            ),
        )
        run_if_needed(
            outputs=[artifacts.source_consistency_pdf],
            dependencies=[artifacts.source_consistency_tsv, *module_paths],
            action=lambda: write_source_consistency_pdf(
                source_rows,
                artifacts.source_consistency_pdf,
            ),
        )
    update_plot_manifest(
        manifest_path,
        artifacts=artifacts,
        metrics=bundle.metrics,
        gene_sets=bundle.gene_sets,
        overlap_rows=bundle.overlap_rows,
        annotation_consistency=bundle.annotation_consistency,
    )
    return artifacts


__all__ = [
    "PipelinePlotArtifacts",
    "parse_pipeline_plot_metrics",
    "parse_step_log",
    "run_pipeline_plots",
    "summarize_pipeline_plots",
]
