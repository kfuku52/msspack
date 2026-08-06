from __future__ import annotations

import time
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from . import __version__
from .annotation_consistency import (
    audit_annotation_consistency,
    run_family_similarity_search,
)
from .annotation_table import build_annotation_table
from .annotation_taxonomy import (
    discover_busco_summary_paths,
    resolve_annotation_taxonomy,
)
from .build_manifest import ManifestRecorder
from .config import MSSPackConfig, load_config
from .ddbj_tools import describe_installation, list_installed
from .execution import module_origin, path_list, run_if_needed
from .functional_annotation import (
    apply_functional_annotations,
    run_cdd_domain_search,
    run_diamond_annotation_search,
    run_pfam_domain_search,
    write_empty_cdd_results,
    write_empty_diamond_results,
    write_empty_pfam_results,
    write_translated_protein_fasta,
)
from .gff_cleanup import (
    drop_duplicate_coordinate_genes,
    fix_gff_semicolons_file,
    trim_gff_to_fasta_bounds,
)
from .padding_tools import write_padding_log_for_gff
from .pipeline_actions import (
    copy_input_fasta,
    copy_input_gff,
    pad_locus_tags,
    render_final_annotation_with_log,
    replace_annotation_products,
    run_apply_padding_to_gff,
    run_fix_gff_to_inframe,
    run_gapjust,
    run_gff2mss,
    run_mss_cds_to_misc,
    run_remove_trailing_ns,
    run_select_one_mrna,
    sort_gff,
    write_gapjust_passthrough,
    write_submission_fasta,
)
from .submission_render import render_header
from .utils import (
    MSSPackError,
    ensure_dir,
    expand_path,
)
from .validation import (
    ValidationArtifacts,
    ValidationOptions,
    run_validation,
    write_validation_not_run_summary,
)


@dataclass
class PipelineOutputs:
    root: Path
    intermediate: Path
    logs: Path
    final: Path
    ann_path: Path
    fasta_path: Path
    manifest_path: Path


@dataclass(frozen=True)
class ModulePaths:
    pipeline: Path
    pipeline_actions: Path
    fasta_steps: Path
    gap_normalization: Path
    gff_cleanup: Path
    coordinate_duplicates: Path
    submission_render: Path
    gff: Path
    gff_adjustments: Path
    gff_inframe: Path
    gff_padding_update: Path
    transcript_selection: Path
    mss_postprocess: Path
    padding_tools: Path
    annotation_table: Path
    annotation_consistency: Path
    annotation_taxonomy: Path
    functional_annotation: Path
    product_names: Path
    mss_converter_core: Path
    mss_converter_features: Path
    mss_converter_gaps: Path
    mss_converter_models: Path
    mss_converter_render: Path


@dataclass(frozen=True)
class PipelineContext:
    config: MSSPackConfig
    config_path: Path
    outputs: PipelineOutputs
    modules: ModulePaths
    manifest: ManifestRecorder

    def run_step(
        self,
        *,
        name: str,
        outputs: object,
        dependencies: object,
        action: Callable[[], object],
    ) -> bool:
        output_paths = path_list(outputs)
        dependency_paths = path_list(dependencies)
        started = time.perf_counter()
        ran = run_if_needed(
            outputs=output_paths,
            dependencies=dependency_paths,
            action=action,
            cache_key={"msspack_version": __version__, "step": name},
        )
        self.manifest.record_stage(
            name=name,
            outputs=output_paths,
            dependencies=dependency_paths,
            ran=ran,
            duration_seconds=time.perf_counter() - started,
        )
        return ran


@dataclass(frozen=True)
class PreparedInputs:
    prepared_fasta: Path
    prepared_gff: Path
    downstream_fasta: Path
    downstream_gff: Path


@dataclass(frozen=True)
class PreparedGff:
    padded_gff: Path
    genes_with_stops: Path


@dataclass(frozen=True)
class AnnotationArtifacts:
    annotation_table: Path
    postprocessed_gff2mss: Path


def _resolve_modules() -> ModulePaths:
    return ModulePaths(
        pipeline=Path(__file__).resolve(),
        pipeline_actions=module_origin("msspack.pipeline_actions"),
        fasta_steps=module_origin("msspack.fasta_steps"),
        gap_normalization=module_origin("msspack.gap_normalization"),
        gff_cleanup=module_origin("msspack.gff_cleanup"),
        coordinate_duplicates=module_origin("msspack.coordinate_duplicates"),
        submission_render=module_origin("msspack.submission_render"),
        gff=module_origin("msspack.gff"),
        gff_adjustments=module_origin("msspack.gff_adjustments"),
        gff_inframe=module_origin("msspack.gff_inframe"),
        gff_padding_update=module_origin("msspack.gff_padding_update"),
        transcript_selection=module_origin("msspack.transcript_selection"),
        mss_postprocess=module_origin("msspack.mss_postprocess"),
        padding_tools=module_origin("msspack.padding_tools"),
        annotation_table=module_origin("msspack.annotation_table"),
        annotation_consistency=module_origin("msspack.annotation_consistency"),
        annotation_taxonomy=module_origin("msspack.annotation_taxonomy"),
        functional_annotation=module_origin("msspack.functional_annotation"),
        product_names=module_origin("msspack.product_names"),
        mss_converter_core=module_origin("msspack.mss_converter.core"),
        mss_converter_features=module_origin("msspack.mss_converter.features"),
        mss_converter_gaps=module_origin("msspack.mss_converter.gaps"),
        mss_converter_models=module_origin("msspack.mss_converter.models"),
        mss_converter_render=module_origin("msspack.mss_converter.render"),
    )


def _build_outputs(config: MSSPackConfig) -> PipelineOutputs:
    output_root = ensure_dir(config.output_dir)
    intermediate = ensure_dir(output_root / "intermediate")
    logs = ensure_dir(output_root / "logs")
    final = ensure_dir(output_root / "final")
    return PipelineOutputs(
        root=output_root,
        intermediate=intermediate,
        logs=logs,
        final=final,
        ann_path=final / f"{config.sample.locus_tag}.ann.txt",
        fasta_path=final / f"{config.sample.locus_tag}.fasta",
        manifest_path=output_root / "build-manifest.json",
    )


def _initialize_pipeline(
    config_file: str | Path,
) -> tuple[MSSPackConfig, Path, PipelineOutputs, ModulePaths, PipelineContext]:
    config_path = Path(config_file).expanduser().resolve()
    config = load_config(config_file)
    outputs = _build_outputs(config)
    modules = _resolve_modules()
    manifest = ManifestRecorder(
        config=config,
        config_path=config_path,
        output_root=outputs.root,
        ann_path=outputs.ann_path,
        fasta_path=outputs.fasta_path,
    )
    return (
        config,
        config_path,
        outputs,
        modules,
        PipelineContext(
            config=config,
            config_path=config_path,
            outputs=outputs,
            modules=modules,
            manifest=manifest,
        ),
    )


def prepare_pipeline_for_busco(config_file: str | Path) -> PipelineOutputs:
    """Prepare normalized FASTA/GFF artifacts required by BUSCO without annotation."""
    _config, _config_path, outputs, _modules, ctx = _initialize_pipeline(config_file)
    prepared_inputs = _prepare_inputs(ctx)
    _prepare_gff(ctx, prepared_inputs)
    return outputs


def _prepare_inputs(ctx: PipelineContext) -> PreparedInputs:
    config = ctx.config
    input_fasta = config.fasta_path
    input_gff = config.gff_path
    if not input_fasta.exists():
        raise MSSPackError(f"Input FASTA not found: {input_fasta}")
    if not input_gff.exists():
        raise MSSPackError(f"Input GFF not found: {input_gff}")

    intermediate = ctx.outputs.intermediate
    logs = ctx.outputs.logs
    prepared_fasta = intermediate / "00.input.genome.fasta"
    copy_fasta_metrics = logs / "00.copy-input-fasta.metrics.json"
    ctx.run_step(
        name="00.copy-input-fasta",
        outputs=[prepared_fasta, logs / "00.copy-input-fasta.log", copy_fasta_metrics],
        dependencies=[input_fasta, ctx.modules.pipeline_actions],
        action=lambda: copy_input_fasta(
            input_path=input_fasta,
            output_path=prepared_fasta,
            log_path=logs / "00.copy-input-fasta.log",
            metrics_path=copy_fasta_metrics,
        ),
    )

    prepared_gff = intermediate / "00.input.annotations.gff"
    copy_gff_metrics = logs / "00.copy-input-gff.metrics.json"
    ctx.run_step(
        name="00.copy-input-gff",
        outputs=[prepared_gff, logs / "00.copy-input-gff.log", copy_gff_metrics],
        dependencies=[input_gff, ctx.modules.pipeline_actions],
        action=lambda: copy_input_gff(
            input_path=input_gff,
            output_path=prepared_gff,
            log_path=logs / "00.copy-input-gff.log",
            metrics_path=copy_gff_metrics,
        ),
    )

    trailing_fasta = intermediate / "01.fasta.trailing-ns-removed.fasta"
    trailing_metrics = logs / "01.remove-trailing-ns.metrics.json"
    ctx.run_step(
        name="01.remove-trailing-ns",
        outputs=[trailing_fasta, logs / "01.remove-trailing-ns.log", trailing_metrics],
        dependencies=[prepared_fasta, ctx.modules.fasta_steps, ctx.modules.pipeline_actions],
        action=lambda: run_remove_trailing_ns(
            input_path=prepared_fasta,
            output_path=trailing_fasta,
            log_path=logs / "01.remove-trailing-ns.log",
            metrics_path=trailing_metrics,
        ),
    )

    gapjust_fasta = intermediate / "02.gap-normalized.genome.fasta"
    gapjust_gff = intermediate / "02.gap-normalized.annotations.gff"
    gapjust_log = logs / "02.gapjust.log"
    gapjust_metrics = logs / "02.gapjust.metrics.json"
    if config.pipeline.run_gapjust:
        ctx.run_step(
            name="02.gapjust",
            outputs=[gapjust_fasta, gapjust_gff, gapjust_log, gapjust_metrics],
            dependencies=[
                trailing_fasta,
                prepared_gff,
                ctx.config_path,
                ctx.modules.pipeline_actions,
                ctx.modules.gap_normalization,
            ],
            action=lambda: run_gapjust(
                input_fasta=trailing_fasta,
                output_fasta=gapjust_fasta,
                input_gff=prepared_gff,
                output_gff=gapjust_gff,
                gap_len=config.pipeline.gapjust_gap_len,
                gap_just_min=config.pipeline.gapjust_min,
                gap_just_max=config.pipeline.gapjust_max,
                log_path=gapjust_log,
                metrics_path=gapjust_metrics,
            ),
        )
    else:
        ctx.run_step(
            name="02.gapjust-skip",
            outputs=[gapjust_fasta, gapjust_gff, gapjust_log, gapjust_metrics],
            dependencies=[trailing_fasta, prepared_gff, ctx.config_path, ctx.modules.pipeline_actions],
            action=lambda: write_gapjust_passthrough(
                input_fasta=trailing_fasta,
                input_gff=prepared_gff,
                output_fasta=gapjust_fasta,
                output_gff=gapjust_gff,
                log_path=gapjust_log,
                metrics_path=gapjust_metrics,
            ),
        )
    downstream_fasta = gapjust_fasta
    downstream_gff = gapjust_gff

    stage_fasta = intermediate / "03.fasta.submission.fasta"
    submission_fasta_metrics = logs / "03.write-submission-fasta.metrics.json"
    ctx.run_step(
        name="03.write-submission-fasta",
        outputs=[stage_fasta, ctx.outputs.fasta_path, logs / "03.write-submission-fasta.log", submission_fasta_metrics],
        dependencies=[downstream_fasta, ctx.modules.fasta_steps, ctx.modules.pipeline_actions],
        action=lambda: write_submission_fasta(
            input_path=downstream_fasta,
            stage_output_path=stage_fasta,
            final_output_path=ctx.outputs.fasta_path,
            log_path=logs / "03.write-submission-fasta.log",
            metrics_path=submission_fasta_metrics,
        ),
    )

    return PreparedInputs(
        prepared_fasta=prepared_fasta,
        prepared_gff=prepared_gff,
        downstream_fasta=downstream_fasta,
        downstream_gff=downstream_gff,
    )


def _prepare_gff(ctx: PipelineContext, inputs: PreparedInputs) -> PreparedGff:
    config = ctx.config
    intermediate = ctx.outputs.intermediate
    logs = ctx.outputs.logs

    semicolon_fixed_gff = intermediate / "04.gff.semicolons-fixed.gff"
    semicolon_metrics = logs / "04.fix-gff-semicolons.metrics.json"
    ctx.run_step(
        name="04.fix-gff-semicolons",
        outputs=[semicolon_fixed_gff, logs / "04.fix-gff-semicolons.log", semicolon_metrics],
        dependencies=[inputs.downstream_gff, ctx.modules.gff_cleanup],
        action=lambda: fix_gff_semicolons_file(
            input_path=inputs.downstream_gff,
            output_path=semicolon_fixed_gff,
            log_path=logs / "04.fix-gff-semicolons.log",
            metrics_path=semicolon_metrics,
        ),
    )

    trimmed_gff = intermediate / "05.gff.coordinates-trimmed.gff"
    trim_metrics = logs / "05.trim-gff-coordinates.metrics.json"
    ctx.run_step(
        name="05.trim-gff-coordinates",
        outputs=[trimmed_gff, logs / "05.trim-gff-coordinates.log", trim_metrics],
        dependencies=[semicolon_fixed_gff, inputs.downstream_fasta, ctx.modules.gff_cleanup],
        action=lambda: trim_gff_to_fasta_bounds(
            gff_path=semicolon_fixed_gff,
            fasta_path=inputs.downstream_fasta,
            output_path=trimmed_gff,
            log_path=logs / "05.trim-gff-coordinates.log",
            metrics_path=trim_metrics,
        ),
    )

    dedup_gff = intermediate / "06.gff.duplicate-genes-removed.gff"
    dedup_removed_gene_ids = logs / "06.drop-duplicate-coordinate-gene.changed-gene-ids.txt"
    dedup_duplicate_map = logs / "06.drop-duplicate-coordinate-gene.duplicate-map.tsv"
    dedup_metrics = logs / "06.drop-duplicate-coordinate-gene.metrics.json"
    ctx.run_step(
        name="06.drop-duplicate-coordinate-gene",
        outputs=[
            dedup_gff,
            logs / "06.drop-duplicate-coordinate-gene.log",
            dedup_removed_gene_ids,
            dedup_duplicate_map,
            dedup_metrics,
        ],
        dependencies=[
            trimmed_gff,
            inputs.downstream_fasta,
            ctx.config_path,
            ctx.modules.gff_cleanup,
            ctx.modules.coordinate_duplicates,
        ],
        action=lambda: drop_duplicate_coordinate_genes(
            input_path=trimmed_gff,
            fasta_path=inputs.downstream_fasta,
            output_path=dedup_gff,
            log_path=logs / "06.drop-duplicate-coordinate-gene.log",
            genetic_code=config.sample.genetic_code,
            selection_policy=config.pipeline.coordinate_duplicate_policy,
            removed_gene_ids_path=dedup_removed_gene_ids,
            duplicate_map_path=dedup_duplicate_map,
            metrics_path=dedup_metrics,
        ),
    )

    one_mrna_gff = intermediate / "07.gff.single-mrna-per-gene.gff"
    transcript_changed_gene_ids = logs / "07.select-one-mrna.changed-gene-ids.txt"
    transcript_metrics = logs / "07.select-one-mrna.metrics.json"
    ctx.run_step(
        name="07.select-one-mrna",
        outputs=[
            one_mrna_gff,
            logs / "07.select-one-mrna.log",
            transcript_changed_gene_ids,
            transcript_metrics,
        ],
        dependencies=[dedup_gff, ctx.modules.transcript_selection, ctx.modules.pipeline_actions],
        action=lambda: run_select_one_mrna(
            input_path=dedup_gff,
            output_path=one_mrna_gff,
            log_path=logs / "07.select-one-mrna.log",
            changed_gene_ids_path=transcript_changed_gene_ids,
            metrics_path=transcript_metrics,
        ),
    )

    padded_ids_gff = intermediate / "08.gff.locus-tags-padded.gff"
    locus_tag_metrics = logs / "08.pad-locus-tags.metrics.json"
    ctx.run_step(
        name="08.pad-locus-tags",
        outputs=[padded_ids_gff, logs / "08.pad-locus-tags.log", locus_tag_metrics],
        dependencies=[one_mrna_gff, ctx.config_path, ctx.modules.pipeline_actions],
        action=lambda: pad_locus_tags(
            one_mrna_gff,
            padded_ids_gff,
            config.sample.locus_tag,
            config.sample.locus_tag_digits,
            log_path=logs / "08.pad-locus-tags.log",
            metrics_path=locus_tag_metrics,
        ),
    )

    frame_fixed_gff = intermediate / "09.gff.inframe-fixed.gff"
    frame_fix_summary_log = intermediate / "09.gff.inframe-fix.summary.log"
    inframe_updated_gene_ids = logs / "09.update-gff-to-inframe.changed-gene-ids.txt"
    inframe_metrics = logs / "09.update-gff-to-inframe.metrics.json"
    ctx.run_step(
        name="09.update-gff-to-inframe",
        outputs=[
            frame_fixed_gff,
            frame_fix_summary_log,
            logs / "09.update-gff-to-inframe.log",
            inframe_updated_gene_ids,
            inframe_metrics,
        ],
        dependencies=[
            padded_ids_gff,
            ctx.modules.gff_adjustments,
            ctx.modules.gff_inframe,
            ctx.modules.pipeline_actions,
        ],
        action=lambda: run_fix_gff_to_inframe(
            input_path=padded_ids_gff,
            output_path=frame_fixed_gff,
            summary_log_path=frame_fix_summary_log,
            command_log_path=logs / "09.update-gff-to-inframe.log",
            updated_gene_ids_path=inframe_updated_gene_ids,
            metrics_path=inframe_metrics,
        ),
    )

    frame_fixed_pad_log = intermediate / "10.gff.padding-analysis.log"
    extract_padding_metrics = logs / "10.extract-cds-for-padding.metrics.json"
    pad_analysis_metrics = logs / "10.pad-from-gff.metrics.json"
    ctx.run_step(
        name="10.pad-from-gff",
        outputs=[
            logs / "10.extract-cds-for-padding.log",
            extract_padding_metrics,
            frame_fixed_pad_log,
            pad_analysis_metrics,
        ],
        dependencies=[
            inputs.downstream_fasta,
            frame_fixed_gff,
            ctx.config_path,
            ctx.modules.pipeline_actions,
            ctx.modules.padding_tools,
        ],
        action=lambda: write_padding_log_for_gff(
            fasta_path=inputs.downstream_fasta,
            gff_path=frame_fixed_gff,
            extract_log_path=logs / "10.extract-cds-for-padding.log",
            log_path=frame_fixed_pad_log,
            genetic_code=config.sample.genetic_code,
            extract_metrics_path=extract_padding_metrics,
            metrics_path=pad_analysis_metrics,
        ),
    )

    padded_gff = intermediate / "11.gff.padding-applied.gff"
    genes_with_stops = intermediate / "11.gff.genes-with-stops.txt"
    updated_genes = intermediate / "11.gff.updated-genes.txt"
    padding_metrics = logs / "11.update-gff-with-padding.metrics.json"
    ctx.run_step(
        name="11.update-gff-with-padding",
        outputs=[
            padded_gff,
            genes_with_stops,
            updated_genes,
            logs / "11.update-gff-with-padding.log",
            padding_metrics,
        ],
        dependencies=[
            frame_fixed_gff,
            frame_fixed_pad_log,
            ctx.modules.gff_adjustments,
            ctx.modules.gff_padding_update,
            ctx.modules.pipeline_actions,
        ],
        action=lambda: run_apply_padding_to_gff(
            gff_path=frame_fixed_gff,
            padding_log_path=frame_fixed_pad_log,
            output_path=padded_gff,
            genes_with_stops_path=genes_with_stops,
            updated_genes_path=updated_genes,
            command_log_path=logs / "11.update-gff-with-padding.log",
            metrics_path=padding_metrics,
        ),
    )

    final_sorted_gff = intermediate / "12.gff.final-sorted.gff"
    sort_metrics = logs / "12.gff3sort.metrics.json"
    ctx.run_step(
        name="12.gff3sort",
        outputs=[final_sorted_gff, logs / "12.gff3sort.log", sort_metrics],
        dependencies=[padded_gff, ctx.config_path, ctx.modules.pipeline_actions, ctx.modules.gff],
        action=lambda: sort_gff(
            padded_gff,
            final_sorted_gff,
            command_override=config.tools.gff3sort,
            log_path=logs / "12.gff3sort.log",
            metrics_path=sort_metrics,
        ),
    )

    return PreparedGff(
        padded_gff=final_sorted_gff,
        genes_with_stops=genes_with_stops,
    )


def _build_annotation_artifacts(
    ctx: PipelineContext,
    inputs: PreparedInputs,
    prepared_gff: PreparedGff,
) -> AnnotationArtifacts:
    config = ctx.config
    intermediate = ctx.outputs.intermediate
    logs = ctx.outputs.logs

    raw_annotation_table = intermediate / "13.annotation-table.raw.tsv"
    annotation_metrics = logs / "13.gff2annotation-table.metrics.json"
    ctx.run_step(
        name="13.gff2annotation-table",
        outputs=[raw_annotation_table, logs / "13.gff2annotation-table.log", annotation_metrics],
        dependencies=[prepared_gff.padded_gff, ctx.config_path, ctx.modules.annotation_table],
        action=lambda: build_annotation_table(
            gff_path=prepared_gff.padded_gff,
            output_path=raw_annotation_table,
            locus_tag_prefix=config.sample.locus_tag,
            log_path=logs / "13.gff2annotation-table.log",
            metrics_path=annotation_metrics,
        ),
    )

    annotation_table = intermediate / "14.annotation-table.normalized.tsv"
    annotation_normalize_metrics = logs / "14.normalize-annotation-products.metrics.json"
    ctx.run_step(
        name="14.normalize-annotation-products",
        outputs=[annotation_table, logs / "14.normalize-annotation-products.log", annotation_normalize_metrics],
        dependencies=[raw_annotation_table, ctx.config_path, ctx.modules.pipeline_actions],
        action=lambda: replace_annotation_products(
            raw_annotation_table,
            annotation_table,
            config.pipeline.replace_product_patterns,
            config.pipeline.replace_product_with,
            log_path=logs / "14.normalize-annotation-products.log",
            metrics_path=annotation_normalize_metrics,
        ),
    )

    if config.functional_annotation.enabled:
        annotation_config = config.functional_annotation
        protein_fasta = intermediate / "14a.functional-annotation.proteins.fasta"
        protein_metrics = logs / "14a.functional-annotation-extract-proteins.metrics.json"
        ctx.run_step(
            name="14a.functional-annotation-extract-proteins",
            outputs=[
                protein_fasta,
                logs / "14a.functional-annotation-extract-proteins.log",
                protein_metrics,
            ],
            dependencies=[
                inputs.downstream_fasta,
                prepared_gff.padded_gff,
                ctx.config_path,
                ctx.modules.functional_annotation,
            ],
            action=lambda: write_translated_protein_fasta(
                fasta_path=inputs.downstream_fasta,
                gff_path=prepared_gff.padded_gff,
                output_path=protein_fasta,
                genetic_code=config.sample.genetic_code,
                log_path=logs / "14a.functional-annotation-extract-proteins.log",
                metrics_path=protein_metrics,
            ),
        )

        taxonomy_context = ctx.outputs.final / "functional-annotation-taxonomy.json"
        taxonomy_metrics = logs / "14a2.functional-annotation-taxonomy.metrics.json"
        busco_summary_paths = (
            discover_busco_summary_paths(ctx.outputs.root)
            if annotation_config.taxonomy.busco_crosscheck
            else []
        )
        ctx.run_step(
            name="14a2.functional-annotation-resolve-taxonomy",
            outputs=[
                taxonomy_context,
                logs / "14a2.functional-annotation-taxonomy.log",
                taxonomy_metrics,
            ],
            dependencies=[
                ctx.config_path,
                ctx.modules.annotation_taxonomy,
                *busco_summary_paths,
            ],
            action=lambda: resolve_annotation_taxonomy(
                scientific_name=config.sample.scientific_name,
                configured_busco_lineage=config.busco.lineage_dataset,
                busco_summary_paths=busco_summary_paths,
                output_path=taxonomy_context,
                log_path=logs / "14a2.functional-annotation-taxonomy.log",
                metrics_path=taxonomy_metrics,
                cache_dir=config.database_dir / "taxonomy",
                config=annotation_config.taxonomy,
            ),
        )

        diamond_hits = intermediate / "14b.functional-annotation.diamond.tsv"
        diamond_metadata = intermediate / "14b.functional-annotation.diamond-metadata.tsv"
        diamond_provenance = logs / "14b.functional-annotation-primary.database.json"
        diamond_metrics = logs / "14b.functional-annotation-primary-search.metrics.json"
        diamond_dependencies: list[Path] = [
            protein_fasta,
            taxonomy_context,
            ctx.config_path,
            ctx.modules.functional_annotation,
            ctx.modules.annotation_taxonomy,
        ]
        if annotation_config.swissprot_fasta:
            diamond_dependencies.append(
                expand_path(annotation_config.swissprot_fasta, config.base_dir)
            )
        if annotation_config.reference_proteins:
            diamond_dependencies.append(
                expand_path(annotation_config.reference_proteins, config.base_dir)
            )
        has_diamond_source = bool(
            annotation_config.swissprot_enabled or annotation_config.reference_proteins.strip()
        )
        ctx.run_step(
            name="14b.functional-annotation-primary-search",
            outputs=[
                diamond_hits,
                diamond_metadata,
                diamond_provenance,
                logs / "14b.functional-annotation-primary-search.log",
                diamond_metrics,
            ],
            dependencies=diamond_dependencies,
            action=(
                lambda: run_diamond_annotation_search(
                    protein_fasta_path=protein_fasta,
                    output_path=diamond_hits,
                    metadata_path=diamond_metadata,
                    provenance_path=diamond_provenance,
                    log_path=logs / "14b.functional-annotation-primary-search.log",
                    metrics_path=diamond_metrics,
                    config=annotation_config,
                    base_dir=config.base_dir,
                    cache_dir=config.database_dir,
                    lock_settings=config.database_lock_settings,
                    taxonomy_context_path=taxonomy_context,
                    taxonomy_cache_dir=(
                        config.database_dir / "taxonomy"
                    ),
                    source_group="primary",
                    step_name="functional-annotation-primary-search",
                )
            )
            if has_diamond_source
            else lambda: write_empty_diamond_results(
                output_path=diamond_hits,
                metadata_path=diamond_metadata,
                provenance_path=diamond_provenance,
                log_path=logs / "14b.functional-annotation-primary-search.log",
                metrics_path=diamond_metrics,
                step_name="functional-annotation-primary-search",
            ),
        )

        uniref90_hits = intermediate / "14c.functional-annotation.uniref90.tsv"
        uniref90_metadata = intermediate / "14c.functional-annotation.uniref90-metadata.tsv"
        uniref90_provenance = logs / "14c.functional-annotation-uniref90.database.json"
        uniref90_metrics = logs / "14c.functional-annotation-uniref90-search.metrics.json"
        uniref90_dependencies = [
            protein_fasta,
            diamond_hits,
            diamond_metadata,
            taxonomy_context,
            ctx.config_path,
            ctx.modules.functional_annotation,
            ctx.modules.annotation_taxonomy,
        ]
        if annotation_config.uniref90_fasta:
            uniref90_dependencies.append(
                expand_path(annotation_config.uniref90_fasta, config.base_dir)
            )
        ctx.run_step(
            name="14c.functional-annotation-uniref90-search",
            outputs=[
                uniref90_hits,
                uniref90_metadata,
                uniref90_provenance,
                logs / "14c.functional-annotation-uniref90-search.log",
                uniref90_metrics,
            ],
            dependencies=uniref90_dependencies,
            action=(
                lambda: run_diamond_annotation_search(
                    protein_fasta_path=protein_fasta,
                    output_path=uniref90_hits,
                    metadata_path=uniref90_metadata,
                    provenance_path=uniref90_provenance,
                    log_path=logs / "14c.functional-annotation-uniref90-search.log",
                    metrics_path=uniref90_metrics,
                    config=annotation_config,
                    base_dir=config.base_dir,
                    cache_dir=config.database_dir,
                    lock_settings=config.database_lock_settings,
                    taxonomy_context_path=taxonomy_context,
                    taxonomy_cache_dir=(
                        config.database_dir / "taxonomy"
                    ),
                    source_group="uniref90",
                    prior_similarity_inputs=((diamond_hits, diamond_metadata),),
                    step_name="functional-annotation-uniref90-search",
                )
            )
            if annotation_config.uniref90_enabled
            else lambda: write_empty_diamond_results(
                output_path=uniref90_hits,
                metadata_path=uniref90_metadata,
                provenance_path=uniref90_provenance,
                log_path=logs / "14c.functional-annotation-uniref90-search.log",
                metrics_path=uniref90_metrics,
                reason="UniRef90 fallback is disabled.",
                step_name="functional-annotation-uniref90-search",
            ),
        )

        pfam_hits = intermediate / "14d.functional-annotation.pfam.domtblout"
        pfam_metadata = intermediate / "14d.functional-annotation.pfam-metadata.tsv"
        pfam_provenance = logs / "14d.functional-annotation-pfam.database.json"
        pfam_metrics = logs / "14d.functional-annotation-pfam-search.metrics.json"
        pfam_dependencies: list[Path] = [
            protein_fasta,
            diamond_hits,
            diamond_metadata,
            uniref90_hits,
            uniref90_metadata,
            ctx.config_path,
            ctx.modules.functional_annotation,
        ]
        if annotation_config.pfam_hmm:
            pfam_dependencies.append(expand_path(annotation_config.pfam_hmm, config.base_dir))
        ctx.run_step(
            name="14d.functional-annotation-pfam-search",
            outputs=[
                pfam_hits,
                pfam_metadata,
                pfam_provenance,
                logs / "14d.functional-annotation-pfam-search.log",
                pfam_metrics,
            ],
            dependencies=pfam_dependencies,
            action=(
                lambda: run_pfam_domain_search(
                    protein_fasta_path=protein_fasta,
                    output_path=pfam_hits,
                    metadata_path=pfam_metadata,
                    provenance_path=pfam_provenance,
                    log_path=logs / "14d.functional-annotation-pfam-search.log",
                    metrics_path=pfam_metrics,
                    config=annotation_config,
                    base_dir=config.base_dir,
                    cache_dir=config.database_dir,
                    lock_settings=config.database_lock_settings,
                    diamond_hits_path=diamond_hits,
                    diamond_metadata_path=diamond_metadata,
                    additional_similarity_inputs=((uniref90_hits, uniref90_metadata),),
                )
            )
            if annotation_config.pfam_enabled
            else lambda: write_empty_pfam_results(
                output_path=pfam_hits,
                metadata_path=pfam_metadata,
                provenance_path=pfam_provenance,
                log_path=logs / "14d.functional-annotation-pfam-search.log",
                metrics_path=pfam_metrics,
            ),
        )

        cdd_hits = intermediate / "14e.functional-annotation.cdd.tsv"
        cdd_metadata = intermediate / "14e.functional-annotation.cdd-metadata.tsv"
        cdd_provenance = logs / "14e.functional-annotation-cdd.database.json"
        cdd_metrics = logs / "14e.functional-annotation-cdd-search.metrics.json"
        cdd_dependencies: list[Path] = [
            protein_fasta,
            diamond_hits,
            diamond_metadata,
            uniref90_hits,
            uniref90_metadata,
            ctx.config_path,
            ctx.modules.functional_annotation,
        ]
        if annotation_config.cdd_database:
            configured_cdd = expand_path(annotation_config.cdd_database, config.base_dir)
            cdd_dependencies.append(
                configured_cdd / "Cdd.aux"
                if configured_cdd.is_dir()
                else Path(str(configured_cdd) + ".aux")
            )
        if annotation_config.cdd_data_dir:
            cdd_dependencies.append(
                expand_path(annotation_config.cdd_data_dir, config.base_dir) / "cddid.tbl"
            )
        ctx.run_step(
            name="14e.functional-annotation-cdd-search",
            outputs=[
                cdd_hits,
                cdd_metadata,
                cdd_provenance,
                logs / "14e.functional-annotation-cdd-search.log",
                cdd_metrics,
            ],
            dependencies=cdd_dependencies,
            action=(
                lambda: run_cdd_domain_search(
                    protein_fasta_path=protein_fasta,
                    output_path=cdd_hits,
                    metadata_path=cdd_metadata,
                    provenance_path=cdd_provenance,
                    log_path=logs / "14e.functional-annotation-cdd-search.log",
                    metrics_path=cdd_metrics,
                    config=annotation_config,
                    base_dir=config.base_dir,
                    cache_dir=config.database_dir,
                    lock_settings=config.database_lock_settings,
                    similarity_inputs=(
                        (diamond_hits, diamond_metadata),
                        (uniref90_hits, uniref90_metadata),
                    ),
                )
            )
            if annotation_config.cdd_enabled
            else lambda: write_empty_cdd_results(
                output_path=cdd_hits,
                metadata_path=cdd_metadata,
                provenance_path=cdd_provenance,
                log_path=logs / "14e.functional-annotation-cdd-search.log",
                metrics_path=cdd_metrics,
            ),
        )

        functional_annotation_table = intermediate / "14f.annotation-table.functional.tsv"
        functional_evidence = ctx.outputs.final / "functional-annotation.tsv"
        functional_domain_comparison = (
            ctx.outputs.final / "functional-domain-search-comparison.tsv"
        )
        name_standardization_summary = (
            ctx.outputs.final / "functional-annotation-name-standardization.tsv"
        )
        functional_metrics = (
            logs
            / "14f.functional-annotation-assign-and-standardize-products.metrics.json"
        )
        ctx.run_step(
            name="14f.functional-annotation-assign-and-standardize-products",
            outputs=[
                functional_annotation_table,
                functional_evidence,
                functional_domain_comparison,
                name_standardization_summary,
                logs / "14f.functional-annotation-assign-and-standardize-products.log",
                functional_metrics,
            ],
            dependencies=[
                annotation_table,
                diamond_hits,
                diamond_metadata,
                uniref90_hits,
                uniref90_metadata,
                pfam_hits,
                pfam_metadata,
                cdd_hits,
                cdd_metadata,
                pfam_metrics,
                logs / "14d.functional-annotation-pfam-search.log",
                cdd_metrics,
                logs / "14e.functional-annotation-cdd-search.log",
                ctx.config_path,
                ctx.modules.functional_annotation,
                ctx.modules.product_names,
                taxonomy_context,
            ],
            action=lambda: apply_functional_annotations(
                annotation_table_path=annotation_table,
                diamond_hits_path=diamond_hits,
                diamond_metadata_path=diamond_metadata,
                uniref90_hits_path=uniref90_hits,
                uniref90_metadata_path=uniref90_metadata,
                pfam_hits_path=pfam_hits,
                pfam_metadata_path=pfam_metadata,
                cdd_hits_path=cdd_hits,
                cdd_metadata_path=cdd_metadata,
                output_path=functional_annotation_table,
                evidence_path=functional_evidence,
                log_path=(
                    logs
                    / "14f.functional-annotation-assign-and-standardize-products.log"
                ),
                metrics_path=functional_metrics,
                config=annotation_config,
                missing_product=config.pipeline.replace_product_with,
                domain_comparison_path=functional_domain_comparison,
                pfam_search_metrics_path=pfam_metrics,
                pfam_search_log_path=logs / "14d.functional-annotation-pfam-search.log",
                cdd_search_metrics_path=cdd_metrics,
                cdd_search_log_path=logs / "14e.functional-annotation-cdd-search.log",
                taxonomy_context_path=taxonomy_context,
                name_standardization_summary_path=name_standardization_summary,
            ),
        )
        annotation_table = functional_annotation_table

        if annotation_config.consistency.enabled:
            family_similarity = (
                intermediate / "14g.functional-annotation.family-similarity.tsv"
            )
            family_provenance = (
                logs / "14g.functional-annotation-family-search.database.json"
            )
            family_metrics = (
                logs / "14g.functional-annotation-family-search.metrics.json"
            )
            ctx.run_step(
                name="14g.functional-annotation-family-search",
                outputs=[
                    family_similarity,
                    family_provenance,
                    logs / "14g.functional-annotation-family-search.log",
                    family_metrics,
                ],
                dependencies=[
                    protein_fasta,
                    ctx.config_path,
                    ctx.modules.annotation_consistency,
                ],
                action=lambda: run_family_similarity_search(
                    protein_fasta_path=protein_fasta,
                    output_path=family_similarity,
                    provenance_path=family_provenance,
                    log_path=logs / "14g.functional-annotation-family-search.log",
                    metrics_path=family_metrics,
                    config=annotation_config,
                ),
            )

            consistent_annotation_table = (
                intermediate / "14h.annotation-table.consistent.tsv"
            )
            consistency_gene_table = (
                ctx.outputs.final / "functional-annotation-consistency.tsv"
            )
            consistency_family_table = (
                ctx.outputs.final / "functional-annotation-families.tsv"
            )
            consistency_pair_table = (
                ctx.outputs.final / "functional-annotation-consistency-pairs.tsv"
            )
            consistency_review_table = (
                ctx.outputs.final / "functional-annotation-conflicts.tsv"
            )
            consistency_summary_table = (
                ctx.outputs.final / "functional-annotation-consistency-summary.tsv"
            )
            consistency_source_pair_table = (
                ctx.outputs.final / "functional-annotation-source-pairs.tsv"
            )
            consistency_metrics = (
                logs / "14h.functional-annotation-consistency-audit.metrics.json"
            )
            ctx.run_step(
                name="14h.functional-annotation-consistency-audit",
                outputs=[
                    consistent_annotation_table,
                    consistency_gene_table,
                    consistency_family_table,
                    consistency_pair_table,
                    consistency_review_table,
                    consistency_summary_table,
                    consistency_source_pair_table,
                    logs / "14h.functional-annotation-consistency-audit.log",
                    consistency_metrics,
                ],
                dependencies=[
                    functional_annotation_table,
                    functional_evidence,
                    family_similarity,
                    ctx.config_path,
                    ctx.modules.annotation_consistency,
                    ctx.modules.product_names,
                ],
                action=lambda: audit_annotation_consistency(
                    annotation_table_path=functional_annotation_table,
                    evidence_path=functional_evidence,
                    similarity_path=family_similarity,
                    output_annotation_table_path=consistent_annotation_table,
                    gene_output_path=consistency_gene_table,
                    family_output_path=consistency_family_table,
                    pair_output_path=consistency_pair_table,
                    review_output_path=consistency_review_table,
                    summary_output_path=consistency_summary_table,
                    source_pair_output_path=consistency_source_pair_table,
                    log_path=logs / "14h.functional-annotation-consistency-audit.log",
                    metrics_path=consistency_metrics,
                    config=annotation_config,
                ),
            )
            annotation_table = consistent_annotation_table

    raw_gff2mss = intermediate / "15.mss.raw.txt"
    gff2mss_metrics = logs / "15.gff2mss.metrics.json"
    ctx.run_step(
        name="15.gff2mss",
        outputs=[raw_gff2mss, logs / "15.gff2mss.log", gff2mss_metrics],
        dependencies=[
            inputs.downstream_fasta,
            prepared_gff.padded_gff,
            annotation_table,
            ctx.config_path,
            ctx.modules.pipeline_actions,
            ctx.modules.mss_converter_core,
            ctx.modules.mss_converter_features,
            ctx.modules.mss_converter_gaps,
            ctx.modules.mss_converter_models,
            ctx.modules.mss_converter_render,
        ],
        action=lambda: run_gff2mss(
            config,
            fasta_path=inputs.downstream_fasta,
            gff_path=prepared_gff.padded_gff,
            ann_table_path=annotation_table,
            output_path=raw_gff2mss,
            log_path=logs / "15.gff2mss.log",
            metrics_path=gff2mss_metrics,
        ),
    )

    postprocessed_gff2mss = intermediate / "16.mss.cds-to-misc.txt"
    misc_converted_gene_ids = logs / "16.mss-cds-to-misc.changed-gene-ids.txt"
    misc_metrics = logs / "16.mss-cds-to-misc.metrics.json"
    ctx.run_step(
        name="16.mss-cds-to-misc",
        outputs=[
            postprocessed_gff2mss,
            logs / "16.mss-cds-to-misc.log",
            misc_converted_gene_ids,
            misc_metrics,
        ],
        dependencies=[raw_gff2mss, prepared_gff.genes_with_stops, ctx.modules.mss_postprocess, ctx.modules.pipeline_actions],
        action=lambda: run_mss_cds_to_misc(
            mss_input_path=raw_gff2mss,
            genes_input_path=prepared_gff.genes_with_stops,
            mss_output_path=postprocessed_gff2mss,
            log_path=logs / "16.mss-cds-to-misc.log",
            converted_gene_ids_path=misc_converted_gene_ids,
            metrics_path=misc_metrics,
            locus_tag_prefix=config.sample.locus_tag,
        ),
    )

    return AnnotationArtifacts(
        annotation_table=annotation_table,
        postprocessed_gff2mss=postprocessed_gff2mss,
    )


def run_pipeline(config_file: str | Path, *, validate: bool = True) -> PipelineOutputs:
    config, config_path, outputs, modules, ctx = _initialize_pipeline(config_file)
    manifest = ctx.manifest
    try:
        prepared_inputs = _prepare_inputs(ctx)
        prepared_gff = _prepare_gff(ctx, prepared_inputs)
        annotation_artifacts = _build_annotation_artifacts(ctx, prepared_inputs, prepared_gff)

        ctx.run_step(
            name="17.render-final-annotation",
            outputs=[
                outputs.ann_path,
                outputs.logs / "17.render-final-annotation.log",
                outputs.logs / "17.render-final-annotation.metrics.json",
            ],
            dependencies=[
                annotation_artifacts.postprocessed_gff2mss,
                config_path,
                modules.pipeline,
                modules.submission_render,
                modules.pipeline_actions,
            ],
            action=lambda: render_final_annotation_with_log(
                header=render_header(config),
                body_path=annotation_artifacts.postprocessed_gff2mss,
                output_path=outputs.ann_path,
                log_path=outputs.logs / "17.render-final-annotation.log",
                metrics_path=outputs.logs / "17.render-final-annotation.metrics.json",
            ),
        )

        validation_options = ValidationOptions.from_config(config)
        validation_artifacts = ValidationArtifacts(
            ann_path=outputs.ann_path,
            fasta_path=outputs.fasta_path,
            logs_dir=outputs.logs,
            output_dir=outputs.final,
        )
        validation_requested = (
            validation_options.run_parser or validation_options.run_transchecker
        )
        installed_tools = (
            list_installed(config.cache_dir) if validate and validation_requested else {}
        )
        tool_details = {
            component: describe_installation(installed_tools[component])
            for component in ("parser", "transchecker")
            if component in installed_tools
        }
        validation_manifest_options = {
            "run_parser": validation_options.run_parser,
            "run_transchecker": validation_options.run_transchecker,
            "parallel": validation_options.parallel,
            "heap": validation_options.heap,
            "java": validation_options.java_cmd,
            "tools": tool_details,
        }
        if validate:
            manifest.set_validation(
                enabled=validation_requested,
                result_paths={
                    "validation_summary": validation_artifacts.validation_summary,
                },
                options=validation_manifest_options,
            )
            validation_outputs = run_validation(
                options=validation_options,
                artifacts=validation_artifacts,
            )
            manifest.set_validation(
                enabled=validation_requested,
                result_paths=validation_outputs,
                options=validation_manifest_options,
            )
        else:
            write_validation_not_run_summary(artifacts=validation_artifacts)
            manifest.set_validation(
                enabled=False,
                result_paths={"validation_summary": validation_artifacts.validation_summary},
                options=validation_manifest_options,
            )

        manifest.mark_completed()
        return outputs
    except Exception as exc:
        manifest.mark_failed(exc)
        raise
    finally:
        manifest.write(outputs.manifest_path)
