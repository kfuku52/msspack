from __future__ import annotations

import argparse
import importlib.resources
import sys
from pathlib import Path

from . import __version__
from .busco import cleanup_busco_cache, run_busco_comparison, summarize_busco_artifacts
from .config import load_config
from .ddbj_tools import DDBJ_LICENSE_URL, install_component, list_installed
from .doctor import doctor_succeeded, render_doctor_report, run_doctor
from .internal_cli import add_internal_parser, handle_internal
from .pipeline import run_pipeline
from .pipeline_plots import run_pipeline_plots, summarize_pipeline_plots
from .report import run_html_report
from .utils import MSSPackError, write_text
from .validation import validate_existing


def _example_config_text() -> str:
    template = importlib.resources.files("msspack").joinpath("templates/msspack.example.toml")
    return (
        "# Replace the placeholder paths and submission metadata for your project.\n"
        + template.read_text(encoding="utf-8")
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="msspack")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    init_parser = subparsers.add_parser("init", help="write a starter config")
    init_parser.add_argument("output", help="path to the new TOML config")
    init_parser.add_argument(
        "--force",
        action="store_true",
        help="overwrite an existing config file",
    )

    doctor_parser = subparsers.add_parser("doctor", help="check runtime dependencies")
    doctor_parser.add_argument("--config", help="optional config TOML")

    tools_parser = subparsers.add_parser("tools", help="manage DDBJ tools")
    tools_subparsers = tools_parser.add_subparsers(dest="tools_command", required=True)

    tools_install = tools_subparsers.add_parser("install", help="download DDBJ tools")
    tools_install.add_argument(
        "components",
        nargs="*",
        default=["parser", "transchecker"],
        help="components to install",
    )
    tools_install.add_argument("--reinstall", action="store_true")
    tools_install.add_argument("--cache-dir")

    tools_list = tools_subparsers.add_parser("list", help="show installed DDBJ tools")
    tools_list.add_argument("--cache-dir")

    pack_parser = subparsers.add_parser("pack", help="build MSS submission files")
    pack_parser.add_argument("--config", required=True, help="config TOML")
    pack_parser.add_argument("--no-validate", action="store_true")

    plot_parser = subparsers.add_parser(
        "plot",
        help="render stage-wise pipeline Sankey and event-count plots from existing build logs",
    )
    plot_parser.add_argument("--config", required=True, help="config TOML")
    plot_parser.add_argument("--force", action="store_true", help="rerender plot outputs")

    report_parser = subparsers.add_parser(
        "report",
        help="render an HTML summary report from an existing build, plots, validation, and BUSCO outputs",
    )
    report_parser.add_argument("--config", required=True, help="config TOML")
    report_parser.add_argument("--force", action="store_true", help="rerender plot outputs before writing the report")

    busco_parser = subparsers.add_parser(
        "busco",
        help="run BUSCO on genome FASTA and GFF-derived CDS inputs/outputs and compare them",
    )
    busco_parser.add_argument("--config", required=True, help="config TOML")
    busco_parser.add_argument("--lineage-dataset", help="override BUSCO lineage dataset")
    busco_parser.add_argument("--threads", type=int, help="override BUSCO CPU count")
    busco_parser.add_argument("--force", action="store_true", help="rerun BUSCO even if cached")
    busco_parser.add_argument(
        "--clean-cache",
        action="store_true",
        help="remove temporary BUSCO staging workspaces under the msspack cache before running",
    )
    busco_parser.add_argument(
        "--genome",
        dest="run_genome",
        action="store_true",
        default=None,
        help="also run BUSCO on genome FASTA",
    )
    busco_parser.add_argument(
        "--no-genome",
        dest="run_genome",
        action="store_false",
        help="skip BUSCO on genome FASTA even if enabled in the config",
    )
    busco_parser.add_argument(
        "--cds",
        dest="run_cds",
        action="store_true",
        default=None,
        help="run BUSCO on GFF-derived CDS FASTA",
    )
    busco_parser.add_argument(
        "--no-cds",
        dest="run_cds",
        action="store_false",
        help="skip BUSCO on GFF-derived CDS FASTA",
    )
    busco_parser.add_argument(
        "--auto-lineage",
        dest="auto_lineage",
        action="store_true",
        default=None,
        help="override config to auto-select lineage from the input FASTA",
    )
    busco_parser.add_argument(
        "--no-auto-lineage",
        dest="auto_lineage",
        action="store_false",
        help="override config to require an explicit lineage dataset",
    )
    busco_parser.add_argument(
        "--auto-lineage-scope",
        choices=["all", "euk", "prok"],
        help="override BUSCO auto-lineage scope",
    )

    validate_parser = subparsers.add_parser(
        "validate", help="run Parser/transChecker on existing files"
    )
    validate_parser.add_argument("--config", help="optional config TOML")
    validate_parser.add_argument("--ann", required=True, help="annotation file")
    validate_parser.add_argument("--fasta", required=True, help="sequence file")

    add_internal_parser(subparsers)

    return parser


def _handle_init(args: argparse.Namespace) -> int:
    output_path = Path(args.output).expanduser().resolve()
    if output_path.exists() and not args.force:
        raise MSSPackError(
            f"Refusing to overwrite existing config: {output_path}. Pass --force to replace it."
        )
    write_text(output_path, _example_config_text())
    print(output_path)
    return 0


def _handle_doctor(args: argparse.Namespace) -> int:
    config = load_config(args.config) if args.config else None
    checks = run_doctor(config)
    print(render_doctor_report(config, checks=checks))
    return 0 if doctor_succeeded(checks) else 1


def _handle_tools(args: argparse.Namespace) -> int:
    if args.tools_command == "install":
        print(
            "DDBJ validation tools are separately licensed; downloading or using them "
            f"is subject to the DDBJ agreement: {DDBJ_LICENSE_URL}",
            file=sys.stderr,
        )
        for component in args.components:
            item = install_component(
                component,
                cache_dir=args.cache_dir,
                reinstall=args.reinstall,
            )
            print(f"{item.component}\t{item.version}\t{item.root}")
        return 0

    if args.tools_command == "list":
        for item in list_installed(args.cache_dir).values():
            print(f"{item.component}\t{item.version}\t{item.root}")
        return 0

    raise MSSPackError(f"Unsupported tools command: {args.tools_command}")


def _handle_pack(args: argparse.Namespace) -> int:
    outputs = run_pipeline(args.config, validate=not args.no_validate)
    print(outputs.ann_path)
    print(outputs.fasta_path)
    return 0


def _handle_plot(args: argparse.Namespace) -> int:
    artifacts = run_pipeline_plots(args.config, force=args.force)
    print(artifacts.root)
    print(artifacts.summary_tsv)
    print(artifacts.gene_flow_svg)
    print(artifacts.gene_flow_pdf)
    print(artifacts.event_counts_svg)
    print(artifacts.event_counts_pdf)
    print(artifacts.overlap_svg)
    print(artifacts.overlap_pdf)
    if artifacts.name_consistency_svg.exists():
        print(artifacts.name_consistency_svg)
        print(artifacts.name_consistency_pdf)
        print(artifacts.source_consistency_svg)
        print(artifacts.source_consistency_pdf)
    for line in summarize_pipeline_plots(artifacts):
        print(line)
    return 0


def _handle_report(args: argparse.Namespace) -> int:
    artifacts = run_html_report(args.config, force=args.force)
    print(artifacts.root)
    print(artifacts.index_html)
    return 0


def _handle_busco(args: argparse.Namespace) -> int:
    cleaned_root = cleanup_busco_cache() if args.clean_cache else None
    outputs = run_busco_comparison(
        args.config,
        lineage_dataset=args.lineage_dataset or "",
        threads=args.threads,
        force=args.force,
        clean_cache=False,
        auto_lineage=args.auto_lineage,
        auto_lineage_scope=args.auto_lineage_scope or "",
        run_genome=args.run_genome,
        run_cds=args.run_cds,
    )
    if cleaned_root is not None:
        print(f"cleaned_staging_cache\t{cleaned_root}")
    print(outputs.root)
    if outputs.cds is not None:
        print(outputs.cds.comparison_tsv)
        print(outputs.cds.comparison_svg)
        print(outputs.cds.comparison_pdf)
    if outputs.genome is not None:
        print(outputs.genome.comparison_tsv)
        print(outputs.genome.comparison_svg)
        print(outputs.genome.comparison_pdf)
    for line in summarize_busco_artifacts(outputs):
        print(line)
    return 0


def _handle_validate(args: argparse.Namespace) -> int:
    outputs = validate_existing(
        config_file=args.config,
        ann_path=args.ann,
        fasta_path=args.fasta,
    )
    for path in outputs.values():
        print(path)
    return 0


def _dispatch(args: argparse.Namespace) -> int:
    handlers = {
        "init": _handle_init,
        "doctor": _handle_doctor,
        "tools": _handle_tools,
        "pack": _handle_pack,
        "plot": _handle_plot,
        "report": _handle_report,
        "busco": _handle_busco,
        "validate": _handle_validate,
        "internal": handle_internal,
    }
    handler = handlers.get(args.command)
    if handler is None:
        raise MSSPackError(f"Unsupported command: {args.command}")
    return handler(args)


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        return _dispatch(args)
    except MSSPackError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    except (KeyError, OSError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
