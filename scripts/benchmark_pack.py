#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import shutil
import statistics
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from msspack.config import load_config  # noqa: E402
from msspack.pipeline import run_pipeline  # noqa: E402


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Benchmark msspack pack runs against a config file.",
    )
    parser.add_argument("--config", required=True, help="Path to msspack TOML config")
    parser.add_argument(
        "--repeats",
        type=int,
        default=3,
        help="Number of timed runs to execute",
    )
    parser.add_argument(
        "--no-validate",
        action="store_true",
        help="Skip Parser/transChecker during the benchmark",
    )
    parser.add_argument(
        "--clean-first",
        action="store_true",
        help="Remove the build directory before the first run",
    )
    parser.add_argument(
        "--clean-between-runs",
        action="store_true",
        help="Remove the build directory before every timed run",
    )
    parser.add_argument(
        "--json-out",
        help="Optional path for a JSON summary",
    )
    return parser


def remove_output_dir(config_path: Path) -> None:
    config = load_config(config_path)
    if config.output_dir.exists():
        shutil.rmtree(config.output_dir)


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    config_path = Path(args.config).expanduser().resolve()

    if args.clean_first or args.clean_between_runs:
        remove_output_dir(config_path)

    runs: list[float] = []
    for run_index in range(args.repeats):
        if run_index > 0 and args.clean_between_runs:
            remove_output_dir(config_path)
        started = time.perf_counter()
        run_pipeline(config_path, validate=not args.no_validate)
        elapsed = time.perf_counter() - started
        runs.append(elapsed)
        print(f"run {run_index + 1}: {elapsed:.3f}s")

    summary = {
        "config": str(config_path),
        "repeats": args.repeats,
        "validate": not args.no_validate,
        "clean_first": args.clean_first,
        "clean_between_runs": args.clean_between_runs,
        "runs_seconds": runs,
        "min_seconds": min(runs),
        "max_seconds": max(runs),
        "mean_seconds": statistics.mean(runs),
        "median_seconds": statistics.median(runs),
    }
    print(json.dumps(summary, indent=2, sort_keys=True))

    if args.json_out:
        output_path = Path(args.json_out).expanduser().resolve()
        output_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"wrote {output_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
