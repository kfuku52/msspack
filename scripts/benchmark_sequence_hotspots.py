#!/usr/bin/env python3
"""Repeatable offline benchmarks on synthetic, non-submittable gene models."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import resource
import statistics
import sys
import tempfile
import time
from collections.abc import Callable
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from msspack.demo import write_demo_dataset  # noqa: E402
from msspack.gap_normalization import normalize_gap_lengths  # noqa: E402
from msspack.gff import child_ids, parse_attributes  # noqa: E402
from msspack.mss_converter import ConversionOptions, convert_gff_to_mss  # noqa: E402
from msspack.pipeline import PipelineOutputs, run_pipeline  # noqa: E402


def prepare_inputs(base: Path, size: int) -> list[str]:
    write_demo_dataset(base)
    sequence = "ATG" + "AAA" * 98 + "TAA" + "N" * 20
    (base / "genome.fa").write_text(
        ">synthetic_NOT_FOR_SUBMISSION\n" + (sequence + "\n") * size, encoding="utf-8",
    )
    attributes: list[str] = []
    rows = ["##gff-version 3"]
    annotations = ["ID\tDescription\tLocus_tag"]
    for index in range(size):
        start, end = index * len(sequence) + 1, index * len(sequence) + 300
        gene_id, transcript_id = f"g{index}", f"t{index}"
        attrs = (
            f"ID={gene_id}",
            f"ID={transcript_id};Parent={gene_id};product=test protein",
            f"ID=e{index};Parent={transcript_id}",
            f"ID=c{index};Parent={transcript_id};product=test protein",
        )
        for feature, phase, attr in zip(
            ("gene", "mRNA", "exon", "CDS"), (".", ".", ".", "0"), attrs, strict=True,
        ):
            rows.append(
                f"synthetic_NOT_FOR_SUBMISSION\t.\t{feature}\t{start}\t{end}\t.\t+\t{phase}\t{attr}"
            )
        attributes.extend(attrs)
        # Include escaped list members and UTF-8 alongside the common ASCII path.
        attributes.append(f"ID=shared{index};Parent=t{index}%2Cextra,t{index};Note=caf%C3%A9")
        annotations.append(f"{transcript_id}\ttest protein\tTEST_{gene_id}")
    (base / "annotation.gff3").write_text("\n".join(rows) + "\n", encoding="utf-8")
    (base / "annotation.tsv").write_text("\n".join(annotations) + "\n", encoding="utf-8")
    config = base / "config.toml"
    config.write_text(
        config.read_text(encoding="utf-8").replace(
            "run_gapjust = false", "run_gapjust = true\ngapjust_min = 0",
        ),
        encoding="utf-8",
    )
    return attributes


def benchmark(case: str, base: Path, size: int, repeats: int) -> dict[str, object]:
    attributes = prepare_inputs(base, size)
    action: Callable[[], object]
    checksum: Callable[[object], str]

    def file_checksum(paths: list[Path]) -> str:
        digest = hashlib.sha256()
        for path in paths:
            digest.update(path.read_bytes())
        return digest.hexdigest()

    if case == "attributes":
        def action() -> object:
            return [parse_attributes(attr) for attr in attributes * 10]

        def checksum(result: object) -> str:
            # Keep escaped Parent list boundaries in the equivalence check.
            assert isinstance(result, list)
            payload = [(attrs, child_ids(attrs.get("Parent"))) for attrs in result]
            return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    elif case == "gap-normalization":
        def action() -> object:
            normalize_gap_lengths(
                fasta_path=base / "genome.fa", output_fasta_path=base / "normalized.fa",
                input_gff_path=base / "annotation.gff3", output_gff_path=base / "normalized.gff3",
                log_path=base / "gapjust.log", gap_len=100,
            )

        def checksum(result: object) -> str:
            return file_checksum([base / "normalized.fa", base / "normalized.gff3"])
    elif case == "conversion":
        options = ConversionOptions(
            fasta_path=base / "genome.fa", gff_path=base / "annotation.gff3",
            annotation_path=base / "annotation.tsv", output_path=base / "output.ann.txt",
            locus_tag_prefix="TEST", organism_name="Msspackia fictitia",
        )

        def action() -> object:
            return convert_gff_to_mss(options)

        def checksum(result: object) -> str:
            return file_checksum([options.output_path])
    else:
        # Each repetition uses a fresh output root; no user outputs are removed.
        config_template = (base / "config.toml").read_text(encoding="utf-8")
        run_number = 0

        def action() -> object:
            nonlocal run_number
            config = base / "config.toml"
            config.write_text(
                config_template.replace("build/MSSPackDemo", f"build/run{run_number}"),
                encoding="utf-8",
            )
            run_number += 1
            return run_pipeline(config, validate=False)

        def checksum(result: object) -> str:
            assert isinstance(result, PipelineOutputs)
            return file_checksum([result.ann_path, result.fasta_path])

    expected = checksum(action())  # Untimed warmup, including imports and filesystem caches.
    elapsed: list[float] = []
    cpu_elapsed: list[float] = []
    for _ in range(repeats):
        started = time.perf_counter()
        cpu_started = time.process_time()
        result = action()
        cpu_elapsed.append(time.process_time() - cpu_started)
        elapsed.append(time.perf_counter() - started)
        if checksum(result) != expected:
            raise RuntimeError("Benchmark output changed between repetitions")
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return {
        "case": case, "genes": size, "gff_records": size * 4, "gaps": size,
        "repeats": repeats, "warmups": 1, "seconds": elapsed,
        "median_seconds": statistics.median(elapsed), "output_sha256": expected,
        "process_cpu_seconds": cpu_elapsed,
        "median_process_cpu_seconds": statistics.median(cpu_elapsed),
        "process_peak_rss_mib": peak_rss / (1024**2 if sys.platform == "darwin" else 1024),
        "python": platform.python_version(), "platform": platform.platform(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case", required=True,
        choices=("attributes", "gap-normalization", "conversion", "pack"),
    )
    parser.add_argument("--size", type=int, default=2000, help="Number of genes and N runs")
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args()
    if args.size <= 0 or args.repeats <= 0:
        parser.error("--size and --repeats must be positive")
    with tempfile.TemporaryDirectory(prefix="msspack-hotspots-") as tmp:
        print(json.dumps(benchmark(args.case, Path(tmp), args.size, args.repeats), indent=2))


if __name__ == "__main__":
    main()
