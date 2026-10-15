#!/usr/bin/env python3
"""Offline sequence, annotation, cache and memory benchmarks on synthetic inputs."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import resource
import statistics
import sys
import tempfile
import time
import tracemalloc
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))


@dataclass
class Workload:
    action: Callable[[], object]
    outputs: Callable[[], dict[str, Path]]
    semantic: Callable[[], object] = lambda: None

    def checksum(self) -> str:
        hashes = {}
        for name, path in self.outputs().items():
            digest = hashlib.sha256()
            with path.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
            hashes[name] = digest.hexdigest()
        payload = {"files": hashes, "semantic": self.semantic()}
        return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def prepare_models(base: Path, genes: int, exons: int) -> None:
    from msspack.demo import write_demo_dataset

    write_demo_dataset(base)
    (base / "genome.fa").write_text(">chr1\n" + "A" * (genes * exons * 120) + "\n")
    with (base / "annotation.gff3").open("w") as gff, (base / "annotation.tsv").open("w") as ann:
        gff.write("##gff-version 3\n")
        ann.write("ID\tDescription\tLocus_tag\n")
        for gene in range(genes):
            start, end = gene * exons * 120 + 1, (gene + 1) * exons * 120
            strand = "-" if gene % 2 else "+"
            gff.write(f"chr1\t.\tgene\t{start}\t{end}\t.\t{strand}\t.\tID=g{gene}\n")
            gff.write(f"chr1\t.\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID=t{gene};Parent=g{gene}\n")
            utr_start = start if strand == "+" else end - 32
            gff.write(f"chr1\t.\tfive_prime_UTR\t{utr_start}\t{utr_start+2}\t.\t{strand}\t.\tParent=t{gene}\n")
            for exon in range(exons):
                left, right = start + exon * 120, start + exon * 120 + 89
                gff.write(f"chr1\t.\texon\t{left}\t{right}\t.\t{strand}\t.\tID=e{gene}_{exon};Parent=t{gene}\n")
                cds_start = left + (3 if strand == "+" and exon == 0 else 0)
                cds_end = right - (3 if strand == "-" and exon == exons - 1 else 0)
                gff.write(f"chr1\t.\tCDS\t{cds_start}\t{cds_end}\t.\t{strand}\t0\tID=c{gene}_{exon};Parent=t{gene}\n")
            ann.write(f"t{gene}\ttest protein\tTEST_g{gene}\n")


def prepare_case(base: Path, args: argparse.Namespace) -> Workload:
    size, case = args.size, args.case
    paths: dict[str, Path] = {}
    values: dict[str, object] = {}
    if case in {"exons", "annotation", "containment", "counts", "pack"}:
        prepare_models(base, size, 1 if case == "annotation" else args.exons)
    if case == "padding":
        from msspack.padding_tools import _process_padding

        sequence = "ATG" + "AAA" * 998 + "TAA"
        inputs = [(str(i), sequence if i % 2 == 0 else sequence[:1500] + "TAG" + sequence[1503:])
                  for i in range(size)]

        def action() -> object:
            values["padding"] = [_process_padding(name, seq, "1") for name, seq in inputs]

        return Workload(action, lambda: {}, lambda: values["padding"])
    if case == "containment":
        from msspack.gff_feature_sync import validate_parent_child_containment

        rows = [line.split("\t") for line in (base / "annotation.gff3").read_text().splitlines()
                if not line.startswith("#")]
        rows[4][4] = "999999999"
        rows[2][4] = str(int(rows[3][4]) + 5)

        def action() -> object:
            values["issues"] = [vars(issue) for issue in validate_parent_child_containment(rows)]

        return Workload(action, lambda: {}, lambda: values["issues"])
    if case == "annotation":
        from msspack.annotation_table import build_annotation_table

        return Workload(lambda: build_annotation_table(
            gff_path=base / "annotation.gff3", output_path=base / "out.tsv",
            locus_tag_prefix="TEST", log_path=base / "log",
        ), lambda: {"table": base / "out.tsv"})
    if case == "exons":
        from msspack.mss_converter import ConversionOptions, convert_gff_to_mss

        options = ConversionOptions(
            fasta_path=base / "genome.fa", gff_path=base / "annotation.gff3",
            annotation_path=base / "annotation.tsv", output_path=base / "out.ann",
            locus_tag_prefix="TEST", organism_name="Msspackia fictitia",
        )
        return Workload(lambda: convert_gff_to_mss(options), lambda: {"annotation": options.output_path})
    if case == "counts":
        from msspack.pipeline_actions import run_apply_padding_to_gff

        (base / "padding.log").write_text("")

        def counts() -> object:
            metrics = json.loads((base / "metrics.json").read_text())
            return {key: metrics[key] for key in ("input_total", "changed_total", "output_total")}

        return Workload(lambda: run_apply_padding_to_gff(
            gff_path=base / "annotation.gff3", padding_log_path=base / "padding.log",
            output_path=base / "out.gff", genes_with_stops_path=base / "stops",
            updated_genes_path=base / "updated", command_log_path=base / "log",
            metrics_path=base / "metrics.json",
        ), lambda: {"gff": base / "out.gff", "stops": base / "stops", "updated": base / "updated"}, counts)
    if case == "consistency":
        from msspack.annotation_consistency import SIMILARITY_FIELDS, audit_annotation_consistency
        from msspack.config_models import (
            FunctionalAnnotationConfig,
            FunctionalAnnotationConsistencyConfig,
        )

        names = ("protein kinase domain-containing protein", "serine threonine protein kinase",
                 "cellulose synthase-like protein B4", "cellulose synthase-like protein H1",
                 "BTB/POZ domain-containing protein", "BTB domain-containing protein")
        with (base / "annotation.tsv").open("w") as ann, (base / "evidence.tsv").open("w") as evidence:
            ann.write("ID\tDescription\tLocus_tag\n")
            evidence.write("ID\tLocus_tag\tassigned_product\tsource\tconfidence\tbitscore\tconsensus_score\n")
            for i in range(size):
                product = names[i % len(names)] + f" marker {i}"
                ann.write(f"t{i}\t{product}\tG{i}\n")
                evidence.write(f"t{i}\tG{i}\t{product}\tpfam\tlow\t200\t0.7\n")
        with (base / "similarity.tsv").open("w") as similarity:
            similarity.write("\t".join(SIMILARITY_FIELDS) + "\n")
            for i in range(size):
                for other in range(i + 1, min(size, i + 11)):
                    similarity.write(f"t{i}\tt{other}\t75\t90\t90\t1e-30\t200\t300\t300\n")
        outputs = {name: base / (name + ".tsv") for name in (
            "output_annotation_table_path", "gene_output_path", "family_output_path", "pair_output_path",
            "review_output_path", "summary_output_path", "source_pair_output_path",
        )}
        config = FunctionalAnnotationConfig(consistency=FunctionalAnnotationConsistencyConfig(enabled=True))
        return Workload(lambda: audit_annotation_consistency(
            annotation_table_path=base / "annotation.tsv", evidence_path=base / "evidence.tsv",
            similarity_path=base / "similarity.tsv", config=config,
            log_path=base / "log", metrics_path=base / "metrics.json", **outputs,
        ), lambda: outputs)
    if case == "descendants":
        from msspack.transcript_selection import select_one_mrna_per_gene

        rows = [f"chr1\t.\tregion\t3\t8\t.\t+\t.\tID=child{i};Parent="
                + (f"child{i-1}" if i else "drop") for i in reversed(range(size))]
        rows += ["chr1\t.\tgene\t1\t9\t.\t+\t.\tID=g",
                 "chr1\t.\tmRNA\t1\t9\t.\t+\t.\tID=keep;Parent=g",
                 "chr1\t.\tCDS\t1\t9\t.\t+\t0\tParent=keep",
                 "chr1\t.\tmRNA\t2\t8\t.\t+\t.\tID=drop;Parent=g",
                 "chr1\t.\tCDS\t3\t8\t.\t+\t0\tParent=drop",
                 f"chr1\t.\tregion\t3\t8\t.\t+\t.\tID=shared;Parent=child{size-1},keep"]
        (base / "in.gff").write_text("\n".join(rows) + "\n")
        return Workload(lambda: select_one_mrna_per_gene(
            input_gff_path=base / "in.gff", output_gff_path=base / "out.gff",
        ), lambda: {"gff": base / "out.gff"})
    if case == "memory":
        from msspack.gap_normalization import normalize_gap_lengths

        with (base / "in.fa").open("w") as fasta:
            for i in range(size):
                sequence = "A" * (1024 * 1024)
                if args.late_edit and i == size - 1:
                    sequence = sequence[:500] + "NNNN" + sequence[504:]
                fasta.write(f">contig{i}\n{sequence}\n")
        return Workload(lambda: normalize_gap_lengths(
            fasta_path=base / "in.fa", output_fasta_path=base / "out.fa",
            log_path=base / "log", gap_len=2,
        ), lambda: {"fasta": base / "out.fa"})
    if case == "cache":
        from msspack.execution import is_up_to_date, run_if_needed

        run_number = 0

        def action() -> object:
            nonlocal run_number
            paths.clear()
            for i in range(size):
                output = base / f"output-{run_number}-{i}"
                run_if_needed(outputs=[output], dependencies=[],
                              action=lambda output=output: output.write_text("same output\n"),
                              cache_key={"step": i})
                paths[str(i)] = output
            run_number += 1

        def cache_status() -> object:
            fresh = [is_up_to_date([path], [], cache_key={"step": int(i)}) for i, path in paths.items()]
            if not all(fresh):
                raise RuntimeError("Benchmark produced a non-reusable cache entry")
            return fresh

        return Workload(action, lambda: paths, cache_status)
    if case == "workflow":
        from msspack.demo import write_demo_dataset
        from msspack.workflow import run_all

        write_demo_dataset(base)

        def action() -> object:
            result = run_all(base / "config.toml", run_busco=False, validate=False)
            paths.update(annotation=result.pipeline.ann_path, fasta=result.pipeline.fasta_path,
                         svg=result.plots.gene_flow_svg)

        return Workload(action, lambda: paths)
    from msspack.pipeline import run_pipeline

    template = (base / "config.toml").read_text()
    run_number = 0

    def action() -> object:
        nonlocal run_number
        config = base / "config.toml"
        config.write_text(template.replace("build/MSSPackDemo", f"build/run{run_number}"))
        run_number += 1
        result = run_pipeline(config, validate=False)
        paths.update(annotation=result.ann_path, fasta=result.fasta_path)

    return Workload(action, lambda: paths)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True, choices=(
        "padding", "containment", "annotation", "exons", "counts", "consistency",
        "descendants", "cache", "workflow", "memory", "pack",
    ))
    parser.add_argument("--size", type=int, default=1000)
    parser.add_argument("--exons", type=int, default=10)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--late-edit", action="store_true", help="Edit the last contig in memory case")
    args = parser.parse_args()
    if min(args.size, args.exons, args.repeats) <= 0:
        parser.error("--size, --exons and --repeats must be positive")
    with tempfile.TemporaryDirectory(prefix="msspack-followup-") as directory:
        workload = prepare_case(Path(directory), args)
        workload.action()
        expected = workload.checksum()
        wall, cpu, allocations = [], [], []
        for _ in range(args.repeats):
            if args.case == "memory":
                tracemalloc.start()
            started, cpu_started = time.perf_counter(), time.process_time()
            workload.action()
            cpu.append(time.process_time() - cpu_started)
            wall.append(time.perf_counter() - started)
            if args.case == "memory":
                allocations.append(tracemalloc.get_traced_memory()[1] / 1024**2)
                tracemalloc.stop()
            if workload.checksum() != expected:
                raise RuntimeError("Output changed between benchmark repetitions")
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        print(json.dumps({
            "case": args.case, "size": args.size, "exons": args.exons, "late_edit": args.late_edit,
            "repeats": args.repeats, "warmups": 1, "wall_seconds": wall, "cpu_seconds": cpu,
            "median_wall_seconds": statistics.median(wall), "median_cpu_seconds": statistics.median(cpu),
            "python_allocation_peaks_mib": allocations,
            "process_peak_rss_mib": rss / (1024**2 if sys.platform == "darwin" else 1024),
            "output_sha256": expected, "python": platform.python_version(),
            "platform": platform.platform(), "hash_seed": os.environ.get("PYTHONHASHSEED"),
        }, indent=2))


if __name__ == "__main__":
    main()
