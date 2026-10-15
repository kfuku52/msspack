# Sequence processing performance review

This review covers the CLI/config entry points, workflow/stage execution and
cache identity, FASTA/GFF parsing, transcript selection and adjustment, MSS
rendering, functional annotation/database wrappers, plots/reports, and the local
test/packaging workflow. Changes target measured Python bottlenecks while
preserving the existing biological and submission contracts.

## Findings and changes

| Area | Finding and action |
| --- | --- |
| Gap normalization | Each feature boundary sorted all gap edits and scanned preceding edits. Build sorted edits and cumulative offsets once per sequence, then use binary search. Coordinate queries fall from O(G log G) to O(log G), with O(G) index storage. |
| MSS rendering | Each exon/CDS segment checked every gap on its contig, including non-overlapping gaps. Index the sorted, disjoint N runs once per contig and query only overlapping runs in O(log G + K). The existing plain-list calling convention remains supported. |
| GFF attributes | Common unescaped ASCII values went through UTF-8 encoding/decoding; list attributes were decoded twice. Bypass encoding for plain ASCII, decode each list member once, and reuse parsed attributes during coordinate deduplication. Percent-escape, UTF-8, duplicate-key, control-character and list-boundary validation remain intact. |
| Pipeline/cache/publication | The demo profile shows repeated implementation fingerprints are significant on tiny inputs. They also detect source/dependency changes during a stage. Preserve those checks and the publication/locking contract. |
| Transcript handling | CDS extraction already groups models by sequence. Descendant-removal rescans and repeated transcript-coverage calculations are further candidates for deep hierarchies or very large exon counts; this benchmark does not establish their cost. Preserve selection, strand, phase and genetic-code policies. |
| External annotation/databases | DIAMOND, HMMER, CDD and BUSCO already have index reuse, query filtering and integrity checks. No large databases or external-tool workloads were provided, so no performance claim is made for these paths. |
| Plots/reports/config/CLI | Inspect their shared readers and cache consumers; exercise them in the offline suite and demo. No separate speedup claim is made for these paths. |
| CI/packaging | Preserve test coverage and delivery checks. No workflow or dependency constraints were changed. |

## Reproducing the measurements

Use a supported development environment from [CONTRIBUTING](../CONTRIBUTING.md).
The benchmark creates and removes only its own temporary synthetic data. It
never uses a user config, genome, database or existing output directory.

```bash
python scripts/benchmark_sequence_hotspots.py --case attributes --size 2000 --repeats 5
python scripts/benchmark_sequence_hotspots.py --case gap-normalization --size 2000 --repeats 5
python scripts/benchmark_sequence_hotspots.py --case conversion --size 2000 --repeats 5
python scripts/benchmark_sequence_hotspots.py --case pack --size 2000 --repeats 5
```

Each case runs in a separate process, performs one untimed warmup, then reports
all wall-time and process-CPU samples, their medians, process peak RSS, and an
output SHA-256. Input generation and output hashing are outside the timed region.
RSS includes imports, input setup and the warmup; it is a whole-process high-water
mark, not an allocation measurement of the optimized function. The RSS metric
uses the Unix `resource` module (macOS/Linux).

The input is one 640,000-base contig with 2,000 complete 300-base coding models,
8,000 GFF rows and 2,000 N runs. Normalization expands 20 Ns to 100.
The attribute case parses 100,000 strings per iteration, including escaped
commas and UTF-8. Conversion includes parsing, indexing, rendering and writing.
The pack case uses a new output root each time, enables gap normalization with
`gapjust_min = 0`, disables external validation and functional annotation, and
compares the final annotation/FASTA bytes. It does not time BUSCO, plots or HTML.
The pack pipeline removes the terminal N run before normalizing its remaining
1,999 intergenic gaps.

For the original implementation, export the fixed baseline and copy the same
benchmark driver into that export:

```bash
set -e
set -o pipefail
baseline_dir=$(mktemp -d)
git archive 00f63a6895b47fdbb647f40497d742168b5540a1 | tar -x -C "$baseline_dir"
cp scripts/benchmark_sequence_hotspots.py "$baseline_dir/scripts/"
python "$baseline_dir/scripts/benchmark_sequence_hotspots.py" --case pack --size 2000 --repeats 5
```

Run each baseline/changed pair consecutively in the same environment. Compare
`output_sha256` before interpreting timing differences. Increase `--size` to
investigate scaling; the original implementation can become slow on many gaps.

## Results: 2026-10-15

Baseline: `00f63a6895b47fdbb647f40497d742168b5540a1` (0.8.12). Changed:
the sequence optimizations in 0.8.13. Environment: Python 3.14.7, Biopython 1.88,
ReportLab 5.0.1, macOS 27.0.1 arm64. Each row uses five measured runs after one
warmup, at the default size above. Wall and CPU times are seconds.

| Case | Median process CPU before → after | CPU reduction | Median wall before → after | Process peak RSS MiB before → after |
| --- | --- | --- | --- | --- |
| Attributes | 0.494 → 0.332 | 32.8% | 1.541 → 1.317 | 201.9 → 198.7 |
| Gap normalization | 1.521 → 0.049 | 96.8% | 5.651 → 0.145 | 50.9 → 50.5 |
| MSS conversion | 0.470 → 0.096 | 79.5% | 1.065 → 0.388 | 59.0 → 60.5 |
| Pack | 3.728 → 1.178 | 68.4% | 25.919 → 3.390 | 60.6 → 71.4 |

The shared machine had substantial scheduling contention. These wall medians
are observations, not reliable speedup ratios; use process CPU for the numerical
improvement claims. Whole-process RSS varied between batches: an earlier pack
pair was 73.6 → 69.9 MiB, whereas the final pair above was 60.6 → 71.4 MiB.
The indexes add O(G) storage; these measurements do not establish a consistent
memory improvement or an index-caused 10.8 MiB increase.

Process-CPU samples, rounded to milliseconds:

| Case | Before | After |
| --- | --- | --- |
| Attributes | 0.452, 0.494, 0.538, 0.546, 0.491 | 0.265, 0.314, 0.350, 0.332, 0.343 |
| Gap normalization | 1.537, 1.492, 1.524, 1.414, 1.521 | 0.042, 0.051, 0.043, 0.049, 0.049 |
| MSS conversion | 0.552, 0.470, 0.455, 0.465, 0.502 | 0.098, 0.094, 0.096, 0.091, 0.106 |
| Pack | 3.728, 3.674, 3.742, 3.821, 3.607 | 1.194, 1.143, 1.153, 1.274, 1.178 |

All four before/after SHA-256 values matched:

```text
attributes         88c192a66d639222c4e37934f385e0311a52776a5f9bdbd830bb0e554ae9f98f
gap-normalization  5b40778b6b546b731c063cab9aef5560c4fb27863512af50862d74c658d6b2a6
conversion         af7c13b0322d6e78776ab3465e618f741cc0432933cec316970c741566908c1e
pack               167f2b79edd1b0f0cd16f25bd9c597a50f5b093cab0ecb4e50f19611789533e3
```

## Validation and limits

Regressions cover feature boundaries before/inside/after expanding and shrinking
gaps, zero-length deletions, unsorted features on multiple contigs, and gap
splitting/flags on both strands. Existing fixture outputs remain unchanged.
The benchmark hashes also include decoded Parent-list boundaries in the
attribute case and exact normalized/MSS/final FASTA bytes in the other cases.

These are synthetic core-pipeline measurements, not a representative real-genome
or external-tool benchmark. Real DDBJ validation, BUSCO and the unbundled
Triphyophyllum regression require their own tools/data. Wall-time measurements
on a shared machine may vary substantially; retain the samples and process-CPU
measurements when comparing results. No reduction in peak memory is claimed.
