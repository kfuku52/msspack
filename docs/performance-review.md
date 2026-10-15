# Processing performance review

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
| Pipeline/cache/publication | Keep implementation checks before and after each stage; in 0.8.14 reuse the verified second fingerprint when recording cache state instead of computing a third. Preserve dependency checks, publication and locking. |
| Transcript handling | CDS extraction already groups models by sequence. In 0.8.14 build transcript coverage once for exon/UTR queries and propagate descendant removals through parent edges. Preserve selection, strand, phase, repeated-ID and shared-parent policies. |
| External annotation/databases | DIAMOND, HMMER, CDD and BUSCO already have index reuse, query filtering and integrity checks. No large databases or external-tool workloads were provided, so no performance claim is made for these paths. |
| Plots/reports/config/CLI | In 0.8.14 the full workflow prepares plots once and renders the report using the loaded configuration. Standalone report commands still refresh their plots. |
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
measurements when comparing results. The 0.8.13 measurements above do not establish
a reduction in peak memory; assembly streaming is measured separately below.

## Additional optimizations in 0.8.14

The follow-up review measured ten further costs against
`7b0c904f143102fee6a937cefe6c32b0fcc59e0b` (0.8.13):

| Area | Change and constraint |
| --- | --- |
| Assembly memory | Stream normalized contigs. With no edits retain the original FASTA bytes; when an edit first occurs, reread and canonically wrap the preceding unchanged contigs once. This bounds sequence storage by contig size, while edit metadata and IDs remain proportional to the input. |
| Internal stop counting | Cache a regular expression per stop-codon set, count only in-frame nonterminal matches, and retain overlapping matches. Genetic codes, ambiguous bases, lowercase input and partial terminal codons keep their previous meaning. |
| GFF containment | Decode attributes once per active row within each validation call, including escaped parent-list boundaries. Do not retain mutable row metadata across calls. |
| Annotation table | Collect metadata while the existing validated transcript-model builder consumes one GFF stream, avoiding a second parse. |
| Exon/UTR rendering | Merge transcript coverage once and use binary search to test containment. Existing strand, gap and feature-retention policies remain in effect. |
| Name consistency | Bound normalization/token/modifier caches to one audit, with 16,384 entries per transformation. Reuse immutable transformations and tier specifications; calculate decisions against the current IDF weights and thresholds. |
| Padding-stage counts | Build padding models from the already parsed rows and reuse output counts for logs and metrics. Keep strict validation of non-CDS rows as well as coding models. |
| Workflow plots | Render HTML after the workflow's existing plot pass; retain plot preparation for a standalone report. |
| Stage fingerprints | Record the implementation identity verified after the action; keep before/after mutation detection. A later source change invalidates the entry instead of relabelling old outputs. |
| Descendant removals | Replace convergence scans with parent-edge propagation. Selection keeps its per-row repeated-ID rules; coordinate cleanup keeps the union of parents across repeated IDs. Surviving shared parents and unreachable cycles remain intact. |

Use the same driver against an exported baseline and the current checkout:

```bash
set -e
set -o pipefail
baseline_dir=$(mktemp -d)
git archive 7b0c904f143102fee6a937cefe6c32b0fcc59e0b | tar -x -C "$baseline_dir"
cp scripts/benchmark_followup.py "$baseline_dir/scripts/"
PYTHONHASHSEED=0 python "$baseline_dir/scripts/benchmark_followup.py" --case pack --size 1000 --exons 10 --repeats 5
PYTHONHASHSEED=0 python scripts/benchmark_followup.py --case pack --size 1000 --exons 10 --repeats 5
```

For each case below, replace the arguments on both commands with the listed case,
size and exon count. The driver generates its own data, warms up once, times five
runs, checks each repetition, and prints wall/CPU samples, RSS and SHA-256. The
memory case also measures Python allocation peaks with `tracemalloc`; its timed
runs therefore include tracing overhead. File hashing and input generation are
outside the timed region. `PYTHONHASHSEED=0` makes set iteration reproducible.

The model workloads alternate plus/minus strands and include UTRs. The pack case
uses fresh outputs for each run; the workflow case measures a cached demo run
including plots and HTML. It hashes the deterministic annotation, FASTA and SVG,
excluding HTML timestamps and temporary paths. Consistency hashes all seven
output tables; containment hashes diagnostics; padding/count/cache cases hash
their relevant results as well as any output files.

### Follow-up measurements: 2026-10-15

Environment: Python 3.14.7, Biopython 1.88, ReportLab 5.0.1, macOS 27.0.1
arm64. All 13 baseline/changed output SHA-256 values matched. Times are seconds.

| Workload and reproduction arguments | Median CPU before → after | Median wall before → after | Process peak RSS MiB before → after |
| --- | --- | --- | --- |
| Padding: `--case padding --size 10000` (3,000 bp each) | 1.484 → 0.500 | 2.737 → 0.551 | 296.4 → 293.8 |
| Containment: `--case containment --size 2000 --exons 10` | 0.414 → 0.178 | 0.520 → 0.303 | 73.7 → 85.6 |
| Table: `--case annotation --size 10000` (one exon each) | 0.384 → 0.187 | 0.466 → 0.212 | 64.7 → 65.2 |
| MSS: `--case exons --size 100 --exons 50` | 0.139 → 0.083 | 0.142 → 0.089 | 63.7 → 60.1 |
| MSS: `--case exons --size 100 --exons 200` | 1.078 → 0.431 | 1.281 → 0.613 | 131.1 → 135.6 |
| Padding stage: `--case counts --size 1000 --exons 10` | 0.428 → 0.269 | 0.846 → 0.444 | 64.8 → 69.5 |
| Names: `--case consistency --size 2000` (19,945 pairs) | 0.946 → 0.485 | 1.289 → 0.584 | 92.5 → 93.1 |
| Removal: `--case descendants --size 1000` (deep, child-first) | 0.142 → 0.004 | 0.175 → 0.005 | 29.2 → 29.2 |
| Cache: `--case cache --size 20` (new tiny stages) | 0.255 → 0.197 | 0.324 → 0.230 | 32.2 → 31.8 |
| Cached demo: `--case workflow` (initial batch) | 0.197 → 0.209 | 0.272 → 0.628 | 55.0 → 49.2 |
| FASTA: `--case memory --size 128` (128 × 1 MiB; no edits) | 1.077 → 1.120 | 1.910 → 1.520 | 203.9 → 42.2 |
| FASTA: `--case memory --size 128 --late-edit` | 3.231 → 3.349 | 4.214 → 5.331 | 213.1 → 46.7 |
| Pack: `--case pack --size 1000 --exons 10` | 2.410 → 2.022 | 2.865 → 2.363 | 111.8 → 108.4 |

The mixed-strand, many-exon pack uses 23,000 GFF rows and a 1.2-million-base
assembly. Its CPU median decreased 16.1%; this is an additional improvement over
the already optimized 0.8.13 baseline. The isolated CPU reductions were 66.3% for
padding, 57.0% for containment, 51.3% for table generation, 40.7–60.1% for
50/200-exon conversion, 37.1% for the padding stage including counts/models,
48.8% for name consistency, 96.9% for the deliberately deep hierarchy, and 22.8%
for tiny fresh cache stages. The deep hierarchy is a scaling case, not a typical
annotation structure. Shared-machine wall times remain noisy.

Median traced Python allocation peaks for 128 × 1 MiB contigs fell from
129.15 to 3.15 MiB with no edits (97.6% less), and from 131.14 to 6.42 MiB
when the last contig was edited (95.1% less). These are function allocation
measurements, distinct from whole-process RSS. Streaming does not bound memory
independently of the largest contig. Rereading an unchanged prefix trades extra
input I/O for lower retained memory; the traced CPU medians were about 4% higher
in these two memory cases. Do not infer a runtime improvement from them.

Containment now retains per-row metadata and parent indexes: RSS in this pair
was 73.7 → 85.6 MiB. The audit caches are bounded but also retain names/tokens
for the duration of an audit. These are time/space tradeoffs. The other small
RSS differences include imports, warmup, allocation high-water marks and normal
measurement variation.

The initial cached-demo pair did not establish a speedup (CPU 0.197 → 0.209 s).
Additional alternating batches were used to investigate that noisy small workload.
Four batches of five runs per version, alternating which version ran first,
gave pooled CPU medians of 0.205 → 0.182 s and wall medians of 0.274 → 0.245 s.
Batch CPU medians were 0.212 → 0.163, 0.191 → 0.183, 0.234 → 0.182,
and 0.179 → 0.201 s. All checksums matched. The pooled reduction is 11.2%, but
the reversed last batch and initial pair show that a small-run speedup remains
sensitive to machine load; the pack and memory results are stronger evidence.
The regression independently checks that the full workflow prepares plots once
and emits the HTML/SVG and completed manifest; standalone reports still refresh
plots.

<details>
<summary>CPU samples and output checksums for the initial pairs</summary>

| Case | Before CPU samples | After CPU samples |
| --- | --- | --- |
| padding | 1.477, 1.513, 1.458, 1.561, 1.484 | 0.501, 0.503, 0.500, 0.498, 0.498 |
| containment | 0.445, 0.443, 0.414, 0.413, 0.381 | 0.155, 0.169, 0.187, 0.191, 0.178 |
| annotation | 0.388, 0.444, 0.317, 0.302, 0.384 | 0.239, 0.186, 0.187, 0.182, 0.196 |
| exons50 | 0.139, 0.137, 0.136, 0.144, 0.141 | 0.083, 0.088, 0.102, 0.079, 0.082 |
| exons200 | 1.250, 1.147, 1.078, 1.076, 1.073 | 0.423, 0.441, 0.431, 0.381, 0.452 |
| counts | 0.356, 0.434, 0.428, 0.396, 0.464 | 0.269, 0.266, 0.289, 0.268, 0.270 |
| consistency | 0.951, 0.932, 0.957, 0.946, 0.914 | 0.425, 0.545, 0.450, 0.496, 0.485 |
| descendants | 0.136, 0.144, 0.135, 0.142, 0.145 | 0.005, 0.004, 0.004, 0.004, 0.004 |
| cache | 0.280, 0.247, 0.255, 0.254, 0.283 | 0.192, 0.161, 0.201, 0.199, 0.197 |
| workflow | 0.150, 0.207, 0.192, 0.197, 0.243 | 0.209, 0.186, 0.233, 0.175, 0.217 |
| memory | 1.057, 1.067, 1.088, 1.077, 1.142 | 1.037, 1.121, 1.097, 1.126, 1.120 |
| memory-late | 3.218, 3.157, 3.231, 3.236, 3.300 | 3.349, 3.352, 3.351, 3.330, 3.249 |
| pack | 2.352, 2.544, 2.451, 2.410, 2.373 | 2.022, 1.840, 1.969, 2.242, 2.138 |

```text
padding        55b561d72424f6329e149e60ab75adb01846c9cf0dd1ec1d4c2030724d1d8169
containment    2ab148ae3d0ea86cc098a884966539a57f7bc76429e807e570d2a66755309e97
annotation     b9aef80467723674255e7cd492f3fbde0c530285213a6f4abc36b125fa36bc77
exons50        64acf1b6bcc72f470a50422a07b284090e1db9c3289a921530d4a741be6f479d
exons200       e7b6fe7ee01850a486320dc1c48a52441401708102bdf62047988b836a2a215a
counts         69ba890444f08e1513f6ee278818f51c06466b1a38d291c1ffd3eac24c96de59
consistency    2d34dea771923b72ee6b4d9a017a922b06dcb96677fa9fe0bcf185bb4206aaa1
descendants    c72a34bb869911f2f6c3741c80f7191b7f6622f29d3297e5f8de10060797f068
cache          620bf13219fd8b664b6e7a1445983c1790c99682997e9c288f780f4fdd103278
workflow       c5de1de8029608cec59aeb2405acc811bbded5c8b0d76f4be6deb8dd4c8e3c4b
memory         9752b893bb7000dc05d911f729fda4f687068fb35c261289e1b50aa2856e8663
memory-late    f2bcd944afb291e7189795cdbfc97913392761acc7842cb41b93858c67c6d6a1
pack           88322cfd93529a2b9d863a692f1469b61a590819872baced5d74bed013bcfd4e
```

</details>

All regressions use offline fixtures or generated data. No expected biological
fixture outputs were regenerated. Real genomes, DDBJ validation and BUSCO remain
outside these measurements and require the tools/data described in CONTRIBUTING.
