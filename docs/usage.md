# Inputs, configuration, and outputs

Start with the [README quick start](../README.md#quick-start). This guide describes
local behavior; successful file generation is not confirmation of DDBJ acceptance.

## Inputs

`inputs.fasta` is a nucleotide FASTA and `inputs.gff` is a nine-column,
tab-separated GFF3. Both accept gzip files with a `.gz` suffix. A separate FASTA
is required even when the GFF3 contains a `##FASTA` section. FASTA identifiers are
the first whitespace-delimited token after `>`; they must be unique and have
nonempty sequences. Every GFF3 sequence ID must match a FASTA identifier exactly.

GFF3 columns are `seqid`, `source`, `type`, `start`, `end`, `score`, `strand`,
`phase`, and `attributes`. Coordinates are one-based and inclusive. Use GFF3 `.`
for unavailable fields, not CSV-style empty cells or `NA`; CDS phase should be
`0`, `1`, or `2` and coding strand `+` or `-`. `ID` and `Parent` attributes express
the feature hierarchy; a gene row is not required for a coding transcript.
Parent lists and percent-escaped IDs follow the rules in
[input integrity](execution-integrity.md). `doctor` checks columns, coordinates,
FASTA IDs, and sequence-ID correspondence, but is not a full biological validator.

The main pipeline selects one mRNA/transcript per parent gene; the ranking is
[described in the README](../README.md#gff3-feature-handling). It also trims trailing
Ns, optionally normalizes N runs, and adjusts gene models. Thus `pack` is
not a sequence-preserving annotation update. Use
[`prepare-update`](annotation-updates.md) for an existing assembly whose bases
and CDS locations must remain unchanged.

## Configuration and overrides

Config files are TOML. Unknown sections/keys, wrong value types, and non-finite
numbers are rejected. Missing optional keys use loader defaults; required keys
must be present. Empty strings are accepted only for settings that allow them,
not as a general missing-value convention. Quote dates and genetic codes as in
the starter file. TOML date literals are not accepted for string-valued dates.

Input, output, local database, and tool-cache paths in the config are resolved
relative to the **config file's directory**, not the shell's working directory.
`~` is expanded; `$VARIABLE` is not shell-expanded inside TOML path strings.
If `project.output_dir` is omitted, outputs go to `build/<project.name>/` beside
the config. CLI file paths such as `--config`, `--ann`, and `--output-dir` are
relative to the working directory.

Database-root precedence, highest first:

1. Nonempty `--db-dir` on `run` or `db status`.
2. Nonempty `MSSPACK_DB_DIR` environment variable.
3. `databases.root`, defaulting to `msspack_db` beside the config.

Relative database-root overrides, including `--db-dir`, are also resolved from
the config directory. An explicit `busco.download_path` independently overrides
`<database-root>/busco`. There is no general environment-to-TOML override mapping.

`run --no-busco` skips BUSCO, and `run`/`pack --no-validate` skip configured DDBJ
validation without editing the config. `run --no-report` skips only HTML; plots
are still generated. `pack` does not run BUSCO, plots, or the HTML report.
`busco` has its own lineage, thread, mode-selection, and auto-lineage overrides;
use `msspack busco --help` for their names and choices. Config validation occurs
before these overrides; see the [unresolved BUSCO issue](documentation-audit.md#unresolved).

The starter config is not an exhaustive schema or a list of omission defaults.
In particular, its BUSCO `threads = 8`, `sample.linkage_evidence = "proximity ligation"`,
product replacement patterns, hold date, and metadata are explicit example values.
Omitting them yields 1 BUSCO thread, `"paired-ends"`, no replacement patterns,
no hold date, or a required-field error as appropriate.

Useful settings omitted from the starter file include:

| Key | Omission default | Meaning |
| --- | --- | --- |
| `pipeline.gapjust_gap_len` | `100` | Normalized internal gap length, in bases |
| `pipeline.gapjust_min`, `pipeline.gapjust_max` | `80`, `120` | Lower bound for expanding shorter gaps / upper bound for shrinking longer gaps, in bases |
| `pipeline.replace_product_with` | `"hypothetical protein"` | Replacement for configured product patterns |
| `pipeline.validate_with_parser`, `pipeline.validate_with_transchecker` | `true`, `true` | Enable validation in `pack`/`run` |
| `tools.java` | `"java"` | Java executable |
| `tools.java_heap` | `"16G"` | Java maximum heap setting |
| `tools.cache_dir` | platform cache | DDBJ tool cache location |
| `tools.gff3sort` | unset | Legacy accepted setting, ignored; sorting always uses the internal sorter |

The accepted keys/types and omission defaults are defined in
[config_loading.py](../src/msspack/config_loading.py); allowed values and ranges
are checked in [config_validation.py](../src/msspack/config_validation.py).
Annotation identity and query/subject coverage thresholds use percentages
(0–100); `pfam_min_domain_coverage`, `near_top_bitscore_ratio`, and
`min_token_score` use fractions/scores in 0–1. These are filtering settings,
not estimated probabilities of a correct biological assignment.

## External tools and validation

The pip install supplies Python dependencies, not Java, BUSCO, DIAMOND, HMMER,
RPS-BLAST, or DDBJ tools. The core demo needs none of these external tools.
Default `run` settings enable CDS BUSCO and DDBJ validation; use the quick-start
skip flags when these tools/data have not been prepared. Functional annotation
is off by default; enabled remote database sources and taxonomy lookup can
require network access. The bundled functional demo uses only a local reference
but still requires DIAMOND.

DDBJ validation requires Java, bash, and separately installed Parser/transChecker.
`tools install` and `tools list` accept `--cache-dir`, but do not read a project
config. For a custom cache, pass the same absolute directory to `tools install`
and set it as `tools.cache_dir` in the TOML. The default is
`~/Library/Caches/msspack` on macOS, or `$XDG_CACHE_HOME/msspack` on Linux when
XDG_CACHE_HOME is absolute (otherwise `~/.cache/msspack`).

`validate --ann FILE --fasta FILE` always requests **both** DDBJ tools, even if
`--config` disables pipeline validation. Its optional config supplies Java/cache
settings; without it, platform cache and `java` with a `16G` heap are used.
Standalone validation writes `logs/` and `validation/` beside the resolved
annotation file, not under `project.output_dir`. It does not upload a submission.

`doctor` is a diagnostic command. In particular its BUSCO check is currently
optional even when `run` would require BUSCO; see the
[audit finding](documentation-audit.md#unresolved). Check the reported details,
not only its exit status.

## Outputs and interpretation

Paths below are relative to `project.output_dir`:

| Path | Contents |
| --- | --- |
| `final/<sample.locus_tag>.ann.txt` | Headerless, five-column MSS annotation: entry, feature, location, qualifier, value; blank cells continue the current entry/feature |
| `final/<sample.locus_tag>.fasta` | MSS nucleotide sequences with `//` after each entry, not plain FASTA for arbitrary downstream tools |
| `final/ddbj-validation-summary.json` | Structured validation status, including a skipped result when validation was disabled |
| `intermediate/` | Numbered prepared FASTA/GFF and annotation tables (`ID`, `Description`, `Locus_tag` for the core annotation tables) |
| `logs/` | Per-stage logs, metrics JSON, and ID/duplicate audit tables |
| `build-manifest.json` | Latest build status, stage execution/cache decisions, timings, and optional orchestration/plot/report metadata |
| `plots/pipeline-flow-summary.tsv` | Columns `metric`, `value`, `unit`, `description` |
| `plots/pipeline-event-counts.tsv` | Columns `metric`, `label`, `count`, `unit` |
| `plots/pipeline-gene-flow.tsv` | Columns `source`, `target`, `count`, `source_label`, `target_label`, `source_stage`, `target_stage` |
| `plots/pipeline-gene-flow.sankey.svg` and `.pdf` | Gene-flow diagram |
| `plots/coordinate-duplicate-gene-models.tsv` | Kept/removed model feature rows for all duplicate removals; figure limits do not truncate this table |
| `report/index.html` | HTML report linked to the other outputs; keep the output directory together when moving it |
| `busco/` | Optional CDS/genome comparison summaries and figures |

Event counts mix genes and transcripts as indicated by `unit`. Events at different
stages can affect the same gene; do not sum them to obtain unique affected genes.
Gene-flow `count` values describe transitions between stages, not independent
biological replicates. A missing optional BUSCO or annotation result is not a zero
score. A skipped validation status is not PASS. See the README's
[validation description](../README.md#example-outputs) for the tool count checks.

For cache inspection, each manifest stage has `ran` and `duration_seconds`;
`ran = false` denotes reuse. The top-level `status` describes the pack operation;
`run.status`, when present, describes the enclosing `run` operation. Later
standalone commands can retain earlier orchestration metadata, so also check
timestamps. Final-file existence alone cannot establish that the latest run passed.

Unchanged reruns reuse stage caches; reruns may publish another complete final
generation. `--force-compute` invalidates analysis caches, retains databases and
published generations, and reruns enabled work. It is not an output cleanup command.
`init` and `demo` refuse to replace their files unless `--force` is passed; demo
force replaces bundled files rather than clearing unrelated files or old results.
`prepare-update` requires a new output directory. See
[execution integrity](execution-integrity.md) for atomic publication, concurrent
writers, failure recovery, and cache invalidation.
