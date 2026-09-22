# Contributing

Run the commands below from the repository root in an activated development
virtual environment. `pyproject.toml` defines dependencies, Ruff, and mypy;
[CI](.github/workflows/ci.yml) tests Python 3.11–3.14, minimum dependencies, and
an installed wheel on macOS. No separate environment lock or test framework is needed.

## Development setup

For a new environment (Python 3.11 or newer):

```bash
python3.11 -m venv .venv
source .venv/bin/activate
python --version
python -c 'import sys; assert sys.version_info >= (3, 11), "Python 3.11+ required"'
python -m pip install -e ".[dev]"
python -m pip check
```

Replace `python3.11` with an available supported interpreter, e.g. `python3.14`.
Check an existing environment's Python version before installing or testing: an
old `.venv` can activate successfully but fail during imports. Preserve it if it
is still needed; create and activate a fresh environment elsewhere:

```bash
dev_env=$(mktemp -d)/venv
python3.14 -m venv "$dev_env"
source "$dev_env/bin/activate"
```

Then run the version and pip commands above. Installation needs package-index access; the core checks below do not need Java,
BUSCO, DIAMOND, external databases, or real genomes.

## Quick core checks

The existing unittest fixtures already provide a smoke test. Run both:

```bash
PYTHONPATH=src python -m unittest discover -s tests -p 'test_integration_pack.py' -v
PYTHONPATH=src python -m unittest discover -s tests -p 'test_demo.py' -v
```

Expect a nonzero test count and `OK`, with no skips. The minimal fixture checks
exact MSS/FASTA output and cache reuse/invalidation; the demo checks biological
metrics, plots, and report generation. They copy data into temporary directories
and clean up automatically. Do not run pack directly on the tracked fixture config.

For an inspectable CLI run, use a temporary copy (keep its printed path for review):

```bash
smoke_dir=$(mktemp -d)
msspack demo --output "$smoke_dir/demo" &&
msspack run --config "$smoke_dir/demo/config.toml" --no-busco --no-validate
```

Success prints output paths and produces final MSS files, a completed build
manifest, plots, and an HTML report under the temporary demo. This does not prove
DDBJ acceptance; fictional demo accessions must never be submitted.

## Checks by change

Run quick core checks for pipeline/output changes, plus relevant tests below.
Each name denotes `tests/test_<name>.py`; for example:

```bash
PYTHONPATH=src python -m unittest discover -s tests -p 'test_config.py' -v
```

| Change | Additional focused tests |
| --- | --- |
| CLI, config, templates | `cli`, `config`, `doctor`; CLI init checks packaged/example template equality |
| GFF, sequence, CDS, transcript selection | `gff`, `gff_adjustments`, `gap_normalization`, `padding_tools`, `select_one_mrna`, `audit_regressions` |
| MSS rendering or annotation-only updates | `mss_converter`, `mss_postprocess`, `header`, `submission_update` |
| Stage/cache/publication behavior | `pipeline`, `pipeline_actions`, `workflow`, `audit_regressions` |
| Functional annotation or databases | `functional_annotation`, `annotation_taxonomy`, `annotation_consistency`, `annotation_table`, `product_names`, `databases`, `database_lock` |
| Plots/report | `pipeline_plots`, `coordinate_duplicate_plots`, `report`, `demo` |
| External tool wrappers | `tools`, `validation`, `busco`, `doctor`; real-tool checks below when needed |
| Helper scripts | `scripts` |
| Documentation/Skill only | Execute changed commands where feasible; check links and Skill frontmatter |

Select across rows when code has several consumers. See
[the test-suite review](docs/test-suite-review.md) for the defects covered, and
[execution integrity](docs/execution-integrity.md) before changing publication or
cache behavior. Keep new functionality in library modules before CLI wiring.

## Full local checks

```bash
python -m compileall -q src tests
python -m ruff check .
python -m mypy src
MSSPACK_RUN_DDBJ_EXTERNAL=0 MSSPACK_RUN_BUSCO_EXTERNAL=0 PYTHONPATH=src python -m unittest discover -s tests -v
```

Expect successful exit codes; full discovery normally skips the three opt-in
external-tool tests. The explicit flags prevent inherited shell settings from
starting downloads. Any other skip needs explanation. Lint/type checks use the
active Python environment, avoiding tools from a different environment on PATH.

## Delivery checks (before push)

Run the full local checks above, then the existing audit and packaging checks:

```bash
python -m pip_audit .
python -m build
python scripts/check_distribution.py
python -m check_wheel_contents dist/*.whl
python -m twine check dist/*
```

The audit queries vulnerability services; build isolation can download build
requirements. All commands must exit zero; distribution verification prints each
verified archive. Use fresh packaging outputs so stale wheels are not checked.
For an existing checkout with valuable build artifacts, perform packaging in a
scratch copy of the current sources instead of deleting them. Build writes
`build/`, `dist/`, and package metadata; these are not source changes.

Version/changelog and clean installation checks are in [RELEASE.md](RELEASE.md).
A failed or unavailable check is a reported gap, never a reason to weaken CI.

## External and real-data regressions

These are separate from the offline suite. Confirm the tools/data and scope first;
do not automatically install tools or download databases for ordinary edits.
The dedicated CI job runs on schedules, manual dispatch, and version tags.

With Java available and the DDBJ tool agreement reviewed, this test downloads the
reviewed tools into temporary directories and runs Parser/transChecker:

```bash
MSSPACK_RUN_DDBJ_EXTERNAL=1 PYTHONPATH=src python -m unittest discover -s tests -p 'test_external_tools.py' -k ExternalDdbjTests -v
```

For BUSCO, follow the environment and lineage used by the external-tools job in
[CI](.github/workflows/ci.yml). `MSSPACK_RUN_BUSCO_EXTERNAL=1` enables its test;
`MSSPACK_BUSCO_LINEAGE`, `MSSPACK_BUSCO_DOWNLOAD_PATH`, and
`MSSPACK_BUSCO_OFFLINE=1` select a prepared local lineage/cache. Offline mode
requires the data to exist; do not interpret missing data as a passing regression.

For main-pipeline changes, also run the established real `Triphyophyllum` pack +
validation regression when its local config/data and tools are available and the
workload is authorized. It is not bundled. Use a copied config with a separate
output directory, preserve biological parameters, and report an unavailable run.
Inspect `build-manifest.json` for stage reuse and validation results.

## Benchmarking and cleanup

Use [benchmark-performance](.agents/skills/benchmark-performance/SKILL.md) for
performance work and `python scripts/benchmark_pack.py --help` for the harness.
Benchmark cleanup flags remove configured analysis outputs; never use them on
valuable results. `python scripts/clean_artifacts.py --dry-run` previews repository
cleanup, including BUSCO downloads. Neither cleanup nor benchmarking is needed
for the normal test loop.
