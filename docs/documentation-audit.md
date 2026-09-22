# Documentation/implementation audit — 0.8.4

The audit started from `main` at `ccd82f64ba7ac673d1a911eae7ee792251c18f1a`
(0.8.3), with no uncommitted changes. It prioritized installation, the bundled
core workflow, public CLI/config semantics, output interpretation, and reruns.
Runtime behavior and dependencies were not changed; the version increment follows
the repository's push policy. Historical changelog entries and example-run figures
were not rewritten to represent a new run.

## Corrected (A: documentation error or omission)

| Location / previous problem | Evidence | Correction |
| --- | --- | --- |
| README installation: interpreter/environment and Git prerequisites were implicit | `pyproject.toml` requires Python >=3.11, Biopython and ReportLab; console entry point is `msspack.cli:main` | Use `python -m pip`; identify Git and separately installed external tools |
| README quick start: `cd msspack-demo` made later `my_submission.toml` commands look in the wrong directory | `init` and `demo` CLI path handling; shell reproduction | Run the demo in a subshell and identify the actual output root |
| README configuration: starter described as the complete schema/defaults | `config_loading.py`, `config_models.py`, `test_config.py`; omitted tools/gap/validation keys and BUSCO default of 1 versus example 8 | Correct the claim, annotate example values in both synchronized templates, link to the usage guide |
| README database command lacked required `--config`; override/path rules absent | `cli.py`, `MSSPackConfig.database_dir`, `database_directory_override`, `test_cli.py` | Complete the command; describe CLI > environment > TOML precedence and config-relative resolution |
| README GFF feature discussion did not distinguish converter isoform support from pipeline selection | `_prepare_gff`, `_selected_transcript_ids`, `test_select_one_mrna.py`, demo's removed-mRNA metric | Describe the main pipeline's one-transcript selection and ranking |
| Main output names, TSV columns/units, standalone validation destinations, and reuse interpretation absent | `pipeline._build_outputs`, `fasta_steps.write_mss_fasta`, `pipeline_plot_render`, `build_manifest`, `ValidationOptions`/`ValidationArtifacts`; demo output and integration tests | Add the usage guide, including MSS `//` terminators, event overlap, `ran`, and skipped validation |
| Legacy `tools.gff3sort` behavior and custom tool cache setup undocumented | `pipeline_actions.sort_gff`, `utils.default_cache_dir`, `cli._handle_tools`, `test_pipeline_actions.py` | Document that the accepted sorter setting is ignored and how installation/config cache paths must agree |

## Findings resolved in 0.8.5

The reproductions below record 0.8.4 behavior. Both B findings were fixed in
0.8.5: project configs now accept both BUSCO modes disabled; explicit `busco`
requires a mode after CLI overrides; doctor requires the BUSCO executable and
database directory whenever a configured run enables BUSCO. `run --no-busco`
still skips these requirements. Real config/CLI regressions now cover both cases,
including an unusable BUSCO database path and explicit `busco --cds` re-enabling.
The 0.8.5 follow-up passed 233 tests (three opt-in external tests skipped),
compileall, Ruff, mypy, and the dependency audit on Python 3.14.7. A temporary
demo with both modes disabled completed doctor/run and produced nine MSS entries
and the HTML report without BUSCO. Real BUSCO execution remains unverified.

### B: BUSCO-disable configuration rejected before the orchestrator can skip it

README says config settings and `--no-*` options determine optional work.
`workflow.run_all` explicitly supports both BUSCO modes being false, and
`test_run_all_skips_busco_when_config_disables_both_modes` asserts that behavior
with a mocked loader. But `config_validation.validate_busco_config` rejects that
config before orchestration or CLI skip flags are applied.

Reproduction: copy the bundled demo into a temporary directory; append to its
`config.toml` (which has no BUSCO table):

```toml
[busco]
run_cds = false
run_genome = false
```

`msspack run --config config.toml --no-busco --no-validate` exits 1 with
`At least one of 'busco.run_cds' or 'busco.run_genome' must be true`.
No external tools are needed to reproduce it. The intended scope of config-level
validation needs resolution; the loader and orchestrator were left unchanged.
The documented demo remains usable with its original config and `--no-busco`.

### B: doctor success does not establish BUSCO readiness for run

README describes `doctor` as checking the completed config and required tools.
`doctor.run_doctor` always marks the BUSCO executable optional, whereas
`workflow._preflight` marks it required when BUSCO will run. Existing doctor tests
assert optional status, so this is a diagnostic-contract question as well as a
potential implementation defect.

Reproduction: append the following to a separate, original demo config:

```toml
[busco]
command = "msspack-audit-nonexistent-busco"
```

`msspack doctor --config config.toml` exits 0 and shows BUSCO as OPTIONAL;
`msspack run --config config.toml --no-validate` exits 1 at preflight with the same
executable MISSING. No BUSCO process or download starts. The usage guide flags the
limitation without redefining doctor success as proof of readiness. Aligning the
checks is deferred because this audit does not change runtime behavior.

No additional C-class contradiction was established in the reviewed scope.
Scientific validity of all naming thresholds, mappings, and representative
real-data figures was not independently re-established from code.

## Execution and limits

Local verification used macOS and a new temporary Python 3.14.7 virtual environment.
`python -m pip install -e ".[dev]"` and `python -m pip check` succeeded. No existing
environment, real dataset, shared database, or generated repository output was cleared.
The exact Git installation command also succeeded in a second isolated environment
(`python -m pip install git+https://github.com/kfuku52/msspack.git`), installing
the then-published 0.8.3; `pip check` passed. `python -m pip_audit .` on the
modified checkout found no known vulnerabilities (it reported cache-deserialization
warnings, then completed successfully).

- Executed the README `init`, `demo`, subshell `run --no-busco --no-validate`
  commands verbatim inside a temporary working directory. The original directory
  and starter file remained available afterward.
- Executed the bundled demo README's core `pack --no-validate`, `plot`, and `report`
  commands verbatim on a temporary demo copy. Inspected nine MSS sequence entries,
  the annotation pair, `not_run` validation summary, completed manifest, HTML,
  SVG/PDF outputs, and TSV headers. Rerunning reused all 19 pack stages.
- Used the demo config as a small-input substitute for generic project examples:
  `db status`, `run --no-busco --no-validate --no-report`, overwrite refusals,
  and the two unresolved BUSCO cases. Confirmed plots still exist with no report
  requested and verified database override precedence with the configuration API.
- Ran `python -m compileall -q src tests`, `python -m ruff check .`,
  `python -m mypy src`, and
  `MSSPACK_RUN_DDBJ_EXTERNAL=0 MSSPACK_RUN_BUSCO_EXTERNAL=0 PYTHONPATH=src python -m unittest discover -s tests -v`.
  All succeeded: 231 tests, three opt-in external tests skipped, 62 modules type-checked.
  Existing tests cover CLI/template equality, configuration, update preservation,
  exact fixture output, cache invalidation, and demo events; no redundant test was added.

Delivery checks also passed in a scratch source copy: `python -m build`,
`python scripts/check_distribution.py`, `python -m check_wheel_contents dist/*.whl`,
and `python -m twine check dist/*`. A third clean environment installed the 0.8.4
wheel with `python -m pip install dist/*.whl`; outside the checkout, `pip check`,
`msspack --version`, and the demo `run --no-busco --no-validate` succeeded.
Imports were confirmed to come from that environment, and its final pair, completed
manifest, and HTML were checked. All four example TOMLs loaded, the packaged and
repository templates matched, and 14 CLI help entry points exited successfully.
The local Markdown check found all 40 referenced targets/anchors.

External DDBJ installation/validation, BUSCO/lineage downloads, DIAMOND functional
annotation, conda solving, and real-genome runs were not executed. The default
full workflow with external tools and DDBJ acceptance remain unverified. The
annotation-only workflow was covered by synthetic existing tests, not by a public
accession mapping or live submission. Native Windows, Python versions other than
3.14, every internal CLI/docstring, full annotation-evidence schemas, and remote
link availability were outside the execution scope. No documentation build system
is configured; local Markdown targets and anchors were checked separately.
