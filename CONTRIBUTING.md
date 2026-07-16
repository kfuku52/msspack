# Contributing

## Development setup

```bash
python3.11 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

Use Python 3.11 or newer. Recreate environments made for older `msspack` releases rather
than reusing stale package and dependency metadata.

## Local checks

```bash
python -m compileall -q src tests
ruff check .
mypy src
pip-audit .
PYTHONPATH=src python -m unittest discover -s tests -v
rm -rf build dist src/msspack.egg-info
python -m build
python scripts/check_distribution.py
check-wheel-contents dist/*.whl
twine check dist/*
```

## Regression checks

- Run the minimal fixture integration test in `tests/test_integration_pack.py`.
- When touching the main pipeline, run a real `Triphyophyllum` pack + validation regression with a local config.
- Inspect `build-manifest.json` for stage reuse and validation outputs when debugging cache behavior.

## Benchmarking

Use the benchmark harness to compare fresh reruns:

```bash
python scripts/benchmark_pack.py --config /path/to/config.toml --repeats 3 --no-validate
```

## Code structure

- `src/msspack/pipeline.py`: orchestration and stage graph
- `src/msspack/pipeline_actions.py`: stage action wrappers
- `src/msspack/mss_converter/`: MSS rendering and conversion logic
- `src/msspack/config_*.py`: config schema, loading, and validation
- `src/msspack/ddbj_tools.py`: DDBJ tool download/install/metadata

Keep new functionality in library modules first, then expose it through CLI wiring only where needed.
