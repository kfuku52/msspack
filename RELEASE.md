# Release process

1. Run local release checks:

```bash
python -m compileall -q src tests
ruff check .
mypy src
PYTHONPATH=src python -m unittest discover -s tests -v
rm -rf build dist src/msspack.egg-info
python -m build --wheel
```

2. Update [`CHANGELOG.md`](CHANGELOG.md).
3. Bump [`src/msspack/__init__.py`](src/msspack/__init__.py).
4. Build a fresh wheel and verify install in a clean virtual environment.
5. Run at least one real MSS regression with validation.
6. Tag and publish the release.
