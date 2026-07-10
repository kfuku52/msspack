# Release process

1. Run local release checks:

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

2. Update [`CHANGELOG.md`](CHANGELOG.md).
3. Bump [`src/msspack/__init__.py`](src/msspack/__init__.py).
4. Build a fresh wheel and verify both the wheel install and the unpacked sdist test suite in clean environments.
5. Review the DDBJ validation-tool agreement, then run at least one real MSS regression with validation.
6. Tag and publish the release.
