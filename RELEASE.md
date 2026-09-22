# Push and release process

## Every GitHub push

1. Update [CHANGELOG.md](CHANGELOG.md) and bump `__version__` in
   [src/msspack/__init__.py](src/msspack/__init__.py), including documentation-only
   pushes. Use a patch increment for compatible maintenance; local-only commits
   need no bump.
2. Run [delivery checks](CONTRIBUTING.md#delivery-checks-before-push) after the bump.
3. Review the diff and outgoing commits, then follow the repository Git workflow
   in [AGENTS.md](AGENTS.md). A push alone does not request a tag or release.

## Tagged release (only when requested)

In addition to the push checks:

1. Verify a fresh wheel install and the unpacked sdist test suite in clean
   Python 3.11+ environments. Run outside the checkout without `PYTHONPATH=src`
   when testing the installed wheel so imports cannot silently use checkout code.
2. Review the DDBJ validation-tool agreement, then run at least one real MSS
   regression with validation as described in CONTRIBUTING. Record unavailable
   tools/data as a release validation gap.
3. Tag and publish the requested release.
