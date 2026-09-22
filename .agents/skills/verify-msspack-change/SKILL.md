---
name: verify-msspack-change
description: Select and run msspack regressions for a code, config, CLI, or contributor-command change. Use existing offline fixtures first; not for benchmarking or publishing releases.
---

# Verify an msspack change

## Inputs and scope

Use the requested change, current diff (including untracked additions), active
Python environment, and any explicitly supplied real-data/tool constraints.
Read [CONTRIBUTING.md](../../../CONTRIBUTING.md) as the command reference; do not
maintain another check list here. Documentation-only changes need their changed
commands checked, not an invented biological regression.

## Procedure

1. Check the interpreter and dependencies using Development setup. An activated
   environment is not proof of Python 3.11+ compatibility. Preserve an old
   environment and use a fresh supported one when necessary.
2. Map changed behavior and its consumers to Checks by change. For pipeline or
   output changes, run both Quick core checks before expanding to focused modules.
   Use unittest discovery from the repository root; confirm tests were collected.
3. Interpret the existing assertions: integration checks exact annotation/FASTA
   output, unchanged-run cache reuse, and content-based invalidation; demo checks
   expected biological events and plot/report artifacts. If a new behavior is not
   covered, add its smallest meaningful regression in the existing unittest suite.
   Do not rewrite expected outputs without explaining the intended semantic change.
4. When checking CLI/documented workflow changes, run the temporary demo command
   from CONTRIBUTING. Inspect its manifest status and printed artifact paths.
   Missing BUSCO resources are expected with the core demo's BUSCO disabled;
   a successful no-validation run does not establish DDBJ acceptance.
5. Run the applicable full/delivery checks for the task phase. Stop on a failed
   command and diagnose it; automated command sequences must stop on failure
   (for example, `set -e`), not return the last successful command's status.

## Results and limits

Return commands, interpreter version, test counts and outcomes, and any relevant
artifact/manifest observations in the task report. No permanent report file is
required. Separate executed checks, static inspection, expected skips, and checks
blocked by tools/data/network. Temporary fixture runs must not modify real data.

If setup or a check fails, retain useful diagnostics and report the actual blocker.
Do not patch runtime behavior to accommodate an unsupported Python environment,
weaken assertions, or count zero tests as success. Use the explicit external-test
controls in CONTRIBUTING only when the required downloads/workload are in scope;
otherwise state that real DDBJ/BUSCO or genome validation remains unverified.
