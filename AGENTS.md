<!-- BEGIN KF AGENT POLICY: source=https://github.com/kfuku52/kf-agent-policy; version=12; sha256=d38271eaba87f3f2c5c92b6c4ea821eca662d38f6f0db388fda5270edb5e9afa -->
# Common agent policy

Repository-specific instructions override these defaults.

- Follow the user's task scope within higher-priority instructions and execution
  permissions. Complete implementation through affected verification and a result
  report; a plan or investigation ends with its requested deliverable. Continue
  authorized work without repeated approval; identify actual blocking boundaries.
- Inspect the worktree and preserve unrelated changes. Refresh remote information
  when needed; do not merge, rebase, or switch branches merely to inspect it.
- Prefer the default branch when starting work without an established branch.
  Preserve an existing task branch; follow explicit user branch instructions.
  Never create or switch branches solely for a commit, push, release, or PR.
- Change or recommend branch protection only when explicitly asked. Honor explicit
  repository-specific direct-push exceptions; otherwise report a rejected push
  without bypassing protection or inventing a branch or PR.
- Prefer the simplest design that meets established requirements. Avoid
  complexity for speculative needs; extract helpers when they improve clarity
  or enable reuse.
- Unpublished implementation details may be redesigned; preserve existing public
  APIs, file formats, and saved-data compatibility unless a breaking change is
  authorized. Update affected producers, consumers, tests, examples, and docs.
- Validate untrusted input at the responsible boundary; rely on established
  invariants internally. Revalidate when mutation or a new trust boundary
  invalidates those guarantees.
- Fix verified root causes; do not hide failures with fallbacks or weaker checks.
  Document unavoidable workarounds and their removal conditions.
- Read relevant docs and run the repository's check entrypoint for the change and
  phase. Verify affected behavior; report checks run and omitted. Repeat or broaden
  successful checks only for new changes, failures, or unresolved concerns.
- For library metadata, require demonstrated incompatibility for exact pins or
  upper bounds; keep reproducibility locks separate.
- When editing READMEs, keep them concise with useful visuals inline; put extended
  guides in linked documentation.
- For GitHub push/release work, use `prepare-github-push` in `.agents/skills/`.
  Local-only commits need no version bump; GitHub pushes require one.
- Optimize bottlenecks in RAM usage or CPU time whenever found. For software
  performance work, use `benchmark-performance` in `.agents/skills/`.
  Performance claims require comparable measurements and equivalent output.
- For GitHub Actions edits, use `optimize-github-actions` in `.agents/skills/`.
  Preserve required coverage; never run untrusted PR code on self-hosted runners.
<!-- END KF AGENT POLICY -->

# Repository Instructions

## Git workflow

- Work directly on `main` by default.
- Do not create or switch to a feature branch unless the user explicitly asks for one.
- When the user asks to commit or push changes, commit and push them directly to `main`.
- After an explicitly requested branch is merged, delete it from both the local repository
  and the remote.

## Start here

- Read [CONTRIBUTING.md](CONTRIBUTING.md) for setup, checks by change, and delivery
  checks. Read the relevant README workflow/configuration section for user behavior.
- Enter through `src/msspack/cli.py` (commands), `workflow.py` (`run`), and
  `pipeline.py` / `pipeline_actions.py` (stage graph/actions). Config definitions,
  loading, and validation live in `config_models.py`, `config_loading.py`, and
  `config_validation.py`; MSS rendering lives in `mss_converter/`.
- Use [verify-msspack-change](.agents/skills/verify-msspack-change/SKILL.md) to
  select and interpret existing regressions. Do not require real genomes for the
  initial feedback loop.

## Contracts and data

- Preserve configured genetic codes, CDS phase/strand handling, transcript and
  duplicate-selection policies, annotation thresholds, and database/lineage choices
  unless the task explicitly changes them. Defaults are documented in
  `examples/msspack.example.toml`; its packaged template must stay synchronized.
- Preserve CLI/config compatibility and MSS annotation/FASTA semantics. Do not
  regenerate expected fixture outputs just to make tests pass. For publication,
  cache identity, and shared database changes, read
  [execution integrity](docs/execution-integrity.md). For annotation-only updates,
  read [the update contract](docs/annotation-updates.md): bases, CDS locations, and
  existing protein IDs must survive accession renaming.
- Do not edit user configs, real input datasets, shared databases, or generated
  `build/`, `dist/`, caches, submission generations, and run outputs as source.
  Species examples are sanitized; keep placeholders and fictional demo metadata.
  Copy fixtures/demo data to temporary directories for manual runs. The broad
  cleanup script also removes BUSCO downloads; it is not a prerequisite for tests.

## Finish

Run the applicable checks in CONTRIBUTING; before push run its delivery checks
and follow RELEASE for version/changelog policy. Report changes, commands and
Python version, passed/failed/skipped checks, and unmet external-data/tool needs.
Distinguish local checks from CI and real DDBJ validation. Do not claim skipped
external tests passed.
