# Changelog

## Unreleased

## 0.7.4 - 2026-08-06

- Record structured DDBJ Parser and transChecker results, including tool versions,
  diagnostics, translated FASTA record counts, and failure state, and show the
  results in the pipeline Sankey and HTML report.
- Use the newest validation summary from either a full pack run or explicit
  `msspack validate`, including its current logs and translated FASTA links.
- Require transChecker AA and nucleotide record counts to match the final annotation's
  CDS count, preserve skipped-tool state after sequential failures, and reject
  contradictory validation summary JSON.
- Embed subsetted Bitstream Vera fonts in generated PDFs so chart text renders
  consistently across Poppler, Ghostscript, and platform PDF viewers.
- Materialize whitespace-free CDD database aliases with file links or copies so
  RPS-BLAST can memory-map databases stored in cloud-synchronized paths.
- Apply Python's safe tar extraction filter when installing DDBJ tools on Python
  3.12+, preserving Python 3.11 support and avoiding the Python 3.14 behavior change.

## 0.7.3 - 2026-07-31

- Map functional-annotation evidence back to source GFF gene IDs when custom
  submission locus tags are used, so pipeline Sankey plots can reconcile gene
  events with annotation and consistency groups.
- Improve Sankey and event-count label placement for large unchanged groups
  and long event descriptions.

## 0.7.2 - 2026-07-30

- Accept ISO 8601 collection-date intervals in addition to single dates, and reject
  ranges whose start is later than their end.
- Calculate intron sizes directly from adjacent CDS coordinates so one-base CDS
  segments are converted without location-string parsing failures.

## 0.7.1 - 2026-07-29

- Treat BUSCO raw output as temporary: remove it after a summary is durably cached,
  retain it when BUSCO or summary parsing fails, and clean up raw output left by
  successful cached runs.
- Stop recording temporary BUSCO raw and short-summary paths in durable summary JSON
  while continuing to read summaries written by earlier releases.

## 0.7.0 - 2026-07-28

- Add genome-browser-style coordinate-duplicate gene-model plots in TSV, SVG, and
  multipage PDF formats, show the first 50 removals by default, and expose a
  configurable plot limit in pipeline reports and manifests.
- Select exact-coordinate duplicate genes from the input GFF3 and reference FASTA
  using translated-CDS validity, internal stops, start/stop completeness, ambiguous
  amino acids, CDS length, intron count, and stable input-order tie breaking.
- Record kept and removed transcripts, selection evidence, splice motifs, decision
  reasons, and low-confidence flags in a backward-compatible duplicate audit map;
  retain legacy `first` and opt-out `keep_all` policies.

## 0.6.0 - 2026-07-24

- Add `msspack run` for the complete BUSCO-to-report workflow, project-local
  `msspack_db` storage with an absolute shared-root override, database readiness
  reporting, and amalgkit-style heartbeat locks for concurrent downloads and index
  builds.
- Harden shared-database concurrency and recovery: validate lock timing relationships,
  serialize online BUSCO auto-lineage selection, require complete lineage markers,
  retain immutable content-addressed CDD versions, validate writable database roots,
  protect broad `--force-compute` targets, and record failed all-stage runs explicitly.
- Add `msspack demo` with a compact pseudonymized genome/GFF fixture, deliberately
  invalid submission metadata, deterministic pipeline-event coverage, and an optional
  offline local-reference DIAMOND configuration.
- Resolve the input scientific name through NCBI Taxonomy, cross-check configured or
  auto-selected BUSCO lineages, and preserve a cached taxonomy provenance artifact.
- Preserve UniProt/UniRef subject TaxIDs and organisms, weight description consensus by
  phylogenetic proximity for plants, animals, fungi, and prokaryotes, and generalize
  lineage-specific names from distant low-identity hits.
- Continue UniRef90 fallback after moderate/low-confidence Swiss-Prot assignments and
  choose the strongest combined taxonomic and sequence evidence.
- Add an optional DIAMOND all-vs-all functional-annotation name-consistency audit with
  stable family IDs, conservative safe-equivalent harmonization, review evidence tables,
  threshold/source plots, Sankey integration, and HTML report links.
- Label family thresholds directly, standardize all figure text at 8 pt and widths at
  3.6 or 7.2 inches, and add a gene-level name-consistency pie below the Sankey.
- Place both BUSCO pies and the name-consistency pie in one Sankey summary row, order
  functional-annotation outcomes by database priority, and use higher-contrast
  stage-colored ribbons.
- Remove name consistency from the Sankey flow itself, retain it as a summary pie, connect
  the Input GFF, CDS boundary adjustment, and functional-annotation columns to their boxed
  pies with dashed guides, and label boundary outcomes as Adjusted or No adjustment.
- Rename the adjusted BUSCO panel to Boundary-adjusted CDS, show each BUSCO sample size,
  and label the transcript outcome Already one mRNA per gene.
- Evaluate the gene-level name-consistency pie at the close-family 70/80 threshold, replace
  the misleading No comparable family category with No annotated close-family peer, and
  retain near-identical conflicts as the high-severity subset.
- Separate BUSCO lineage gene-set size from the number of CDS input sequences in Sankey
  panels, and label the comparison tier Close family peer.
- Automatically resolve name conflicts without a manual-review requirement: propagate a
  unique higher-priority product only across direct near-identical pairs and retain other
  independently supported family variation.
- Standardize every functional product against DDBJ/NCBI/EMBL-EBI naming conventions
  before homolog-family consistency analysis, preserve proposed and standardized names in
  the evidence table, generalize common cross-domain descriptions, and report all naming
  actions and residual warnings.
- Add opt-in AHRD-inspired product annotation using DIAMOND against Swiss-Prot or a
  close-reference proteome, with conservative Pfam domain fallback
- Add sequential UniRef90 DIAMOND fallback, including disk-saving taxon-scoped downloads
- Add CDD RPS-BLAST/rpsbproc fallback and a durable Pfam/CDD timing-and-yield comparison
- Cache and checksum official annotation database downloads, preserve existing products,
  and write per-transcript functional evidence and provenance artifacts
- Check annotation runtimes with `msspack doctor` and document offline/local database use
- Add functional-annotation outcomes and their exact gene flows to the pipeline Sankey
- Limit Pfam fallback to similarity-unassigned proteins and parallelize HMMER query shards
- Show Swiss-Prot, UniRef90, Pfam, and CDD outcomes separately in the pipeline Sankey

## 0.5.0 - 2026-07-22

- Add stage-aligned BUSCO composition summaries to pipeline gene-flow Sankey plots
- Improve plot sizing, spacing, label alignment, and omission of zero-count Sankey branches
- Clarify pipeline-stage terminology and optional BUSCO-to-plot workflow documentation

## 0.4.0 - 2026-07-16

- Require Python 3.11+ and remove the Python 3.10 `tomli` compatibility dependency
- Use package-index-safe README image URLs and verify all README assets in source distributions

## 0.3.0 - 2026-07-10

- Include test fixtures, examples, scripts, and documentation assets in source distributions
- Refuse to overwrite existing configs from `msspack init` unless `--force` is supplied
- Require explicit DDBJ tool installation after showing the separate license and platform requirements
- Test declared minimum development dependencies and harden GitHub Actions permissions and action pins
- Modernize Python 3.10 type syntax, cache-directory handling, and command-log quoting

## 0.2.0 - 2026-07-10

- Preserve non-mRNA annotations, comments, shared-parent features, and embedded FASTA sections during GFF processing
- Correct partial-CDS padding coordinates and retain models without explicit exon records
- Resolve products from CDS attributes and prevent ambiguous locus-tag suffix matching
- Replace mtime-only reuse with content-addressed cache state and explicit command-option cache keys
- Add atomic output publication, strict configuration/GFF validation, and safer cleanup/doctor behavior
- Pin reviewed DDBJ tool archive checksums and verify downloads before execution
- Include third-party notices in distributions and add distribution/dependency checks to CI
- Require Python 3.10+ and Biopython 1.87+ to avoid known vulnerable dependency versions

## 0.1.0 - 2026-05-15

- Internalized the MSS converter into `src/msspack/mss_converter/`
- Removed external runtime dependencies on `gff2mss`, `gffread`, `cdskit`, and `gff3sort.pl`
- Added built-in GFF sorting, padding, and annotation-table generation
- Added sequential descriptive stage names and `build-manifest.json`
- Added repo-local integration fixtures, CI, and static analysis hooks
