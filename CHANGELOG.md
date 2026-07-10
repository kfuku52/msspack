# Changelog

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
