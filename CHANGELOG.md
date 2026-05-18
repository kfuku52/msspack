# Changelog

## 0.1.0 - 2026-05-15

- Internalized the MSS converter into `src/msspack/mss_converter/`
- Removed external runtime dependencies on `gff2mss`, `gffread`, `cdskit`, and `gff3sort.pl`
- Added built-in GFF sorting, padding, and annotation-table generation
- Added sequential descriptive stage names and `build-manifest.json`
- Added repo-local integration fixtures, CI, and static analysis hooks
