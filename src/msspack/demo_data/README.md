# MSSPack demonstration data — not for submission

This small dataset exercises the main `msspack` workflow with realistic sequence
and gene-model structures. Its organism name, sequence IDs, locus IDs, submitter
details, and database accessions are fictional.

**Do not submit the generated MSS files.** In particular,
`Msspackia fictitia`, `PRJDB_MSSPACK_TEST_ONLY`,
`SAMD_MSSPACK_TEST_ONLY`, and `DRR_MSSPACK_TEST_ONLY` are deliberate
non-production values that must never be replaced implicitly.

The nucleotide and protein sequences are short excerpts from biological data
whose contributor authorized redistribution for this test fixture. Original
sequence and feature identifiers have been removed. Sequence similarity can
still reveal biological provenance, so this is pseudonymized test data rather
than anonymous sequence data.

## Files

- `genome.fa` and `annotation.gff3`: ten compact gene models on nine test contigs.
- `config.toml`: offline core-pipeline configuration.
- `config.functional.toml`: optional local-reference annotation configuration.
- `reference.faa`: five small reference proteins for the optional DIAMOND run.
- `expected-summary.json`: expected pipeline-event counts.

The fixture includes coordinate duplication, multiple mRNAs, frame correction,
CDS boundary adjustment, an unresolved stop-containing CDS, and conversion to
`misc_feature`.

## Run

From the directory containing these files:

```bash
msspack pack --config config.toml --no-validate
msspack plot --config config.toml
msspack report --config config.toml
```

The optional functional-annotation example requires DIAMOND but does not
download UniProt, UniRef90, Pfam, CDD, or taxonomy data:

```bash
msspack pack --config config.functional.toml --no-validate
msspack plot --config config.functional.toml
```
