# Annotation-only updates of existing entries

For an existing assembly whose sequence must stay unchanged, use already prepared
MSS files with `prepare-update`. This command does not rerun trimming, gap
normalization, CDS adjustment, or protein-ID assignment. Supply a tab-separated
mapping verified against the public records, with this exact header and column order:

```text
entry	accession	submitter_seqid
scaffold1	BAAGII010000001	scaffold1
```

The row above illustrates the format; provide the complete verified mapping for
your own dataset. Do not infer accession correspondence from row order alone.

```bash
msspack prepare-update --ann original.ann.txt --fasta original.fasta \
  --mapping entries.tsv --output-dir update
msspack validate --ann update/update.ann.txt --fasta update/update.fasta
```

Install the DDBJ tools before the `validate` step. For custom Java/cache settings,
pass `--config your_config.toml`; explicit `validate` requests both tools even if
pipeline validation is disabled. Its logs and results are written beside the
annotation file; see [validation behavior](usage.md#external-tools-and-validation).

The command writes a new directory with `update.ann.txt`, `update.fasta`, and
`update-manifest.json`. It renames annotation entries and FASTA headers together,
sets explicit original `submitter_seqid` values, and preserves sequence bases,
CDS locations, locus tags, and existing `protein_id` qualifiers. The manifest
records per-sequence lengths/SHA-256 values and feature counts. Mapping rows must
be unique and match both input files exactly; conflicting source identifiers and
existing output directories are rejected. Input files are never overwritten.

Before updating, expand any `@@[entry]@@` values in `COMMON` into individual
entries using the original entry names. These macros are rejected because leaving
them in `COMMON` would change their meaning after accession renaming. Move a
COMMON `submitter_seqid` qualifier into each entry's source feature as well.

Exon and, by default, UTR blocks are removed only when a same-locus mRNA represents
their location on the same strand. Unrepresented blocks are retained and counted
for review. Use `--retain-utr-features` to keep independent UTR annotations.
The command does not verify the mapping against DDBJ, infer protein-ID inheritance,
or run official validation. Verify the mapping and sequence identity against the
published assembly, review the manifest, and run Parser/transChecker before use.
For annotation-only updates, arrange email/SFTP delivery with DDBJ-Update/MSS rather
than assuming a new MSS form is required. No files are uploaded by this command.
