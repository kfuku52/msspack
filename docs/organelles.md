# Organelle submissions

Configure nuclear and organelle inputs together, then generate separate submission
files. DDBJ requires separate applications for a nuclear draft genome (WGS) and
complete organelle sequences (MISC). Mitochondrial and plastid entries can share a
MISC application when the contact, data type and hold date agree. See the
[DDBJ MSS instructions](https://www.ddbj.nig.ac.jp/ddbj/mss-e.html).

## External annotation workflow

Run and review annotations with an external tool such as Chloe or PMGA first.
Provide the resulting valid annotated GenBank, or FASTA plus GFF, to msspack.
msspack's organelle stage imports the supplied models, converts them to MSS,
checks identifiers and format, and produces validation results and review files.

Annotation tools, reference databases and tool-specific output repair scripts
remain outside msspack. Its organelle stage does not install or launch annotators,
predict genes, choose new model boundaries, infer RNA editing, rotate circular
sequences or fill in missing biological product annotations. Resolve these with
the external annotation workflow before importing the result. If annotation is
absent, msspack reports that absence instead of creating a model.

## Inputs

Keep the existing nuclear `[inputs]` and metadata. Add annotated GenBank inputs:

```toml
[[organelles]]
name = "mt"
organelle = "mitochondrion"
genbank = "mitochondrion.gb"
# Per-record values override GenBank topology. Do not infer circularity from a name.
topology = { ctg1 = "circular", ctg2 = "linear" }
completeness = { ctg1 = "complete", ctg2 = "partial" }
# Optional externally calculated functional classes and reference-gene assessment.
# assessment_json = "mitochondrion.assessment.json"

[[organelles]]
name = "pt"
organelle = "plastid:chloroplast"
genbank = "plastid.gb"
completeness = { ctg1 = "complete" }
```

Alternatively use `fasta` and `gff` instead of `genbank`, with an explicit
`genetic_code` and a `topology` mapping for every sequence. GFF CDS coordinates,
strand and phase are preserved. Transcript-linked CDS segments are joined in
transcription order; ambiguous multi-parent, trans-spliced or fuzzy GFF locations
require explicit annotated GenBank locations. Unknown feature types/attributes,
missing parents and inconsistent phases fail rather than discard annotations.
FASTA without GFF is supported as a **source-only candidate**, clearly flagged as
missing biological annotation in the report. msspack does not predict organelle
genes; prepare and review those annotations with an appropriate annotator first.

`genetic_code` is optional for GenBank when each CDS already specifies its
translation table. A configured code must agree with each supplied CDS table.
Missing completeness is reported as unconfirmed; it is never promoted to a
complete genome. `complete` describes the individual sequence, including one
component of a multipart assembly, and produces “complete sequence” in its title.

## Commands and outputs

`msspack run --config config.toml` and `msspack pack --config config.toml` prepare
both submissions when `[[organelles]]` is configured. Nuclear output paths and
processing remain compatible. The organelle stage runs no nuclear cleanup,
frame adjustment, isoform selection, product reassignment or BUSCO.

To add organelles while preserving an already completed nuclear run, use a copied
config with a **new output directory**, then run only:

```bash
msspack organelles --config config.toml --nuclear-annotation /path/to/completed/nuclear.ann.txt
```

This command does not read or recompute the nuclear FASTA/GFF.
`--nuclear-annotation` checks for locus-tag collisions against
the existing MSS file without changing it. `run`/`pack` perform this check against
the nuclear pair they prepared. For independent organelle-only work the option
can be omitted. Under the configured output directory, `organelles/final/` contains:

| File | Contents |
| --- | --- |
| `<prefix>.organelles.ann.txt` / `.fasta` | One MISC annotation/sequence pair, with all organelle entries |
| `summary.json` | Input hashes, per-entry sequence hashes, lengths, topology, completeness, counts and warnings |
| `feature-audit.json` | Input/output qualifiers, locations and extracted feature sequence hashes |
| `ddbj-validation-summary.json` | Actual Parser/transChecker status |
| `index.html`, `sequences.svg`, `sequences.pdf` | Self-contained review and sequence-length plot |
| `organelle-<name>-gene-flow.sankey.svg` / `.pdf`, `.tsv` | Supplied CDS/RNA features retained in MSS, one flow per organelle input |
| `organelle-<name>-functions.svg` / `.pdf` | Externally supplied functional classes, counted per CDS feature |
| `organelle-<name>-gene-content.svg` / `.pdf`, `.tsv` | Reference-gene model recovery, with complete/fragmented/not-detected/not-assessed states |
| `organelle-<name>-assessment.json` | Imported assessment provenance and derived counts; no annotation changes |

With nuclear and organelle inputs, the main report's **Pipeline gene flow** section
stacks the nuclear chart, then each organelle chart in configuration order.
Organelle flows aggregate all records of that input (including non-circular
records), and count CDS/RNA features without counting gene container rows twice.
Each panel has its own scale. They show supplied annotations and retained MSS
features; nuclear cleanup, annotation searches and BUSCO are not run for these
panels. Source-only inputs show an explicit annotation-missing message.
`msspack plot` and `report` can render existing prepared organelle results without
repeating conversion, validation or scientific analysis; changed inputs require
organelle preparation again.

## External functional and gene-content assessments

An optional `assessment_json` path imports external results and adds functional
class and reference-gene panels below that organelle's flow. msspack does not
run BLAST, predict genes, install a reference set, or reassign products here.
The source annotation and per-record sequence SHA-256 values must match the
current inputs. The supplied locus tags must exist in the conversion audit.
Changing the external result updates report panels without rerunning annotation.

The Sankey's **Swiss-Prot annotation** column counts CDS features with an accepted
external Swiss-Prot hit, no qualifying hit, or no assessment. RNA features take a
separate **not applicable** path. Existing `gene`/`product` qualifiers and curated
functional classes do not establish Swiss-Prot hits. This column does not change
products, gene models or MSS content. Without external Swiss-Prot evidence, CDS
are explicitly **not assessed**, not counted as no-hit results.

To supply this evidence, add an optional `swissprot` object to the same JSON:

```json
{
  "swissprot": {
    "database": {"name": "Swiss-Prot", "version": "<release>",
      "fasta_sha256": "<SHA-256 of reference FASTA>"},
    "method": "External alignment tool/version and accepted-hit definition",
    "thresholds": {"evalue": 1e-10, "min_identity_percent": 35,
      "min_query_coverage_percent": 70, "min_subject_coverage_percent": 70},
    "results": [{"entry": "pt_ctg1", "locus_tag": "Example_pt000001",
      "status": "hit", "accession": "P12345",
      "evidence": "Accepted alignment identity, coverage, E-value and score"}]
  }
}
```

Supply exactly one result per supplied CDS locus, including IR copies. Allowed
states are `hit`, `no_hit` and `not_evaluated`. Only a `hit` has an accession;
all states require evidence. Use `not_evaluated` explicitly for unsearched CDS.
Database release, FASTA SHA-256, method and thresholds are required. Native
`organelle-<name>-swissprot.tsv` exports the evidence; the assessment JSON retains
additional alignment/provenance fields. The nuclear Sankey groups genes by
annotation source priority, whereas these organelle counts are independent
Swiss-Prot hit status per CDS feature, including externally named proteins.

Without this input, CDS are shown as **Not classified** and reference-gene
recovery as **Not assessed**, with no numeric score. An imported gene-content
assessment must identify its reference set/version, reference accessions,
definition, scope and per-gene evidence. msspack counts each distinct marker
once, even when an IR contains two copies. A marker with only partial/pseudo
CDS models cannot be marked complete. Any `not_evaluated` marker makes the
overall percentage unavailable rather than silently changing the denominator.

Minimal format (replace example hashes with actual SHA-256 values):

```json
{
  "schema_version": 1,
  "name": "pt",
  "organelle": "plastid:chloroplast",
  "input_sha256": ["<SHA-256 of supplied GenBank, or sorted FASTA/GFF hashes>"],
  "sequence_sha256": {"ctg1": "<SHA-256 of the upper/lower-case-preserved record sequence>"},
  "method": "External tool, version and parameters; no model changes",
  "functional_classification": {
    "assignments": [{"entry": "pt_ctg1", "locus_tag": "Example_pt000001",
      "category": "Photosystem I", "evidence": "External supplied gene/product evidence"}]
  },
  "gene_content": {
    "scope": "supplied_annotations",
    "marker_set": {"id": "example_reference_panel", "version": "1",
      "definition": "Explicitly defined expected protein-gene set",
      "references": ["REFERENCE_ACCESSION.VERSION"]},
    "markers": [{"gene": "psaA", "status": "complete",
      "matches": [{"entry": "pt_ctg1", "locus_tag": "Example_pt000001"}],
      "evidence": "External match coverage and model-completeness evidence"}]
  },
  "limitations": ["Reference-dependent gene recovery is not genome assembly completeness."]
}
```

Allowed marker states are `complete`, `fragmented`, `not_detected` and
`not_evaluated`; allowed scopes are `supplied_annotations` and
`assembly_and_annotation`. The complete-model percentage is complete markers /
all markers, and is available only when every marker was assessed. Biological
gene loss, missing annotation, RNA editing, reference bias and unresolved genome
structure require separate interpretation. This is not a BUSCO score.

Record names are prefixed with the configured organelle `name`, preventing
mitochondrial/plastid ID collisions. GenBank coding/RNA features retain their
locations, qualifiers, genetic code and existing locus tags. Absent tags are
assigned within genes; adjacent joined intervals are merged without changing the
extracted sequence. DDBJ represents gene information as qualifiers on CDS/RNA
features, so GenBank `gene` rows are projected onto those features and retained
in the audit. A gene without a CDS/RNA is retained as `misc_feature`.

Ordinary supplied translations are retained in the audit and recomputed by
transChecker. Translations accompanying `/exception` are preserved in the MSS
file, as required by Parser. The official transChecker omits these exceptions
from its FASTA output: the summary counts them separately and verifies exact
entry/location coverage of every remaining CDS in both AA and nucleotide outputs.
RNA editing and trans-splicing annotations require biological review, not
automatic correction of genomic bases. Unsupported remote locations stop with an
explicit error.

MISC is the **application category**, not a valid `DATATYPE` qualifier value in
the annotation file. The generated MISC pair omits the WGS `DATATYPE`, WGS keywords,
`submitter_seqid`, and the nuclear assembly's structured comment.

Each candidate is checked with the configured official tools before atomic
publication. A failed validation keeps the preceding complete MISC generation
available and retains the failed candidate/logs. `--no-validate` marks the report
as unvalidated; file creation does not establish DDBJ acceptance. No command
uploads a submission or sends messages.

The organelle path recognizes one specific Parser warning: JP0082 on a CDS that
explicitly declares `trans_splicing`. It independently parses and sorts the DNA
intervals and accepts warning-only exit 2 only when every interval is disjoint,
every diagnostic has that code, and Parser finished normally. All original warning
counts/messages remain in the report for curator review. True overlaps, other
warnings and errors still fail; nuclear/explicit validation remains strict.
