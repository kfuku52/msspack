# Input and output integrity

GFF3 parent lists are split before percent decoding, so an ID containing `%2C`
remains one identifier. Coordinate duplicate removal preserves features with a
surviving parent. Transcript extraction is shared by padding and annotation,
accepts child-first input, and applies the initial CDS phase in transcription
order after splicing. A `gene` row is not required for a coding transcript.
Parentless transcripts participate in annotation tables and pipeline metrics as
gene models identified by their transcript IDs. Their generated locus tags use
that ID with the configured prefix, and padding uses the transcript's strand.
Isoform selection uses the declared transcript span even when CDS/UTR rows appear
first. Without exon rows, frame and padding adjustments retain the combined
CDS/UTR span. MSS conversion joins repeated CDS IDs directly under a gene or
without a parent; anonymous gene-direct CDS rows are grouped by their parent.
Distinct gene-direct CDS IDs remain separate coding units throughout extraction,
frame/padding adjustment, functional annotation, BUSCO, and MSS rendering. A
single gene-direct coding unit retains its parent ID in extracted FASTA and
annotation tables; multiple units use their CDS IDs. Gene-level plots count the
parent once, prefer the source order existing/reference/Swiss-Prot/UniRef90/Pfam/CDD,
and show the most actionable consistency status (review, resolved, unannotated,
no close-family peer, consistent). Per-CDS evidence remains in the TSV files.
Stop-containing genes retain the existing gene-wide CDS-to-misc conversion policy;
converted feature counts and distinct converted gene counts are recorded separately.

Gap normalization shifts coordinates after a changed N run. Boundaries within
the run retain their offset where possible and are clipped to the last retained
N when the run shrinks. Removing a complete annotated feature with a zero-length
gap is an error. In MSS conversion, an exon/CDS segment entirely within an N run
retains its uncertain location; the configured gap policy applies to the whole
feature whenever any segment overlaps a gap.
The same policy applies to gene-direct and parentless CDS features.

Supported semicolon repairs run before strict GFF parsing and gap normalization.
Their changes remain in `04.fix-gff-semicolons.log` and its metrics JSON, even
though that stage now runs earlier. Empty sequences after terminal-N removal
are errors; they are not silently removed or emitted with an invalid source.
BUSCO input CDS extraction also uses the repaired GFF with the original genome
coordinates. Its `proteins` mode translates with the configured genetic code,
strand, and initial CDS phase. Locus-tag padding changes identifier attributes
and references without changing sequence IDs, free-text metadata, or embedded FASTA.

## Submission publication

The sample locus-tag prefix is also used in output filenames. Configuration rejects
POSIX/Windows path separators and colons before creating output directories.

`pack` builds candidate final files in `.msspack-work/final/` under the output
root. Enabled validation runs against these candidates. Only successful candidates
are published to `final/`, including their validation summary. Failed runs keep
their logs and candidate files for diagnosis, and mark the build manifest failed;
the previously published annotation/FASTA pair remains available.

On Linux/macOS, `final/` points to a complete directory under
`.msspack-generations/`. The pointer is replaced atomically. Consumers requiring a
stable snapshot across several file reads should resolve `final/` once and read
both files from that resolved directory. Existing generations remain available;
they are not automatically deleted while another reader could still use them.

Upgrading a legacy real `final/` directory requires a directory-to-pointer
migration. Windows uses a complete-directory rename instead of a symlink. These
transitions can briefly leave `final/` absent, but never publish half of a pair.
On an ordinary failure the old directory is restored. If interrupted between
renames, `.msspack-previous-final/` retains the old pair and the next publication
recovers it. This is a filesystem limitation of replacing a nonempty directory;
subsequent POSIX pointer updates do not have that gap.

Project-changing commands hold an OS-backed `.msspack-output.lock`. Another
writer to the same resolved output directory fails with an actionable error.
Nested stages in a single operation share the lock, and the OS releases it after
a process exits. Run independent projects with different output directories;
database sharing continues to use separate database locks.

## Cache identity and database integrity

Stage cache identity includes all Python sources in the installed msspack package,
in addition to declared input/output content hashes and settings. Changing a
transitively imported helper invalidates caches without requiring a version bump.
Pipeline external commands and validation Java commands include executable content
fingerprints, so replacing a binary at the same pathname invalidates its cache.
BUSCO manifests record content hashes of configuration, source inputs, and compared
stage inputs. Plots, reports, and taxonomy inference exclude results whose inputs
have changed. A `run --no-busco` marks retained BUSCO results skipped and excludes
them from current outputs without deleting the historical summaries.

Materialized databases and content-addressed objects are checked against recorded
digests before reuse. DIAMOND and Pfam indexes include builder identity and content
hashes. A different builder gets a separate index location; corrupt files are
rebuilt under the existing database locks. This does not automatically select a
new upstream database release: changing the configured source is still explicit.
CDD support-file versions are checked against their expected content hashes and
provenance before reuse; corrupt versions are rebuilt under the version lock.

Use `run --force-compute` when deliberately recomputing analyses. Downloaded
databases and published generations are retained. The atomic publication contract
applies to the final submission files; analysis logs, plots and the top-level
build manifest describe the latest attempted operation.
