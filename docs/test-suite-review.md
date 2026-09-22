# Test suite review (2026-09-22)

Reviewed all 34 original test modules against the question: what realistic defect
would removing this check allow? Test counts and coverage percentages were not
retention criteria. Production behavior is unchanged.

## Removed or consolidated

| Area | Decision and remaining protection |
| --- | --- |
| Plotting | Remove the synthetic full-pipeline plotting fixture and its duplicate helpers. The real demo now exercises the plot/report CLI, manifest paths, generated files, obsolete-output cleanup, and expected biological event counts. Keep legacy log parsing, functional-annotation flows, validation labels, and large-flow bounds. |
| Chart styling | Remove exact label offsets, font widths, page dimensions, annotation display order, and CSS snapshots from pipeline, BUSCO, and consistency plots. These constrained layout changes without establishing readable or accurate charts. Keep semantic labels, counts, bounds, and embedded-font protection. |
| CLI/template | Consolidate init creation, refusal to overwrite, forced replacement, and example/template synchronization through actual CLI output. Remove the private-template test module. Use only consumed fields in remaining boundary stubs. |
| CLI mocks | Exercise real database override handling and doctor report rendering. Plot/report wiring moves to the demo; transcript-selection wiring moves to the reversed-order tie regression. Keep argument translation tests for expensive external workflows. |
| Configuration | Replace 13 separately constructed invalid-config scenarios with one table using the valid minimal fixture and field-specific errors. Reuse that fixture in manifest, report, and explicit-validation setup. Distinct validation rules remain represented. |
| Cache | Remove the extra fingerprint helper scenario: the action test covers change detection and the integration regression changes content without changing size or mtime. Fold unchanged-rerun reuse into the existing minimal fixture run. |
| Execution | Replace sleep-based overlap detection with a barrier. Drop the sequential worker-count duplicate; sequential validation failure still verifies that the next tool does not run. |
| Validation | Remove the no-op job test and path/default structure checks. Use actual job scheduling with only external tools stubbed. Merge persisted warning/failure checks into the stronger sequential failure scenario. |
| Conversion | Remove direct source-qualifier rendering duplicated by the exact expected annotation fixture. Fold phase normalization into shared-parent indexing; mixed-feature conversion covers ordinary parent traversal. Remove a Biopython start-codon-table lookup check. |
| Installation | Combine incomplete-root recovery with installation metadata verification. Keep corruption, checksum, archive traversal, missing-tool, Java command, and whitespace-path checks. |
| BUSCO | Remove the extra transcriptome parsing setup and cache-root concatenation test. Summary consumers already parse transcriptome fixtures. Cleanup now checks the real calculated path before the temporary-directory context removes it. |
| Input/selection | Consolidate three invalid FASTA inputs with the same no-partial-output assertion. Remove the weaker same-order transcript tie scenario; keep reversed order, shared parents, orphan features, and embedded FASTA cases. |
| Vacuous checks | Remove a heartbeat test whose nondecreasing-mtime assertion passed without any heartbeat, regex matches against hard-coded demo strings, and a POSIX shlex wrapper round-trip. |

## Reviewed and retained

These scenarios have distinct failure consequences; superficial similarities do
not make them duplicates.

| Modules (`tests/test_*.py`) | Concrete failures still targeted |
| --- | --- |
| `annotation_consistency`, `annotation_taxonomy`, `annotation_table`, `product_names` | Wrong name harmonization, lineage weighting, CDS product fallback, or destructive product-name normalization. |
| `functional_annotation` | Active database replacement, mixed concurrent database versions, wrong evidence selection, malformed external metadata, skipped or repeated search queries. Local files and downloaded databases use different paths. |
| `database_lock`, `databases`, `busco` | Unsafe lock recovery, lost mutual exclusion, accepting incomplete databases, wrong lineage selection, losing failure diagnostics or persistent summaries. |
| `gff`, `gff_adjustments`, `gap_normalization`, `padding_tools` | Wrong strand/phase/coordinates, lost parents or features, stale terminal codons, or misaligned sequence and annotation. |
| `mss_converter`, `mss_postprocess`, `header`, `submission_update` | Invalid submission structure, unintended CDS conversion, changed sequence bases/protein IDs, accession mismatches, or overwritten output. |
| `coordinate_duplicate_plots` | Truncating exported data when only the display should be limited, ignoring configured limits, or breaking older duplicate maps. |
| `doctor`, `tools`, `scripts`, `workflow` | Missed unusable inputs/tools, unsafe extraction or deletion, damaged shared databases, wrong execution order, or missing failure records. |
| `report`, `pipeline_plots`, `demo`, `integration_pack` | Stale validation status, unescaped HTML, incorrect metrics, invalid output wiring, or unsafe cache reuse. |
| `external_tools` | Real DDBJ/BUSCO compatibility cannot be established by local stubs. Keep the opt-in external checks unchanged. |
| `config`, `cli`, `pipeline`, `pipeline_actions`, `select_one_mrna`, `utils` | Invalid user settings, wrong dispatch/options, cached failed work, invalid FASTA publication, wrong transcript selection, or ignored platform cache locations. |

The normal suite changes from 257 to 216 test methods; table cases are still
executed as subtests. No runtime improvement is claimed from these single runs.
The three external-tool tests remain opt-in and are skipped in normal discovery.
