# Decision log

## D001 — DerSimonian-Laird rather than REML

**Decision:** implement DerSimonian-Laird as the default random-effects estimator.

**Reason:** it has an analytically testable closed form and avoids adding a specialised meta-analysis dependency. **Rejected alternative:** REML, which generally has better small-sample properties but needs iterative fitting and stronger convergence handling. A future implementation should add REML as an ablation, not silently replace committed results.

## D002 — CPM pseudobulk with paired log differences

**Decision:** sum counts per donor and use paired log2(CPM+1) differences.

**Reason:** the estimand is transparent and works for synthetic fixtures and generic count matrices. **Rejected alternative:** fitting edgeR/DESeq2-style negative-binomial models, which would add an R runtime or heavy Python dependency and complicate a standalone CLI. Real-data users should treat CPM inference as a conservative baseline, not a final publication model.

## D003 — Exact stimulus in the meta-analysis key

**Decision:** meta-analyse only identical stimulus labels inside a shared axis.

**Reason:** merging TCR, IFN and TNF effects would produce a biologically uninterpretable headline. **Rejected alternative:** pooling all “activated” conditions.

## D004 — Study-balanced positive fractions

**Decision:** average per-study cell-positive fractions rather than pooling cells.

**Reason:** prevents a large study from dominating the footprint. **Rejected alternative:** pooled-cell fraction, retained only as a possible future sensitivity analysis.

## D005 — YAML rule engine for evidence class

**Decision:** thresholds are loaded from versioned YAML and evaluated generically.

**Reason:** supports auditability and sensitivity without hidden code constants.

## D006 — Real registry entries are excluded by default

**Decision:** candidate real studies ship with `verified: false` and `default: false` until accession contents, conditions, donor pairing, licence and DOI are checked together.

**Reason:** a correct accession with an incorrect design claim is still fabricated provenance. The sandbox does not execute source downloads.

## D007 — Synthetic tumour reference

**Decision:** synthetic fixtures include a same-study tumour-reference CPM for testing `selectivity_erosion`.

**Reason:** this tests formula and provenance guards without claiming real tumour safety evidence.

## D008 — Offline lint/type tool bundle

**Decision:** install Ruff 0.16.0 and mypy 2.3.0 from user-provided Linux/Python 3.13 wheels, and Pyright 1.1.411 from a user-provided npm tarball. Preserve exact versions in the verification record.

**Reason:** the configured package index could not resolve these tools, but a local artifact bundle permits real execution without weakening the gates. `mypy --strict` remains authoritative; Pyright is supplemental basic-mode checking because scientific Python library stubs are incomplete.

## D009 — Optimised sensitivity

**Decision:** calculate donor effects and meta-analysis once per dataset, then recompute threshold-dependent fractions and ranking scores for each sweep point.

**Reason:** detection thresholds do not alter the pseudobulk LFC. Re-running all statistics for every grid point is redundant and made the suite unnecessarily slow.

## D010 — Planted RNA-only artefact in the synthetic benchmark

**Decision:** synthetic activated T cells include a strong CD3E RNA shift without matching ADT change.

**Reason:** a benchmark in which every positive is obvious and every negative is clean cannot demonstrate the value of protein corroboration or abstention. The artefact is explicitly synthetic and technical; it is not a biological claim that CD3E is inflammation-induced. The protein-aware score improves ranking but deliberately leaves one low-confidence holdout error, so increasing abstention has measurable benefit.

## D011 — Four frozen endpoints for the real benchmark

**Decision:** use IFN-beta 18 h (GSE157857), IFN-gamma 6 h and LPS 6 h (GSE178429), and anti-CD3/CD28 24 h (GSE140244) as the only primary real endpoints.

**Reason:** the initial maximum across all available times created unequal opportunities and could reward noise. **Rejected alternative:** selecting each target's best observed time after inspection. Earlier and later time points remain a secondary time-course output.

## D012 — Lineage-constrained 95% lower-confidence score

**Decision:** score each target by the maximum `effect - 1.96 × SE` across its pre-declared expressing lineages and the four frozen endpoints.

**Reason:** a liability detector should surface any credible covered induction, but a three-donor noisy point estimate should not dominate a stable effect. **Rejected alternatives:** unconstrained maximum across all cell types and raw fold change without uncertainty penalty.

## D013 — Coverage abstention rather than absent-target negatives

**Decision:** require a pre-declared lineage to be present and a maximum positive-cell fraction of at least 5% in either state. Otherwise return abstention.

**Reason:** epithelial or plasma targets absent from PBMC cohorts are not evidence of non-inducibility. **Rejected alternative:** assign zero scores, which made the first benchmark artificially easy.

## D014 — Leave-target-out cell annotation

**Decision:** when `CD3E`, `CD79A` or `MS4A1` is evaluated in GSE178429, remove that gene from its lineage marker module before assigning cells.

**Reason:** a target must not directly select the cells in which it is benchmarked. **Rejected alternative:** excluding these controls entirely; leave-target-out retains hard lineage controls while removing direct feature leakage.

## D015 — Single-study heterogeneity is non-estimable

**Decision:** return `NaN`/JSON `null` for I², Q and tau² when one study contributes to an exact endpoint.

**Reason:** zero means observed homogeneity, not absence of information. Single-study rows may meet B/C rules from donor evidence but can never meet class A.

## D016 — Real holdout is diagnostic

**Decision:** label the current real target holdout `DIAGNOSTIC_NOT_CONFIRMATORY`.

**Reason:** the deterministic split was preserved, but the constrained endpoints and score were finalised after an earlier unconstrained real run was inspected. **Rejected alternative:** describing it as untouched. GSE96583 and a paired tissue cohort are reserved for confirmation under frozen v0.3 rules.

## D017 — Do not reconstruct GSE96583 gene order

**Decision:** do not analyse the available GSE96583 matrices until `GSE96583_batch2.genes.tsv.gz` is present and checksum-recorded.

**Reason:** matrix row order is provenance, not an inferable convenience. A plausible annotation reconstructed from another release could silently corrupt every target effect.

## D018 — Real protein benchmark unavailable

**Decision:** report GSE157857 ADT concordance only for seven exploratory targets and mark the control benchmark protein ablation unavailable.

**Reason:** receptor and ligand names are not interchangeable, and none of the ADT features exactly matches a benchmark control target.

## D019 — Open GSE96583 exactly once under frozen v0.3 rules

**Decision:** analyse GSE96583 only after the authoritative `GSE96583_batch2.genes.tsv.gz` file is present and checksum-recorded. Keep the original four-endpoint benchmark unchanged and write GSE96583 to separate confirmation artefacts.

**Reason:** this preserves the previous primary result and prevents an external cohort from silently becoming another tuning endpoint. The lineage map, 95% lower-confidence score, 5% observability threshold, deterministic target split and frozen classification cutoff are reused without adjustment. **Rejected alternative:** add GSE96583 to the maximum-across-endpoints primary score after seeing its effects.

## D020 — Treat the GSE96583 target holdout as inconclusive

**Decision:** report successful ranking of observable controls but do not call the preassigned target holdout confirmed.

**Reason:** `IL2RA`, `PDCD1LG2` and `TNFRSF9` all fall below the frozen 5% observability threshold in GSE96583. The cohort therefore contains no observable positive holdout target. **Rejected alternatives:** lower the primary threshold after inspection, count abstentions as errors, or claim confirmation from the remaining negatives alone.

## D021 — Resolve duplicated GSE96583 barcodes condition-specifically

**Decision:** reverse the deposited `-11` row-name suffix only in stimulated metadata rows, then require exact set equality and uniqueness against the stimulated matrix barcodes.

**Reason:** 313 raw 10x barcode sequences occur in both libraries. The metadata keeps them unique by appending one extra `1` to the second occurrence. **Rejected alternatives:** merge metadata by raw barcode without condition, drop all duplicated sequences, or infer sample identity from expression.

## D022 — Keep cross-time IFN-I heterogeneity descriptive

**Decision:** compare GSE96583 IFN-beta 6 h with GSE157857 IFN-beta 18 h as a labelled sensitivity analysis, not as exact-endpoint replication.

**Reason:** time point and disease context differ. Two-study random-effects and `I²` values can expose disagreement but cannot establish homogeneity for a common endpoint. **Rejected alternative:** pool them as if they were exchangeable technical replicates.

## D023 — Treat tissue validation as a separate frozen extension

**Decision:** preserve all v0.4 results unchanged and write GSE134809/GSE228421 outputs to `public_v0_5_0`.

**Reason:** the tissue cohorts add biological coverage but should not retroactively improve the diagnostic blood benchmark. **Rejected alternative:** merge new endpoints into v0.4 target scores after seeing tissue effects.

## D024 — Use nine Crohn pairs as primary and eleven as sensitivity

**Decision:** exclude rp6 and rp16 from the primary contrast exactly as declared by the source analysis, but process both in a committed sensitivity.

**Reason:** this respects published QC while exposing the effect of the exclusion. **Rejected alternatives:** silently drop the pairs or override the source QC without sensitivity reporting.

## D025 — Disable epithelial claims in GSE134809

**Decision:** retain epithelial cells only in sample-QC accounting and exclude that lineage from Crohn pseudobulk inference.

**Reason:** the source dissociation protocol was not designed to recover epithelium reliably. **Rejected alternative:** use sparse epithelial recovery to improve lineage coverage.

## D026 — Repair skin annotation before final benchmark and discard provisional outputs

**Decision:** replace `SDC1/XBP1` plasma scoring with a multi-marker plasma gate and add `KRT1/KRT10` to epithelium after marker-only QC identified keratinocyte contamination.

**Reason:** the diagnosis used canonical markers, not target effects or labels. The configuration status records the repair. **Rejected alternatives:** keep implausible plasma counts, tune markers using control performance, or report the provisional result.

## D027 — Process each psoriasis sample in a fresh process

**Decision:** extract only ten baseline matrix/features pairs and process one sample per child process, persisting compact pseudobulk statistics.

**Reason:** raw matrices contain millions of barcode columns and caused sparse-memory fragmentation in a monolithic run. **Rejected alternatives:** commit filtered derivative matrices, use cells as a statistical unit, or require excessive memory.

## D028 — Keep context-specific control failures

**Decision:** retain the positive score for negative control `CD79A` in psoriasis B cells and the negative score for positive control `CD69` in Crohn NK cells.

**Reason:** changing tiers or contexts after inspection would make the benchmark circular. These failures are informative about control-label scope. **Rejected alternative:** move them to tier 2 after seeing the result.
