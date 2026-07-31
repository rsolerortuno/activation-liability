# Implementation report — v0.5.0

## What was built

Version 0.5.0 retains the frozen v0.4 blood/sorted-cell benchmark and external GSE96583 replication, then adds a separate paired tissue-inflammation extension:

- GSE134809 Crohn ileum: nine primary involved/uninvolved patient pairs; two publication-excluded pairs retained only as sensitivity;
- GSE228421 psoriasis skin: five baseline lesional/non-lesional patient pairs; post-treatment visits excluded;
- condition-blind broad-lineage annotation with marker modules disjoint from all benchmark targets;
- compact donor pseudobulk sufficient statistics instead of a persistent cell-level target table;
- paired log-CPM effects, positive-cell fractions, footprint expansion and target scoring;
- donor-fixed-effect negative-binomial GLM sensitivity with log-library offsets;
- leave-one-donor-out robustness, observability-threshold sensitivity and Crohn 9-pair versus 11-pair sensitivity;
- source checksums, result hashes, claims contract and standalone HTML report;
- a memory-isolated psoriasis builder that processes one very large raw MatrixMarket sample per child process.

The final CLI adds `alia tissue-validate`. The reproducible release helper is `scripts/build_tissue_release_assets.py`.

## Frozen analysis design

The tissue score is the maximum one-sided 95% lower confidence bound across the two predeclared tissue endpoints and target-relevant broad lineages. A target abstains if its lineage is not covered or if its maximum resting/activated positive-cell fraction is below 5%.

Crohn epithelial claims are explicitly disabled because the source study's dissociation design did not support epithelial recovery. Psoriasis uses only baseline lesional versus baseline non-lesional samples. The two additional Crohn patients are not used in the primary score.

## What broke and how it was fixed

### Large Drive files could not be materialised directly

The connector truncated downloads above approximately 100 MB. The user-supplied split parts were reconstructed in order. Every part was checked against its recorded hash and each complete TAR was checked against its source MD5 before extraction.

### Psoriasis deposited raw barcode universes, not filtered matrices

Individual matrices contained millions of barcode columns and tens of millions of non-zero entries. Reading ten samples in one Python process caused memory fragmentation. The final builder reads only baseline matrix/features members and processes one sample per fresh child process, persisting only donor-level sufficient statistics.

### Initial marker-only skin annotation misclassified keratinocytes as plasma cells

The first broad classifier used `SDC1/XBP1`. In non-lesional P1, 94% of provisional plasma-labelled cells expressed `KRT10`, while fewer than 1% expressed canonical plasma markers `MZB1`, `JCHAIN` or `DERL3`. Those provisional tissue results were discarded. Before generating the final benchmark, the classifier was repaired with a multi-marker plasma gate and `KRT1/KRT10` epithelial markers. No benchmark target, control label or target effect was used in the repair.

### The tissue control benchmark was not perfectly separated

This was retained rather than tuned away. `CD79A`, a tier-1 negative control, increased in psoriasis B cells. `CD69`, a tier-1 positive control, decreased in Crohn NK cells. These failures reduce AUROC but demonstrate that the benchmark remains context-sensitive and adversarial.

## Final real-data results

### Primary tissue controls

- observable tier-1 controls: 15 (8 positives, 7 negatives);
- AUROC: 0.821429;
- average precision: 0.887351;
- exact one-sided Mann-Whitney p-value: 0.020047;
- expanded tier-1+2 AUROC: 0.761364;
- expanded tier-1+2 average precision: 0.904126.

### Added coverage

Four targets that were not observable in the v0.4 blood/sorted-cell benchmark became observable in tissue: `ERBB2`, `FAP`, `TNFRSF17` and `VCAM1`.

Eight targets showed a positive broad-lineage footprint expansion at the 10% positivity threshold. The largest expansion was three lineages.

### Count-model sensitivity

Of 572 attempted paired negative-binomial fits, all converged; 69 emitted perfect-separation warnings and were excluded from concordance summaries. Among 503 stable fits:

- Spearman correlation with paired log-CPM effects: 0.856452;
- direction agreement: 0.916501.

### Patient robustness

- 27 of 30 observable target drivers remained directionally stable after removing any one donor;
- Crohn primary-nine versus all-eleven effects had Spearman correlation 0.939467 and direction agreement 0.924342.

## Engineering verification

The release contains 48 tests before final packaging, strict Ruff/mypy/Pyright gates, README truthfulness checks for synthetic, primary real, external confirmation and tissue JSON, and configured non-I/O coverage above 85%. The final verification record is regenerated after packaging.

## What remains unsupported

- clinical toxicity prediction;
- surface-protein induction for the tissue targets;
- a real same-study tumour comparator and real `selectivity_erosion`;
- true delivery, exposure or therapeutic-index modelling;
- epithelial conclusions from GSE134809;
- pristine positive-target holdout confirmation.

## Three highest-risk assumptions

1. Broad condition-blind marker modules remain sufficiently stable under disease-associated tissue remodelling.
2. Maximum-LCB scoring across two tissue contexts remains an acceptable ranking rule despite context-specific control failures.
3. RNA footprint expansion is useful for target-safety triage even when matched target-level surface protein is unavailable.
