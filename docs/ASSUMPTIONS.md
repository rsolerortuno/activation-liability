# Assumption ledger

| ID | Assumption | Why needed | Failure mode | Test or mitigation | Status |
|---|---|---|---|---|---|
| A01 | Donor identifiers connect resting and activated samples within study. | Paired inference. | Mispaired samples create false effects. | Manifest and runtime pairing checks; incomplete pairs excluded and counted. | Tested |
| A02 | Pseudobulk CPM is comparable between paired conditions in one study. | Effect estimation. | Global RNA-content shifts can bias compositional normalisation. | Report library sizes; support alternative size factors as future ablation. | Partly tested |
| A03 | At least three donor pairs give a minimally estimable variance. | Standard error. | Rare cell types yield unstable effects. | Hard abstention below three pairs. | Tested |
| A04 | Broad lineage labels remain comparable across conditions within study. | Cell-type-resolved contrasts. | Activation changes marker abundance or boundaries. | Condition-blind broad modules, abstention and leave-target-out marker evaluation. | Partially tested in GSE178429 |
| A05 | Study effects address the same stimulus axis and cell type. | Meta-analysis. | Pooling unlike stimuli obscures biology. | Group key includes stimulus; exact matching required. | Tested |
| A06 | DerSimonian-Laird is adequate for the small synthetic benchmark. | Dependency-light random effects. | Tau² is biased with few studies. | Report study count and use abstention; REML remains a documented alternative. | Tested analytically |
| A07 | Raw count greater than zero is a meaningful default detection event. | Positive-cell fractions. | Ambient RNA and depth alter detection. | Mandatory threshold sensitivity. | Tested synthetically |
| A08 | Matched ADT features correspond to the named target epitope. | RNA-protein corroboration. | Antibody quality or isoforms create discordance. | Deposited feature names are retained; discordance is reported. | Tested only for seven exploratory GSE157857 targets |
| A09 | Tumour reference is same-study and platform-compatible. | Selectivity erosion. | Cross-study batch masquerades as selectivity. | Reject incompatible provenance. | Tested |
| A10 | Synthetic batch effects are sufficiently adversarial. | Ablation validation. | Invalid cross-study method may appear competitive. | Test expects within-study AUROC to exceed cross-study AUROC. | Tested |
| A11 | Literature-derived controls are independent of consumed data. | Non-circular benchmark. | Labels leak from the same expression studies. | Citation/accession disjointness test. | Tested structurally |
| A12 | Absence of protein evidence cannot support a protein-level claim. | Claim truthfulness. | RNA inducibility presented as surface abundance. | Claims contract marks protein claim unsupported. | Tested |
| A13 | Public stimulation datasets represent only sampled contexts. | Interpretation. | Unmeasured inflammatory contexts are mistaken for safety. | Coverage fields and abstention; explicit limitation. | Enforced |
| A14 | Equal study weighting for positive fractions is preferable to cell weighting. | Avoid giant-study domination. | Small noisy studies receive equal influence. | Report per-study fractions; alternative weighting is an ablation candidate. | Documented |
| A15 | GSE96583 gene rows follow the deposited `batch2.genes` order for both matrices. | Assign gene identities without reconstruction. | A mismatched order silently relabels every expression row. | Require exactly 35,635 gene rows and equality with both MatrixMarket row dimensions; pin SHA-256. | Tested |
| A16 | GSE96583 metadata's extra stimulated suffix is a serialization device for barcodes shared across libraries. | Reconcile 313 repeated raw 10x sequences. | Incorrect suffix reversal swaps condition metadata or duplicates cells. | Reverse only within the stimulated sample; require uniqueness and exact equality to each matrix barcode set. | Tested |
| A17 | Demuxlet singlet donor calls are adequate pseudobulk identifiers in GSE96583. | Pair control and IFN-beta cells across eight donors. | Donor misassignment biases paired effects. | Retain singlets only; require every donor × broad-lineage pair in both conditions; record counts. | Structurally tested; dependent on deposited calls |
| A18 | SLE-PBMC replication is informative about ranking robustness but not healthy-donor generalisation. | Use the available untouched IFN-I cohort honestly. | Disease-state expression is mistaken for normal baseline liability. | Label disease context in every output; keep cohort external and prohibit healthy-donor or protein claims. | Enforced |
| A19 | Abstention of all positive holdout targets at 5% is more honest than lowering the threshold after inspection. | Preserve the frozen confirmation protocol. | Post-hoc threshold tuning creates false confirmation. | Keep 5% result primary; report 1–10% only as sensitivity and mark holdout inconclusive. | Enforced |

## A023 — Broad tissue marker modules remain usable under inflammation

**Needed because:** the deposited tissue matrices do not provide one harmonised broad annotation suitable for both cohorts.

**Failure mode:** activation or tissue remodelling changes marker expression and moves cells between broad lineages.

**Coverage:** marker modules are condition-blind, disjoint from benchmark targets, require score margins and allow `Unknown`; the skin marker-only failure was repaired before final target effects were generated. Residual misclassification remains possible.

## A024 — Source-publication Crohn exclusions are appropriate for the primary set

**Needed because:** rp6 had low cell recovery and rp16 had unusually similar paired profiles in the source analysis.

**Failure mode:** exclusion could remove biologically legitimate variation or improve results.

**Coverage:** all eleven pairs are analysed as a named sensitivity; effect correlation and direction agreement are committed.

## A025 — Maximum-LCB tissue scoring is not dominated by opportunistic context selection

**Needed because:** liability can appear in one lineage or inflammatory context and not another.

**Failure mode:** taking the maximum across two diseases inflates scores.

**Coverage:** endpoints and target-relevant lineages were frozen before final effects; per-study metrics, threshold sensitivity and donor robustness are reported. No further context is added after inspection.

## A026 — Method-of-moments NB dispersion is adequate as a sensitivity check

**Needed because:** a count-model comparison was required without introducing a larger empirical-Bayes framework.

**Failure mode:** dispersion can be unstable with five or nine pairs and many zeros.

**Coverage:** the NB model is secondary, fit accounting is explicit, perfect-separation warnings are excluded from concordance summaries, and the primary paired estimator is retained.
