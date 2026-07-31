# Methods

## Scope

`activation-liability` detects a specific evidence gap: a target can appear selective when normal tissue is represented only by resting cells, yet become detectable on normal cells after inflammation or immune activation. The software is a liability detector with abstention. It does not predict toxicity, therapeutic index, delivery, or clinical outcome.

## Statistical unit and estimand

For single-cell RNA data, cells are observations used to estimate donor-level aggregates; donors are the statistical units. For each study, target, cell type, stimulus and donor, counts are summed separately in resting and activated conditions. Library-size normalisation is performed as counts per million (CPM), and the paired donor effect is

\[
d_i = \log_2(\mathrm{CPM}_{i,activated}+c)-\log_2(\mathrm{CPM}_{i,resting}+c),
\]

where `c=1` is the default pseudocount. The within-study effect is the mean paired difference

\[
\hat\theta_s = n_s^{-1}\sum_i d_i,
\]

with standard error `sd(d)/sqrt(n)`. A study requires at least three complete donor pairs. No primary analysis compares absolute expression across studies.

## Random-effects meta-analysis

Study effects for the same `(target, cell type, stimulus)` are combined using a DerSimonian-Laird random-effects model. With inverse-variance fixed-effect weights `w_s=1/v_s`, Cochran's statistic is

\[
Q=\sum_s w_s(\hat\theta_s-\hat\theta_F)^2.
\]

The between-study variance is

\[
\tau^2=\max\left(0,\frac{Q-(k-1)}{\sum w_s-\frac{\sum w_s^2}{\sum w_s}}\right).
\]

Random-effects weights are `1/(v_s+tau²)`. The meta-analytic effect, standard error and 95% interval are derived from these weights. Heterogeneity is

\[
I^2=\max(0,(Q-(k-1))/Q)\times100\%.
\]

A single study is reported without a heterogeneity claim.

## Multiple testing

Within-study paired tests use a two-sided one-sample t-test on donor log-fold changes. Benjamini-Hochberg correction is applied across all tested target/cell-type/stimulus combinations in one audit run. The effect estimate and confidence interval remain primary; adjusted p-values are supporting evidence.

## Metrics

### `inducibility_lfc`

The random-effects meta-analytic log2 fold change for activated versus resting normal cells, with standard error, 95% confidence interval, number of studies, donor pairs and I².

### `positive_fraction_resting` and `positive_fraction_activated`

For a target and cell type, the fraction of cells whose raw count exceeds `detection_count_threshold` (default `0`). Fractions are calculated within study and condition and then averaged with equal study weight, avoiding domination by the largest cell dataset.

### `footprint_expansion`

For target `g`, let `P_rest(g,c)` and `P_act(g,c)` denote the study-balanced positive fractions for normal cell type `c`. With positivity fraction cutoff `q` (default `0.10`):

\[
\mathrm{footprint\_expansion}(g)=\sum_c 1[P_{act}(g,c)\ge q]-\sum_c1[P_{rest}(g,c)\ge q].
\]

The unit is **number of normal cell types**. Negative values are possible and are retained.

### `selectivity_erosion`

When a tumour reference from the same study and compatible platform is available, selectivity is

\[
S_x(g)=\log_2(\mathrm{CPM}_{tumour}(g)+1)-\max_c\log_2(\mathrm{CPM}_{normal,x}(g,c)+1),
\]

for normal state `x`. Selectivity erosion is

\[
E(g)=S_{resting}(g)-S_{activated}(g).
\]

The unit is log2 fold-selectivity lost. A positive value means activation reduced apparent tumour selectivity. Cross-study tumour-normal comparisons are unsupported.

### `protein_concordance`

For matched ADT data, RNA and ADT donor-paired effects are computed on the same study/cell-type/stimulus contrast. `CONCORDANT` requires the same direction and both effects at or above the configured LFC cutoff. `DISCORDANT` means opposite direction or RNA induction without protein induction. `UNAVAILABLE` means no matched protein measurement. The output includes RNA and ADT effects; the category is not a substitute for them.

### `evidence_class`

Evidence class is deterministic and configured in `config/evidence_classes.yaml`.

- `A`: at least two studies, at least six total donor pairs, I² no greater than 50%, RNA LFC at least 1, activated positive fraction at least 0.10, and concordant protein evidence.
- `B`: at least one study, at least four donor pairs, I² no greater than 75%, RNA LFC at least 0.75 and activated positive fraction at least 0.10. Protein may be unavailable; claims remain RNA-level.
- `C`: at least one study, at least three donor pairs, I² no greater than 90%, RNA LFC at least 0.5 and activated positive fraction at least 0.05.
- `INSUFFICIENT`: any hard abstention rule fires, or no class threshold is met.

Hard abstention rules are: target undetected in both states; no covered stimulus relevant to the target's observed lineage; I² above 90%; fewer than three paired donors; or a requested surface-protein claim without protein corroboration.

## Threshold sensitivity

The `sensitivity` command sweeps detection count thresholds, cell-positive-fraction cutoffs and inducibility LFC cutoffs. For each setting it ranks targets by an evidence-weighted score and reports Spearman rank correlation and top-k overlap against the default. Statistical effects are computed once because changing a cell-detection threshold does not change the paired pseudobulk LFC. If median rank correlation is below 0.80, the ranking is labelled unstable.

## Benchmark

Controls are stored in `data/controls/controls.yaml`; labels are based on literature evidence and are never derived from registry expression datasets. Tier 1 contains conservative controls used for primary metrics. Thirty percent of tier-1 positives are assigned to a deterministic holdout before benchmark metric calculation. No parameter tuning uses holdout labels.

Primary metrics are AUROC, average precision, precision at 10 and a coverage-versus-accuracy curve based on increasingly strict evidence thresholds. Required ablations compare within-study versus invalid cross-study contrasts, donor pseudobulk versus cell-level testing, RNA-only versus protein-corroborated scoring, and detection thresholds.

## Deliberately invalid ablations

Cross-study absolute contrasts and cell-level hypothesis tests are included only under `ablations`. They are named `INVALID_FOR_INFERENCE` in result artefacts. They must never populate primary audit fields or claims.

## Real public-data benchmark protocol (v0.4.0 release; frozen v0.3 primary rules)

### Frozen endpoints

The primary public-data benchmark uses exactly four endpoints:

| Axis | Study | Endpoint |
|---|---|---|
| IFN-I | GSE157857 | IFN-beta, 18 h, HTO-resolved broad myeloid cells |
| IFN-II | GSE178429 | IFN-gamma, 6 h, broad PBMC lineages |
| TLR | GSE178429 | LPS, 6 h, broad PBMC lineages |
| Lymphocyte activation | GSE140244 | anti-CD3/CD28, 24 h, sorted CD4-memory T cells |

Earlier GSE178429 time points and the remaining GSE140244 time course are secondary. They cannot replace the primary endpoint after results are inspected.

### Lineage-constrained uncertainty-adjusted score

For target \(t\), endpoint \(e\) and pre-declared expressing lineage \(l\), the paired donor effect and standard error are \(\hat\beta_{tel}\) and \(SE_{tel}\). The conservative endpoint value is

\[
LCB_{tel} = \hat\beta_{tel} - 1.96\,SE_{tel}.
\]

The real benchmark score is

\[
S_t = \max_{(e,l) \in \mathcal{C}_t} LCB_{tel},
\]

where \(\mathcal{C}_t\) contains only the four frozen endpoints and the target's lineages declared in `config/real_benchmark.yaml`. This maximum represents the liability-detector question “is there at least one covered normal lineage and inflammatory axis with credible induction?” It is not a global average response.

A target is observable only when at least one allowed row has

\[
\max(f^{rest}_{tel}, f^{act}_{tel}) \ge 0.05.
\]

Otherwise the benchmark abstains. Targets assigned to epithelial or plasma-cell lineages abstain when those lineages are absent; they are not assigned a favourable zero score.

### Annotation leakage protection

GSE178429 broad lineages are assigned with condition-blind marker modules. When the evaluated target is itself a module marker (`CD3E`, `CD79A` or `MS4A1`), that gene is removed from the marker set before cells are assigned for that target. This leave-target-out procedure prevents the target's own counts from directly determining its evaluation subset.

### Heterogeneity

For one contributing study, the effect and its within-study standard error are reported, but \(I^2\), \(Q\) and \(\tau^2\) are non-estimable. They are not replaced by zero. Evidence classes B and C may be assigned from one study if donor, effect, fraction and relevance thresholds pass; class A still requires at least two studies and concordant protein evidence.

### Donor robustness

For each target's selected endpoint, the benchmark also reports:

- fraction of paired donors with positive log2-CPM difference;
- one-sided paired Wilcoxon p-value;
- minimum and maximum mean effect after removing each donor once;
- whether every leave-one-donor-out mean remains positive.

These checks do not replace the primary paired effect. They expose three-donor instability and outlier dependence.

### Real holdout and external confirmation status

The target split remains deterministic, but the v0.3 primary endpoint and score rules were finalised after an earlier unconstrained real benchmark was inspected. Therefore the original real holdout remains diagnostic.

GSE96583 was subsequently opened once under the unchanged v0.3 lineage map, score definition, `z = 1.96`, 5% observability threshold and deterministic target split. Its only confirmation endpoint is donor-paired IFN-beta at 6 h. The external-cohort ranking is evaluated on all observable controls and the frozen v0.3 classification cutoff is applied without refitting.

The preassigned positive holdout targets are `IL2RA`, `PDCD1LG2` and `TNFRSF9`. All three are below the frozen 5% observability threshold in GSE96583, so the target-holdout result is `INCONCLUSIVE_NO_OBSERVABLE_POSITIVE_HOLDOUT_TARGETS`. This is not converted into a negative call. Thresholds from 1% to 10% are reported only as sensitivity analysis and cannot replace the primary rule.

GSE96583 therefore provides partial external confirmation of the observable ranking, not full confirmation of the target holdout. An independent paired tissue cohort is still required.

### GSE96583 barcode and lineage safeguards

The two deposited 10x libraries share 313 raw barcode sequences. The metadata serialisation keeps row names unique by appending one extra `1` to duplicated stimulated barcodes (`-1` becomes `-11`). The adapter reverses this suffix only inside the stimulated condition and requires exact one-to-one equality between matrix barcodes and condition-filtered metadata before reading counts.

Donor identity and singlet status come from deposited demuxlet metadata. Sample identity, not expression, defines control versus IFN-beta. Broad PBMC lineages are assigned with the same condition-blind marker modules used elsewhere; `CD3E`, `CD79A` and `MS4A1` are re-annotated with the evaluated target removed from its marker module.

### IFN-I cross-time replication sensitivity

GSE96583 at 6 h and GSE157857 at 18 h are compared only as a descriptive IFN-I sensitivity analysis in broad myeloid cells. Score/effect correlation, direction agreement and two-study random-effects summaries are reported. Because time point and donor disease context differ, these rows do not resolve exact-endpoint heterogeneity and cannot be used to claim low `I²`.

### Protein benchmark status

The GSE157857 ADT panel contains `CCR2`, `CD1C`, `FCER1A`, `FCGR1A`, `MRC1`, `TNFSF18` and `TNFSF9`. None is a benchmark control target. RNA/ADT concordance is therefore exploratory for these seven targets, and no protein-aware improvement metric is computed for the control benchmark.

## Paired tissue-inflammation extension (v0.5.0)

### Frozen cohorts and contrasts

The tissue extension is separate from the v0.4 PBMC/sorted-cell benchmark.

- GSE134809: involved versus uninvolved ileum within patient. The primary set contains nine pairs. Two source-publication QC exclusions are included only in sensitivity. Epithelial claims are disabled.
- GSE228421: baseline lesional versus baseline non-lesional skin within each of five patients. Post-treatment visits are excluded.

### Condition-blind broad lineage annotation

Each sample is annotated independently of condition using fixed broad-lineage marker modules. The union of marker genes is asserted to be disjoint from all benchmark targets. Cells abstain as `Unknown` when the top score is below the configured minimum or the margin over the second score is too small. Plasma assignment additionally requires at least two detected plasma markers.

The annotation is a coarse stratification mechanism, not a biological reference atlas. Target-level claims always include the broad lineage and cohort that drove them.

### Compact pseudobulk sufficient statistics

For each donor, condition, broad lineage and target, the adapter stores only:

- target UMI sum `y`;
- total library sum `L`;
- positive-cell count;
- retained cell count.

Raw cells are never treated as statistical replicates and are not persisted in the release.

### Primary tissue estimator

For donor `i`, target `g` and lineage `c`:

\[
d_{igc}=\log_2(10^6 y_{igc,A}/L_{ic,A}+1)-
         \log_2(10^6 y_{igc,R}/L_{ic,R}+1).
\]

The within-cohort effect is `mean(d)`, with standard error `sd(d)/sqrt(n)`. The reported 95% interval is `effect ± 1.96 × SE`. At least three complete pairs are required. Benjamini-Hochberg correction is applied across all tissue target-lineage tests.

Heterogeneity statistics are not reported for a single study. The two diseases are not pooled as if they represented one exchangeable biological endpoint.

### Tissue target score

For each control target, only predeclared relevant lineages are eligible. The score is

\[
S_g=\max_{e,c}\left(\hat\theta_{gec}-1.96\,SE_{gec}\right),
\]

where `e` is one of the two frozen tissue contexts and `c` is an eligible covered lineage. A target abstains when no eligible lineage is covered or its maximum resting/activated positive fraction is below 5%.

### Tissue footprint expansion

For each target, study-specific positive fractions are first collapsed by taking the maximum within each broad lineage. At positivity threshold `p=0.10`:

\[
F_R=\sum_c I(f_{gc,R}\ge p),\qquad
F_A=\sum_c I(f_{gc,A}\ge p),\qquad
\Delta F=F_A-F_R.
\]

The same broad lineage is counted once even when represented in both tissue studies.

### Negative-binomial sensitivity

For each target-lineage-study endpoint, a negative-binomial GLM is fitted to pseudobulk counts:

\[
\log E[y_{ij}]=\log L_{ij}+\beta_0+\beta_A A_{ij}+\gamma_i,
\]

where `A` indicates inflamed/lesional condition and `gamma_i` are donor fixed effects. Dispersion is estimated by a method-of-moments rule after scaling counts to the median library size. Fits with perfect-separation warnings are reported but excluded from estimator-concordance summaries.

The NB model is a sensitivity analysis. The paired log-CPM estimator remains primary because the NB dispersion estimate is intentionally simple and the number of pairs is small.
