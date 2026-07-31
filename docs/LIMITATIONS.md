# Limitations

This document intentionally precedes results.

1. **Acute stimulation is not lifelong exposure.** Public experiments usually measure hours or days of activation. They cannot reproduce chronic inflammation, repeated dosing, tissue repair, developmental expression, or cumulative exposure.
2. **RNA is not surface protein.** Transcript induction can fail to produce accessible surface protein because of translation, trafficking, shedding, internalisation or epitope masking. ADT helps but is itself antibody- and panel-dependent.
3. **There is no true delivery or toxicity ground truth.** The benchmark labels activation inducibility, not clinical toxicity. The tool does not model biodistribution, target density thresholds, payload release, Fc biology, on-target pharmacology or organ reserve.
4. **Stimulus coverage is incomplete.** IFN, TLR, TNF, TCR and tissue-damage panels do not span the diversity, timing and combinations of clinical inflammation.
5. **Public stimulation panels do not span clinical inflammatory diversity.** A negative audit means “not induced in covered contexts,” never “safe.”
6. **Cell-type annotation can move under activation.** Resting and activated cells may be labelled differently, creating apparent gains or losses.
7. **Pseudobulk normalisation is compositional.** Strong global shifts can affect CPM-based fold changes.
8. **Random-effects estimates are fragile with few studies.** I² is imprecise when `k` is small and should not be overinterpreted.
9. **Protein corroboration is sparse and target-panel limited.** Missing ADT is not biological absence.
10. **The synthetic benchmark remains an implementation test.** A public real-data benchmark is now computed, but it is small, coverage-limited and diagnostic rather than confirmatory.


## Additional limitations of the v0.4 real benchmark

11. **The real target holdout is not pristine confirmation.** The constrained protocol was finalised after inspection of an earlier unconstrained real result. The split is retained for diagnostics only.
12. **Coverage is predominantly blood and sorted CD4 T cells.** Four tier-1 negative controls abstain because epithelial and plasma-cell lineages are absent. There is no evidence about tissue-resident accessibility from these cohorts.
13. **Protein overlap is zero for benchmark controls.** The available ADT panel supports exploratory RNA/ADT checks for seven other targets but cannot validate the control ranking at protein level.
14. **Exact-endpoint heterogeneity is not estimable.** Each endpoint currently has one contributing study. Leave-one-endpoint-out analysis measures dependence on axes, not between-study heterogeneity.
15. **Three-donor studies have coarse inference.** For GSE157857 and IFN-gamma in GSE178429, the smallest attainable one-sided paired Wilcoxon p-value is 0.125.
16. **Gene mapping for GSE140244 is indirect.** Stable Ensembl identifiers are mapped using the deposited GSE157857 GRCh38 feature table. This is audited, but a release-matched annotation source would be preferable.
17. **GSE96583 is an SLE-patient cohort, not a healthy-donor replication.** It independently tests the frozen computational rules but disease-associated baseline state may alter both observability and inducibility.
18. **The GSE96583 target holdout is not evaluable at the frozen threshold.** All three preassigned positive holdout targets abstain at 5% observability. The perfect observable-control ranking therefore does not constitute target-holdout confirmation.
19. **GSE96583 has no matched surface-protein layer.** Its confirmation is RNA-only and cannot resolve RNA/protein discordance.
20. **The two IFN-I studies differ in time point and context.** GSE96583 is 6 h in SLE PBMCs; GSE157857 is 18 h in another donor context and contributes only broad myeloid cells. Cross-time `I²` is descriptive.
21. **Tissue coverage remains limited to Crohn ileum and psoriasis skin.** The v0.5 paired extension does not span clinical inflammatory diversity, other organs or chronic exposure.
22. **No real same-study tumour comparator is included.** Real `selectivity_erosion` therefore remains unsupported.
23. **Perfect separation can be fragile.** The primary tier-1 result contains nine observable positives and six observable negatives; the GSE96583 replication contains six positives and five negatives. Control-set selection uncertainty is not captured by exact Mann–Whitney p-values.

## Additional limitations of the v0.5 tissue extension

24. **Broad annotation is heuristic.** Marker modules are condition-blind and leakage-controlled, but they do not replace source-author or reference-mapped cell labels.
25. **Psoriasis has only five pairs.** Confidence intervals and leave-one-patient-out results remain coarse.
26. **The primary interval uses a normal multiplier.** `1.96 × SE` can be optimistic for five or nine pairs; small-sample hierarchical intervals are future work.
27. **Maximum-across-context scoring can inflate ranking.** Only two frozen tissue contexts are used, but the score still selects the strongest eligible context.
28. **Crohn epithelial liability is unsupported.** The source dissociation design precludes epithelial inference, even though epithelial rows appear in raw QC summaries.
29. **Tissue protein corroboration remains unavailable.** Both paired tissue cohorts are RNA-only.
30. **Control labels are context-limited.** `CD79A` and `CD69` violate their expected labels in specific tissue lineages, demonstrating that literature control sets are not universal ground truth.
31. **The NB model is deliberately simple.** Method-of-moments dispersion and donor fixed effects provide sensitivity, not publication-grade empirical-Bayes shrinkage.
32. **Tissue composition and state remain entangled.** Cell-type resolution reduces but cannot fully eliminate inflammation-associated shifts, doublets or transitional states.
