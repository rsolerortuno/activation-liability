# Phase gates and adversarial review

This log records each build phase, the gate result, failures found during adversarial review,
and the decision taken before continuing. Synthetic success never promotes a real-data claim.

## Phase 0 — Design

**Gate:** passed.

**Completed:** formulas, estimands, evidence classes, abstention conditions, assumptions and
limitations were written before implementation in `docs/METHODS.md`, `docs/ASSUMPTIONS.md`
and `docs/LIMITATIONS.md`.

**Adversarial review:** an early design would have pooled every activated condition. That would
make a TCR-stimulation effect and an interferon effect share one estimate. The key was changed to
`target × cell type × exact stimulus`; the rejected pooling is recorded in Decision D003.

**Residual risk:** DerSimonian-Laird and I² are unstable when only a few studies contribute.

## Phase 1 — Skeleton

**Gate:** passed with an environment qualification.

**Completed:** Python `src/` package, Typer CLI, Pydantic manifest schema, GitHub Actions,
configuration, licence, citation metadata and initial tests.

**Adversarial review:** build isolation could not resolve `setuptools` from the restricted package
index. Editable installation succeeded with the already installed backend using
`--no-build-isolation`. No lockfile or false clean-install claim was created.

**Residual risk:** build isolation and a resolved lockfile remain unavailable in the restricted package environment. Ruff, mypy and Pyright are now locally verified from the offline bundle.

## Phase 2 — Synthetic fixtures

**Gate:** passed.

**Completed:** deterministic donor-, study-, cell-type- and condition-resolved count fixtures with
planted inducible targets, lineage targets, study batch effects, donor variation, dropout, ADT and
a same-study tumour reference.

**Adversarial review:** the first fixture was too clean and produced perfect ranking. A documented
RNA-only CD3E technical artefact without matching ADT induction was added so protein
corroboration and abstention have a falsifiable role. It is not a biological CD3E claim.

**Residual risk:** simulation realism is bounded by the mechanisms deliberately planted.

## Phase 3 — Statistics core

**Gate:** passed.

**Completed:** donor pseudobulk, complete-pair checks, paired log2-CPM differences, standard
errors, Benjamini-Hochberg correction, DerSimonian-Laird random-effects estimates, 95% confidence
intervals and I².

**Adversarial review:** an undefined cell-level p-value contaminated every Benjamini-Hochberg
adjustment with `NaN`. Correction now operates on finite p-values and restores missing values in
place. Analytic tests cover the paired effect, null, weighting and heterogeneity cases.

**Residual risk:** CPM-based inference is a transparent baseline, not a replacement for a
negative-binomial real-data sensitivity analysis.

## Phase 4 — Audit layer

**Gate:** passed.

**Completed:** positive fractions, footprint expansion, same-study selectivity erosion,
RNA/ADT concordance, YAML-driven evidence classes, hard abstention and generated `claims.json`.

**Adversarial review:** a protein-level claim could otherwise inherit an RNA-only call. The claim
contract now permits protein corroboration only for rows labelled `CONCORDANT` and explicitly
lists RNA-to-protein inference as unsupported.

**Residual risk:** evidence thresholds are deterministic but not calibrated to a real clinical
error rate.

## Phase 5 — Benchmark

**Gate:** passed for synthetic data only.

**Completed:** literature-derived tiered controls, deterministic 30% tier-1 positive holdout,
AUROC, average precision, precision@10, coverage-versus-accuracy and all required ablations.

**Adversarial review:** control citations were audited independently from registry data. Twenty-five
citations or rationales that were generic, mismatched or too weak were replaced. A test asserts that
registry accessions do not appear in control provenance. The invalid cell-level analysis produces
66 significant calls versus 19 with donor pseudobulk; the invalid cross-study ranking degrades.

**Residual risk:** control membership remains context-dependent, particularly for tier-2 targets.

## Phase 6 — Real-data ingestion

**Gate:** passed for three frozen primary cohorts and one independent external confirmation cohort.

**Completed:** schema-validated candidate manifests, checksum-pinned downloads, mocked fetch tests, a prioritized catalog, a resumable bootstrapper, real adapters, conversion QC and donor-paired analyses for GSE157857, GSE178429, GSE140244 and GSE96583. The primary benchmark remains isolated from the confirmation cells and outputs.

**Adversarial review:** the first unconstrained real benchmark reached AUROC 1.0 partly through absent-lineage negatives, repeated time opportunities and marker leakage. It was rejected. The replacement froze four endpoints, pre-declared target lineages, used a one-sided 95% lower-confidence-bound score, applied observability abstention and reannotated marker targets leave-target-out. GSE96583 was then opened without retuning. Its observable ranking replicated, but the frozen cutoff recovered only half of observable positives and every preassigned positive holdout target abstained. The result is recorded as partial confirmation with an inconclusive holdout rather than promoted to full confirmation.

A second ingestion review found 313 raw 10x barcode sequences shared across the two GSE96583 libraries. The deposited metadata disambiguates stimulated duplicates with a condition-specific suffix. The adapter reverses that encoding only within the stimulated library and requires exact one-to-one barcode equality before conversion.

**Residual risk:** no control-overlapping protein panel has been executed. GSE96583 uses SLE PBMCs and a 6 h IFN-beta endpoint, so its 18 h cross-study comparison is sensitivity-only. The tissue extension covers only Crohn ileum and psoriasis skin, and the preassigned positive holdout still lacks coverage at the frozen threshold.

## Phase 7 — Documentation and packaging

**Gate:** passed subject to the verification status in `docs/VERIFICATION.md`.

**Completed:** README, methods, assumptions, decisions, limitations, data guide, runbook,
implementation report, critical evaluation, truthfulness tests, generated synthetic artefacts and
release archive.

**Adversarial review:** README numbers are parsed from marked spans and compared with committed
JSON. Stale or copied headline values fail the test. Generated caches, editable-install metadata
and bytecode are excluded from the release archive.

**Residual risk:** software and public-data gates pass and the observable ranking has an independent frozen-rule replication. However, the positive target holdout is inconclusive, and paired tissue coverage plus protein-overlapping controls are still required before a 9/10 claim is defensible.

## v0.5 tissue extension gate

- Both large TARs reconstructed and checksum-verified: PASS.
- Exact patient-condition pair maps verified: PASS.
- Condition-blind markers disjoint from controls: PASS.
- Marker-only skin annotation QC and repair recorded before final effects: PASS.
- Donor-pseudobulk tissue effects generated: PASS.
- Negative-binomial sensitivity generated with fit accounting: PASS.
- Crohn source-QC exclusion sensitivity generated: PASS.
- Tissue benchmark retains context-specific control failures: PASS.
- Protein-level tissue corroboration: UNAVAILABLE, explicitly unsupported.
