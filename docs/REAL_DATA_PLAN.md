# Real-data validation plan

## Goal and score boundary

A 9/10 methods-and-portfolio grade requires a preregistered real benchmark, not more synthetic
polish. The minimum defensible package is: one donor-paired RNA+ADT discovery study, at least two
independent donor-paired RNA studies spanning different stimuli, one tissue-inflammation
validation, a sealed target holdout, and an analysis report that preserves study-level estimates.
Clinical toxicity prediction remains explicitly out of scope.

## Recommended sequence

### Stage 1 — primary discovery and protein corroboration

**GSE157857** is the primary study. GEO reports PBMCs from three distinct donors cultured for 18
hours with or without IFN-beta, followed by CD1c dendritic-cell/monocyte enrichment, TotalSeq-A and
10x profiling. RNA, ADT and HTO libraries are present. Donor identities must be recovered from HTO
metadata before pseudobulk aggregation.

**GSE156473** is retained only as an exploratory activation/protein dataset. It must not contribute
to inferential confidence until independent donor replication is demonstrated from metadata.

### Stage 2 — independent stimulus validation

**GSE178429** supplies explicit donor-labelled control, IFN-gamma, LPS and PMA/ionomycin samples at
one and six hours. The primary plan excludes GolgiPlug samples or treats them as a separate stratum.

**GSE96583** supplies an independent IFN-beta PBMC response. The apparent study-level two-library
design is usable only if demuxlet metadata provide cell-to-donor assignments; otherwise it is an
exploratory cell-level atlas and not a valid pseudobulk replication.

### Stage 3 — lineage kinetics and inflamed tissue

**GSE140244** contains a 24-donor anti-CD3/CD28 time course in sorted CD4 memory T cells. Processed
counts and metadata are public; raw reads require dbGaP phs002259. Time contrasts must be
pre-registered before target outcomes are viewed.

**GSE134809** and **GSE228421** have now been executed as the v0.5 paired tissue extension.
Exact Crohn pairing was resolved from the source-author sample index; psoriasis uses five baseline
lesional/non-lesional pairs. The primary Crohn set contains nine pairs and two publication-QC
exclusions are committed as sensitivity. These cohorts test generalization but do not replace
RNA+ADT corroboration.

## Preregistered analysis decisions

1. Freeze target controls and the 30% tier-1 positive holdout before downloading expression files.
2. Use only within-study, within-donor contrasts. No cross-study absolute-expression inference.
3. Aggregate single-cell counts by donor, cell type, condition and stimulus.
4. Keep stimulus and time point in the meta-analysis key; do not pool IFN-I, IFN-II, TLR and T-cell
   activation into one effect.
5. Retain paired donor log-CPM as the transparent primary estimator and fit a donor-fixed-effect
   negative-binomial pseudobulk model as a publication-oriented sensitivity analysis.
6. Require at least three paired donors for an inferential study-level estimate. Two-donor results
   are descriptive; one-donor results force abstention.
7. Require ADT direction agreement for protein-level claims. RNA-only induction can support only an
   RNA-level liability claim.
8. Report target-level leave-one-study-out stability, I-squared, threshold sensitivity and the full
   abstention curve.
9. Do not tune thresholds on the sealed holdout.
10. Preserve every downloaded file checksum and every conversion log.

## What can realistically raise the grade

- **8.0–8.5:** all software gates pass plus GSE157857 and GSE178429 run end to end, with honest
  failures and donor pairing demonstrated.
- **8.5–9.0:** add an independent IFN-I replication, a T-cell activation time course, a sealed real
  holdout and negative-binomial pseudobulk sensitivity.
- **9.0–9.5:** paired tissue inflammation, donor robustness and a frozen release artefact are now
  present; the remaining differentiator is benchmark-target protein corroboration in more than one
  study and another untouched tissue replication.
- **10:** not supportable from public observational/stimulation data alone. It would require
  prospective external validation linked to target accessibility, exposure and clinically relevant
  safety outcomes.

## Access and API requirements

No API key is required for public GEO supplementary files. A dbGaP authorization is required only
for controlled raw reads such as phs002259; the public processed GSE140244 matrices are sufficient
for the planned bulk validation. An NCBI API key is optional for high-rate metadata queries and is
not needed by the committed downloader.
