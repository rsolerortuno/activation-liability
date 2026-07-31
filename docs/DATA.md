# Data contract

## Input table

The generic long-form CSV reader expects:

- `study`, `donor`, `cell_id`, `cell_type`, `condition`, `stimulus`, `target`
- `rna_count`, `library_size`
- optional `adt_count`, `adt_library_size`, `tumour_reference_cpm`

Conditions must be exactly `resting` or `activated`. A donor is eligible only when both conditions exist for the same study, cell type and stimulus.

## Registry

Each YAML file in `data/registry` is validated independently. Only entries with both `verified: true` and `default: true` are fetched by default. The shipped candidate entries are deliberately excluded.

## Cache

Fetchers write to a user-supplied cache directory. Files are stored under a SHA-256 content address and accompanied by metadata containing source URL, retrieval time and checksum. No downloaded data are committed.


## Executed public-data layer (v0.4.0)

The committed real-data artefacts were generated from four downloaded public cohorts:

- `GSE157857`: HTO-resolved resting versus 18 h IFN-beta CITE-seq, restricted to broad myeloid cells;
- `GSE178429`: donor-paired PBMC IFN-gamma and LPS contrasts, with the 6 h endpoints used in the frozen primary benchmark;
- `GSE140244`: sorted CD4-memory anti-CD3/CD28 time course, with 24 h used as the frozen primary endpoint;
- `GSE96583`: eight-donor control versus 6 h IFN-beta PBMC scRNA-seq, used only as the frozen external confirmation cohort.

Raw GEO files are not committed. `results/real/public_v0_4_0/real_cells.csv.gz` contains the target-restricted primary analysis table. `gse96583_confirmation_cells.csv.gz` is kept separate so the unopened cohort cannot silently alter the original primary benchmark. Neither file replaces the source deposits.

`source_checksums.json` records SHA-256 hashes and sizes for all nine local input files. `ingestion_qc.json` and `gse96583_confirmation_ingestion_qc.json` record donor counts, condition alignment, exclusions and annotation-leakage controls. `gse96583_confirmation_benchmark.json` records the unchanged-rule external evaluation.

For GSE96583, `GSM2560248` is the batch-2 control and `GSM2560249` is the batch-2 IFN-beta 6 h sample. The deposited gene-order file has 35,635 rows and is required to interpret both MatrixMarket files. The metadata contains eight demuxlet donor identities and singlet/doublet status. A condition-specific serialisation suffix distinguishes 313 barcode sequences shared between the two libraries; the adapter verifies the corrected mapping exactly.


## Paired tissue extension (v0.5.0)

### GSE134809

The complete `GSE134809_RAW.tar` is checksum-verified locally but not committed. Exact Crohn patient/condition mapping comes from the source authors' `input/tables/sample_index.csv` at commit `d1ff0f9099a9017552d7b9b8d582bcc5a9314ae2`.

The primary analysis includes nine involved/uninvolved ileum pairs: rp5, rp7, rp8, rp10, rp11, rp12, rp13, rp14 and rp15. rp6 and rp16 are sensitivity-only. PBMC samples are excluded. Epithelial inference is disabled.

### GSE228421

The complete raw TAR is checksum-verified locally but not committed. Only ten baseline V1 members are read: one lesional and one non-lesional matrix for each of P1–P5. Day-3 and day-14 post-risankizumab samples are excluded from the primary endpoint.

Because the deposited matrices include millions of raw barcode columns, each baseline sample is processed independently and only donor-level sufficient statistics are persisted.

### Committed tissue artefacts

`results/real/public_v0_5_0/` contains:

- `tissue_pseudobulk.csv.gz`: compact donor/condition/lineage target sums;
- `tissue_audit_rows.csv`: paired primary effects;
- `tissue_negative_binomial.csv`: count-model sensitivity;
- `tissue_target_scores.csv`: coverage, abstention and drivers;
- `tissue_footprint.csv`: broad-lineage footprint counts;
- `tissue_donor_robustness.csv`: leave-one-donor-out checks;
- `crohn_all11_sensitivity_audit_rows.csv`: source-QC exclusion sensitivity;
- `tissue_ingestion_qc.json`: sample-level QC and lineage counts;
- `tissue_source_checksums.json`: data and frozen-analysis input hashes;
- `tissue_benchmark.json`, `tissue_claims.json`, `tissue_execution_summary.json`;
- `tissue_result_manifest.json` and `tissue_report.html`.
