# Real-data runbook — v0.4.0

## 1. Install and verify

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
ruff format --check .
ruff check .
mypy --strict src/activation_liability
pyright
coverage run -m pytest -q
coverage report -m
```

The verified offline bundle uses Ruff 0.16.0, mypy 2.3.0 and Pyright 1.1.411.

## 2. Required directory layout

```text
activation-liability-real-data/
├── 01_DISCOVERY_PRIMARY_RNA_ADT/
│   └── GSE157857_IFNB_CITEseq_PRIMARY/
│       └── GSE157857_RAW.tar
├── 02_VALIDATION_IFN_TNF_TLR/
│   ├── GSE178429_MULTI_STIMULUS_scRNAseq/
│   │   ├── GSE178429_PBMCs_stim_scRNAseq_cellMeta.txt.gz
│   │   ├── GSE178429_PBMCs_stim_scRNAseq_counts.txt.gz
│   │   └── GSE178429_PBMCs_stim_scRNAseq_geneNames.txt.gz
│   └── GSE96583_IFNB_demuxlet_scRNAseq/
│       ├── GSE96583_RAW.tar
│       ├── GSE96583_batch2.total.tsne.df.tsv.gz
│       └── GSE96583_batch2.genes.tsv.gz
└── 04_SORTED_BULK_TIMECOURSES/
    └── GSE140244_CD4_ACTIVATION_TIMECOURSE/
        ├── GSE140244_rnaseq_gene_counts.txt.gz
        └── GSE140244_rnaseq_meta_data.txt.gz
```

Never reconstruct a missing gene-order file or infer donor/condition labels from expression.
The GSE96583 adapter verifies the deposited 35,635-row order against both matrices and fails
closed if barcode reconciliation is not exact.

## 3. Execute the frozen benchmark and independent confirmation

```bash
alia real-validate \
  --data-root /path/to/activation-liability-real-data \
  --output results/real/public_v0_4_0 \
  --controls data/controls/controls.yaml \
  --rules config/evidence_classes.yaml \
  --benchmark-config config/real_benchmark.yaml

python scripts/build_real_release_assets.py \
  --data-root /path/to/activation-liability-real-data \
  --output results/real/public_v0_4_0
```

The command keeps the frozen primary benchmark and GSE96583 confirmation separate. GSE96583
cannot alter the original primary cells, scores, cutoff or target split.

Expected key files:

```text
results/real/public_v0_4_0/
├── ingestion_qc.json
├── real_cells.csv.gz
├── audit.json
├── target_summary.csv
├── claims.json
├── benchmark.json
├── target_scores.csv
├── donor_robustness.csv
├── timecourse_validation.csv
├── gse96583_confirmation_ingestion_qc.json
├── gse96583_confirmation_cells.csv.gz
├── gse96583_confirmation_audit.json
├── gse96583_confirmation_claims.json
├── gse96583_confirmation_benchmark.json
├── gse96583_confirmation_target_scores.csv
├── gse96583_confirmation_donor_robustness.csv
├── source_checksums.json
├── result_manifest.json
└── report.html
```

## 4. Inspect mandatory guardrails

Before quoting performance, confirm:

- primary `holdout.integrity_status` is `DIAGNOSTIC_NOT_CONFIRMATORY`;
- confirmation `status` is `PARTIAL_EXTERNAL_CONFIRMATION_WITH_TARGET_HOLDOUT_ABSTENTION`;
- confirmation `target_holdout.status` is `INCONCLUSIVE_NO_OBSERVABLE_POSITIVE_HOLDOUT_TARGETS`;
- confirmation uses the frozen cutoff recorded in its protocol;
- `protein_benchmark.status` remains `UNAVAILABLE_NO_CONTROL_TARGET_OVERLAP`;
- exact-endpoint heterogeneity remains non-estimable where only one study contributes;
- uncovered targets abstain instead of being scored as zero;
- the primary endpoint list exactly matches `config/real_benchmark.yaml`.

## 5. Interpret GSE96583 correctly

GSE96583 is an external SLE-PBMC IFN-beta 6 h cohort with eight demuxlet-resolved donors. It can
support no-retuning ranking and classification checks, but not healthy-donor generalisation or
surface-protein claims.

The observable control ranking replicated. The frozen cutoff recovered half of observable positives
with no false positives. The three preassigned positive holdout targets were not observable at 5%,
so the holdout is inconclusive—not a negative result and not confirmation. Lower-threshold results
are sensitivity analyses only.

## 6. Add paired tissue validation

A paired inflamed/uninflamed tissue cohort must preserve patient pairing and broad tissue cell types.
Register the accession, DOI, licence, exact checksums, pairing fields and exclusions before inspecting
target-level results. Tissue endpoints must remain separate and must not retune the PBMC benchmark.

## 7. Promotion and publication checklist

- Validate every registry YAML file.
- Archive source URLs, sizes, SHA-256 and licence snapshots.
- Run from a clean checkout on a second machine.
- Compare generated artefacts byte-for-byte with the committed manifest.
- Report weakened or inconclusive confirmation unchanged.
- Do not claim toxicity prediction, therapeutic window or protein corroboration without matching evidence.

## Paired tissue extension — v0.5.0

### Recommended memory-isolated build

```bash
PYTHONPATH=src python scripts/build_tissue_release_assets.py \
  --crohn-tar /data/GSE134809_RAW.tar \
  --psoriasis-tar /data/GSE228421_RAW.tar \
  --psoriasis-directory /cache/GSE228421_baseline \
  --controls data/controls/controls.yaml \
  --config config/tissue_extension.yaml \
  --baseline-benchmark results/real/public_v0_4_0/benchmark.json \
  --output results/real/public_v0_5_0
```

When the ten psoriasis baseline members have already been extracted and compact pseudobulk/QC have already been generated, the CLI can be used directly:

```bash
alia tissue-validate \
  --crohn-tar /data/GSE134809_RAW.tar \
  --psoriasis-directory /cache/GSE228421_baseline \
  --psoriasis-pseudobulk /cache/gse228421_pseudobulk.csv.gz \
  --psoriasis-qc /cache/gse228421_qc.json \
  --output results/real/public_v0_5_0
```

Expected headline output is recorded in `tissue_execution_summary.json`. Do not compare absolute counts between the two studies or enable Crohn epithelial claims.
