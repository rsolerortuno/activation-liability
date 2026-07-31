# Verification record — v0.5.0

## Scope

This record covers the frozen v0.4 primary/confirmation results plus the v0.5 paired tissue extension. Raw GEO files are not committed.

## Source archive verification

### GSE134809

- reconstructed size: 607,467,520 bytes;
- MD5: `7fa49ce56aa3aec763a9458bb3a422e6`;
- SHA-256: `2f55271ec3cdb37f9680588b2cbd8d491536694af7a585062f075ad7aa437a6c`;
- 31 valid 10x sample groups in the TAR;
- exact patient/condition mapping sourced from the authors' `sample_index.csv` at commit `d1ff0f9099a9017552d7b9b8d582bcc5a9314ae2`.

### GSE228421

- reconstructed size: 2,651,084,800 bytes;
- MD5: `2050b5505809a67b070835b175c46c25`;
- SHA-256: `15e242a9479fb0d9e5eb07f0d0b0e205d191c9d2dd3fd40962fdd3a882b84fa2`;
- 28 of 28 split parts matched their individual recorded hashes;
- exactly ten baseline V1 matrix/features members were used;
- V2 day-3 and V3 day-14 lesional samples were excluded from the primary contrast.

## Tissue ingestion accounting

### GSE134809 primary

- 9 complete patient pairs / 18 samples;
- 78,138 QC-passing cells;
- 71,431 retained broad-lineage cells;
- epithelial output excluded from inference.

### GSE134809 sensitivity

- 11 complete patient pairs / 22 samples;
- 87,749 QC-passing cells;
- 79,142 retained broad-lineage cells.

### GSE228421 primary

- 5 complete baseline pairs / 10 samples;
- 48,327,070 raw barcode columns;
- 13,709,517 non-empty barcodes;
- 120,222 QC-passing cells;
- 118,405 retained broad-lineage cells.

## Final benchmark checks

- tier-1 AUROC: 0.821429;
- tier-1 average precision: 0.887351;
- observable tier-1 controls: 15;
- NB/primary effect Spearman: 0.856452 across 503 stable fits;
- NB/primary direction agreement: 0.916501;
- Crohn 9-pair/11-pair effect Spearman: 0.939467;
- leave-one-donor-out stable drivers: 27 of 30.

## Software commands

```bash
ruff format --check .
ruff check .
mypy --strict src/activation_liability
pyright
python -m compileall -q src tests
pytest -q
coverage erase
coverage run -m pytest -q
coverage report -m
alia validate-manifest data/registry
```

Pre-package verification on 2026-07-31:

- Python 3.13.5 runtime against a Python 3.11+ package contract;
- Ruff 0.16.0: PASS;
- mypy 2.3.0 strict: PASS across 21 source modules;
- Pyright 1.1.411: 0 errors, 0 warnings;
- pytest 9.0.2: 50 core v0.5.0 tests passed before public-release validators were added;
- configured non-I/O coverage: 85.3801%;
- manifest validation: 3 manifests valid, 0 verified/default entries enabled.

The same gates are repeated after extracting the release ZIP.

## Real tissue rebuild command

```bash
PYTHONPATH=src python scripts/build_tissue_release_assets.py \
  --crohn-tar /data/GSE134809_RAW.tar \
  --psoriasis-tar /data/GSE228421_RAW.tar \
  --psoriasis-directory /cache/GSE228421_baseline \
  --output results/real/public_v0_5_0
```

The psoriasis builder uses one child process per sample to release sparse-matrix memory between samples.

## Truthfulness and exclusion checks

- README numeric claims are resolved against committed JSON in tests.
- Marker genes are asserted to be disjoint from all benchmark controls.
- No raw TAR, MatrixMarket, FASTQ, BAM, H5AD, environment, cache or downloaded lint tool is allowed into the release ZIP.
- `tissue_claims.json` explicitly lists unsupported clinical, protein and epithelial-Crohn claims.


## Clean release stage

A separate 120-file v0.5.0 release tree passed Ruff, mypy, Pyright, compileall and all 50 core tests. The subsequent public repository adds two packaging/download validator tests, for 52 total; see `PUBLIC_RELEASE_VERIFICATION.md`.
The stage contained zero TAR, MatrixMarket, FASTQ, BAM or H5AD inputs and excluded caches,
bytecode, virtual environments, downloaded type/lint tools and package metadata.
