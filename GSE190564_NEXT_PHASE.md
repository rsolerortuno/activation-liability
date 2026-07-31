# GSE190564 next phase — paired colitis CITE-seq

## Why this dataset

GSE190564 is the highest-value next addition because it combines single-cell gene expression and
CITE-seq surface-protein measurements in colon and blood. The GEO design includes active
ulcerative colitis (`UC_I`), paired non-inflamed areas from the same endoscopy (`UC_NI`), checkpoint
colitis contexts and healthy controls. Pools include PBMC, epithelium-enriched and stromal/immune
fractions.

This dataset can test whether RNA liabilities detected by `alia` are corroborated at the protein
level and whether tissue inflammation expands surface-target coverage.

## Official files

- GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190564
- Processed data: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190564/suppl/GSE190564_processed_data.tar.gz
- Series metadata: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190564/soft/GSE190564_family.soft.gz

## Drive destination

https://drive.google.com/drive/folders/11hz8VhGV2bcSGp3dWmeQlufIsIeZQAzF

Expected contents:

```text
GSE190564_UC_PAIRED_CITESEQ/
├── GSE190564_processed_data.tar.gz
├── GSE190564_family.soft.gz
├── download_manifest.json
├── archive_inventory.tsv
├── pool_inventory.tsv
└── split/                         # optional 95 MB parts for portability
```

## Frozen analysis questions

Before looking at target outcomes, the next phase should freeze:

1. Which UC patients have both `UC_I` and `UC_NI` tissue.
2. Which hashtag/pool identifiers map cells back to patients and conditions.
3. Which GEX and ADT barcodes match exactly.
4. Which antibodies overlap benchmark targets and which are isotypes or controls.
5. Which epithelial and stromal/immune contrasts satisfy minimum donor and cell coverage.
6. How ADT background is normalised and how an antibody is declared observable.

## Primary endpoint

The primary endpoint should be a within-patient contrast:

```text
UC inflamed colon − UC non-inflamed colon
```

RNA and ADT should be analysed separately, then compared at target × lineage level. Cross-patient
or cross-pool absolute-expression comparisons must not replace paired inference.

## Promotion criteria

GSE190564 may support a stronger evidence class only when:

- patient/condition demultiplexing is unambiguous;
- GEX and ADT cell identifiers align;
- at least three paired patients contribute to a target-lineage estimate;
- antibody identity and controls are documented;
- results reproduce under donor leave-one-out and alternative ADT normalisation;
- the current v0.5 scoring and target labels are not retuned after inspection.

## What it cannot prove

Even matched ADT does not establish in-vivo drug delivery, receptor occupancy, internalisation or
clinical toxicity. Those require separate pharmacology and safety evidence.

## Notebook

Use `notebooks/download_GSE190564_to_drive.ipynb`. It downloads resumably, validates the archive,
computes SHA-256, inventories GEX/ADT members and optionally creates 95 MB split parts.
