# Public data and Google Drive index

This document separates source code, derived results and raw public inputs.

## Public root

- Google Drive root: https://drive.google.com/drive/folders/1f3K-MzEQDsUFmMIb5cnFxH47_G-HYPKN
- Results: https://drive.google.com/drive/folders/1AQjydxf1Z_Y0bdJ27IQc8Ca1gwVlHPcy
- Releases and manifests: https://drive.google.com/drive/folders/1vWfWckJ-euzDClLFYG-hpmDGjyGL0yyw

## What belongs in GitHub

Commit these items:

- source code, tests and configuration;
- compact target-restricted or pseudobulk derived results;
- reports, plots and machine-readable claims;
- source URLs, checksums and ingestion QC;
- download/validation notebooks.

Do not commit multi-gigabyte GEO archives, FASTQ/BAM files, virtual environments or credentials.

## Executed public cohorts

| Accession | Context | Modality | Role | Drive location |
|---|---|---|---|---|
| GSE157857 | IFN-beta PBMC/myeloid | RNA + ADT + HTO | Primary discovery | `01_DISCOVERY_PRIMARY_RNA_ADT/` |
| GSE178429 | IFN-gamma/LPS PBMC | RNA | Multi-stimulus validation | `02_VALIDATION_IFN_TNF_TLR/` |
| GSE96583 | IFN-beta PBMC | RNA | Frozen external replication | `02_VALIDATION_IFN_TNF_TLR/` |
| GSE140244 | CD4-memory activation | bulk RNA | Time-course validation | `04_SORTED_BULK_TIMECOURSES/` |
| GSE134809 | paired Crohn ileum | RNA | Paired tissue extension | `03_TISSUE_INFLAMMATION/` |
| GSE228421 | paired psoriasis skin | RNA | Paired tissue extension | `03_TISSUE_INFLAMMATION/` |

## Next phase: GSE190564

- Drive destination: https://drive.google.com/drive/folders/11hz8VhGV2bcSGp3dWmeQlufIsIeZQAzF
- GEO record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190564
- Processed archive: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190564/suppl/GSE190564_processed_data.tar.gz
- SOFT metadata: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190564/soft/GSE190564_family.soft.gz
- Expected processed archive size: approximately 11.1 GB.

GSE190564 contains colon and PBMC pools with GEX, ADT, TCR and BCR data. The target use case is
paired ulcerative-colitis inflamed versus non-inflamed colon, separated into epithelium and
stromal/immune fractions where metadata support it. No GSE190564 result is part of v0.5.0.

## Existing technical protein reference outside the core benchmark

The connected Drive also contains a Xenium renal-cell-carcinoma dataset with RNA and a 27-protein
panel including PD-1, PD-L1, LAG-3 and HLA-DR. It may support technical RNA/protein concordance
checks, but it is tumour tissue and does not replace paired normal inflamed tissue.

## Sharing recommendation

Use **Anyone with the link → Viewer** for the public root. Grant **Editor** only to named
collaborators. Keep a checksum manifest in every release folder so accidental or malicious changes
are detectable.
