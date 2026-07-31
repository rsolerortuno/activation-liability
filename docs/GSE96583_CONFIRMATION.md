# GSE96583 external confirmation

## Inputs

```text
GSE96583_RAW.tar
GSE96583_batch2.total.tsne.df.tsv.gz
GSE96583_batch2.genes.tsv.gz
```

All three files are SHA-256 pinned in
`results/real/public_v0_4_0/source_checksums.json`.

## Integrity checks

- 35,635 deposited gene rows match both MatrixMarket row dimensions.
- Sample accessions define control and IFN-beta 6 h conditions.
- Eight demuxlet donor identities and singlet calls are retained.
- 313 raw 10x barcode sequences occur in both libraries; condition-specific metadata suffixes are reconciled and checked one-to-one.
- Every retained donor × broad-lineage combination has both conditions.

## Frozen protocol

The cohort was analysed with the v0.3 lineage map, 5% observability rule, one-sided 95% LCB score,
control labels, target split and classification cutoff unchanged. Primary and confirmation datasets are
kept separate.

## Result

The 11 observable tier-1 controls ranked with AUROC and average precision of 1.0. At the frozen
classification cutoff, balanced accuracy was 0.75, sensitivity 0.50 and specificity 1.00. Score
correlation with the primary benchmark was 0.8411.

The preassigned positive holdout targets `IL2RA`, `PDCD1LG2` and `TNFRSF9` all failed the frozen 5%
observability requirement. The target holdout is therefore inconclusive. Lower-threshold results are
retained only as sensitivity analyses.

## Boundaries

GSE96583 contains SLE-patient PBMCs, has no matched surface-protein measurement and uses IFN-beta at
6 h rather than the primary cohort's 18 h endpoint. It supports a no-retuning ranking replication,
not healthy-donor, protein-level or exact-time confirmation.
