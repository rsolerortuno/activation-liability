"""Conservative adapters for the real public validation cohorts.

The adapters deliberately emit only benchmark targets and broad cell classes. They do not
perform cross-study normalization. Each emitted contrast remains paired within donor and study.
"""

from __future__ import annotations

import gzip
import tarfile
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import IO

import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from scipy.io import mmread

BROAD_MARKERS: dict[str, tuple[str, ...]] = {
    "T_cell": ("CD3D", "CD3E", "TRAC", "CD247"),
    "NK_cell": ("NKG7", "GNLY", "KLRD1", "PRF1", "CTSW"),
    "B_cell": ("MS4A1", "CD79A", "CD79B", "CD37", "CD22"),
    "Myeloid": ("LST1", "TYROBP", "FCER1G", "LILRB1", "CTSS", "AIF1"),
    "Platelet": ("PPBP", "PF4", "NRGN", "RGS18"),
}

ADT_TARGETS: dict[str, str] = {
    "ADT_GITRL-CGAGAACCGACCAAA": "TNFSF18",
    "ADT_CD1c-GAGCTACTTCACTCG": "CD1C",
    "ADT_CD64-AAGTATGCCCTACGA": "FCGR1A",
    "ADT_FCERIA-CTCGTTTCCGTATCG": "FCER1A",
    "ADT_CCR2-GAGTTCCCTTACCTG": "CCR2",
    "ADT_CD206-TCAGAACGTCTAACT": "MRC1",
    "ADT_4-1BBL-ATTCGCCTTACGCAA": "TNFSF9",
}

GSE178429_CONTRASTS: tuple[tuple[str, str, str, str], ...] = (
    ("IFN_1h", "Control_1h", "IFN_II_1h", "GSE178429_IFN_1h"),
    ("IFN_6h", "Control_6h", "IFN_II_6h", "GSE178429_IFN_6h"),
    ("LPS_1h", "Control_1h", "TLR_1h", "GSE178429_LPS_1h"),
    ("LPS_6h", "Control_6h", "TLR_6h", "GSE178429_LPS_6h"),
)


@dataclass(frozen=True)
class RealBuildResult:
    """A standardised long table plus machine-readable quality-control summaries."""

    cells: pd.DataFrame
    qc: dict[str, object]


def load_control_targets(path: Path) -> list[str]:
    """Load externally curated benchmark target names."""

    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    return sorted(str(entry["target"]) for entry in payload["controls"])


def _read_gzip_lines(handle: IO[bytes]) -> list[str]:
    with gzip.GzipFile(fileobj=handle) as stream:
        return [line.decode("utf-8").rstrip("\n") for line in stream]


def _read_gzip_matrix(handle: IO[bytes]) -> sparse.csr_matrix:
    with gzip.GzipFile(fileobj=handle) as stream:
        return sparse.csr_matrix(mmread(stream))


def _tar_member_lines(archive: tarfile.TarFile, member: str) -> list[str]:
    handle = archive.extractfile(member)
    if handle is None:
        raise ValueError(f"unable to extract {member}")
    return _read_gzip_lines(handle)


def _tar_member_matrix(archive: tarfile.TarFile, member: str) -> sparse.csr_matrix:
    handle = archive.extractfile(member)
    if handle is None:
        raise ValueError(f"unable to extract {member}")
    return _read_gzip_matrix(handle)


def _find_member(names: Iterable[str], *, sample: str, token: str) -> str:
    matches = [name for name in names if sample in name and token in name]
    if len(matches) != 1:
        raise ValueError(f"expected one {sample}/{token} member, found {matches}")
    return matches[0]


def _first_gene_index(genes: list[str]) -> dict[str, int]:
    index: dict[str, int] = {}
    for position, gene in enumerate(genes):
        index.setdefault(gene, position)
    return index


def broad_lineage_labels(
    counts: sparse.csr_matrix,
    genes: list[str],
    *,
    min_score: float = 0.08,
    min_margin: float = 0.02,
    excluded_marker: str | None = None,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Assign broad PBMC lineages using condition-blind, non-response marker modules.

    The classifier abstains when marker evidence is weak or tied. It intentionally avoids
    activation markers and does not claim subtypes such as classical versus non-classical
    monocytes.
    """

    if counts.shape[0] != len(genes):
        raise ValueError("gene list length does not match matrix rows")
    library = np.asarray(counts.sum(axis=0)).ravel().astype(float)
    gene_index = _first_gene_index(genes)
    labels = list(BROAD_MARKERS)
    scores: list[np.ndarray] = []
    for label in labels:
        marker_genes = [gene for gene in BROAD_MARKERS[label] if gene != excluded_marker]
        indices = [gene_index[gene] for gene in marker_genes if gene in gene_index]
        if not indices:
            scores.append(np.zeros(counts.shape[1], dtype=float))
            continue
        selected = counts[indices, :].astype(float)
        scaled = selected.multiply(10_000.0 / np.maximum(library, 1.0))
        score = np.asarray(np.log1p(scaled.toarray()).mean(axis=0)).ravel()
        scores.append(score)
    score_matrix = np.vstack(scores)
    order = np.argsort(score_matrix, axis=0)
    top_index = order[-1, :]
    second_index = order[-2, :]
    column_index = np.arange(counts.shape[1])
    top_score = score_matrix[top_index, column_index]
    margin = top_score - score_matrix[second_index, column_index]
    assigned = np.asarray(labels, dtype=object)[top_index]
    assigned[(top_score < min_score) | (margin < min_margin)] = "Unknown"
    score_frame = pd.DataFrame(score_matrix.T, columns=labels)
    score_frame["top_score"] = top_score
    score_frame["margin"] = margin
    return assigned, score_frame


def assign_hash_donors(
    hto_counts: sparse.csr_matrix,
    *,
    min_top_count: int = 20,
    min_top_to_second_ratio: float = 3.0,
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    """Conservatively demultiplex three hashtag libraries without using RNA expression."""

    dense = hto_counts.toarray().T.astype(float)
    if dense.shape[1] < 2:
        raise ValueError("HTO demultiplexing requires at least two hashtag features")
    sorted_counts = np.sort(dense, axis=1)
    top = sorted_counts[:, -1]
    second = sorted_counts[:, -2]
    ratio = (top + 1.0) / (second + 1.0)
    singlet = (top >= min_top_count) & (ratio >= min_top_to_second_ratio)
    donor = np.asarray([f"HTO{index + 1}" for index in dense.argmax(axis=1)], dtype=object)
    diagnostics = pd.DataFrame(
        {
            "top_hto_count": top,
            "second_hto_count": second,
            "top_to_second_ratio": ratio,
            "hto_singlet": singlet,
            "donor": donor,
        }
    )
    return donor, singlet, diagnostics


def _long_cell_table(
    *,
    matrix: sparse.csr_matrix,
    genes: list[str],
    targets: list[str],
    metadata: pd.DataFrame,
    study: str,
    stimulus: str,
    adt_matrix: sparse.csr_matrix | None = None,
    adt_features: list[str] | None = None,
) -> pd.DataFrame:
    """Expand selected sparse gene rows to the tool's standard long input contract."""

    if matrix.shape[1] != len(metadata):
        raise ValueError("metadata rows do not match matrix columns")
    gene_index = _first_gene_index(genes)
    available_targets = [target for target in targets if target in gene_index]
    if not available_targets:
        raise ValueError("none of the requested targets is present")
    target_indices = [gene_index[target] for target in available_targets]
    rna = matrix[target_indices, :].T.toarray().astype(float)
    library = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
    n_cells = len(metadata)
    n_targets = len(available_targets)
    output = pd.DataFrame(
        {
            "study": np.repeat(study, n_cells * n_targets),
            "donor": np.repeat(metadata["donor"].astype(str).to_numpy(), n_targets),
            "cell_id": np.repeat(metadata["cell_id"].astype(str).to_numpy(), n_targets),
            "cell_type": np.repeat(metadata["cell_type"].astype(str).to_numpy(), n_targets),
            "condition": np.repeat(metadata["condition"].astype(str).to_numpy(), n_targets),
            "stimulus": np.repeat(stimulus, n_cells * n_targets),
            "target": np.tile(np.asarray(available_targets, dtype=object), n_cells),
            "rna_count": rna.ravel(order="C"),
            "library_size": np.repeat(library, n_targets),
        }
    )
    output["adt_count"] = np.nan
    output["adt_library_size"] = np.nan
    if adt_matrix is not None and adt_features is not None:
        adt_library = np.asarray(adt_matrix.sum(axis=0)).ravel().astype(float)
        feature_index = _first_gene_index(adt_features)
        for feature, target in ADT_TARGETS.items():
            if feature not in feature_index or target not in available_targets:
                continue
            cell_values = np.asarray(adt_matrix[feature_index[feature], :].toarray()).ravel()
            target_mask = output["target"].eq(target).to_numpy()
            output.loc[target_mask, "adt_count"] = cell_values
            output.loc[target_mask, "adt_library_size"] = adt_library
    return output


def build_gse157857(
    raw_tar: Path,
    *,
    targets: list[str],
) -> RealBuildResult:
    """Build the paired IFN-beta myeloid CITE-seq cohort with HTO donor recovery."""

    sample_map = {
        "resting": {"RNA": "GSM4776909", "ADT": "GSM4776907", "HTO": "GSM4776908"},
        "activated": {"RNA": "GSM4776912", "ADT": "GSM4776910", "HTO": "GSM4776911"},
    }
    frames: list[pd.DataFrame] = []
    condition_qc: dict[str, object] = {}
    with tarfile.open(raw_tar) as archive:
        names = archive.getnames()
        for condition, samples in sample_map.items():
            rna_matrix = _tar_member_matrix(
                archive,
                _find_member(names, sample=samples["RNA"], token="matrix.mtx.gz"),
            )
            rna_features = _tar_member_lines(
                archive,
                _find_member(names, sample=samples["RNA"], token="features.tsv.gz"),
            )
            genes = [entry.split("\t")[1] for entry in rna_features]
            rna_barcodes = _tar_member_lines(
                archive,
                _find_member(names, sample=samples["RNA"], token="barcodes.tsv.gz"),
            )
            adt_matrix = _tar_member_matrix(
                archive,
                _find_member(names, sample=samples["ADT"], token="matrix.mtx.gz"),
            )
            adt_features = _tar_member_lines(
                archive,
                _find_member(names, sample=samples["ADT"], token="features.tsv.gz"),
            )
            hto_matrix = _tar_member_matrix(
                archive,
                _find_member(names, sample=samples["HTO"], token="matrix.mtx.gz"),
            )
            labels, annotation_scores = broad_lineage_labels(rna_matrix, genes)
            donor, singlet, hto_qc = assign_hash_donors(hto_matrix)
            library = np.asarray(rna_matrix.sum(axis=0)).ravel()
            detected = np.asarray((rna_matrix > 0).sum(axis=0)).ravel()
            keep = singlet & (library >= 500) & (detected >= 200) & (labels == "Myeloid")
            selected = np.flatnonzero(keep)
            metadata = pd.DataFrame(
                {
                    "donor": donor[selected],
                    "cell_id": [
                        f"GSE157857_{condition}_{rna_barcodes[index]}" for index in selected
                    ],
                    "cell_type": "Myeloid",
                    "condition": condition,
                }
            )
            expanded_targets = sorted(set(targets) | set(ADT_TARGETS.values()))
            frames.append(
                _long_cell_table(
                    matrix=rna_matrix[:, selected],
                    genes=genes,
                    targets=expanded_targets,
                    metadata=metadata,
                    study="GSE157857",
                    stimulus="IFN_I_18h",
                    adt_matrix=adt_matrix[:, selected],
                    adt_features=adt_features,
                )
            )
            condition_qc[condition] = {
                "n_filtered_matrix_cells": int(rna_matrix.shape[1]),
                "n_hto_singlets": int(singlet.sum()),
                "n_primary_myeloid_cells": int(keep.sum()),
                "donor_counts": metadata["donor"].value_counts().sort_index().to_dict(),
                "broad_annotation_counts": pd.Series(labels).value_counts().to_dict(),
                "median_annotation_margin_primary": float(
                    annotation_scores.loc[selected, "margin"].median()
                ),
                "median_top_hto_count_primary": float(
                    hto_qc.loc[selected, "top_hto_count"].median()
                ),
                "median_hto_ratio_primary": float(
                    hto_qc.loc[selected, "top_to_second_ratio"].median()
                ),
            }
    cells = pd.concat(frames, ignore_index=True)
    return RealBuildResult(
        cells=cells,
        qc={
            "accession": "GSE157857",
            "primary_population": "broad Myeloid only",
            "stimulus": "IFN-beta, 18 h",
            "conditions": condition_qc,
            "protein_targets": sorted(ADT_TARGETS.values()),
        },
    )


def build_gse178429(
    counts_path: Path,
    genes_path: Path,
    metadata_path: Path,
    *,
    targets: list[str],
) -> RealBuildResult:
    """Build paired IFN-gamma and LPS PBMC contrasts, excluding GolgiPlug and PMA."""

    metadata = pd.read_csv(metadata_path, sep="\t")
    genes = pd.read_csv(genes_path, header=None)[0].astype(str).tolist()
    with gzip.open(counts_path, "rb") as stream:
        matrix = sparse.csr_matrix(mmread(stream))
    if matrix.shape[1] != len(metadata):
        raise ValueError("GSE178429 metadata does not align with matrix columns")
    labels, scores = broad_lineage_labels(matrix, genes)
    metadata = metadata.copy()
    metadata["cell_type"] = labels
    metadata["cell_id"] = metadata["cellBarcode"].astype(str)
    metadata["donor"] = metadata["Donor"].astype(str)
    quality_mask = metadata["nCount_RNA"].ge(500) & metadata["nFeature_RNA"].ge(200)
    supported_lineages = ["T_cell", "NK_cell", "B_cell", "Myeloid"]
    qc_mask = quality_mask & metadata["cell_type"].isin(supported_lineages)
    all_markers = {gene for markers in BROAD_MARKERS.values() for gene in markers}
    leakage_targets = sorted(set(targets) & all_markers)
    non_leakage_targets = sorted(set(targets) - set(leakage_targets))
    leave_one_target_out_labels: dict[str, np.ndarray] = {}
    leave_one_target_out_scores: dict[str, pd.DataFrame] = {}
    for target in leakage_targets:
        target_labels, target_scores = broad_lineage_labels(
            matrix,
            genes,
            excluded_marker=target,
        )
        leave_one_target_out_labels[target] = target_labels
        leave_one_target_out_scores[target] = target_scores

    frames: list[pd.DataFrame] = []
    contrast_qc: dict[str, object] = {}
    for active_label, control_label, stimulus, study in GSE178429_CONTRASTS:
        condition_mask = metadata["Condition"].isin([active_label, control_label])
        selected_mask = qc_mask & condition_mask
        selected = np.flatnonzero(selected_mask.to_numpy())
        contrast_meta = metadata.loc[selected].copy()
        contrast_meta["condition"] = np.where(
            contrast_meta["Condition"].eq(active_label), "activated", "resting"
        )
        contrast_meta["cell_id"] = study + "_" + contrast_meta["cell_id"].astype(str)
        if non_leakage_targets:
            frames.append(
                _long_cell_table(
                    matrix=matrix[:, selected],
                    genes=genes,
                    targets=non_leakage_targets,
                    metadata=contrast_meta[["donor", "cell_id", "cell_type", "condition"]],
                    study=study,
                    stimulus=stimulus,
                )
            )
        for target in leakage_targets:
            target_metadata = metadata.copy()
            target_metadata["cell_type"] = leave_one_target_out_labels[target]
            target_mask = (
                quality_mask
                & target_metadata["cell_type"].isin(supported_lineages)
                & condition_mask
            )
            target_selected = np.flatnonzero(target_mask.to_numpy())
            target_meta = target_metadata.loc[target_selected].copy()
            target_meta["condition"] = np.where(
                target_meta["Condition"].eq(active_label), "activated", "resting"
            )
            target_meta["cell_id"] = study + "_" + target_meta["cell_id"].astype(str)
            frames.append(
                _long_cell_table(
                    matrix=matrix[:, target_selected],
                    genes=genes,
                    targets=[target],
                    metadata=target_meta[["donor", "cell_id", "cell_type", "condition"]],
                    study=study,
                    stimulus=stimulus,
                )
            )
        counts = (
            contrast_meta.groupby(["donor", "cell_type", "condition"], observed=True)
            .size()
            .rename("n_cells")
            .reset_index()
        )
        paired = counts.pivot_table(
            index=["donor", "cell_type"], columns="condition", values="n_cells"
        ).dropna()
        contrast_qc[study] = {
            "active_label": active_label,
            "control_label": control_label,
            "n_cells": int(len(contrast_meta)),
            "n_complete_donor_celltype_pairs": int(len(paired)),
            "cell_type_counts": contrast_meta["cell_type"].value_counts().to_dict(),
        }
    selected_primary = qc_mask & metadata["Condition"].isin(
        [entry for pair in GSE178429_CONTRASTS for entry in pair[:2]]
    )
    return RealBuildResult(
        cells=pd.concat(frames, ignore_index=True),
        qc={
            "accession": "GSE178429",
            "included_stimuli": ["IFN-gamma", "LPS"],
            "excluded_primary_conditions": [
                "ControlGolgiPlug_6h",
                "IFNGolgiPlug_6h",
                "LPSGolgiPlug_6h",
                "PMAGolgiPlug_6h",
                "PMA_1h",
                "PMA_6h",
            ],
            "broad_annotation_counts_all_cells": metadata["cell_type"].value_counts().to_dict(),
            "median_annotation_margin_primary": float(
                scores.loc[selected_primary.to_numpy(), "margin"].median()
            ),
            "annotation_leakage_control": {
                "method": "leave-target-out marker modules",
                "targets": leakage_targets,
                "median_margins": {
                    target: float(
                        leave_one_target_out_scores[target]
                        .loc[selected_primary.to_numpy(), "margin"]
                        .median()
                    )
                    for target in leakage_targets
                },
            },
            "contrasts": contrast_qc,
        },
    )


def _gse96583_condition_metadata(
    metadata: pd.DataFrame,
    barcodes: list[str],
    *,
    deposited_condition: str,
) -> tuple[pd.DataFrame, int]:
    """Align one GSE96583 library to deposited metadata without barcode guessing.

    The control and stimulated 10x libraries share 313 raw barcode sequences. The deposited
    metadata preserves uniqueness by appending one extra ``1`` to the stimulated duplicate
    row names (``-1`` becomes ``-11``). This adapter reverses that documented serialization
    artefact only inside the stimulated condition and then requires exact one-to-one coverage.
    """

    selected = metadata[metadata["stim"].eq(deposited_condition)].copy()
    selected["raw_barcode"] = selected.index.astype(str)
    duplicate_suffix_rows = 0
    if deposited_condition == "stim":
        duplicate_mask = selected["raw_barcode"].str.endswith("-11")
        duplicate_suffix_rows = int(duplicate_mask.sum())
        selected.loc[duplicate_mask, "raw_barcode"] = selected.loc[
            duplicate_mask, "raw_barcode"
        ].str.replace(r"-11$", "-1", regex=True)
    if selected["raw_barcode"].duplicated().any():
        raise ValueError(f"GSE96583 {deposited_condition} metadata barcodes are not unique")
    expected = set(barcodes)
    observed = set(selected["raw_barcode"].astype(str))
    if expected != observed:
        missing = sorted(expected - observed)[:5]
        extra = sorted(observed - expected)[:5]
        raise ValueError(
            "GSE96583 metadata/barcode mismatch: "
            f"missing={missing}, extra={extra}, condition={deposited_condition}"
        )
    aligned = selected.set_index("raw_barcode").loc[barcodes].copy()
    if not aligned["stim"].eq(deposited_condition).all():
        raise ValueError("GSE96583 condition labels do not match library identity")
    return aligned, duplicate_suffix_rows


def build_gse96583(
    raw_tar: Path,
    genes_path: Path,
    metadata_path: Path,
    *,
    targets: list[str],
) -> RealBuildResult:
    """Build the unopened eight-donor IFN-beta replication under frozen v0.3 rules.

    Sample identity defines condition; deposited demuxlet metadata supplies donor and singlet
    status. Broad lineage modules are condition-blind and benchmark marker targets use the same
    leave-target-out annotation guardrail as GSE178429.
    """

    with gzip.open(genes_path, "rt", encoding="utf-8") as stream:
        fields = [line.rstrip("\n").split("\t") for line in stream]
    if not fields or any(len(entry) < 2 for entry in fields):
        raise ValueError("GSE96583 gene-order file must contain Ensembl ID and symbol columns")
    genes = [entry[1] for entry in fields]
    metadata_all = pd.read_csv(metadata_path, sep="\t", index_col=0)
    required_metadata = {"ind", "stim", "multiplets"}
    missing_metadata = required_metadata - set(metadata_all.columns)
    if missing_metadata:
        raise ValueError(f"GSE96583 metadata missing columns: {sorted(missing_metadata)}")

    sample_map = {
        "resting": ("GSM2560248", "ctrl"),
        "activated": ("GSM2560249", "stim"),
    }
    supported_lineages = ["T_cell", "NK_cell", "B_cell", "Myeloid"]
    all_markers = {gene for markers in BROAD_MARKERS.values() for gene in markers}
    leakage_targets = sorted(set(targets) & all_markers)
    non_leakage_targets = sorted(set(targets) - set(leakage_targets))
    frames: list[pd.DataFrame] = []
    condition_qc: dict[str, object] = {}
    cell_metadata_frames: list[pd.DataFrame] = []

    with tarfile.open(raw_tar) as archive:
        names = archive.getnames()
        for condition, (sample, deposited_condition) in sample_map.items():
            matrix = _tar_member_matrix(
                archive,
                _find_member(names, sample=sample, token="mtx.gz"),
            )
            if matrix.shape[0] != len(genes):
                raise ValueError("GSE96583 gene-order length does not match matrix rows")
            barcodes = _tar_member_lines(
                archive,
                _find_member(names, sample=sample, token="barcodes.tsv.gz"),
            )
            if matrix.shape[1] != len(barcodes):
                raise ValueError("GSE96583 barcode length does not match matrix columns")
            metadata, repaired_duplicates = _gse96583_condition_metadata(
                metadata_all,
                barcodes,
                deposited_condition=deposited_condition,
            )
            labels, annotation_scores = broad_lineage_labels(matrix, genes)
            library = np.asarray(matrix.sum(axis=0)).ravel()
            detected = np.asarray((matrix > 0).sum(axis=0)).ravel()
            quality_mask = (
                metadata["multiplets"].eq("singlet").to_numpy()
                & (library >= 500)
                & (detected >= 200)
            )
            base_keep = quality_mask & np.isin(labels, supported_lineages)
            selected = np.flatnonzero(base_keep)
            base_metadata = pd.DataFrame(
                {
                    "donor": metadata.iloc[selected]["ind"].astype(str).to_numpy(),
                    "cell_id": [f"GSE96583_{condition}_{barcodes[index]}" for index in selected],
                    "cell_type": labels[selected],
                    "condition": condition,
                }
            )
            cell_metadata_frames.append(base_metadata)
            if non_leakage_targets:
                frames.append(
                    _long_cell_table(
                        matrix=matrix[:, selected],
                        genes=genes,
                        targets=non_leakage_targets,
                        metadata=base_metadata,
                        study="GSE96583",
                        stimulus="IFN_I_6h_confirmation",
                    )
                )

            leave_target_out_margins: dict[str, float] = {}
            for target in leakage_targets:
                target_labels, target_scores = broad_lineage_labels(
                    matrix,
                    genes,
                    excluded_marker=target,
                )
                target_keep = quality_mask & np.isin(target_labels, supported_lineages)
                target_selected = np.flatnonzero(target_keep)
                target_metadata = pd.DataFrame(
                    {
                        "donor": metadata.iloc[target_selected]["ind"].astype(str).to_numpy(),
                        "cell_id": [
                            f"GSE96583_{condition}_{barcodes[index]}" for index in target_selected
                        ],
                        "cell_type": target_labels[target_selected],
                        "condition": condition,
                    }
                )
                frames.append(
                    _long_cell_table(
                        matrix=matrix[:, target_selected],
                        genes=genes,
                        targets=[target],
                        metadata=target_metadata,
                        study="GSE96583",
                        stimulus="IFN_I_6h_confirmation",
                    )
                )
                leave_target_out_margins[target] = float(
                    target_scores.loc[target_selected, "margin"].median()
                )

            condition_qc[condition] = {
                "sample": sample,
                "deposited_condition": deposited_condition,
                "n_matrix_cells": int(matrix.shape[1]),
                "n_metadata_rows_aligned": int(len(metadata)),
                "n_serialized_duplicate_barcodes_repaired": repaired_duplicates,
                "n_demuxlet_singlets": int(metadata["multiplets"].eq("singlet").sum()),
                "n_primary_cells": int(base_keep.sum()),
                "donor_counts": base_metadata["donor"].value_counts().sort_index().to_dict(),
                "broad_annotation_counts": pd.Series(labels).value_counts().to_dict(),
                "median_annotation_margin_primary": float(
                    annotation_scores.loc[selected, "margin"].median()
                ),
                "leave_target_out_median_margins": leave_target_out_margins,
            }

    cells = pd.concat(frames, ignore_index=True)
    cell_metadata = pd.concat(cell_metadata_frames, ignore_index=True)
    counts = (
        cell_metadata.groupby(["donor", "cell_type", "condition"], observed=True)
        .size()
        .rename("n_cells")
        .reset_index()
    )
    complete_pairs = counts.pivot_table(
        index=["donor", "cell_type"],
        columns="condition",
        values="n_cells",
    ).dropna(subset=["resting", "activated"])
    donors = sorted(cell_metadata["donor"].astype(str).unique())
    return RealBuildResult(
        cells=cells,
        qc={
            "accession": "GSE96583",
            "role": "frozen independent IFN-I replication",
            "disease_context": "SLE patient PBMCs",
            "stimulus": "IFN-beta, 6 h",
            "n_donors": len(donors),
            "donors": donors,
            "conditions": condition_qc,
            "n_complete_donor_celltype_pairs": int(len(complete_pairs)),
            "annotation_leakage_control": {
                "method": "leave-target-out marker modules",
                "targets": leakage_targets,
            },
            "condition_source": "GSM library identity; metadata alignment verified exactly",
        },
    )


def ensembl_symbol_map_from_gse157857(raw_tar: Path) -> dict[str, str]:
    """Reuse the deposited GRCh38 feature table as a local Ensembl-to-symbol mapping."""

    with tarfile.open(raw_tar) as archive:
        names = archive.getnames()
        feature_member = _find_member(names, sample="GSM4776909", token="features.tsv.gz")
        features = _tar_member_lines(archive, feature_member)
    result: dict[str, str] = {}
    for entry in features:
        fields = entry.split("\t")
        if len(fields) >= 2:
            result.setdefault(fields[0].split(".")[0], fields[1])
    return result


def build_gse140244(
    counts_path: Path,
    metadata_path: Path,
    *,
    targets: list[str],
    ensembl_to_symbol: dict[str, str],
    time_points: tuple[int, ...] = (2, 4, 8, 12, 24, 48, 72),
) -> RealBuildResult:
    """Build donor-paired CD4-memory activation contrasts against each donor's baseline."""

    metadata = pd.read_csv(metadata_path, sep="\t")
    counts = pd.read_csv(counts_path, sep="\t")
    sample_columns = [column for column in counts.columns if column != "GENE_ID"]
    if sample_columns != metadata["ExpressionMatrix_SampleID"].astype(str).tolist():
        raise ValueError("GSE140244 metadata order does not match count columns")
    symbols = counts["GENE_ID"].astype(str).str.split(".").str[0].map(ensembl_to_symbol)
    counts = counts.assign(symbol=symbols)
    selected_counts = counts[counts["symbol"].isin(targets)].copy()
    selected_counts = selected_counts.groupby("symbol", observed=True)[sample_columns].sum()
    library_sizes = counts[sample_columns].sum(axis=0)
    frames: list[pd.DataFrame] = []
    contrast_qc: dict[str, object] = {}
    for time_point in time_points:
        selected_meta = metadata[metadata["Time_point"].isin([0, time_point])].copy()
        donors_by_time = selected_meta.groupby("Donor_ID")["Time_point"].agg(set)
        expected_times = {0, time_point}
        complete_donors = donors_by_time[donors_by_time.map(expected_times.__eq__)].index
        selected_meta = selected_meta[selected_meta["Donor_ID"].isin(complete_donors)].copy()
        selected_meta["condition"] = np.where(
            selected_meta["Time_point"].eq(0), "resting", "activated"
        )
        study = f"GSE140244_{time_point}h"
        stimulus = f"lymphocyte_activation_{time_point}h"
        rows: list[dict[str, object]] = []
        for sample in selected_meta.itertuples(index=False):
            sample_id = str(sample.ExpressionMatrix_SampleID)
            for target in targets:
                count = (
                    float(selected_counts.at[target, sample_id])
                    if target in selected_counts.index
                    else 0.0
                )
                rows.append(
                    {
                        "study": study,
                        "donor": str(sample.Donor_ID),
                        "cell_id": sample_id,
                        "cell_type": "T_cell",
                        "condition": str(sample.condition),
                        "stimulus": stimulus,
                        "target": target,
                        "rna_count": count,
                        "library_size": float(library_sizes[sample_id]),
                        "adt_count": np.nan,
                        "adt_library_size": np.nan,
                    }
                )
        frame = pd.DataFrame.from_records(rows)
        frames.append(frame)
        contrast_qc[study] = {
            "time_point_hours": time_point,
            "n_complete_donors": int(len(complete_donors)),
            "n_libraries": int(len(selected_meta)),
            "replicate_counts": selected_meta["Replicate"].value_counts().to_dict(),
        }
    return RealBuildResult(
        cells=pd.concat(frames, ignore_index=True),
        qc={
            "accession": "GSE140244",
            "cell_type": "sorted CD4 memory T cells",
            "baseline": "0 h",
            "contrasts": contrast_qc,
            "mapping_source": "GSE157857 deposited GRCh38 feature table",
        },
    )


def combine_real_cohorts(results: Iterable[RealBuildResult]) -> RealBuildResult:
    """Combine already-standardised cohorts without cross-study normalization."""

    materialised = list(results)
    if not materialised:
        raise ValueError("at least one cohort is required")
    return RealBuildResult(
        cells=pd.concat([result.cells for result in materialised], ignore_index=True),
        qc={"cohorts": [result.qc for result in materialised]},
    )
