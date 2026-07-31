from __future__ import annotations

import gzip
import io
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

from activation_liability.real_data import (
    _gse96583_condition_metadata,
    broad_lineage_labels,
    build_gse96583,
)


def test_leave_target_out_marker_removes_direct_annotation_feature() -> None:
    genes = ["CD3E", "TRAC", "MS4A1", "CD79B"]
    counts = sparse.csr_matrix(
        np.asarray(
            [
                [100, 0],
                [20, 0],
                [0, 100],
                [0, 20],
            ]
        )
    )
    labels, scores = broad_lineage_labels(counts, genes, min_score=0.0, min_margin=0.0)
    excluded_labels, excluded_scores = broad_lineage_labels(
        counts,
        genes,
        min_score=0.0,
        min_margin=0.0,
        excluded_marker="CD3E",
    )
    assert labels.tolist() == ["T_cell", "B_cell"]
    assert excluded_labels.tolist() == ["T_cell", "B_cell"]
    assert excluded_scores.loc[0, "T_cell"] < scores.loc[0, "T_cell"]


def test_gse96583_duplicate_barcode_serialization_is_condition_specific() -> None:
    metadata = pd.DataFrame(
        {
            "stim": ["ctrl", "stim"],
            "ind": ["D1", "D2"],
            "multiplets": ["singlet", "singlet"],
        },
        index=["AAAC-1", "AAAC-11"],
    )
    control, control_repairs = _gse96583_condition_metadata(
        metadata,
        ["AAAC-1"],
        deposited_condition="ctrl",
    )
    stimulated, stimulated_repairs = _gse96583_condition_metadata(
        metadata,
        ["AAAC-1"],
        deposited_condition="stim",
    )
    assert control.index.tolist() == ["AAAC-1"]
    assert stimulated.index.tolist() == ["AAAC-1"]
    assert control_repairs == 0
    assert stimulated_repairs == 1
    assert stimulated.iloc[0]["ind"] == "D2"


def _gzip_matrix(matrix: sparse.csr_matrix) -> bytes:
    buffer = io.BytesIO()
    mmwrite(buffer, matrix)
    return gzip.compress(buffer.getvalue(), mtime=0)


def _add_tar_bytes(archive: tarfile.TarFile, name: str, payload: bytes) -> None:
    info = tarfile.TarInfo(name)
    info.size = len(payload)
    info.mtime = 0
    archive.addfile(info, io.BytesIO(payload))


def test_build_gse96583_validates_gene_order_and_donor_pairs(tmp_path: Path) -> None:
    genes = ["CD3D", "CD3E", "TRAC", "CD247", "CD69"] + [f"GENE{index}" for index in range(205)]
    n_cells = 8
    base = np.full((len(genes), n_cells), 3, dtype=int)
    for marker in ("CD3D", "CD3E", "TRAC", "CD247"):
        base[genes.index(marker), :] = 60
    base[genes.index("CD69"), :] = 10
    control = sparse.csr_matrix(base)
    stimulated_array = base.copy()
    stimulated_array[genes.index("CD69"), :] = 80
    stimulated = sparse.csr_matrix(stimulated_array)
    barcodes = [f"BC{index:02d}-1" for index in range(n_cells)]

    raw_tar = tmp_path / "GSE96583_RAW.tar"
    with tarfile.open(raw_tar, "w") as archive:
        _add_tar_bytes(archive, "GSM2560248_2.1.mtx.gz", _gzip_matrix(control))
        _add_tar_bytes(
            archive,
            "GSM2560248_barcodes.tsv.gz",
            gzip.compress(("\n".join(barcodes) + "\n").encode(), mtime=0),
        )
        _add_tar_bytes(archive, "GSM2560249_2.2.mtx.gz", _gzip_matrix(stimulated))
        _add_tar_bytes(
            archive,
            "GSM2560249_barcodes.tsv.gz",
            gzip.compress(("\n".join(barcodes) + "\n").encode(), mtime=0),
        )

    genes_path = tmp_path / "GSE96583_batch2.genes.tsv.gz"
    with gzip.open(genes_path, "wt") as stream:
        for index, gene in enumerate(genes):
            stream.write(f"ENSG{index:011d}\t{gene}\n")

    donors = [f"D{index // 2 + 1}" for index in range(n_cells)]
    metadata = pd.DataFrame(
        {
            "ind": donors + donors,
            "stim": ["ctrl"] * n_cells + ["stim"] * n_cells,
            "multiplets": ["singlet"] * (2 * n_cells),
        },
        index=barcodes + [barcode + "1" for barcode in barcodes],
    )
    metadata_path = tmp_path / "GSE96583_batch2.total.tsne.df.tsv.gz"
    metadata.to_csv(metadata_path, sep="\t", compression="gzip")

    result = build_gse96583(
        raw_tar,
        genes_path,
        metadata_path,
        targets=["CD3E", "CD69"],
    )
    assert set(result.cells["condition"]) == {"resting", "activated"}
    assert set(result.cells["target"]) == {"CD3E", "CD69"}
    assert set(result.cells["cell_type"]) == {"T_cell"}
    assert result.qc["n_donors"] == 4
    conditions = result.qc["conditions"]
    assert isinstance(conditions, dict)
    assert conditions["resting"]["n_serialized_duplicate_barcodes_repaired"] == 0
    assert conditions["activated"]["n_serialized_duplicate_barcodes_repaired"] == n_cells
