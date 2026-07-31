from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

from activation_liability.benchmark import load_controls
from activation_liability.tissue import (
    PSORIASIS_BASELINE_MAP,
    build_gse228421_tissue_directory,
    condition_blind_lineage_labels,
    negative_binomial_sensitivity,
    paired_tissue_effects,
    ranking_metrics,
    score_tissue_controls,
    tissue_footprint,
)


def test_tissue_annotation_markers_do_not_overlap_controls(root: Path) -> None:
    import yaml

    config = yaml.safe_load((root / "config/tissue_extension.yaml").read_text())
    controls = load_controls(root / "data/controls/controls.yaml")
    markers = {str(gene) for values in config["annotation"]["markers"].values() for gene in values}
    assert markers.isdisjoint(set(controls["target"].astype(str)))


def test_condition_blind_lineage_labels_uses_multi_marker_gate() -> None:
    genes = ["CD3D", "TRAC", "MZB1", "JCHAIN", "KRT1", "KRT10"]
    counts = sparse.csr_matrix(
        np.asarray(
            [
                [50, 0, 0],
                [40, 0, 0],
                [0, 60, 0],
                [0, 0, 0],
                [0, 0, 70],
                [0, 0, 60],
            ]
        )
    )
    labels, _ = condition_blind_lineage_labels(
        counts,
        genes,
        {
            "T_cell": ["CD3D", "TRAC"],
            "Plasma_cell": ["MZB1", "JCHAIN"],
            "Epithelial": ["KRT1", "KRT10"],
        },
        min_score=0.0,
        min_margin=0.0,
        minimum_positive_markers={"Plasma_cell": 2},
    )
    assert labels.tolist() == ["T_cell", "Unknown", "Epithelial"]


def _pb_row(
    donor: str,
    condition: str,
    *,
    target: str = "X",
    cell_type: str = "T_cell",
    count: int,
) -> dict[str, object]:
    return {
        "study": "S1",
        "stimulus": "activation",
        "tissue": "blood",
        "donor": donor,
        "condition": condition,
        "cell_type": cell_type,
        "target": target,
        "value_sum": count,
        "library_sum": 10_000,
        "positive_count": min(count, 100),
        "n_cells": 100,
    }


def test_paired_tissue_effects_recovers_known_paired_increase() -> None:
    rows: list[dict[str, object]] = []
    for donor, resting in zip(("D1", "D2", "D3", "D4"), (10, 20, 30, 40), strict=True):
        rows.append(_pb_row(donor, "resting", count=resting))
        rows.append(_pb_row(donor, "activated", count=resting * 4))
    result = paired_tissue_effects(pd.DataFrame(rows)).iloc[0]
    assert result["n_pairs"] == 4
    assert result["inducibility_lfc"] > 1.5
    assert result["ci_low"] > 0
    assert result["fraction_donors_positive"] == 1.0


def test_tissue_scoring_respects_predeclared_lineages_and_observability() -> None:
    controls = pd.DataFrame(
        [
            {"target": "POS", "label": "positive", "tier": 1, "y_true": 1},
            {"target": "NEG", "label": "negative", "tier": 1, "y_true": 0},
            {"target": "MISS", "label": "negative", "tier": 1, "y_true": 0},
        ]
    )
    audit = pd.DataFrame(
        [
            {
                "study": "S",
                "target": "POS",
                "cell_type": "T_cell",
                "stimulus": "A",
                "inducibility_lfc": 2.0,
                "inducibility_se": 0.25,
                "positive_fraction_resting": 0.1,
                "positive_fraction_activated": 0.6,
                "total_pairs": 5,
            },
            {
                "study": "S",
                "target": "POS",
                "cell_type": "Myeloid",
                "stimulus": "A",
                "inducibility_lfc": 20.0,
                "inducibility_se": 0.01,
                "positive_fraction_resting": 0.1,
                "positive_fraction_activated": 0.9,
                "total_pairs": 5,
            },
            {
                "study": "S",
                "target": "NEG",
                "cell_type": "B_cell",
                "stimulus": "A",
                "inducibility_lfc": 0.1,
                "inducibility_se": 0.1,
                "positive_fraction_resting": 0.5,
                "positive_fraction_activated": 0.5,
                "total_pairs": 5,
            },
            {
                "study": "S",
                "target": "MISS",
                "cell_type": "Epithelial",
                "stimulus": "A",
                "inducibility_lfc": 1.0,
                "inducibility_se": 0.1,
                "positive_fraction_resting": 0.0,
                "positive_fraction_activated": 0.01,
                "total_pairs": 5,
            },
        ]
    )
    config = {
        "score": {"z_value": 1.96, "observability_fraction_cutoff": 0.05},
        "target_lineages": {
            "POS": ["T_cell"],
            "NEG": ["B_cell"],
            "MISS": ["Epithelial"],
        },
    }
    scored = score_tissue_controls(audit, controls, config).set_index("target")
    assert scored.loc["POS", "driver_cell_type"] == "T_cell"
    assert scored.loc["POS", "score"] == 2.0 - 1.96 * 0.25
    assert scored.loc["MISS", "abstention_reason"] == (
        "TARGET_NOT_OBSERVABLE_AT_COVERAGE_THRESHOLD"
    )
    metrics = ranking_metrics(scored.reset_index())
    assert metrics["status"] == "COMPUTED"
    assert metrics["auroc"] == 1.0


def test_tissue_footprint_counts_unique_lineages_across_studies() -> None:
    audit = pd.DataFrame(
        [
            {
                "study": "S1",
                "target": "X",
                "cell_type": "T_cell",
                "positive_fraction_resting": 0.05,
                "positive_fraction_activated": 0.20,
            },
            {
                "study": "S2",
                "target": "X",
                "cell_type": "T_cell",
                "positive_fraction_resting": 0.01,
                "positive_fraction_activated": 0.30,
            },
            {
                "study": "S2",
                "target": "X",
                "cell_type": "Myeloid",
                "positive_fraction_resting": 0.20,
                "positive_fraction_activated": 0.25,
            },
        ]
    )
    result = tissue_footprint(audit, cutoff=0.10).iloc[0]
    assert result["tissue_footprint_resting"] == 1
    assert result["tissue_footprint_activated"] == 2
    assert result["tissue_footprint_expansion"] == 1


def test_negative_binomial_sensitivity_recovers_direction() -> None:
    rows: list[dict[str, object]] = []
    for index, donor in enumerate(("D1", "D2", "D3", "D4", "D5"), start=1):
        rows.append(_pb_row(donor, "resting", count=10 * index))
        rows.append(_pb_row(donor, "activated", count=30 * index))
    result = negative_binomial_sensitivity(pd.DataFrame(rows)).iloc[0]
    assert bool(result["nb_converged"])
    assert result["nb_lfc"] > 0


def _write_10x_pair(matrix_path: Path, feature_path: Path, activated: bool) -> None:
    genes = ["CD3D", "TRAC", "CD37", "CD22", "X"]
    matrix = np.asarray(
        [
            [50, 50, 50],
            [40, 40, 40],
            [0, 0, 0],
            [0, 0, 0],
            [30 if activated else 5] * 3,
        ],
        dtype=int,
    )
    with gzip.open(matrix_path, "wb") as stream:
        mmwrite(stream, sparse.csr_matrix(matrix))
    with gzip.open(feature_path, "wt", encoding="utf-8") as stream:
        for index, gene in enumerate(genes):
            stream.write(f"ENSG{index:05d}\t{gene}\tGene Expression\n")


def test_build_psoriasis_directory_validates_all_five_pairs(tmp_path: Path) -> None:
    for prefix, (_, condition) in PSORIASIS_BASELINE_MAP.items():
        _write_10x_pair(
            tmp_path / f"GSM000_{prefix}.matrix.mtx.gz",
            tmp_path / f"GSM000_{prefix}.features.tsv.gz",
            activated=condition == "activated",
        )
    config = {
        "endpoints": [
            {
                "study": "GSE228421",
                "primary_patients": ["P1", "P2", "P3", "P4", "P5"],
                "excluded_lineages": [],
            }
        ],
        "annotation": {
            "markers": {"T_cell": ["CD3D", "TRAC"], "B_cell": ["CD37", "CD22"]},
            "minimum_positive_markers": {},
        },
        "qc": {
            "minimum_library_size": 1,
            "minimum_detected_genes": 1,
            "maximum_mitochondrial_fraction": 1.0,
            "maximum_haemoglobin_fraction": 1.0,
            "annotation_min_score": 0.0,
            "annotation_min_margin": 0.0,
        },
    }
    result = build_gse228421_tissue_directory(tmp_path, targets=["X"], config=config)
    assert result.qc["included_donors"] == ["P1", "P2", "P3", "P4", "P5"]
    assert set(result.pseudobulk["condition"]) == {"resting", "activated"}
    assert set(result.pseudobulk["donor"]) == {"P1", "P2", "P3", "P4", "P5"}
    effects = paired_tissue_effects(result.pseudobulk)
    assert effects.iloc[0]["inducibility_lfc"] > 1.0
