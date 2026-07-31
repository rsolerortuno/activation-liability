"""End-to-end activation-liability audit."""

from __future__ import annotations

from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd

from activation_liability.claims import build_claims
from activation_liability.evidence import classify_evidence, load_evidence_rules
from activation_liability.io import versions
from activation_liability.metrics import (
    add_protein_concordance,
    footprint_by_target,
    is_relevant_stimulus,
    positive_fractions,
    protein_meta,
    selectivity_erosion,
)
from activation_liability.statistics import meta_analyse, paired_effects

EVIDENCE_WEIGHT = {"A": 3.0, "B": 2.0, "C": 1.0, "INSUFFICIENT": 0.0}


def validate_input(cells: pd.DataFrame) -> None:
    """Reject malformed or scientifically unsafe generic inputs."""

    required = {
        "study",
        "donor",
        "cell_id",
        "cell_type",
        "condition",
        "stimulus",
        "target",
        "rna_count",
        "library_size",
    }
    missing = required - set(cells.columns)
    if missing:
        raise ValueError(f"missing input columns: {sorted(missing)}")
    conditions = set(cells["condition"].dropna().unique())
    if not conditions.issubset({"resting", "activated"}) or not conditions:
        raise ValueError("condition must contain only resting and activated")
    if (cells["rna_count"] < 0).any() or (cells["library_size"] <= 0).any():
        raise ValueError("counts must be non-negative and library sizes positive")


def _target_summary(
    rows: pd.DataFrame,
    cells: pd.DataFrame,
    fraction_cutoff: float,
) -> pd.DataFrame:
    footprints = footprint_by_target(rows, positive_fraction_cutoff=fraction_cutoff)
    erosion = selectivity_erosion(cells)
    ordered = rows.copy()
    ordered["class_weight"] = ordered["evidence_class"].map(EVIDENCE_WEIGHT).fillna(0.0)
    ordered["protein_adjustment"] = np.select(
        [
            ordered["protein_concordance"] == "CONCORDANT",
            ordered["protein_concordance"] == "DISCORDANT",
        ],
        [0.35, -0.75],
        default=0.0,
    )
    ordered["row_score"] = (
        ordered["inducibility_lfc"].clip(lower=-4.0, upper=6.0)
        + 0.40 * ordered["class_weight"]
        + ordered["protein_adjustment"]
    )
    best_index = ordered.groupby("target", observed=True)["row_score"].idxmax()
    summary = ordered.loc[
        best_index,
        [
            "target",
            "cell_type",
            "stimulus",
            "inducibility_lfc",
            "ci_low",
            "ci_high",
            "i2",
            "positive_fraction_resting",
            "positive_fraction_activated",
            "protein_concordance",
            "evidence_class",
            "abstention_reasons",
            "protein_adjustment",
            "row_score",
        ],
    ].rename(columns={"cell_type": "driver_cell_type", "stimulus": "driver_stimulus"})
    summary = summary.merge(footprints, on="target", how="left").merge(
        erosion, on="target", how="left"
    )
    summary["rna_only_score"] = summary["inducibility_lfc"] + 0.55 * summary["footprint_expansion"]
    summary["score"] = summary["row_score"] + 0.55 * summary["footprint_expansion"]
    return summary.sort_values(["score", "target"], ascending=[False, True]).reset_index(drop=True)


def run_audit(
    cells: pd.DataFrame,
    *,
    rules_path: Path,
    detection_count_threshold: int = 0,
    positive_fraction_cutoff: float = 0.10,
    protein_lfc_cutoff: float = 0.50,
) -> dict[str, Any]:
    """Run paired pseudobulk statistics, metrics, evidence classes and claims."""

    validate_input(cells)
    within = paired_effects(cells)
    if within.empty:
        raise ValueError("no target-cell-stimulus group has at least three complete donor pairs")
    meta = meta_analyse(within).rename(
        columns={
            "estimate": "inducibility_lfc",
            "standard_error": "inducibility_se",
        }
    )
    fractions = positive_fractions(cells, count_threshold=detection_count_threshold)
    rows = meta.merge(fractions, on=["target", "cell_type", "stimulus"], how="left")
    protein = protein_meta(cells)
    rows = rows.merge(protein, on=["target", "cell_type", "stimulus"], how="left")
    rows = add_protein_concordance(rows, lfc_cutoff=protein_lfc_cutoff)
    rows["relevant_stimulus"] = [
        is_relevant_stimulus(str(cell_type), str(stimulus))
        for cell_type, stimulus in zip(rows["cell_type"], rows["stimulus"], strict=True)
    ]
    rules = load_evidence_rules(rules_path)
    classes: list[str] = []
    reasons: list[list[str]] = []
    for record in rows.to_dict(orient="records"):
        evidence_class, abstention_reasons = classify_evidence(record, rules)
        classes.append(evidence_class)
        reasons.append(abstention_reasons)
    rows["evidence_class"] = classes
    rows["abstention_reasons"] = reasons
    rows = rows.sort_values(["target", "cell_type", "stimulus"]).reset_index(drop=True)
    summary = _target_summary(rows, cells, positive_fraction_cutoff)
    claims = build_claims(rows, summary)
    return {
        "schema_version": 1,
        "status": (
            "SYNTHETIC" if cells["study"].astype(str).str.startswith("SYNTH_").all() else "COMPUTED"
        ),
        "parameters": {
            "detection_count_threshold": detection_count_threshold,
            "positive_fraction_cutoff": positive_fraction_cutoff,
            "protein_lfc_cutoff": protein_lfc_cutoff,
            "meta_analysis": "DerSimonian-Laird random effects",
            "statistical_unit": "donor pseudobulk",
        },
        "versions": versions(),
        "audit_rows": _json_records(rows),
        "target_summary": _json_records(summary),
        "claims": claims,
    }


def _json_records(frame: pd.DataFrame) -> list[dict[str, Any]]:
    """Convert a frame to strict JSON-compatible records."""

    clean = frame.replace({np.nan: None, np.inf: None, -np.inf: None})
    return cast(list[dict[str, Any]], clean.to_dict(orient="records"))
