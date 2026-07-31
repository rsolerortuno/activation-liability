"""Threshold sensitivity without redundant re-fitting."""

from __future__ import annotations

from itertools import product
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from activation_liability.audit import run_audit
from activation_liability.metrics import footprint_by_target, positive_fractions


def _rank_scores(
    base_rows: pd.DataFrame,
    fractions: pd.DataFrame,
    *,
    fraction_cutoff: float,
    lfc_cutoff: float,
) -> pd.DataFrame:
    rows = base_rows.drop(
        columns=["positive_fraction_resting", "positive_fraction_activated"]
    ).merge(fractions, on=["target", "cell_type", "stimulus"], how="left")
    footprints = footprint_by_target(rows, positive_fraction_cutoff=fraction_cutoff)
    rows["eligible_lfc"] = rows["inducibility_lfc"] >= lfc_cutoff
    rows["protein_bonus"] = np.select(
        [
            rows["protein_concordance"] == "CONCORDANT",
            rows["protein_concordance"] == "DISCORDANT",
        ],
        [0.35, -0.35],
        default=0.0,
    )
    rows["row_score"] = np.where(
        rows["eligible_lfc"],
        rows["inducibility_lfc"] - lfc_cutoff + rows["protein_bonus"],
        -2.0 + rows["inducibility_lfc"],
    )
    best = rows.groupby("target", observed=True)["row_score"].max().reset_index()
    best = best.merge(footprints[["target", "footprint_expansion"]], on="target", how="left")
    best["score"] = best["row_score"] + 0.55 * best["footprint_expansion"]
    return best.sort_values(["score", "target"], ascending=[False, True]).reset_index(drop=True)


def run_sensitivity(
    cells: pd.DataFrame,
    *,
    rules_path: Path,
    detection_thresholds: tuple[int, ...] = (0, 1, 2),
    fraction_cutoffs: tuple[float, ...] = (0.05, 0.10, 0.20),
    lfc_cutoffs: tuple[float, ...] = (0.50, 0.75, 1.00),
    base_audit: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Sweep threshold-dependent metrics while fitting effects only once."""

    base = base_audit or run_audit(cells, rules_path=rules_path)
    base_rows = pd.DataFrame(base["audit_rows"])
    all_detection_thresholds = set(detection_thresholds) | {0}
    fractions_by_detection = {
        threshold: positive_fractions(cells, count_threshold=threshold)
        for threshold in sorted(all_detection_thresholds)
    }
    default = _rank_scores(
        base_rows,
        fractions_by_detection[0],
        fraction_cutoff=0.10,
        lfc_cutoff=0.50,
    )
    default_rank = default.set_index("target")["score"].rank(ascending=False, method="average")
    default_top = set(default.head(10)["target"])
    settings: list[dict[str, Any]] = []
    for detection, fraction, lfc in product(
        detection_thresholds,
        fraction_cutoffs,
        lfc_cutoffs,
    ):
        ranked = _rank_scores(
            base_rows,
            fractions_by_detection[detection],
            fraction_cutoff=fraction,
            lfc_cutoff=lfc,
        )
        rank = ranked.set_index("target")["score"].rank(ascending=False, method="average")
        aligned = default_rank.index.intersection(rank.index)
        correlation = float(spearmanr(default_rank.loc[aligned], rank.loc[aligned]).statistic)
        top_overlap = len(default_top & set(ranked.head(10)["target"])) / max(len(default_top), 1)
        settings.append(
            {
                "detection_count_threshold": detection,
                "positive_fraction_cutoff": fraction,
                "lfc_cutoff": lfc,
                "spearman_vs_default": correlation,
                "top10_overlap_vs_default": top_overlap,
                "top10": ranked.head(10)["target"].tolist(),
            }
        )
    median_spearman = float(np.median([entry["spearman_vs_default"] for entry in settings]))
    return {
        "schema_version": 1,
        "default": {
            "detection_count_threshold": 0,
            "positive_fraction_cutoff": 0.10,
            "lfc_cutoff": 0.50,
        },
        "n_settings": len(settings),
        "median_spearman": median_spearman,
        "ranking_stability": "STABLE" if median_spearman >= 0.80 else "UNSTABLE",
        "settings": settings,
    }
