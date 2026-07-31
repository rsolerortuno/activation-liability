"""Synthetic/real benchmark with holdout and required invalid ablations."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import yaml
from sklearn.metrics import average_precision_score, balanced_accuracy_score, roc_auc_score

from activation_liability.audit import run_audit
from activation_liability.sensitivity import run_sensitivity
from activation_liability.statistics import invalid_cell_level_tests, paired_effects


def load_controls(path: Path) -> pd.DataFrame:
    """Load benchmark controls and validate basic independence fields."""

    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    controls = pd.DataFrame(payload["controls"])
    required = {"target", "label", "tier", "citation", "rationale"}
    if required - set(controls.columns):
        raise ValueError("control file is missing required fields")
    if controls["target"].duplicated().any():
        raise ValueError("control targets must be unique")
    controls["y_true"] = (controls["label"] == "positive").astype(int)
    return controls


def deterministic_holdout(controls: pd.DataFrame) -> tuple[set[str], set[str]]:
    """Hold out 30% of tier-1 positives before threshold selection."""

    positives = sorted(controls.loc[(controls["tier"] == 1) & (controls["y_true"] == 1), "target"])
    ranked = sorted(positives, key=lambda target: hashlib.sha256(target.encode()).hexdigest())
    n_holdout = max(1, int(np.ceil(0.30 * len(ranked))))
    holdout = set(ranked[:n_holdout])
    train = set(positives) - holdout
    return train, holdout


def classification_metrics(labels: np.ndarray, scores: np.ndarray) -> dict[str, float]:
    """Compute ranking metrics, retaining finite values only."""

    finite = np.isfinite(scores)
    labels = labels[finite]
    scores = scores[finite]
    if len(np.unique(labels)) < 2:
        return {
            "auroc": float("nan"),
            "average_precision": float("nan"),
            "precision_at_10": float("nan"),
        }
    order = np.argsort(scores)[::-1]
    k = min(10, len(order))
    return {
        "auroc": float(roc_auc_score(labels, scores)),
        "average_precision": float(average_precision_score(labels, scores)),
        "precision_at_10": float(labels[order[:k]].mean()),
    }


def _choose_cutoff(train: pd.DataFrame) -> float:
    unique = np.unique(train["score"].to_numpy())
    if len(unique) == 1:
        return float(unique[0])
    candidates = (unique[:-1] + unique[1:]) / 2.0
    best_cutoff = float(candidates[0])
    best_accuracy = -1.0
    for cutoff in candidates:
        predicted = (train["score"] >= cutoff).astype(int)
        accuracy = balanced_accuracy_score(train["y_true"], predicted)
        if accuracy > best_accuracy:
            best_accuracy = float(accuracy)
            best_cutoff = float(cutoff)
    return best_cutoff


def _coverage_accuracy(frame: pd.DataFrame, cutoff: float) -> list[dict[str, float]]:
    confidence = np.abs(frame["score"].to_numpy() - cutoff)
    predicted = (frame["score"].to_numpy() >= cutoff).astype(int)
    correct = predicted == frame["y_true"].to_numpy()
    rows: list[dict[str, float]] = []
    for coverage in (1.0, 0.75, 0.50, 0.25):
        keep_n = max(1, int(np.ceil(coverage * len(frame))))
        keep = np.argsort(confidence)[::-1][:keep_n]
        rows.append(
            {
                "coverage": keep_n / len(frame),
                "accuracy": float(correct[keep].mean()),
                "abstention": 1.0 - keep_n / len(frame),
            }
        )
    return rows


def invalid_cross_study_scores(cells: pd.DataFrame) -> pd.DataFrame:
    """Deliberately invalid absolute contrast across two different studies."""

    data = cells.copy()
    data["log_expression"] = np.log2(data["rna_count"] / data["library_size"] * 1_000_000.0 + 1.0)
    study_names = sorted(data["study"].unique())
    if len(study_names) < 2:
        raise ValueError("cross-study ablation requires at least two studies")
    activated = data[(data["study"] == study_names[0]) & (data["condition"] == "activated")]
    resting = data[(data["study"] == study_names[1]) & (data["condition"] == "resting")]
    active_mean = activated.groupby(["target", "cell_type"], observed=True)["log_expression"].mean()
    rest_mean = resting.groupby(["target", "cell_type"], observed=True)["log_expression"].mean()
    joined = (
        active_mean.rename("activated").to_frame().join(rest_mean.rename("resting"), how="inner")
    )
    joined["difference"] = joined["activated"] - joined["resting"]
    return (
        joined.reset_index()
        .groupby("target", observed=True)["difference"]
        .max()
        .reset_index(name="score")
    )


def run_benchmark(
    cells: pd.DataFrame,
    *,
    controls_path: Path,
    rules_path: Path,
) -> dict[str, Any]:
    """Run primary tier-1 holdout metrics and all required ablations."""

    controls = load_controls(controls_path)
    audit = run_audit(cells, rules_path=rules_path)
    summary = pd.DataFrame(audit["target_summary"])
    scored = controls.merge(
        summary[["target", "score", "rna_only_score", "protein_concordance"]],
        on="target",
        how="left",
    )
    tier1 = scored[scored["tier"] == 1].copy()
    train_positive, holdout_positive = deterministic_holdout(controls)
    negatives = set(tier1.loc[tier1["y_true"] == 0, "target"])
    train = tier1[tier1["target"].isin(train_positive | negatives)].copy()
    holdout = tier1[tier1["target"].isin(holdout_positive | negatives)].copy()
    cutoff = _choose_cutoff(train)

    train_metrics = classification_metrics(train["y_true"].to_numpy(), train["score"].to_numpy())
    holdout_metrics = classification_metrics(
        holdout["y_true"].to_numpy(), holdout["score"].to_numpy()
    )
    coverage_curve = _coverage_accuracy(holdout, cutoff)

    cross = invalid_cross_study_scores(cells)
    cross_eval = tier1.merge(cross, on="target", how="left", suffixes=("_within", "_cross"))
    cross_metrics = classification_metrics(
        cross_eval["y_true"].to_numpy(),
        cross_eval["score_cross"].to_numpy(),
    )

    donor_tests = paired_effects(cells)
    cell_tests = invalid_cell_level_tests(cells)
    donor_significant = int(
        ((donor_tests["q_value"] < 0.05) & (donor_tests["effect"] > 0.50)).sum()
    )
    cell_significant = int(((cell_tests["q_value"] < 0.05) & (cell_tests["effect"] > 0.50)).sum())

    rna_only = tier1["rna_only_score"].to_numpy(dtype=float)
    protein_score = tier1["score"].to_numpy(dtype=float)
    rna_metrics = classification_metrics(tier1["y_true"].to_numpy(), rna_only)
    protein_metrics = classification_metrics(tier1["y_true"].to_numpy(), protein_score)

    sensitivity = run_sensitivity(cells, rules_path=rules_path, base_audit=audit)
    return {
        "schema_version": 1,
        "status": audit["status"],
        "holdout": {
            "rule": (
                "30% of tier-1 positives, deterministic SHA-256 ordering; "
                "never used to select cutoff"
            ),
            "train_positive_targets": sorted(train_positive),
            "holdout_positive_targets": sorted(holdout_positive),
            "negative_targets_used_in_both": sorted(negatives),
            "selected_cutoff_from_train": cutoff,
            "train_metrics": train_metrics,
            "holdout_metrics": holdout_metrics,
            "coverage_vs_accuracy": coverage_curve,
        },
        "all_tier1_metrics": classification_metrics(
            tier1["y_true"].to_numpy(), tier1["score"].to_numpy()
        ),
        "ablations": {
            "within_study_vs_cross_study": {
                "primary_within_study": classification_metrics(
                    tier1["y_true"].to_numpy(), tier1["score"].to_numpy()
                ),
                "INVALID_FOR_INFERENCE_cross_study": cross_metrics,
            },
            "pseudobulk_vs_cell_level": {
                "pseudobulk_significant_calls": donor_significant,
                "INVALID_FOR_INFERENCE_cell_level_significant_calls": cell_significant,
                "inflation_ratio": cell_significant / max(donor_significant, 1),
            },
            "rna_only_vs_protein_corroborated": {
                "rna_only": rna_metrics,
                "protein_corroborated_score": protein_metrics,
            },
            "detection_threshold_sweep": sensitivity,
        },
        "target_scores": scored[
            [
                "target",
                "label",
                "tier",
                "rna_only_score",
                "score",
                "protein_concordance",
            ]
        ]
        .sort_values("score", ascending=False)
        .to_dict(orient="records"),
        "claims": audit["claims"],
    }
