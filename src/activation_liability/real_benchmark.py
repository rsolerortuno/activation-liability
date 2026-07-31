"""Adversarially constrained benchmark for the public real-data validation layer."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd
import yaml
from scipy import stats

from activation_liability.benchmark import deterministic_holdout, load_controls
from activation_liability.statistics import pseudobulk, random_effects_meta


def load_real_benchmark_config(path: Path) -> dict[str, Any]:
    """Load and minimally validate the pre-registered real benchmark configuration."""

    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("real benchmark configuration must be a mapping")
    required = {"version", "score", "endpoints", "target_lineages"}
    if required - set(payload):
        raise ValueError("real benchmark configuration is missing required sections")
    endpoints = payload["endpoints"]
    if not isinstance(endpoints, list) or not endpoints:
        raise ValueError("real benchmark requires at least one endpoint")
    stimuli = [str(endpoint["stimulus"]) for endpoint in endpoints]
    if len(stimuli) != len(set(stimuli)):
        raise ValueError("real benchmark endpoint stimuli must be unique")
    return cast(dict[str, Any], payload)


def _json_records(frame: pd.DataFrame) -> list[dict[str, Any]]:
    clean = frame.replace({np.nan: None, np.inf: None, -np.inf: None})
    return cast(list[dict[str, Any]], clean.to_dict(orient="records"))


def _ranking_metrics(frame: pd.DataFrame) -> dict[str, Any]:
    use = frame.dropna(subset=["score"]).copy()
    labels = use["y_true"].to_numpy(dtype=int)
    scores = use["score"].to_numpy(dtype=float)
    n_positive = int(labels.sum())
    n_negative = int(len(labels) - n_positive)
    result: dict[str, Any] = {
        "n_targets": int(len(use)),
        "n_positive": n_positive,
        "n_negative": n_negative,
        "targets": sorted(use["target"].astype(str).tolist()),
    }
    if n_positive == 0 or n_negative == 0:
        result.update(
            {
                "auroc": None,
                "average_precision": None,
                "precision_at_10": None,
                "exact_mann_whitney_p_one_sided": None,
                "perfect_separation_observed": None,
            }
        )
        return result
    positive = scores[labels == 1]
    negative = scores[labels == 0]
    pairwise = (positive[:, None] > negative[None, :]).mean() + 0.5 * (
        positive[:, None] == negative[None, :]
    ).mean()
    order = np.argsort(scores)[::-1]
    k = min(10, len(order))
    ranked_labels = labels[order]
    precision_at_10 = float(ranked_labels[:k].mean())
    from sklearn.metrics import average_precision_score

    exact = stats.mannwhitneyu(
        positive,
        negative,
        alternative="greater",
        method="exact",
    )
    result.update(
        {
            "auroc": float(pairwise),
            "average_precision": float(average_precision_score(labels, scores)),
            "precision_at_10": precision_at_10,
            "exact_mann_whitney_p_one_sided": float(exact.pvalue),
            "perfect_separation_observed": bool(positive.min() > negative.max()),
            "interpretation_guardrail": (
                "The exact p-value conditions on this small curated control set and does not "
                "measure dataset-selection or external-validity uncertainty."
            ),
        }
    )
    return result


def _choose_cutoff(train: pd.DataFrame) -> float:
    from sklearn.metrics import balanced_accuracy_score

    unique = np.unique(train["score"].dropna().to_numpy(dtype=float))
    if len(unique) == 0:
        return float("nan")
    if len(unique) == 1:
        return float(unique[0])
    candidates = (unique[:-1] + unique[1:]) / 2.0
    best = float(candidates[0])
    best_accuracy = -1.0
    for cutoff in candidates:
        predicted = (train["score"].to_numpy(dtype=float) >= cutoff).astype(int)
        accuracy = float(balanced_accuracy_score(train["y_true"], predicted))
        if accuracy > best_accuracy:
            best_accuracy = accuracy
            best = float(cutoff)
    return best


def _coverage_accuracy(frame: pd.DataFrame, cutoff: float) -> list[dict[str, float]]:
    if frame.empty or not np.isfinite(cutoff):
        return []
    confidence = np.abs(frame["score"].to_numpy(dtype=float) - cutoff)
    predicted = (frame["score"].to_numpy(dtype=float) >= cutoff).astype(int)
    correct = predicted == frame["y_true"].to_numpy(dtype=int)
    rows: list[dict[str, float]] = []
    for requested in (1.0, 0.75, 0.50, 0.25):
        keep_n = max(1, int(np.ceil(requested * len(frame))))
        keep = np.argsort(confidence)[::-1][:keep_n]
        rows.append(
            {
                "coverage": keep_n / len(frame),
                "accuracy": float(correct[keep].mean()),
                "abstention": 1.0 - keep_n / len(frame),
            }
        )
    return rows


def score_real_controls(
    audit_rows: pd.DataFrame,
    controls: pd.DataFrame,
    config: dict[str, Any],
    *,
    excluded_stimuli: Iterable[str] = (),
) -> pd.DataFrame:
    """Score controls only in pre-declared lineages and endpoints using a 95% LCB."""

    excluded = set(excluded_stimuli)
    endpoints = [
        endpoint for endpoint in config["endpoints"] if str(endpoint["stimulus"]) not in excluded
    ]
    endpoint_stimuli = {str(endpoint["stimulus"]) for endpoint in endpoints}
    lineage_map = cast(dict[str, list[str]], config["target_lineages"])
    z_value = float(config["score"]["z_value"])
    observability_cutoff = float(config["score"]["observability_fraction_cutoff"])
    observed_lineages = set(
        audit_rows.loc[audit_rows["stimulus"].isin(endpoint_stimuli), "cell_type"].astype(str)
    )
    rows: list[dict[str, Any]] = []
    for control in controls.itertuples(index=False):
        target = str(control.target)
        expected = [str(value) for value in lineage_map.get(target, [])]
        covered = sorted(set(expected) & observed_lineages)
        candidates = audit_rows[
            audit_rows["target"].eq(target)
            & audit_rows["stimulus"].isin(endpoint_stimuli)
            & audit_rows["cell_type"].isin(expected)
        ].copy()
        candidates["lcb95"] = (
            candidates["inducibility_lfc"] - z_value * candidates["inducibility_se"]
        )
        if candidates.empty:
            max_fraction = float("nan")
            best: pd.Series[Any] | None = None
        else:
            max_fraction = float(
                candidates[["positive_fraction_resting", "positive_fraction_activated"]]
                .max(axis=1)
                .max()
            )
            best = candidates.loc[candidates["lcb95"].idxmax()]
        if not expected:
            reason = "NO_PREDECLARED_LINEAGE"
        elif not covered:
            reason = "LINEAGE_NOT_COVERED"
        elif not np.isfinite(max_fraction) or max_fraction < observability_cutoff:
            reason = "TARGET_NOT_OBSERVABLE_AT_COVERAGE_THRESHOLD"
        else:
            reason = None
        score = float(best["lcb95"]) if best is not None and reason is None else float("nan")
        rows.append(
            {
                "target": target,
                "label": str(control.label),
                "tier": int(control.tier),
                "y_true": int(control.y_true),
                "expected_lineages": expected,
                "covered_lineages": covered,
                "observable_fraction_max": max_fraction,
                "observable": reason is None,
                "abstention_reason": reason,
                "score": score,
                "driver_stimulus": None if best is None else str(best["stimulus"]),
                "driver_cell_type": None if best is None else str(best["cell_type"]),
                "driver_effect": None if best is None else float(best["inducibility_lfc"]),
                "driver_standard_error": None if best is None else float(best["inducibility_se"]),
                "driver_total_pairs": None if best is None else int(best["total_pairs"]),
                "driver_positive_fraction_resting": None
                if best is None
                else float(best["positive_fraction_resting"]),
                "driver_positive_fraction_activated": None
                if best is None
                else float(best["positive_fraction_activated"]),
            }
        )
    return pd.DataFrame.from_records(rows)


def _holdout_evaluation(scored: pd.DataFrame, controls: pd.DataFrame) -> dict[str, Any]:
    tier1 = scored[(scored["tier"] == 1) & scored["observable"]].copy()
    train_positive, holdout_positive = deterministic_holdout(controls)
    negative_targets = set(tier1.loc[tier1["y_true"] == 0, "target"])
    train = tier1[tier1["target"].isin(train_positive | negative_targets)].copy()
    holdout = tier1[tier1["target"].isin(holdout_positive | negative_targets)].copy()
    cutoff = _choose_cutoff(train)
    return {
        "integrity_status": "DIAGNOSTIC_NOT_CONFIRMATORY",
        "integrity_reason": (
            "The deterministic target split was preserved, but the constrained real endpoint "
            "and scoring rules were finalised after inspecting an earlier unconstrained real "
            "benchmark. This holdout must not be described as untouched confirmation."
        ),
        "future_confirmatory_requirement": (
            "Freeze the v0.3 rules and evaluate the unopened GSE96583 IFN-I replication plus "
            "an independent paired tissue cohort without further tuning."
        ),
        "rule": (
            "30% of tier-1 positives by deterministic SHA-256 ordering; all observable "
            "tier-1 negatives are reused because the negative set is small."
        ),
        "train_positive_targets": sorted(train_positive),
        "holdout_positive_targets": sorted(holdout_positive),
        "negative_targets_used_in_both": sorted(negative_targets),
        "unobservable_tier1_targets": sorted(
            scored.loc[(scored["tier"] == 1) & ~scored["observable"], "target"].astype(str)
        ),
        "selected_cutoff_from_train": cutoff,
        "train_metrics": _ranking_metrics(train),
        "holdout_metrics": _ranking_metrics(holdout),
        "coverage_vs_accuracy": _coverage_accuracy(holdout, cutoff),
    }


def _driver_donor_robustness(
    cells: pd.DataFrame,
    scored: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    endpoint_study = {
        str(endpoint["stimulus"]): str(endpoint["study"]) for endpoint in config["endpoints"]
    }
    bulk = pseudobulk(cells)
    rows: list[dict[str, Any]] = []
    for record in scored[scored["observable"]].itertuples(index=False):
        stimulus = str(record.driver_stimulus)
        study = endpoint_study[stimulus]
        selected = bulk[
            bulk["study"].eq(study)
            & bulk["target"].eq(record.target)
            & bulk["cell_type"].eq(record.driver_cell_type)
            & bulk["stimulus"].eq(stimulus)
        ]
        wide = selected.pivot(index="donor", columns="condition", values="log2_cpm")
        complete = wide.dropna(subset=["resting", "activated"])
        differences = (complete["activated"] - complete["resting"]).to_numpy(dtype=float)
        if len(differences) == 0:
            continue
        if np.allclose(differences, 0.0):
            wilcoxon_p = 1.0
        else:
            wilcoxon_p = float(
                stats.wilcoxon(differences, alternative="greater", method="auto").pvalue
            )
        loo_means = [
            float(np.delete(differences, index).mean())
            for index in range(len(differences))
            if len(differences) > 1
        ]
        rows.append(
            {
                "target": str(record.target),
                "label": str(record.label),
                "study": study,
                "stimulus": stimulus,
                "cell_type": str(record.driver_cell_type),
                "n_paired_donors": int(len(differences)),
                "mean_log2_cpm_difference": float(differences.mean()),
                "median_log2_cpm_difference": float(np.median(differences)),
                "fraction_donors_positive": float((differences > 0).mean()),
                "wilcoxon_p_one_sided": wilcoxon_p,
                "leave_one_donor_out_min_mean": min(loo_means) if loo_means else None,
                "leave_one_donor_out_max_mean": max(loo_means) if loo_means else None,
                "direction_stable_leave_one_donor_out": bool(loo_means and min(loo_means) > 0),
            }
        )
    return pd.DataFrame.from_records(rows)


def _timecourse_summary(audit_rows: pd.DataFrame, config: dict[str, Any]) -> pd.DataFrame:
    stimuli = set(config["secondary_timecourse"]["stimuli"])
    target_lineages = cast(dict[str, list[str]], config["target_lineages"])
    use = audit_rows[audit_rows["stimulus"].isin(stimuli)].copy()
    use = use[use["cell_type"].eq("T_cell")]
    use["time_hours"] = use["stimulus"].str.extract(r"_(\d+)h$")[0].astype(int)
    use["lcb95"] = (
        use["inducibility_lfc"] - float(config["score"]["z_value"]) * use["inducibility_se"]
    )
    use = use[["T_cell" in target_lineages.get(str(target), []) for target in use["target"]]].copy()
    peak = use.loc[use.groupby("target", observed=True)["lcb95"].idxmax()].copy()
    peak = peak.rename(
        columns={
            "time_hours": "peak_time_hours",
            "inducibility_lfc": "peak_effect",
            "inducibility_se": "peak_standard_error",
            "lcb95": "peak_lcb95",
        }
    )
    at_24 = use[use["time_hours"].eq(24)][
        ["target", "inducibility_lfc", "inducibility_se", "lcb95"]
    ].rename(
        columns={
            "inducibility_lfc": "effect_24h",
            "inducibility_se": "standard_error_24h",
            "lcb95": "lcb95_24h",
        }
    )
    return peak[
        ["target", "peak_time_hours", "peak_effect", "peak_standard_error", "peak_lcb95"]
    ].merge(at_24, on="target", how="left")


def _classification_metrics(frame: pd.DataFrame, cutoff: float) -> dict[str, Any]:
    """Evaluate a frozen score cutoff without selecting anything on the new cohort."""

    use = frame.dropna(subset=["score"]).copy()
    if use.empty or not np.isfinite(cutoff):
        return {
            "n_targets": int(len(use)),
            "cutoff": cutoff,
            "accuracy": None,
            "balanced_accuracy": None,
            "sensitivity": None,
            "specificity": None,
            "positive_predictive_value": None,
            "confusion": None,
        }
    truth = use["y_true"].to_numpy(dtype=int)
    predicted = (use["score"].to_numpy(dtype=float) >= cutoff).astype(int)
    true_positive = int(((truth == 1) & (predicted == 1)).sum())
    true_negative = int(((truth == 0) & (predicted == 0)).sum())
    false_positive = int(((truth == 0) & (predicted == 1)).sum())
    false_negative = int(((truth == 1) & (predicted == 0)).sum())
    sensitivity = (
        true_positive / (true_positive + false_negative) if true_positive + false_negative else None
    )
    specificity = (
        true_negative / (true_negative + false_positive) if true_negative + false_positive else None
    )
    precision = (
        true_positive / (true_positive + false_positive) if true_positive + false_positive else None
    )
    balanced = (
        (sensitivity + specificity) / 2.0
        if sensitivity is not None and specificity is not None
        else None
    )
    return {
        "n_targets": int(len(use)),
        "cutoff": cutoff,
        "accuracy": float((truth == predicted).mean()),
        "balanced_accuracy": balanced,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "positive_predictive_value": precision,
        "confusion": {
            "true_positive": true_positive,
            "true_negative": true_negative,
            "false_positive": false_positive,
            "false_negative": false_negative,
        },
    }


def _score_agreement(
    baseline_scores: pd.DataFrame,
    confirmation_scores: pd.DataFrame,
    *,
    cutoff: float,
) -> dict[str, Any]:
    """Compare frozen v0.3 scores with the unopened cohort without fitting a mapping."""

    left = baseline_scores[baseline_scores["observable"]].copy()
    right = confirmation_scores[confirmation_scores["observable"]].copy()
    merged = left.merge(
        right,
        on=["target", "label", "tier", "y_true"],
        suffixes=("_v0_3", "_gse96583"),
    )
    if len(merged) >= 2:
        correlation = stats.spearmanr(
            merged["score_v0_3"].to_numpy(dtype=float),
            merged["score_gse96583"].to_numpy(dtype=float),
        )
        spearman = float(correlation.statistic)
        p_value = float(correlation.pvalue)
    else:
        spearman = float("nan")
        p_value = float("nan")
    if merged.empty:
        classification_agreement = float("nan")
    else:
        baseline_call = merged["score_v0_3"].to_numpy(dtype=float) >= cutoff
        confirmation_call = merged["score_gse96583"].to_numpy(dtype=float) >= cutoff
        classification_agreement = float((baseline_call == confirmation_call).mean())
    return {
        "n_common_observable_targets": int(len(merged)),
        "spearman_score_correlation": spearman,
        "spearman_p_value": p_value,
        "frozen_cutoff_call_agreement": classification_agreement,
        "rows": _json_records(
            merged[
                [
                    "target",
                    "label",
                    "tier",
                    "score_v0_3",
                    "score_gse96583",
                    "driver_stimulus_v0_3",
                    "driver_cell_type_v0_3",
                    "driver_cell_type_gse96583",
                ]
            ]
        ),
    }


def _ifn_i_cross_time_replication(audit_rows: pd.DataFrame) -> dict[str, Any]:
    """Summarise IFN-I agreement across 6 h and 18 h as a labelled sensitivity analysis."""

    baseline = audit_rows[
        audit_rows["stimulus"].eq("IFN_I_18h") & audit_rows["cell_type"].eq("Myeloid")
    ][["target", "inducibility_lfc", "inducibility_se"]].rename(
        columns={
            "inducibility_lfc": "effect_18h",
            "inducibility_se": "se_18h",
        }
    )
    confirmation = audit_rows[
        audit_rows["stimulus"].eq("IFN_I_6h_confirmation") & audit_rows["cell_type"].eq("Myeloid")
    ][["target", "inducibility_lfc", "inducibility_se"]].rename(
        columns={
            "inducibility_lfc": "effect_6h",
            "inducibility_se": "se_6h",
        }
    )
    merged = baseline.merge(confirmation, on="target", how="inner")
    if len(merged) >= 2:
        correlation = stats.spearmanr(
            merged["effect_18h"].to_numpy(dtype=float),
            merged["effect_6h"].to_numpy(dtype=float),
        )
        spearman = float(correlation.statistic)
        p_value = float(correlation.pvalue)
    else:
        spearman = float("nan")
        p_value = float("nan")
    rows: list[dict[str, Any]] = []
    for record in merged.itertuples(index=False):
        meta = random_effects_meta(
            np.asarray([record.effect_18h, record.effect_6h], dtype=float),
            np.asarray([record.se_18h, record.se_6h], dtype=float),
        )
        rows.append(
            {
                "target": str(record.target),
                "effect_6h": float(record.effect_6h),
                "effect_18h": float(record.effect_18h),
                "same_direction": bool(np.sign(record.effect_6h) == np.sign(record.effect_18h)),
                "cross_time_random_effects_estimate": float(meta["estimate"]),
                "cross_time_i2": float(meta["i2"]),
            }
        )
    direction_agreement = (
        float(np.mean([row["same_direction"] for row in rows])) if rows else float("nan")
    )
    return {
        "status": "SENSITIVITY_ONLY_DIFFERENT_TIMEPOINTS_AND_DISEASE_CONTEXTS",
        "comparison": "GSE96583 IFN-beta 6 h versus GSE157857 IFN-beta 18 h, broad Myeloid",
        "n_shared_targets": int(len(merged)),
        "effect_spearman": spearman,
        "effect_spearman_p_value": p_value,
        "direction_agreement": direction_agreement,
        "heterogeneity_warning": (
            "I-squared values use only two studies at different time points and are descriptive, "
            "not evidence that exact-endpoint heterogeneity has been resolved."
        ),
        "rows": rows,
    }


def run_gse96583_confirmation(
    cells: pd.DataFrame,
    audit_rows: pd.DataFrame,
    *,
    baseline_benchmark: dict[str, Any],
    baseline_audit_rows: pd.DataFrame,
    controls_path: Path,
    config_path: Path,
) -> dict[str, Any]:
    """Open the frozen GSE96583 cohort once and evaluate unchanged v0.3 rules."""

    controls = load_controls(controls_path)
    frozen = load_real_benchmark_config(config_path)
    confirmation_config = dict(frozen)
    confirmation_config["endpoints"] = [
        {
            "study": "GSE96583",
            "stimulus": "IFN_I_6h_confirmation",
            "axis": "IFN-I",
            "time_hours": 6,
            "role": "external-confirmation",
        }
    ]
    scored = score_real_controls(audit_rows, controls, confirmation_config)
    tier1 = scored[(scored["tier"] == 1) & scored["observable"]].copy()
    expanded = scored[scored["observable"]].copy()
    frozen_cutoff = float(baseline_benchmark["holdout"]["selected_cutoff_from_train"])
    train_positive, holdout_positive = deterministic_holdout(controls)
    del train_positive
    negative_targets = set(tier1.loc[tier1["y_true"] == 0, "target"].astype(str))
    holdout = tier1[tier1["target"].isin(holdout_positive | negative_targets)].copy()
    observable_holdout_positives = sorted(
        holdout.loc[holdout["y_true"] == 1, "target"].astype(str).tolist()
    )
    if observable_holdout_positives:
        holdout_status = "EVALUABLE"
    else:
        holdout_status = "INCONCLUSIVE_NO_OBSERVABLE_POSITIVE_HOLDOUT_TARGETS"

    sensitivity_rows: list[dict[str, Any]] = []
    for fraction_cutoff in (0.01, 0.025, 0.05, 0.10):
        sensitivity_config = dict(confirmation_config)
        sensitivity_config["score"] = dict(frozen["score"])
        sensitivity_config["score"]["observability_fraction_cutoff"] = fraction_cutoff
        sensitivity_scored = score_real_controls(audit_rows, controls, sensitivity_config)
        sensitivity_tier1 = sensitivity_scored[
            (sensitivity_scored["tier"] == 1) & sensitivity_scored["observable"]
        ]
        sensitivity_holdout = sensitivity_tier1[
            sensitivity_tier1["target"].isin(holdout_positive | negative_targets)
        ]
        sensitivity_rows.append(
            {
                "observability_fraction_cutoff": fraction_cutoff,
                "observable_tier1_controls": int(len(sensitivity_tier1)),
                "observable_holdout_positives": sorted(
                    sensitivity_holdout.loc[sensitivity_holdout["y_true"] == 1, "target"].astype(
                        str
                    )
                ),
                "tier1_metrics": _ranking_metrics(sensitivity_tier1),
                "frozen_cutoff_classification": _classification_metrics(
                    sensitivity_tier1,
                    frozen_cutoff,
                ),
            }
        )

    robustness = _driver_donor_robustness(cells, scored, confirmation_config)
    replication_columns = [
        "target",
        "cell_type",
        "stimulus",
        "inducibility_lfc",
        "inducibility_se",
    ]
    combined_audit = pd.concat(
        [
            baseline_audit_rows[replication_columns],
            audit_rows[replication_columns],
        ],
        ignore_index=True,
    )
    baseline_scores = pd.DataFrame(baseline_benchmark["target_scores"])
    return {
        "schema_version": 1,
        "status": "PARTIAL_EXTERNAL_CONFIRMATION_WITH_TARGET_HOLDOUT_ABSTENTION",
        "protocol": {
            "rule_source": "v0.3 frozen endpoint-lineage-observability-score specification",
            "endpoint": confirmation_config["endpoints"][0],
            "score_definition": str(frozen["score"]["definition"]),
            "observability_fraction_cutoff": float(
                frozen["score"]["observability_fraction_cutoff"]
            ),
            "frozen_classification_cutoff": frozen_cutoff,
            "no_tuning_statement": (
                "No endpoint, lineage mapping, score definition, z value, observability cutoff "
                "or target split was changed after GSE96583 was opened."
            ),
        },
        "coverage": {
            "all_controls": int(len(scored)),
            "observable_controls": int(scored["observable"].sum()),
            "tier1_controls": int((scored["tier"] == 1).sum()),
            "observable_tier1_controls": int(len(tier1)),
            "abstained_targets": scored.loc[
                ~scored["observable"],
                ["target", "label", "tier", "expected_lineages", "abstention_reason"],
            ].to_dict(orient="records"),
        },
        "tier1_ranking_metrics": _ranking_metrics(tier1),
        "expanded_tier1_tier2_metrics": _ranking_metrics(expanded),
        "frozen_cutoff_classification_tier1": _classification_metrics(tier1, frozen_cutoff),
        "target_holdout": {
            "status": holdout_status,
            "preassigned_positive_targets": sorted(holdout_positive),
            "observable_positive_targets": observable_holdout_positives,
            "observable_negative_targets": sorted(negative_targets),
            "ranking_metrics": _ranking_metrics(holdout),
            "frozen_cutoff_classification": _classification_metrics(holdout, frozen_cutoff),
            "interpretation": (
                "The external cohort is informative for observable controls, but it cannot "
                "confirm the preassigned positive target holdout when all held-out positives "
                "abstain at the frozen 5% observability rule."
            ),
        },
        "agreement_with_v0_3": _score_agreement(
            baseline_scores,
            scored,
            cutoff=frozen_cutoff,
        ),
        "ifn_i_cross_time_replication": _ifn_i_cross_time_replication(combined_audit),
        "observability_threshold_sensitivity": sensitivity_rows,
        "donor_robustness": {
            "definition": (
                "One-sided paired Wilcoxon and leave-one-donor-out means for each target's "
                "GSE96583 driver lineage."
            ),
            "rows": _json_records(robustness),
        },
        "target_scores": _json_records(scored.sort_values("score", ascending=False)),
        "claims": {
            "permitted": [
                (
                    "GSE96583 independently ranks the observable tier-1 controls under the "
                    "unchanged v0.3 rules."
                ),
                "Eight donor-resolved within-study IFN-beta contrasts were analysed.",
            ],
            "conditional": [
                (
                    "The confirmation applies to observable PBMC lineages in SLE patient "
                    "samples after 6 h IFN-beta stimulation."
                ),
                (
                    "Perfect ranking separation does not imply that the frozen target holdout "
                    "was confirmed because its positive targets abstained."
                ),
            ],
            "unsupported": [
                "The complete target holdout is externally confirmed.",
                "The result generalises to healthy donors or tissue-resident lineages.",
                "The cohort provides protein-level corroboration.",
                "The tool predicts clinical toxicity.",
            ],
        },
    }


def run_real_benchmark(
    cells: pd.DataFrame,
    audit_rows: pd.DataFrame,
    *,
    controls_path: Path,
    config_path: Path,
) -> dict[str, Any]:
    """Run the frozen public-data benchmark with explicit coverage and stress tests."""

    controls = load_controls(controls_path)
    config = load_real_benchmark_config(config_path)
    scored = score_real_controls(audit_rows, controls, config)
    primary = scored[(scored["tier"] == 1) & scored["observable"]].copy()
    expanded = scored[scored["observable"]].copy()
    endpoint_ablations: dict[str, Any] = {}
    for endpoint in config["endpoints"]:
        stimulus = str(endpoint["stimulus"])
        ablated = score_real_controls(
            audit_rows,
            controls,
            config,
            excluded_stimuli=[stimulus],
        )
        same_targets = ablated[ablated["target"].isin(primary["target"])].copy()
        endpoint_ablations[stimulus] = {
            "omitted_axis": str(endpoint["axis"]),
            "metrics_on_remaining_observable_primary_targets": _ranking_metrics(
                same_targets[same_targets["observable"]]
            ),
            "coverage_of_frozen_primary_set": float(same_targets["observable"].mean()),
            "newly_abstained_targets": sorted(
                same_targets.loc[~same_targets["observable"], "target"].astype(str)
            ),
        }
    robustness = _driver_donor_robustness(cells, scored, config)
    timecourse = _timecourse_summary(audit_rows, config)
    control_targets = set(controls["target"].astype(str))
    adt_targets = set(cells.loc[cells["adt_count"].notna(), "target"].astype(str))
    protein_overlap = sorted(control_targets & adt_targets)
    return {
        "schema_version": 1,
        "status": "COMPUTED_REAL_PUBLIC_DATA_WITH_COVERAGE_ABSTENTION",
        "pre_registration": {
            "config_version": int(config["version"]),
            "endpoint_count": len(config["endpoints"]),
            "score_definition": str(config["score"]["definition"]),
            "endpoints": config["endpoints"],
            "lineage_mapping_basis": str(config["lineage_mapping_basis"]),
            "annotation_leakage_control": config["annotation_leakage_control"],
        },
        "coverage": {
            "all_controls": int(len(scored)),
            "observable_controls": int(scored["observable"].sum()),
            "tier1_controls": int((scored["tier"] == 1).sum()),
            "observable_tier1_controls": int(((scored["tier"] == 1) & scored["observable"]).sum()),
            "abstained_targets": scored.loc[
                ~scored["observable"],
                [
                    "target",
                    "label",
                    "tier",
                    "expected_lineages",
                    "abstention_reason",
                ],
            ].to_dict(orient="records"),
        },
        "primary_tier1_metrics": _ranking_metrics(primary),
        "expanded_tier1_tier2_metrics": _ranking_metrics(expanded),
        "holdout": _holdout_evaluation(scored, controls),
        "endpoint_leave_one_out": endpoint_ablations,
        "donor_robustness": {
            "definition": (
                "One-sided paired Wilcoxon and leave-one-donor-out means for each target's "
                "pre-registered driver endpoint."
            ),
            "rows": _json_records(robustness),
        },
        "timecourse_validation": {
            "status": "SECONDARY_NOT_USED_TO_SELECT_PRIMARY_ENDPOINT",
            "primary_endpoint_hours": 24,
            "rows": _json_records(timecourse),
        },
        "protein_benchmark": {
            "status": ("AVAILABLE" if protein_overlap else "UNAVAILABLE_NO_CONTROL_TARGET_OVERLAP"),
            "control_targets_with_matched_adt": protein_overlap,
            "exploratory_adt_targets": sorted(adt_targets - control_targets),
            "claim": (
                "No RNA-versus-protein performance improvement is claimed for the control "
                "benchmark when overlap is empty."
            ),
        },
        "heterogeneity": {
            "status": "NOT_ESTIMABLE_FOR_EXACT_ENDPOINTS",
            "reason": (
                "Each exact target-cell-stimulus endpoint currently has one contributing "
                "study. Leave-one-endpoint-out tests dependence but is not an I-squared estimate."
            ),
        },
        "target_scores": _json_records(scored.sort_values("score", ascending=False)),
        "claims": {
            "permitted": [
                (
                    "The evaluated public cohorts separate the observable tier-1 controls "
                    "under the frozen endpoint and lineage rules."
                ),
                (
                    "Within-study donor-paired inducibility is estimated for the reported "
                    "cohorts and endpoints."
                ),
            ],
            "conditional": [
                (
                    "Performance applies only to observable targets in PBMC, myeloid and "
                    "CD4-memory T-cell coverage."
                ),
                (
                    "Perfect separation is an observed property of a small curated control "
                    "set, not proof of generalisation."
                ),
            ],
            "unsupported": [
                "The tool predicts clinical toxicity.",
                "The real benchmark validates epithelial or tissue-resident lineages.",
                "The control benchmark has protein-level corroboration.",
                "Between-study heterogeneity is low.",
                "An abstained or non-induced target is safe.",
            ],
        },
    }
