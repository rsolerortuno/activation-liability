"""Pseudobulk, paired testing, multiplicity and random-effects statistics."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd
from scipy import stats

GROUP_KEYS = ["study", "target", "cell_type", "stimulus"]


def benjamini_hochberg(p_values: Iterable[float]) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values preserving input order."""

    values = np.asarray(list(p_values), dtype=float)
    if values.size == 0:
        return values
    output = np.full(values.shape, np.nan, dtype=float)
    finite_mask = np.isfinite(values)
    finite_values = values[finite_mask]
    if finite_values.size == 0:
        return output
    order = np.argsort(finite_values)
    ranked = finite_values[order]
    adjusted = ranked * finite_values.size / np.arange(1, finite_values.size + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    finite_output = np.empty_like(adjusted)
    finite_output[order] = np.clip(adjusted, 0.0, 1.0)
    output[finite_mask] = finite_output
    return output


def pseudobulk(
    cells: pd.DataFrame,
    *,
    value_column: str = "rna_count",
    library_column: str = "library_size",
    pseudocount: float = 1.0,
) -> pd.DataFrame:
    """Aggregate cells to donor × cell type × condition pseudobulk CPM."""

    required = set(GROUP_KEYS + ["donor", "condition", value_column, library_column])
    missing = required - set(cells.columns)
    if missing:
        raise ValueError(f"missing input columns: {sorted(missing)}")
    use = cells.dropna(subset=[value_column, library_column]).copy()
    grouped = (
        use.groupby(GROUP_KEYS + ["donor", "condition"], observed=True, sort=True)
        .agg(
            value_sum=(value_column, "sum"),
            library_sum=(library_column, "sum"),
            n_cells=("cell_id", "nunique"),
        )
        .reset_index()
    )
    grouped["cpm"] = grouped["value_sum"] / grouped["library_sum"] * 1_000_000.0
    grouped["log2_cpm"] = np.log2(grouped["cpm"] + pseudocount)
    return grouped


def paired_effects(
    cells: pd.DataFrame,
    *,
    value_column: str = "rna_count",
    library_column: str = "library_size",
    min_pairs: int = 3,
) -> pd.DataFrame:
    """Compute within-study paired donor effects."""

    bulk = pseudobulk(
        cells,
        value_column=value_column,
        library_column=library_column,
    )
    rows: list[dict[str, object]] = []
    for key, group in bulk.groupby(GROUP_KEYS, observed=True, sort=True):
        wide = group.pivot(index="donor", columns="condition", values="log2_cpm")
        complete = (
            wide.dropna(subset=["resting", "activated"], how="any")
            if {
                "resting",
                "activated",
            }.issubset(wide.columns)
            else pd.DataFrame()
        )
        n_pairs = int(len(complete))
        if n_pairs < min_pairs:
            continue
        differences = complete["activated"].to_numpy() - complete["resting"].to_numpy()
        effect = float(np.mean(differences))
        if n_pairs > 1:
            standard_error = float(np.std(differences, ddof=1) / np.sqrt(n_pairs))
        else:
            standard_error = float("nan")
        if np.isclose(np.std(differences, ddof=1), 0.0):
            p_value = 0.0 if not np.isclose(effect, 0.0) else 1.0
            standard_error = max(standard_error, 1e-8)
        else:
            p_value = float(stats.ttest_1samp(differences, popmean=0.0).pvalue)
        row = dict(zip(GROUP_KEYS, key, strict=True))
        row.update(
            {
                "effect": effect,
                "standard_error": max(standard_error, 1e-8),
                "p_value": p_value,
                "n_pairs": n_pairs,
            }
        )
        rows.append(row)
    result = pd.DataFrame.from_records(rows)
    if not result.empty:
        result["q_value"] = benjamini_hochberg(result["p_value"])
    return result


def random_effects_meta(effects: np.ndarray, standard_errors: np.ndarray) -> dict[str, float]:
    """DerSimonian-Laird random-effects estimate for one group."""

    effects = np.asarray(effects, dtype=float)
    standard_errors = np.asarray(standard_errors, dtype=float)
    if effects.size == 0 or effects.size != standard_errors.size:
        raise ValueError("effects and standard_errors must be non-empty and equally sized")
    variances = np.maximum(standard_errors**2, 1e-12)
    fixed_weights = 1.0 / variances
    fixed_mean = float(np.sum(fixed_weights * effects) / np.sum(fixed_weights))
    q = float(np.sum(fixed_weights * (effects - fixed_mean) ** 2))
    k = effects.size
    if k == 1:
        tau2 = float("nan")
        q = float("nan")
        i2 = float("nan")
        weighting_tau2 = 0.0
    else:
        c_value = float(np.sum(fixed_weights) - np.sum(fixed_weights**2) / np.sum(fixed_weights))
        tau2 = max(0.0, (q - (k - 1)) / c_value) if c_value > 0 else 0.0
        i2 = max(0.0, (q - (k - 1)) / q) * 100.0 if q > 0 else 0.0
        weighting_tau2 = tau2
    random_weights = 1.0 / (variances + weighting_tau2)
    estimate = float(np.sum(random_weights * effects) / np.sum(random_weights))
    standard_error = float(np.sqrt(1.0 / np.sum(random_weights)))
    return {
        "estimate": estimate,
        "standard_error": standard_error,
        "ci_low": estimate - 1.96 * standard_error,
        "ci_high": estimate + 1.96 * standard_error,
        "tau2": tau2,
        "q": q,
        "i2": i2,
        "n_studies": float(k),
    }


def meta_analyse(within_study: pd.DataFrame) -> pd.DataFrame:
    """Meta-analyse only identical target × cell type × stimulus contrasts."""

    keys = ["target", "cell_type", "stimulus"]
    rows: list[dict[str, object]] = []
    for key, group in within_study.groupby(keys, observed=True, sort=True):
        result = random_effects_meta(
            group["effect"].to_numpy(),
            group["standard_error"].to_numpy(),
        )
        row = dict(zip(keys, key, strict=True))
        row.update(result)
        row["total_pairs"] = int(group["n_pairs"].sum())
        row["min_q_value"] = float(group["q_value"].min())
        rows.append(row)
    return pd.DataFrame.from_records(rows)


def invalid_cell_level_tests(cells: pd.DataFrame) -> pd.DataFrame:
    """Deliberately invalid independent-cell tests for the pseudoreplication ablation."""

    data = cells.copy()
    data["log_expression"] = np.log2(data["rna_count"] / data["library_size"] * 1_000_000.0 + 1.0)
    rows: list[dict[str, object]] = []
    for key, group in data.groupby(GROUP_KEYS, observed=True, sort=True):
        resting = group.loc[group["condition"] == "resting", "log_expression"].to_numpy()
        activated = group.loc[group["condition"] == "activated", "log_expression"].to_numpy()
        if len(resting) < 2 or len(activated) < 2:
            continue
        test = stats.ttest_ind(activated, resting, equal_var=False)
        row = dict(zip(GROUP_KEYS, key, strict=True))
        row.update(
            {
                "effect": float(np.mean(activated) - np.mean(resting)),
                "p_value": float(test.pvalue),
                "n_cells": int(len(resting) + len(activated)),
            }
        )
        rows.append(row)
    result = pd.DataFrame.from_records(rows)
    if not result.empty:
        result["q_value"] = benjamini_hochberg(result["p_value"])
    return result
