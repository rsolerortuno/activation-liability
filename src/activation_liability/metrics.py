"""Audit metrics derived from cells and paired effects."""

from __future__ import annotations

import numpy as np
import pandas as pd

from activation_liability.statistics import meta_analyse, paired_effects, pseudobulk

RELEVANT_STIMULUS: dict[str, tuple[str, ...]] = {
    "T_cell": ("lymphocyte_activation", "IFN_I", "IFN_II"),
    "NK_cell": ("lymphocyte_activation", "IFN_I", "IFN_II"),
    "B_cell": ("lymphocyte_activation", "IFN_I", "IFN_II", "TLR"),
    "Plasma_cell": ("lymphocyte_activation", "IFN_I", "IFN_II", "TLR"),
    "Myeloid": ("IFN_I", "IFN_II", "TLR", "TNF"),
    "Endothelial": ("IFN_I", "IFN_II", "TNF"),
    "Fibroblast": ("tissue_damage", "IFN_I", "IFN_II", "TNF"),
    "Epithelial": ("IFN_I", "IFN_II", "TNF", "tissue_damage"),
}


def is_relevant_stimulus(cell_type: str, stimulus: str) -> bool:
    """Return whether a pre-declared stimulus axis is biologically covered for a lineage."""

    return any(stimulus.startswith(prefix) for prefix in RELEVANT_STIMULUS.get(cell_type, ()))


def positive_fractions(cells: pd.DataFrame, *, count_threshold: int = 0) -> pd.DataFrame:
    """Calculate study-balanced positive-cell fractions per contrast."""

    use = cells.copy()
    use["positive"] = use["rna_count"] > count_threshold
    per_study = (
        use.groupby(
            ["study", "target", "cell_type", "stimulus", "condition"],
            observed=True,
            sort=True,
        )["positive"]
        .mean()
        .reset_index(name="fraction")
    )
    balanced = (
        per_study.groupby(
            ["target", "cell_type", "stimulus", "condition"],
            observed=True,
            sort=True,
        )["fraction"]
        .mean()
        .reset_index()
    )
    wide = balanced.pivot(
        index=["target", "cell_type", "stimulus"],
        columns="condition",
        values="fraction",
    ).reset_index()
    for condition in ("resting", "activated"):
        if condition not in wide.columns:
            wide[condition] = 0.0
    return wide.rename(
        columns={
            "resting": "positive_fraction_resting",
            "activated": "positive_fraction_activated",
        }
    )


def protein_meta(cells: pd.DataFrame) -> pd.DataFrame:
    """Compute matched ADT paired effects and meta-analysis where available."""

    if "adt_count" not in cells or cells["adt_count"].notna().sum() == 0:
        return pd.DataFrame(
            columns=["target", "cell_type", "stimulus", "adt_estimate", "adt_n_studies"]
        )
    effects = paired_effects(
        cells,
        value_column="adt_count",
        library_column="adt_library_size",
    )
    if effects.empty:
        return pd.DataFrame(
            columns=["target", "cell_type", "stimulus", "adt_estimate", "adt_n_studies"]
        )
    meta = meta_analyse(effects)
    return meta[["target", "cell_type", "stimulus", "estimate", "n_studies"]].rename(
        columns={"estimate": "adt_estimate", "n_studies": "adt_n_studies"}
    )


def add_protein_concordance(frame: pd.DataFrame, *, lfc_cutoff: float = 0.5) -> pd.DataFrame:
    """Assign transparent RNA-versus-ADT concordance categories."""

    result = frame.copy()
    categories: list[str] = []
    for row in result.itertuples(index=False):
        adt = row.adt_estimate
        rna = row.inducibility_lfc
        if pd.isna(adt):
            categories.append("UNAVAILABLE")
        elif np.sign(adt) != np.sign(rna) and not (np.isclose(adt, 0.0) or np.isclose(rna, 0.0)):
            categories.append("DISCORDANT")
        elif rna >= lfc_cutoff and adt >= lfc_cutoff:
            categories.append("CONCORDANT")
        elif rna >= lfc_cutoff and adt < lfc_cutoff:
            categories.append("DISCORDANT")
        else:
            categories.append("CONCORDANT")
    result["protein_concordance"] = categories
    return result


def footprint_by_target(
    audit_rows: pd.DataFrame,
    *,
    positive_fraction_cutoff: float = 0.10,
) -> pd.DataFrame:
    """Count the normal-cell-type footprint before and after activation."""

    collapsed = (
        audit_rows.groupby(["target", "cell_type"], observed=True, sort=True)[
            ["positive_fraction_resting", "positive_fraction_activated"]
        ]
        .max()
        .reset_index()
    )
    collapsed["resting_crosses"] = (
        collapsed["positive_fraction_resting"] >= positive_fraction_cutoff
    )
    collapsed["activated_crosses"] = (
        collapsed["positive_fraction_activated"] >= positive_fraction_cutoff
    )
    result = (
        collapsed.groupby("target", observed=True, sort=True)
        .agg(
            footprint_resting=("resting_crosses", "sum"),
            footprint_activated=("activated_crosses", "sum"),
        )
        .reset_index()
    )
    result["footprint_expansion"] = result["footprint_activated"] - result["footprint_resting"]
    return result


def selectivity_erosion(cells: pd.DataFrame) -> pd.DataFrame:
    """Estimate same-study tumour-versus-normal selectivity loss."""

    if "tumour_reference_cpm" not in cells:
        return pd.DataFrame(columns=["target", "selectivity_erosion"])
    bulk = pseudobulk(cells)
    normal = (
        bulk.groupby(["study", "target", "cell_type", "condition"], observed=True)["cpm"]
        .mean()
        .reset_index()
    )
    max_normal = (
        normal.groupby(["study", "target", "condition"], observed=True)["cpm"]
        .max()
        .unstack("condition")
    )
    tumour = (
        cells.groupby(["study", "target"], observed=True)["tumour_reference_cpm"]
        .first()
        .rename("tumour_cpm")
    )
    joined = max_normal.join(tumour, how="inner").dropna()
    if not {"resting", "activated"}.issubset(joined.columns):
        return pd.DataFrame(columns=["target", "selectivity_erosion"])
    joined["selectivity_resting"] = np.log2(joined["tumour_cpm"] + 1.0) - np.log2(
        joined["resting"] + 1.0
    )
    joined["selectivity_activated"] = np.log2(joined["tumour_cpm"] + 1.0) - np.log2(
        joined["activated"] + 1.0
    )
    joined["selectivity_erosion"] = joined["selectivity_resting"] - joined["selectivity_activated"]
    return (
        joined.reset_index()
        .groupby("target", observed=True)[
            ["selectivity_resting", "selectivity_activated", "selectivity_erosion"]
        ]
        .mean()
        .reset_index()
    )
