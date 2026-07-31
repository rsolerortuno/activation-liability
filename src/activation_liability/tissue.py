"""Paired tissue-inflammation extension for activation-liability.

This module deliberately operates on donor-level sufficient statistics rather than a long
cell-by-target table. Raw cells are used only inside one sample at a time for QC, condition-blind
broad-lineage annotation, target sums and positive-cell counts. The persisted statistical unit is
therefore donor pseudobulk.
"""

from __future__ import annotations

import gzip
import math
import tarfile
import warnings
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import IO, Any, cast

import numpy as np
import pandas as pd
import yaml
from scipy import sparse, stats
from scipy.io import mmread
from sklearn.metrics import average_precision_score, roc_auc_score
from statsmodels.genmod.families import NegativeBinomial
from statsmodels.genmod.generalized_linear_model import GLM
from statsmodels.tools.sm_exceptions import PerfectSeparationWarning

from activation_liability.statistics import benjamini_hochberg


@dataclass(frozen=True)
class TissueBuildResult:
    """Compact donor-level data and QC produced from a tissue cohort."""

    pseudobulk: pd.DataFrame
    qc: dict[str, object]


CROHN_SAMPLE_MAP: dict[str, tuple[str, str, str]] = {
    "69": ("rp5", "activated", "GSM3972009"),
    "68": ("rp5", "resting", "GSM3972010"),
    "122": ("rp6", "activated", "GSM3972011"),
    "123": ("rp6", "resting", "GSM3972012"),
    "128": ("rp7", "activated", "GSM3972013"),
    "129": ("rp7", "resting", "GSM3972014"),
    "135": ("rp8", "resting", "GSM3972015"),
    "138": ("rp8", "activated", "GSM3972016"),
    "158": ("rp10", "activated", "GSM3972017"),
    "159": ("rp10", "resting", "GSM3972018"),
    "180": ("rp11", "resting", "GSM3972019"),
    "181": ("rp11", "activated", "GSM3972020"),
    "186": ("rp12", "resting", "GSM3972021"),
    "187": ("rp12", "activated", "GSM3972022"),
    "189": ("rp13", "resting", "GSM3972023"),
    "190": ("rp13", "activated", "GSM3972024"),
    "192": ("rp14", "resting", "GSM3972025"),
    "193": ("rp14", "activated", "GSM3972026"),
    "195": ("rp15", "resting", "GSM3972027"),
    "196": ("rp15", "activated", "GSM3972028"),
    "208": ("rp16", "resting", "GSM3972029"),
    "209": ("rp16", "activated", "GSM3972030"),
}

PSORIASIS_BASELINE_MAP: dict[str, tuple[str, str]] = {
    "P1_V1_L": ("P1", "activated"),
    "P1_V1_NL": ("P1", "resting"),
    "P2_V1_L": ("P2", "activated"),
    "P2_V1_NL": ("P2", "resting"),
    "P3_V1_L": ("P3", "activated"),
    "P3_V1_NL": ("P3", "resting"),
    "P4_V1_L": ("P4", "activated"),
    "P4_V1_NL": ("P4", "resting"),
    "P5_V1_L": ("P5", "activated"),
    "P5_V1_NL": ("P5", "resting"),
}


def load_tissue_config(path: Path) -> dict[str, Any]:
    """Load and validate the separately versioned tissue-extension configuration."""

    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("tissue configuration must be a mapping")
    required = {"version", "score", "qc", "endpoints", "annotation", "target_lineages"}
    if required - set(payload):
        raise ValueError("tissue configuration is missing required sections")
    marker_map = payload["annotation"].get("markers")
    if not isinstance(marker_map, dict) or len(marker_map) < 2:
        raise ValueError("tissue annotation requires at least two marker modules")
    return cast(dict[str, Any], payload)


def _read_gzip_lines(handle: IO[bytes]) -> list[str]:
    with gzip.GzipFile(fileobj=handle) as stream:
        return [line.decode("utf-8").rstrip("\n") for line in stream]


def _read_gzip_matrix(handle: IO[bytes]) -> sparse.csr_matrix:
    with gzip.GzipFile(fileobj=handle) as stream:
        return sparse.csr_matrix(mmread(stream))


def _extract_lines(archive: tarfile.TarFile, member: str) -> list[str]:
    handle = archive.extractfile(member)
    if handle is None:
        raise ValueError(f"unable to extract {member}")
    return _read_gzip_lines(handle)


def _extract_matrix(archive: tarfile.TarFile, member: str) -> sparse.csr_matrix:
    handle = archive.extractfile(member)
    if handle is None:
        raise ValueError(f"unable to extract {member}")
    return _read_gzip_matrix(handle)


def _gene_symbols(feature_lines: list[str]) -> list[str]:
    symbols: list[str] = []
    for line in feature_lines:
        fields = line.split("\t")
        symbols.append(fields[1] if len(fields) > 1 else fields[0])
    return symbols


def _first_index(values: list[str]) -> dict[str, int]:
    result: dict[str, int] = {}
    for index, value in enumerate(values):
        result.setdefault(value, index)
    return result


def _fraction_for_prefixes(
    counts: sparse.csr_matrix,
    genes: list[str],
    library: np.ndarray,
    prefixes: tuple[str, ...],
) -> np.ndarray:
    indices = [index for index, gene in enumerate(genes) if gene.startswith(prefixes)]
    if not indices:
        return np.zeros(counts.shape[1], dtype=float)
    total = np.asarray(counts[indices, :].sum(axis=0)).ravel().astype(float)
    return np.asarray(total / np.maximum(library, 1.0), dtype=float)


def condition_blind_lineage_labels(
    counts: sparse.csr_matrix,
    genes: list[str],
    marker_map: Mapping[str, Iterable[str]],
    *,
    min_score: float,
    min_margin: float,
    minimum_positive_markers: Mapping[str, int] | None = None,
) -> tuple[np.ndarray, pd.DataFrame]:
    """Assign broad tissue lineages with fixed marker modules and abstention."""

    if counts.shape[0] != len(genes):
        raise ValueError("gene list length does not match count rows")
    if len(marker_map) < 2:
        raise ValueError("at least two marker modules are required")
    library = np.asarray(counts.sum(axis=0)).ravel().astype(float)
    index = _first_index(genes)
    labels = list(marker_map)
    score_rows: list[np.ndarray] = []
    marker_coverage: dict[str, int] = {}
    for label, marker_values in marker_map.items():
        marker_genes = [str(gene) for gene in marker_values]
        indices = [index[gene] for gene in marker_genes if gene in index]
        marker_coverage[label] = len(indices)
        if not indices:
            score_rows.append(np.zeros(counts.shape[1], dtype=float))
            continue
        selected = counts[indices, :].astype(float)
        scaled = selected.multiply(10_000.0 / np.maximum(library, 1.0))
        score_rows.append(np.asarray(np.log1p(scaled.toarray()).mean(axis=0)).ravel())
    score_matrix = np.vstack(score_rows)
    order = np.argsort(score_matrix, axis=0)
    columns = np.arange(counts.shape[1])
    top_index = order[-1, :]
    second_index = order[-2, :]
    top_score = score_matrix[top_index, columns]
    margin = top_score - score_matrix[second_index, columns]
    assigned = np.asarray(labels, dtype=object)[top_index]
    assigned[(top_score < min_score) | (margin < min_margin)] = "Unknown"
    requirements = minimum_positive_markers or {}
    for label, required_count in requirements.items():
        marker_genes = [str(gene) for gene in marker_map.get(label, ())]
        indices = [index[gene] for gene in marker_genes if gene in index]
        if not indices:
            assigned[assigned == label] = "Unknown"
            continue
        positive_markers = np.asarray((counts[indices, :] > 0).sum(axis=0)).ravel()
        assigned[(assigned == label) & (positive_markers < int(required_count))] = "Unknown"
    diagnostics = pd.DataFrame(score_matrix.T, columns=labels)
    diagnostics["top_score"] = top_score
    diagnostics["margin"] = margin
    diagnostics.attrs["marker_coverage"] = marker_coverage
    return assigned, diagnostics


def _aggregate_sample(
    *,
    archive: tarfile.TarFile,
    prefix: str,
    feature_token: str,
    donor: str,
    condition: str,
    study: str,
    stimulus: str,
    tissue: str,
    targets: list[str],
    marker_map: Mapping[str, Iterable[str]],
    minimum_positive_markers: Mapping[str, int],
    qc_parameters: Mapping[str, float],
    excluded_lineages: set[str],
) -> tuple[pd.DataFrame, dict[str, object]]:
    names = archive.getnames()
    matrix_member = next(
        (name for name in names if prefix in name and name.endswith("matrix.mtx.gz")), None
    )
    feature_member = next(
        (name for name in names if prefix in name and name.endswith(feature_token)), None
    )
    if matrix_member is None or feature_member is None:
        raise ValueError(f"incomplete 10x members for {prefix}")
    matrix = _extract_matrix(archive, matrix_member)
    genes = _gene_symbols(_extract_lines(archive, feature_member))
    if matrix.shape[0] != len(genes):
        raise ValueError(f"10x feature dimensions do not match for {prefix}")

    library = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
    detected = np.asarray((matrix > 0).sum(axis=0)).ravel().astype(int)
    mitochondrial_fraction = _fraction_for_prefixes(matrix, genes, library, ("MT-", "mt-"))
    haemoglobin_fraction = _fraction_for_prefixes(
        matrix,
        genes,
        library,
        ("HBA", "HBB", "HBD", "HBG"),
    )
    keep = (
        (library >= float(qc_parameters["minimum_library_size"]))
        & (detected >= int(qc_parameters["minimum_detected_genes"]))
        & (mitochondrial_fraction <= float(qc_parameters["maximum_mitochondrial_fraction"]))
        & (haemoglobin_fraction <= float(qc_parameters["maximum_haemoglobin_fraction"]))
    )
    selected = np.flatnonzero(keep)
    if selected.size == 0:
        raise ValueError(f"no cells pass QC for {prefix}")
    qc_matrix = matrix[:, selected]
    qc_library = library[selected]
    labels, annotation = condition_blind_lineage_labels(
        qc_matrix,
        genes,
        marker_map,
        min_score=float(qc_parameters["annotation_min_score"]),
        min_margin=float(qc_parameters["annotation_min_margin"]),
        minimum_positive_markers=minimum_positive_markers,
    )
    label_keep = (labels != "Unknown") & ~np.isin(labels, sorted(excluded_lineages))
    gene_index = _first_index(genes)
    available_targets = [target for target in targets if target in gene_index]
    target_indices = [gene_index[target] for target in available_targets]
    target_matrix = qc_matrix[target_indices, :]
    rows: list[dict[str, object]] = []
    for cell_type in sorted(set(labels[label_keep].astype(str))):
        cells = np.flatnonzero(label_keep & (labels == cell_type))
        if cells.size == 0:
            continue
        selected_targets = target_matrix[:, cells]
        value_sum = np.asarray(selected_targets.sum(axis=1)).ravel().astype(float)
        positive_count = np.asarray((selected_targets > 0).sum(axis=1)).ravel().astype(int)
        library_sum = float(qc_library[cells].sum())
        for target, count_sum, n_positive in zip(
            available_targets,
            value_sum,
            positive_count,
            strict=True,
        ):
            rows.append(
                {
                    "study": study,
                    "stimulus": stimulus,
                    "tissue": tissue,
                    "donor": donor,
                    "condition": condition,
                    "cell_type": cell_type,
                    "target": target,
                    "value_sum": float(count_sum),
                    "library_sum": library_sum,
                    "positive_count": int(n_positive),
                    "n_cells": int(cells.size),
                }
            )
    counts_by_label = pd.Series(labels).value_counts(dropna=False).sort_index()
    sample_qc: dict[str, object] = {
        "sample": prefix,
        "donor": donor,
        "condition": condition,
        "raw_barcodes": int(matrix.shape[1]),
        "qc_pass_cells": int(selected.size),
        "retained_annotated_cells": int(label_keep.sum()),
        "median_library_qc": float(np.median(qc_library)),
        "median_detected_genes_qc": float(np.median(detected[selected])),
        "median_mitochondrial_fraction_qc": float(np.median(mitochondrial_fraction[selected])),
        "median_haemoglobin_fraction_qc": float(np.median(haemoglobin_fraction[selected])),
        "lineage_counts": {str(key): int(value) for key, value in counts_by_label.items()},
        "marker_coverage": annotation.attrs.get("marker_coverage", {}),
        "available_targets": len(available_targets),
    }
    return pd.DataFrame.from_records(rows), sample_qc


def build_gse134809_tissue(
    raw_tar: Path,
    *,
    targets: list[str],
    config: Mapping[str, Any],
    include_sensitivity_patients: bool = False,
) -> TissueBuildResult:
    """Build paired Crohn ileum pseudobulk, excluding published QC failures by default."""

    endpoint = next(item for item in config["endpoints"] if str(item["study"]) == "GSE134809")
    primary = {str(value) for value in endpoint["primary_patients"]}
    sensitivity = {str(value) for value in endpoint.get("sensitivity_only_patients", [])}
    allowed = primary | sensitivity if include_sensitivity_patients else primary
    marker_map = cast(Mapping[str, Iterable[str]], config["annotation"]["markers"])
    minimum_positive_markers = cast(
        Mapping[str, int], config["annotation"].get("minimum_positive_markers", {})
    )
    qc_parameters = cast(Mapping[str, float], config["qc"])
    excluded = {str(value) for value in endpoint.get("excluded_lineages", [])}
    frames: list[pd.DataFrame] = []
    sample_qc: list[dict[str, object]] = []
    with tarfile.open(raw_tar, mode="r") as archive:
        for sample_number, (donor, condition, gsm) in CROHN_SAMPLE_MAP.items():
            if donor not in allowed:
                continue
            frame, qc = _aggregate_sample(
                archive=archive,
                prefix=f"{gsm}_{sample_number}_",
                feature_token="genes.tsv.gz",
                donor=donor,
                condition=condition,
                study="GSE134809",
                stimulus="tissue_inflammation_crohns",
                tissue="ileum",
                targets=targets,
                marker_map=marker_map,
                minimum_positive_markers=minimum_positive_markers,
                qc_parameters=qc_parameters,
                excluded_lineages=excluded,
            )
            frames.append(frame)
            sample_qc.append(qc)
    combined = pd.concat(frames, ignore_index=True)
    _validate_complete_pairs(combined, expected_donors=allowed)
    return TissueBuildResult(
        pseudobulk=combined.sort_values(
            ["study", "target", "cell_type", "donor", "condition"]
        ).reset_index(drop=True),
        qc={
            "accession": "GSE134809",
            "analysis_set": "all_11_sensitivity" if include_sensitivity_patients else "primary_9",
            "included_donors": sorted(allowed),
            "published_qc_exclusions": sorted(sensitivity),
            "sample_qc": sample_qc,
        },
    )


def build_gse228421_tissue(
    raw_tar: Path,
    *,
    targets: list[str],
    config: Mapping[str, Any],
) -> TissueBuildResult:
    """Build five paired baseline psoriasis contrasts; post-treatment visits are excluded."""

    endpoint = next(item for item in config["endpoints"] if str(item["study"]) == "GSE228421")
    allowed = {str(value) for value in endpoint["primary_patients"]}
    marker_map = cast(Mapping[str, Iterable[str]], config["annotation"]["markers"])
    minimum_positive_markers = cast(
        Mapping[str, int], config["annotation"].get("minimum_positive_markers", {})
    )
    qc_parameters = cast(Mapping[str, float], config["qc"])
    excluded = {str(value) for value in endpoint.get("excluded_lineages", [])}
    frames: list[pd.DataFrame] = []
    sample_qc: list[dict[str, object]] = []
    with tarfile.open(raw_tar, mode="r") as archive:
        names = archive.getnames()
        for prefix, (donor, condition) in PSORIASIS_BASELINE_MAP.items():
            matching = [name for name in names if prefix in name and name.endswith("matrix.mtx.gz")]
            if len(matching) != 1:
                raise ValueError(f"expected one baseline matrix for {prefix}, found {matching}")
            gsm_prefix = matching[0].removesuffix("matrix.mtx.gz")
            frame, qc = _aggregate_sample(
                archive=archive,
                prefix=gsm_prefix,
                feature_token="features.tsv.gz",
                donor=donor,
                condition=condition,
                study="GSE228421",
                stimulus="tissue_inflammation_psoriasis",
                tissue="skin",
                targets=targets,
                marker_map=marker_map,
                minimum_positive_markers=minimum_positive_markers,
                qc_parameters=qc_parameters,
                excluded_lineages=excluded,
            )
            frames.append(frame)
            sample_qc.append(qc)
    combined = pd.concat(frames, ignore_index=True)
    _validate_complete_pairs(combined, expected_donors=allowed)
    return TissueBuildResult(
        pseudobulk=combined.sort_values(
            ["study", "target", "cell_type", "donor", "condition"]
        ).reset_index(drop=True),
        qc={
            "accession": "GSE228421",
            "analysis_set": "baseline_V1_paired_only",
            "included_donors": sorted(allowed),
            "excluded_visits": ["V2_day3_lesional", "V3_day14_lesional"],
            "sample_qc": sample_qc,
        },
    )


def _validate_complete_pairs(frame: pd.DataFrame, *, expected_donors: set[str]) -> None:
    observed = set(frame["donor"].astype(str))
    if observed != expected_donors:
        raise ValueError(
            f"donor set mismatch: expected {sorted(expected_donors)}, found {sorted(observed)}"
        )
    donor_conditions = frame.groupby("donor", observed=True)["condition"].agg(lambda x: set(x))
    incomplete = [
        donor for donor, values in donor_conditions.items() if values != {"resting", "activated"}
    ]
    if incomplete:
        raise ValueError(f"incomplete paired donors: {incomplete}")


def paired_tissue_effects(pseudobulk: pd.DataFrame, *, min_pairs: int = 3) -> pd.DataFrame:
    """Compute paired log2-CPM effects from compact donor-level tissue data."""

    required = {
        "study",
        "target",
        "cell_type",
        "stimulus",
        "donor",
        "condition",
        "value_sum",
        "library_sum",
        "positive_count",
        "n_cells",
    }
    missing = required - set(pseudobulk)
    if missing:
        raise ValueError(f"missing tissue pseudobulk columns: {sorted(missing)}")
    data = pseudobulk.copy()
    data["cpm"] = data["value_sum"] / data["library_sum"] * 1_000_000.0
    data["log2_cpm"] = np.log2(data["cpm"] + 1.0)
    keys = ["study", "target", "cell_type", "stimulus"]
    rows: list[dict[str, object]] = []
    for key, group in data.groupby(keys, observed=True, sort=True):
        wide = group.pivot(index="donor", columns="condition", values="log2_cpm")
        if not {"resting", "activated"}.issubset(wide.columns):
            continue
        complete = wide.dropna(subset=["resting", "activated"])
        if len(complete) < min_pairs:
            continue
        differences = (complete["activated"] - complete["resting"]).to_numpy(dtype=float)
        effect = float(differences.mean())
        sd = float(np.std(differences, ddof=1))
        standard_error = max(sd / math.sqrt(len(differences)), 1e-8)
        p_value = (
            1.0
            if np.isclose(sd, 0.0) and np.isclose(effect, 0.0)
            else 0.0
            if np.isclose(sd, 0.0)
            else float(stats.ttest_1samp(differences, popmean=0.0).pvalue)
        )
        positive = group.assign(fraction=group["positive_count"] / group["n_cells"])
        fractions = positive.groupby("condition", observed=True)["fraction"].mean()
        row = dict(zip(keys, key, strict=True))
        row.update(
            {
                "inducibility_lfc": effect,
                "inducibility_se": standard_error,
                "ci_low": effect - 1.96 * standard_error,
                "ci_high": effect + 1.96 * standard_error,
                "p_value": p_value,
                "n_pairs": int(len(complete)),
                "total_pairs": int(len(complete)),
                "positive_fraction_resting": float(fractions.get("resting", 0.0)),
                "positive_fraction_activated": float(fractions.get("activated", 0.0)),
                "fraction_donors_positive": float((differences > 0).mean()),
                "i2": float("nan"),
                "tau2": float("nan"),
                "n_studies": 1,
            }
        )
        rows.append(row)
    result = pd.DataFrame.from_records(rows)
    if not result.empty:
        result["q_value"] = benjamini_hochberg(result["p_value"])
    return result


def _mom_alpha(counts: np.ndarray, libraries: np.ndarray) -> float:
    median_library = float(np.median(libraries))
    scaled = counts / np.maximum(libraries, 1.0) * median_library
    mean = float(np.mean(scaled))
    variance = float(np.var(scaled, ddof=1)) if len(scaled) > 1 else mean
    if mean <= 0:
        return 1e-8
    return max((variance - mean) / (mean**2), 1e-8)


def negative_binomial_sensitivity(
    pseudobulk: pd.DataFrame,
    *,
    min_pairs: int = 3,
) -> pd.DataFrame:
    """Fit paired NB GLMs with donor fixed effects and a log-library offset."""

    keys = ["study", "target", "cell_type", "stimulus"]
    rows: list[dict[str, object]] = []
    for key, group in pseudobulk.groupby(keys, observed=True, sort=True):
        condition_counts = group.groupby("donor", observed=True)["condition"].nunique()
        donors = sorted(str(value) for value in condition_counts[condition_counts == 2].index)
        use = group[group["donor"].astype(str).isin(donors)].copy()
        if len(donors) < min_pairs or use["value_sum"].sum() <= 0:
            continue
        use["activated"] = (use["condition"] == "activated").astype(float)
        donor_dummies = pd.get_dummies(use["donor"].astype(str), prefix="donor", drop_first=True)
        design = pd.concat(
            [
                pd.Series(1.0, index=use.index, name="intercept"),
                use[["activated"]],
                donor_dummies.astype(float),
            ],
            axis=1,
        )
        alpha = _mom_alpha(
            use["value_sum"].to_numpy(dtype=float),
            use["library_sum"].to_numpy(dtype=float),
        )
        perfect_separation = False
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always", PerfectSeparationWarning)
                fit = GLM(
                    use["value_sum"].to_numpy(dtype=float),
                    design.to_numpy(dtype=float),
                    family=NegativeBinomial(alpha=alpha),
                    offset=np.log(use["library_sum"].to_numpy(dtype=float)),
                ).fit(maxiter=200, disp=0)
            perfect_separation = any(
                issubclass(item.category, PerfectSeparationWarning) for item in caught
            )
            coefficient = float(fit.params[1])
            standard_error = float(fit.bse[1])
            p_value = float(fit.pvalues[1])
            converged = bool(getattr(fit, "converged", True))
        except (ValueError, np.linalg.LinAlgError, FloatingPointError):
            coefficient = float("nan")
            standard_error = float("nan")
            p_value = float("nan")
            converged = False
        row = dict(zip(keys, key, strict=True))
        row.update(
            {
                "nb_lfc": coefficient / math.log(2.0),
                "nb_standard_error": standard_error / math.log(2.0),
                "nb_p_value": p_value,
                "nb_alpha": alpha,
                "nb_converged": converged,
                "nb_perfect_separation_warning": perfect_separation,
                "n_pairs": len(donors),
            }
        )
        rows.append(row)
    result = pd.DataFrame.from_records(rows)
    if not result.empty:
        result["nb_q_value"] = benjamini_hochberg(result["nb_p_value"])
    return result


def tissue_footprint(audit_rows: pd.DataFrame, *, cutoff: float = 0.10) -> pd.DataFrame:
    """Count broad tissue lineages positive before and during inflammation."""

    collapsed = (
        audit_rows.groupby(["target", "cell_type"], observed=True, sort=True)[
            ["positive_fraction_resting", "positive_fraction_activated"]
        ]
        .max()
        .reset_index()
    )
    collapsed["resting_crosses"] = collapsed["positive_fraction_resting"] >= cutoff
    collapsed["activated_crosses"] = collapsed["positive_fraction_activated"] >= cutoff
    footprint = (
        collapsed.groupby("target", observed=True, sort=True)
        .agg(
            tissue_footprint_resting=("resting_crosses", "sum"),
            tissue_footprint_activated=("activated_crosses", "sum"),
        )
        .reset_index()
    )
    footprint["tissue_footprint_expansion"] = (
        footprint["tissue_footprint_activated"] - footprint["tissue_footprint_resting"]
    )
    return footprint


def score_tissue_controls(
    audit_rows: pd.DataFrame,
    controls: pd.DataFrame,
    config: Mapping[str, Any],
) -> pd.DataFrame:
    """Score tissue controls using only pre-declared lineages and a one-sided 95% LCB."""

    lineage_map = cast(Mapping[str, Iterable[str]], config["target_lineages"])
    z_value = float(config["score"]["z_value"])
    observability_cutoff = float(config["score"]["observability_fraction_cutoff"])
    observed_lineages = set(audit_rows["cell_type"].astype(str))
    rows: list[dict[str, object]] = []
    for control in controls.itertuples(index=False):
        target = str(control.target)
        expected = [str(value) for value in lineage_map.get(target, [])]
        covered = sorted(set(expected) & observed_lineages)
        candidates = audit_rows[
            audit_rows["target"].eq(target) & audit_rows["cell_type"].isin(expected)
        ].copy()
        candidates["lcb95"] = (
            candidates["inducibility_lfc"] - z_value * candidates["inducibility_se"]
        )
        max_fraction = (
            float(
                candidates[["positive_fraction_resting", "positive_fraction_activated"]]
                .max(axis=1)
                .max()
            )
            if not candidates.empty
            else float("nan")
        )
        if not expected:
            reason = "NO_PREDECLARED_LINEAGE"
        elif not covered:
            reason = "LINEAGE_NOT_COVERED"
        elif not np.isfinite(max_fraction) or max_fraction < observability_cutoff:
            reason = "TARGET_NOT_OBSERVABLE_AT_COVERAGE_THRESHOLD"
        else:
            reason = None
        best = candidates.loc[candidates["lcb95"].idxmax()] if not candidates.empty else None
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
                "driver_study": None if best is None else str(best["study"]),
                "driver_stimulus": None if best is None else str(best["stimulus"]),
                "driver_cell_type": None if best is None else str(best["cell_type"]),
                "driver_effect": None if best is None else float(best["inducibility_lfc"]),
                "driver_standard_error": None if best is None else float(best["inducibility_se"]),
                "driver_total_pairs": None if best is None else int(best["total_pairs"]),
            }
        )
    return pd.DataFrame.from_records(rows)


def ranking_metrics(scored: pd.DataFrame) -> dict[str, object]:
    """Return ranking metrics only when both classes are represented."""

    use = scored[scored["observable"]].dropna(subset=["score"]).copy()
    if use.empty or use["y_true"].nunique() < 2:
        return {
            "status": "INSUFFICIENT_CLASS_COVERAGE",
            "n": int(len(use)),
            "n_positive": int((use["y_true"] == 1).sum()),
            "n_negative": int((use["y_true"] == 0).sum()),
            "auroc": None,
            "average_precision": None,
        }
    return {
        "status": "COMPUTED",
        "n": int(len(use)),
        "n_positive": int((use["y_true"] == 1).sum()),
        "n_negative": int((use["y_true"] == 0).sum()),
        "auroc": float(roc_auc_score(use["y_true"], use["score"])),
        "average_precision": float(average_precision_score(use["y_true"], use["score"])),
    }


def _aggregate_sample_files(
    *,
    matrix_path: Path,
    feature_path: Path,
    donor: str,
    condition: str,
    study: str,
    stimulus: str,
    tissue: str,
    targets: list[str],
    marker_map: Mapping[str, Iterable[str]],
    minimum_positive_markers: Mapping[str, int],
    qc_parameters: Mapping[str, float],
    excluded_lineages: set[str],
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Aggregate one extracted 10x matrix without reading unused raw barcode strings."""

    with gzip.open(matrix_path, "rb") as matrix_stream:
        matrix = sparse.csr_matrix(mmread(matrix_stream))
    with gzip.open(feature_path, "rt", encoding="utf-8") as feature_stream:
        genes = _gene_symbols([line.rstrip("\n") for line in feature_stream])
    if matrix.shape[0] != len(genes):
        raise ValueError(f"10x feature dimensions do not match for {matrix_path.name}")

    library = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
    detected = np.asarray((matrix > 0).sum(axis=0)).ravel().astype(int)
    mitochondrial_fraction = _fraction_for_prefixes(matrix, genes, library, ("MT-", "mt-"))
    haemoglobin_fraction = _fraction_for_prefixes(
        matrix,
        genes,
        library,
        ("HBA", "HBB", "HBD", "HBG"),
    )
    keep = (
        (library >= float(qc_parameters["minimum_library_size"]))
        & (detected >= int(qc_parameters["minimum_detected_genes"]))
        & (mitochondrial_fraction <= float(qc_parameters["maximum_mitochondrial_fraction"]))
        & (haemoglobin_fraction <= float(qc_parameters["maximum_haemoglobin_fraction"]))
    )
    selected = np.flatnonzero(keep)
    if selected.size == 0:
        raise ValueError(f"no cells pass QC for {matrix_path.name}")
    qc_matrix = matrix[:, selected]
    qc_library = library[selected]
    labels, annotation = condition_blind_lineage_labels(
        qc_matrix,
        genes,
        marker_map,
        min_score=float(qc_parameters["annotation_min_score"]),
        min_margin=float(qc_parameters["annotation_min_margin"]),
        minimum_positive_markers=minimum_positive_markers,
    )
    label_keep = (labels != "Unknown") & ~np.isin(labels, sorted(excluded_lineages))
    gene_index = _first_index(genes)
    available_targets = [target for target in targets if target in gene_index]
    target_matrix = qc_matrix[[gene_index[target] for target in available_targets], :]
    rows: list[dict[str, object]] = []
    for cell_type in sorted(set(labels[label_keep].astype(str))):
        cells = np.flatnonzero(label_keep & (labels == cell_type))
        selected_targets = target_matrix[:, cells]
        value_sum = np.asarray(selected_targets.sum(axis=1)).ravel().astype(float)
        positive_count = np.asarray((selected_targets > 0).sum(axis=1)).ravel().astype(int)
        library_sum = float(qc_library[cells].sum())
        for target, count_sum, n_positive in zip(
            available_targets,
            value_sum,
            positive_count,
            strict=True,
        ):
            rows.append(
                {
                    "study": study,
                    "stimulus": stimulus,
                    "tissue": tissue,
                    "donor": donor,
                    "condition": condition,
                    "cell_type": cell_type,
                    "target": target,
                    "value_sum": float(count_sum),
                    "library_sum": library_sum,
                    "positive_count": int(n_positive),
                    "n_cells": int(cells.size),
                }
            )
    counts_by_label = pd.Series(labels).value_counts(dropna=False).sort_index()
    sample_qc: dict[str, object] = {
        "sample": matrix_path.name.removesuffix(".matrix.mtx.gz"),
        "donor": donor,
        "condition": condition,
        "raw_barcodes": int(matrix.shape[1]),
        "nonempty_barcodes": int((library > 0).sum()),
        "qc_pass_cells": int(selected.size),
        "retained_annotated_cells": int(label_keep.sum()),
        "median_library_qc": float(np.median(qc_library)),
        "median_detected_genes_qc": float(np.median(detected[selected])),
        "median_mitochondrial_fraction_qc": float(np.median(mitochondrial_fraction[selected])),
        "median_haemoglobin_fraction_qc": float(np.median(haemoglobin_fraction[selected])),
        "lineage_counts": {str(key): int(value) for key, value in counts_by_label.items()},
        "marker_coverage": annotation.attrs.get("marker_coverage", {}),
        "available_targets": len(available_targets),
    }
    return pd.DataFrame.from_records(rows), sample_qc


def build_gse228421_tissue_directory(
    directory: Path,
    *,
    targets: list[str],
    config: Mapping[str, Any],
) -> TissueBuildResult:
    """Build psoriasis from extracted baseline members for efficient reproducible execution."""

    endpoint = next(item for item in config["endpoints"] if str(item["study"]) == "GSE228421")
    allowed = {str(value) for value in endpoint["primary_patients"]}
    marker_map = cast(Mapping[str, Iterable[str]], config["annotation"]["markers"])
    minimum_positive_markers = cast(
        Mapping[str, int], config["annotation"].get("minimum_positive_markers", {})
    )
    qc_parameters = cast(Mapping[str, float], config["qc"])
    excluded = {str(value) for value in endpoint.get("excluded_lineages", [])}
    frames: list[pd.DataFrame] = []
    sample_qc: list[dict[str, object]] = []
    for prefix, (donor, condition) in PSORIASIS_BASELINE_MAP.items():
        matrices = sorted(directory.glob(f"*_{prefix}.matrix.mtx.gz"))
        features = sorted(directory.glob(f"*_{prefix}.features.tsv.gz"))
        if len(matrices) != 1 or len(features) != 1:
            raise ValueError(
                f"expected one extracted matrix/features pair for {prefix}; "
                f"found {matrices} and {features}"
            )
        frame, qc = _aggregate_sample_files(
            matrix_path=matrices[0],
            feature_path=features[0],
            donor=donor,
            condition=condition,
            study="GSE228421",
            stimulus="tissue_inflammation_psoriasis",
            tissue="skin",
            targets=targets,
            marker_map=marker_map,
            minimum_positive_markers=minimum_positive_markers,
            qc_parameters=qc_parameters,
            excluded_lineages=excluded,
        )
        frames.append(frame)
        sample_qc.append(qc)
    combined = pd.concat(frames, ignore_index=True)
    _validate_complete_pairs(combined, expected_donors=allowed)
    return TissueBuildResult(
        pseudobulk=combined.sort_values(
            ["study", "target", "cell_type", "donor", "condition"]
        ).reset_index(drop=True),
        qc={
            "accession": "GSE228421",
            "analysis_set": "baseline_V1_paired_only",
            "input_mode": "extracted_members",
            "included_donors": sorted(allowed),
            "excluded_visits": ["V2_day3_lesional", "V3_day14_lesional"],
            "sample_qc": sample_qc,
        },
    )
