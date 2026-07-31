"""End-to-end paired tissue-inflammation extension."""

from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd
from scipy import stats

from activation_liability.benchmark import load_controls
from activation_liability.io import write_json
from activation_liability.tissue import (
    TissueBuildResult,
    build_gse134809_tissue,
    build_gse228421_tissue_directory,
    load_tissue_config,
    negative_binomial_sensitivity,
    paired_tissue_effects,
    ranking_metrics,
    score_tissue_controls,
    tissue_footprint,
)


def _canonicalize(value: Any, *, digits: int = 12) -> Any:
    if isinstance(value, float):
        if not math.isfinite(value):
            return None
        return round(value, digits)
    if isinstance(value, dict):
        return {str(key): _canonicalize(item, digits=digits) for key, item in value.items()}
    if isinstance(value, list):
        return [_canonicalize(item, digits=digits) for item in value]
    return value


def _records(frame: pd.DataFrame) -> list[dict[str, Any]]:
    clean = frame.replace({np.nan: None, np.inf: None, -np.inf: None})
    return cast(list[dict[str, Any]], clean.to_dict(orient="records"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(16 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _source_checksums(
    crohn_tar: Path,
    psoriasis_directory: Path,
    *,
    controls_path: Path,
    tissue_config_path: Path,
    baseline_benchmark_path: Path,
) -> dict[str, Any]:
    psoriasis_files = sorted(
        list(psoriasis_directory.glob("*_V1_L.matrix.mtx.gz"))
        + list(psoriasis_directory.glob("*_V1_NL.matrix.mtx.gz"))
        + list(psoriasis_directory.glob("*_V1_L.features.tsv.gz"))
        + list(psoriasis_directory.glob("*_V1_NL.features.tsv.gz"))
    )
    return {
        "GSE134809_RAW.tar": {
            "size": crohn_tar.stat().st_size,
            "sha256": _sha256(crohn_tar),
        },
        "GSE228421_baseline_members": [
            {"name": path.name, "size": path.stat().st_size, "sha256": _sha256(path)}
            for path in psoriasis_files
        ],
        "analysis_inputs": [
            {
                "name": controls_path.name,
                "role": "externally curated benchmark controls",
                "size": controls_path.stat().st_size,
                "sha256": _sha256(controls_path),
            },
            {
                "name": tissue_config_path.name,
                "role": "frozen tissue analysis configuration",
                "size": tissue_config_path.stat().st_size,
                "sha256": _sha256(tissue_config_path),
            },
            {
                "name": baseline_benchmark_path.name,
                "role": "v0.4 coverage reference only",
                "size": baseline_benchmark_path.stat().st_size,
                "sha256": _sha256(baseline_benchmark_path),
            },
        ],
    }


def _tier_metrics(scored: pd.DataFrame, *, tier_max: int) -> dict[str, object]:
    use = scored[(scored["tier"] <= tier_max) & scored["observable"]].copy()
    metrics = ranking_metrics(use)
    if metrics["status"] == "COMPUTED":
        positives = use.loc[use["y_true"] == 1, "score"].to_numpy(dtype=float)
        negatives = use.loc[use["y_true"] == 0, "score"].to_numpy(dtype=float)
        if len(positives) and len(negatives):
            metrics["mann_whitney_p_one_sided"] = float(
                stats.mannwhitneyu(
                    positives, negatives, alternative="greater", method="exact"
                ).pvalue
            )
    metrics["targets"] = sorted(use["target"].astype(str))
    return metrics


def _donor_robustness(pseudobulk: pd.DataFrame, scored: pd.DataFrame) -> pd.DataFrame:
    data = pseudobulk.copy()
    data["log2_cpm"] = np.log2(data["value_sum"] / data["library_sum"] * 1_000_000.0 + 1.0)
    rows: list[dict[str, object]] = []
    for item in scored[scored["observable"]].itertuples(index=False):
        selected = data[
            data["study"].eq(item.driver_study)
            & data["target"].eq(item.target)
            & data["cell_type"].eq(item.driver_cell_type)
            & data["stimulus"].eq(item.driver_stimulus)
        ]
        wide = selected.pivot(index="donor", columns="condition", values="log2_cpm")
        if not {"resting", "activated"}.issubset(wide.columns):
            continue
        differences = (
            wide.dropna(subset=["resting", "activated"])["activated"]
            - wide.dropna(subset=["resting", "activated"])["resting"]
        ).to_numpy(dtype=float)
        if differences.size == 0:
            continue
        loo = [
            float(np.delete(differences, index).mean())
            for index in range(len(differences))
            if len(differences) > 1
        ]
        if np.allclose(differences, 0.0):
            wilcoxon = 1.0
        else:
            wilcoxon = float(
                stats.wilcoxon(differences, alternative="greater", method="auto").pvalue
            )
        rows.append(
            {
                "target": str(item.target),
                "label": str(item.label),
                "study": str(item.driver_study),
                "cell_type": str(item.driver_cell_type),
                "n_pairs": int(len(differences)),
                "mean_difference": float(differences.mean()),
                "median_difference": float(np.median(differences)),
                "fraction_donors_positive": float((differences > 0).mean()),
                "wilcoxon_p_one_sided": wilcoxon,
                "leave_one_donor_out_min_mean": min(loo) if loo else None,
                "leave_one_donor_out_max_mean": max(loo) if loo else None,
                "direction_stable_leave_one_donor_out": bool(loo and min(loo) > 0),
            }
        )
    return pd.DataFrame.from_records(rows)


def _nb_concordance(audit: pd.DataFrame, nb: pd.DataFrame) -> dict[str, object]:
    merged = audit.merge(
        nb,
        on=["study", "target", "cell_type", "stimulus", "n_pairs"],
        how="left",
    )
    valid = merged[
        merged["nb_converged"].eq(True)
        & merged["nb_perfect_separation_warning"].eq(False)
        & merged["nb_lfc"].notna()
    ].copy()
    if len(valid) < 3:
        metrics: dict[str, object] = {
            "status": "INSUFFICIENT_STABLE_FITS",
            "n_valid": int(len(valid)),
        }
    else:
        metrics = {
            "status": "COMPUTED",
            "n_valid": int(len(valid)),
            "spearman_effect_correlation": float(
                stats.spearmanr(valid["inducibility_lfc"], valid["nb_lfc"]).statistic
            ),
            "direction_agreement": float(
                (np.sign(valid["inducibility_lfc"]) == np.sign(valid["nb_lfc"])).mean()
            ),
            "both_estimators_positive_fraction": float(
                ((valid["inducibility_lfc"] > 0) & (valid["nb_lfc"] > 0)).mean()
            ),
        }
    metrics["fit_accounting"] = {
        "attempted": int(len(nb)),
        "converged": int(nb["nb_converged"].sum()) if not nb.empty else 0,
        "perfect_separation_warnings": int(nb["nb_perfect_separation_warning"].sum())
        if not nb.empty
        else 0,
    }
    return metrics


def _crohn_exclusion_sensitivity(
    primary: pd.DataFrame,
    all_eleven: pd.DataFrame,
) -> dict[str, object]:
    keys = ["target", "cell_type", "stimulus"]
    joined = primary.merge(
        all_eleven,
        on=keys,
        suffixes=("_primary9", "_all11"),
        how="inner",
    )
    return {
        "n_matched_effects": int(len(joined)),
        "spearman_effect_correlation": float(
            stats.spearmanr(
                joined["inducibility_lfc_primary9"],
                joined["inducibility_lfc_all11"],
            ).statistic
        ),
        "direction_agreement": float(
            (
                np.sign(joined["inducibility_lfc_primary9"])
                == np.sign(joined["inducibility_lfc_all11"])
            ).mean()
        ),
        "median_absolute_effect_change": float(
            np.median(
                np.abs(joined["inducibility_lfc_primary9"] - joined["inducibility_lfc_all11"])
            )
        ),
        "maximum_absolute_effect_change": float(
            np.max(np.abs(joined["inducibility_lfc_primary9"] - joined["inducibility_lfc_all11"]))
        ),
    }


def _threshold_sensitivity(
    audit: pd.DataFrame,
    controls: pd.DataFrame,
    config: dict[str, Any],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for observability in (0.01, 0.05, 0.10):
        varied = copy.deepcopy(config)
        varied["score"]["observability_fraction_cutoff"] = observability
        scored = score_tissue_controls(audit, controls, varied)
        metrics = _tier_metrics(scored, tier_max=1)
        rows.append(
            {
                "observability_fraction_cutoff": observability,
                "observable_tier1": int(((scored["tier"] == 1) & scored["observable"]).sum()),
                "metrics": metrics,
            }
        )
    return rows


def _coverage_gain(scored: pd.DataFrame, baseline_benchmark: dict[str, Any]) -> dict[str, object]:
    baseline = pd.DataFrame(baseline_benchmark["target_scores"])
    compare = baseline[["target", "observable", "abstention_reason"]].merge(
        scored[["target", "observable", "abstention_reason"]],
        on="target",
        suffixes=("_v0_4", "_tissue"),
    )
    gained = compare[~compare["observable_v0_4"] & compare["observable_tissue"]]
    lost = compare[compare["observable_v0_4"] & ~compare["observable_tissue"]]
    return {
        "newly_observable_targets": sorted(gained["target"].astype(str)),
        "newly_observable_count": int(len(gained)),
        "baseline_observable_but_not_tissue_observable": sorted(lost["target"].astype(str)),
        "note": (
            "Tissue coverage is an extension, not a replacement; loss here means a target was "
            "not observable in the tissue cohorts, not that baseline evidence was invalidated."
        ),
    }


def _lineage_composition(qc: list[dict[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for cohort in qc:
        for sample in cast(list[dict[str, object]], cohort["sample_qc"]):
            counts = cast(dict[str, int], sample["lineage_counts"])
            denominator = sum(counts.values())
            for lineage, count in sorted(counts.items()):
                rows.append(
                    {
                        "accession": str(cohort["accession"]),
                        "donor": str(sample["donor"]),
                        "condition": str(sample["condition"]),
                        "lineage": lineage,
                        "n_cells": int(count),
                        "fraction": count / denominator if denominator else 0.0,
                    }
                )
    return rows


def _write_html_report(payload: dict[str, Any], path: Path) -> None:
    primary = cast(dict[str, Any], payload["tissue_benchmark"]["tier1_metrics"])
    expanded = cast(dict[str, Any], payload["tissue_benchmark"]["tier1_tier2_metrics"])
    footprint = cast(dict[str, Any], payload["footprint_summary"])
    gained = cast(dict[str, Any], payload["coverage_gain"])["newly_observable_targets"]
    gained_text = ", ".join(cast(list[str], gained)) or "none"
    lines = [
        "<!doctype html>",
        '<html lang="en"><head><meta charset="utf-8">',
        "<title>activation-liability v0.5.0</title>",
        "<style>",
        "body{font-family:Arial,sans-serif;max-width:1050px;",
        "margin:40px auto;line-height:1.5}",
        "table{border-collapse:collapse}",
        "th,td{border:1px solid #bbb;padding:7px}",
        "code{background:#eee;padding:2px 4px}",
        "</style></head><body>",
        "<h1>activation-liability v0.5.0 — paired tissue extension</h1>",
        f"<p><strong>Status:</strong> {payload['status']}</p>",
        "<p>This report detects inflammation-associated expansion of ",
        "normal-tissue target expression. It does not predict toxicity.</p>",
        "<h2>Primary tissue benchmark</h2>",
        "<table><tr><th>Set</th><th>N</th><th>Positive</th>",
        "<th>Negative</th><th>AUROC</th><th>Average precision</th></tr>",
        "<tr><td>Tier 1</td>",
        f"<td>{primary.get('n')}</td>",
        f"<td>{primary.get('n_positive')}</td>",
        f"<td>{primary.get('n_negative')}</td>",
        f"<td>{primary.get('auroc')}</td>",
        f"<td>{primary.get('average_precision')}</td></tr>",
        "<tr><td>Tier 1+2</td>",
        f"<td>{expanded.get('n')}</td>",
        f"<td>{expanded.get('n_positive')}</td>",
        f"<td>{expanded.get('n_negative')}</td>",
        f"<td>{expanded.get('auroc')}</td>",
        f"<td>{expanded.get('average_precision')}</td></tr></table>",
        "<h2>Coverage and footprint</h2>",
        f"<p>Newly observable relative to v0.4.0: {gained_text}.</p>",
        "<p>Targets with positive footprint expansion: ",
        f"{footprint.get('targets_with_positive_expansion')}.</p>",
        "<h2>Guardrails</h2><ul>",
        "<li>Crohn primary analysis excludes two patients ",
        "pre-specified by the publication's QC.</li>",
        "<li>Psoriasis uses baseline lesional versus baseline ",
        "non-lesional pairs only.</li>",
        "<li>RNA evidence is not surface-protein corroboration.</li>",
        "<li>Negative-binomial results are sensitivity analyses, ",
        "not a replacement primary estimator.</li>",
        "</ul></body></html>",
    ]
    path.write_text("".join(lines), encoding="utf-8")


def _write_result_manifest(output_root: Path) -> None:
    files = [
        {
            "file_name": path.name,
            "size": path.stat().st_size,
            "sha256": _sha256(path),
        }
        for path in sorted(output_root.iterdir())
        if path.is_file() and path.name != "tissue_result_manifest.json"
    ]
    write_json(
        output_root / "tissue_result_manifest.json",
        {
            "schema_version": 1,
            "software_version": "0.5.0",
            "files": files,
        },
    )


def run_tissue_pipeline(
    *,
    crohn_tar: Path,
    psoriasis_directory: Path,
    psoriasis_pseudobulk_path: Path | None = None,
    psoriasis_qc_path: Path | None = None,
    controls_path: Path,
    tissue_config_path: Path,
    baseline_benchmark_path: Path,
    output_root: Path,
) -> dict[str, Any]:
    """Execute the frozen v0.5 paired tissue extension and all declared sensitivities."""

    config = load_tissue_config(tissue_config_path)
    controls = load_controls(controls_path)
    targets = sorted(controls["target"].astype(str))

    crohn_primary = build_gse134809_tissue(
        crohn_tar,
        targets=targets,
        config=config,
    )
    crohn_all = build_gse134809_tissue(
        crohn_tar,
        targets=targets,
        config=config,
        include_sensitivity_patients=True,
    )
    if psoriasis_pseudobulk_path is not None or psoriasis_qc_path is not None:
        if psoriasis_pseudobulk_path is None or psoriasis_qc_path is None:
            raise ValueError("psoriasis pseudobulk and QC paths must be supplied together")
        psoriasis = TissueBuildResult(
            pseudobulk=pd.read_csv(psoriasis_pseudobulk_path),
            qc=cast(
                dict[str, object],
                json.loads(psoriasis_qc_path.read_text(encoding="utf-8")),
            ),
        )
    else:
        psoriasis = build_gse228421_tissue_directory(
            psoriasis_directory,
            targets=targets,
            config=config,
        )
    primary_pb = pd.concat(
        [crohn_primary.pseudobulk, psoriasis.pseudobulk], ignore_index=True
    ).sort_values(["study", "target", "cell_type", "donor", "condition"])
    audit = paired_tissue_effects(primary_pb)
    crohn_all_audit = paired_tissue_effects(crohn_all.pseudobulk)
    nb = negative_binomial_sensitivity(primary_pb)
    scored = score_tissue_controls(audit, controls, config)
    footprint = tissue_footprint(
        audit,
        cutoff=float(config["score"]["positive_fraction_cutoff"]),
    )
    robustness = _donor_robustness(primary_pb, scored)
    baseline_benchmark = json.loads(baseline_benchmark_path.read_text(encoding="utf-8"))

    per_study: dict[str, object] = {}
    for study in sorted(audit["study"].unique()):
        study_scored = score_tissue_controls(audit[audit["study"] == study], controls, config)
        per_study[str(study)] = {
            "tier1_metrics": _tier_metrics(study_scored, tier_max=1),
            "observable_tier1": int(
                ((study_scored["tier"] == 1) & study_scored["observable"]).sum()
            ),
        }

    positive_expansion = footprint[footprint["tissue_footprint_expansion"] > 0]
    payload: dict[str, Any] = {
        "schema_version": 1,
        "status": "COMPUTED_REAL_PAIRED_TISSUE_EXTENSION",
        "method_version": "tissue-v1",
        "pre_registration": {
            "config": str(tissue_config_path),
            "status": config["status"],
            "crohn_pairing_source": (
                "effiken/martin_et_al_cell_2019 input/tables/sample_index.csv "
                "at commit d1ff0f9099a9017552d7b9b8d582bcc5a9314ae2"
            ),
            "annotation_qc_repair": (
                "The initial skin marker-only QC showed keratinocytes misclassified as plasma "
                "through SDC1/XBP1. Those markers were replaced by multi-marker plasma gates and "
                "KRT1/KRT10 before the final tissue benchmark was generated. No benchmark target "
                "or control label was used in that repair."
            ),
        },
        "cohorts": {
            "GSE134809": {
                "primary_pairs": 9,
                "sensitivity_pairs": 11,
                "tissue": "ileum",
                "contrast": "involved versus uninvolved within patient",
                "epithelial_claims": "UNSUPPORTED_EXCLUDED_BY_ORIGINAL_DISSOCIATION_DESIGN",
            },
            "GSE228421": {
                "primary_pairs": 5,
                "tissue": "skin",
                "contrast": "baseline lesional versus baseline non-lesional within patient",
                "post_treatment_visits": "EXCLUDED_FROM_PRIMARY",
            },
        },
        "tissue_benchmark": {
            "tier1_metrics": _tier_metrics(scored[scored["tier"] == 1], tier_max=1),
            "tier1_tier2_metrics": _tier_metrics(scored, tier_max=2),
            "per_study": per_study,
            "target_scores": _records(scored),
        },
        "coverage_gain": _coverage_gain(scored, baseline_benchmark),
        "footprint_summary": {
            "targets_with_positive_expansion": int(len(positive_expansion)),
            "maximum_expansion": int(footprint["tissue_footprint_expansion"].max()),
            "rows": _records(footprint),
        },
        "crohn_published_exclusion_sensitivity": _crohn_exclusion_sensitivity(
            audit[audit["study"] == "GSE134809"], crohn_all_audit
        ),
        "negative_binomial_sensitivity": _nb_concordance(audit, nb),
        "donor_robustness": {
            "n_rows": int(len(robustness)),
            "direction_stable_count": int(robustness["direction_stable_leave_one_donor_out"].sum())
            if not robustness.empty
            else 0,
            "rows": _records(robustness),
        },
        "threshold_sensitivity": _threshold_sensitivity(audit, controls, config),
        "annotation_composition": _lineage_composition([crohn_primary.qc, psoriasis.qc]),
        "protein_corroboration": {
            "status": "UNAVAILABLE_FOR_TISSUE_TARGETS",
            "reason": "Both tissue cohorts provide RNA counts only.",
        },
        "claims": {
            "permitted": [
                (
                    "The analysed paired cohorts contain normal-cell lineages whose RNA "
                    "footprint expands during tissue inflammation."
                ),
                (
                    "The tissue extension increases biological lineage coverage relative "
                    "to the PBMC-only benchmark."
                ),
            ],
            "conditional": [
                (
                    "Target-level tissue inducibility claims require the named cohort, "
                    "broad lineage, paired-donor coverage and confidence interval reported "
                    "in tissue_audit_rows.csv."
                ),
            ],
            "unsupported": [
                "Clinical toxicity prediction",
                "Surface-protein induction for tissue targets",
                "Delivery, exposure or therapeutic-index prediction",
                "Epithelial conclusions from GSE134809",
            ],
        },
    }
    payload = cast(dict[str, Any], _canonicalize(payload))

    output_root.mkdir(parents=True, exist_ok=True)
    primary_pb.to_csv(
        output_root / "tissue_pseudobulk.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    audit.to_csv(output_root / "tissue_audit_rows.csv", index=False)
    nb.to_csv(output_root / "tissue_negative_binomial.csv", index=False)
    scored.to_csv(output_root / "tissue_target_scores.csv", index=False)
    footprint.to_csv(output_root / "tissue_footprint.csv", index=False)
    robustness.to_csv(output_root / "tissue_donor_robustness.csv", index=False)
    crohn_all_audit.to_csv(output_root / "crohn_all11_sensitivity_audit_rows.csv", index=False)
    write_json(
        output_root / "tissue_ingestion_qc.json",
        {
            "GSE134809_primary": crohn_primary.qc,
            "GSE134809_all11": crohn_all.qc,
            "GSE228421": psoriasis.qc,
        },
    )
    write_json(output_root / "tissue_benchmark.json", payload)
    write_json(output_root / "tissue_claims.json", payload["claims"])
    write_json(
        output_root / "tissue_source_checksums.json",
        _source_checksums(
            crohn_tar,
            psoriasis_directory,
            controls_path=controls_path,
            tissue_config_path=tissue_config_path,
            baseline_benchmark_path=baseline_benchmark_path,
        ),
    )
    write_json(
        output_root / "tissue_execution_summary.json",
        {
            "schema_version": 1,
            "software_version": "0.5.0",
            "status": payload["status"],
            "cohorts": payload["cohorts"],
            "tier1_metrics": payload["tissue_benchmark"]["tier1_metrics"],
            "tier1_tier2_metrics": payload["tissue_benchmark"]["tier1_tier2_metrics"],
            "coverage_gain": payload["coverage_gain"],
            "footprint_summary": {
                "targets_with_positive_expansion": payload["footprint_summary"][
                    "targets_with_positive_expansion"
                ],
                "maximum_expansion": payload["footprint_summary"]["maximum_expansion"],
            },
            "negative_binomial_sensitivity": payload["negative_binomial_sensitivity"],
            "crohn_exclusion_sensitivity": payload["crohn_published_exclusion_sensitivity"],
            "protein_corroboration": payload["protein_corroboration"],
        },
    )
    _write_html_report(payload, output_root / "tissue_report.html")
    _write_result_manifest(output_root)
    return payload
