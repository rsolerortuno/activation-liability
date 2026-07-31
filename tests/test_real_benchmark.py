from __future__ import annotations

from pathlib import Path

import pandas as pd

from activation_liability.benchmark import load_controls
from activation_liability.real_benchmark import (
    load_real_benchmark_config,
    run_gse96583_confirmation,
    score_real_controls,
)


def _row(
    target: str,
    cell_type: str,
    stimulus: str,
    effect: float,
    se: float,
    fraction: float,
) -> dict[str, object]:
    return {
        "target": target,
        "cell_type": cell_type,
        "stimulus": stimulus,
        "inducibility_lfc": effect,
        "inducibility_se": se,
        "positive_fraction_resting": fraction / 2,
        "positive_fraction_activated": fraction,
        "total_pairs": 4,
    }


def test_real_scoring_uses_predeclared_lineage_and_lcb(root: Path) -> None:
    controls = load_controls(root / "data/controls/controls.yaml")
    config = load_real_benchmark_config(root / "config/real_benchmark.yaml")
    rows = pd.DataFrame(
        [
            _row("IL2RA", "T_cell", "lymphocyte_activation_24h", 2.0, 0.25, 0.8),
            _row("IL2RA", "Myeloid", "IFN_I_18h", 10.0, 0.01, 0.9),
            _row("CD19", "B_cell", "IFN_II_6h", 0.1, 0.10, 0.8),
            _row("CLDN18", "Myeloid", "IFN_I_18h", 8.0, 0.01, 0.9),
        ]
    )
    scored = score_real_controls(rows, controls, config).set_index("target")
    assert scored.loc["IL2RA", "driver_cell_type"] == "T_cell"
    assert scored.loc["IL2RA", "score"] == 2.0 - 1.96 * 0.25
    assert bool(scored.loc["CD19", "observable"])
    assert scored.loc["CLDN18", "abstention_reason"] == "LINEAGE_NOT_COVERED"


def test_real_scoring_abstains_below_observability(root: Path) -> None:
    controls = load_controls(root / "data/controls/controls.yaml")
    config = load_real_benchmark_config(root / "config/real_benchmark.yaml")
    rows = pd.DataFrame([_row("FAP", "Myeloid", "IFN_I_18h", 3.0, 0.1, 0.01)])
    scored = score_real_controls(rows, controls, config).set_index("target")
    assert not bool(scored.loc["FAP", "observable"])
    assert scored.loc["FAP", "abstention_reason"] == "TARGET_NOT_OBSERVABLE_AT_COVERAGE_THRESHOLD"


def test_run_real_benchmark_end_to_end(tmp_path: Path) -> None:
    import yaml

    from activation_liability.real_benchmark import run_real_benchmark

    controls_path = tmp_path / "controls.yaml"
    controls_path.write_text(
        yaml.safe_dump(
            {
                "controls": [
                    {
                        "target": "IL2RA",
                        "label": "positive",
                        "tier": 1,
                        "citation": "external-a",
                        "rationale": "positive control",
                    },
                    {
                        "target": "CD274",
                        "label": "positive",
                        "tier": 1,
                        "citation": "external-b",
                        "rationale": "positive control",
                    },
                    {
                        "target": "CD19",
                        "label": "negative",
                        "tier": 1,
                        "citation": "external-c",
                        "rationale": "negative control",
                    },
                    {
                        "target": "CD4",
                        "label": "negative",
                        "tier": 1,
                        "citation": "external-d",
                        "rationale": "negative control",
                    },
                ]
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    config_path = tmp_path / "real.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "version": 1,
                "score": {
                    "definition": "test LCB",
                    "z_value": 1.96,
                    "observability_fraction_cutoff": 0.05,
                },
                "endpoints": [
                    {
                        "study": "GSE157857",
                        "stimulus": "IFN_I_18h",
                        "axis": "IFN-I",
                        "time_hours": 18,
                        "role": "primary",
                    },
                    {
                        "study": "GSE178429_IFN_6h",
                        "stimulus": "IFN_II_6h",
                        "axis": "IFN-II",
                        "time_hours": 6,
                        "role": "primary",
                    },
                    {
                        "study": "GSE178429_LPS_6h",
                        "stimulus": "TLR_6h",
                        "axis": "TLR",
                        "time_hours": 6,
                        "role": "primary",
                    },
                    {
                        "study": "GSE140244_24h",
                        "stimulus": "lymphocyte_activation_24h",
                        "axis": "activation",
                        "time_hours": 24,
                        "role": "primary",
                    },
                ],
                "secondary_timecourse": {
                    "stimuli": [
                        "lymphocyte_activation_2h",
                        "lymphocyte_activation_24h",
                    ]
                },
                "annotation_leakage_control": {
                    "method": "leave-target-out",
                    "marker_targets": [],
                },
                "lineage_mapping_basis": "external test mapping",
                "target_lineages": {
                    "IL2RA": ["T_cell"],
                    "CD274": ["Myeloid"],
                    "CD19": ["B_cell"],
                    "CD4": ["T_cell"],
                },
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    audit_rows = pd.DataFrame(
        [
            _row("IL2RA", "T_cell", "lymphocyte_activation_2h", 1.0, 0.1, 0.8),
            _row("IL2RA", "T_cell", "lymphocyte_activation_24h", 3.0, 0.1, 0.9),
            _row("CD274", "Myeloid", "IFN_I_18h", 2.5, 0.1, 0.8),
            _row("CD19", "B_cell", "IFN_II_6h", 0.1, 0.1, 0.8),
            _row("CD4", "T_cell", "lymphocyte_activation_2h", 0.0, 0.1, 0.9),
            _row("CD4", "T_cell", "lymphocyte_activation_24h", -0.2, 0.1, 0.9),
        ]
    )
    cells_rows: list[dict[str, object]] = []
    drivers = [
        ("IL2RA", "GSE140244_24h", "lymphocyte_activation_24h", "T_cell", 10, 80),
        ("CD274", "GSE157857", "IFN_I_18h", "Myeloid", 5, 50),
        ("CD19", "GSE178429_IFN_6h", "IFN_II_6h", "B_cell", 50, 50),
        ("CD4", "GSE140244_24h", "lymphocyte_activation_24h", "T_cell", 80, 70),
    ]
    for target, study, stimulus, cell_type, resting, activated in drivers:
        for donor in ("D1", "D2", "D3", "D4"):
            for condition, count in (("resting", resting), ("activated", activated)):
                cells_rows.append(
                    {
                        "study": study,
                        "donor": donor,
                        "cell_id": f"{study}_{target}_{donor}_{condition}",
                        "cell_type": cell_type,
                        "condition": condition,
                        "stimulus": stimulus,
                        "target": target,
                        "rna_count": count,
                        "library_size": 1000,
                        "adt_count": float("nan"),
                        "adt_library_size": float("nan"),
                    }
                )
    payload = run_real_benchmark(
        pd.DataFrame(cells_rows),
        audit_rows,
        controls_path=controls_path,
        config_path=config_path,
    )
    assert payload["primary_tier1_metrics"]["auroc"] == 1.0
    assert payload["protein_benchmark"]["status"] == "UNAVAILABLE_NO_CONTROL_TARGET_OVERLAP"
    assert payload["holdout"]["holdout_metrics"]["n_positive"] == 1
    assert payload["donor_robustness"]["rows"]
    assert payload["timecourse_validation"]["rows"]


def test_gse96583_confirmation_uses_frozen_cutoff_and_rules(tmp_path: Path) -> None:
    import yaml

    controls_path = tmp_path / "controls.yaml"
    controls_path.write_text(
        yaml.safe_dump(
            {
                "controls": [
                    {
                        "target": "IL2RA",
                        "label": "positive",
                        "tier": 1,
                        "citation": "external-a",
                        "rationale": "positive control",
                    },
                    {
                        "target": "CD274",
                        "label": "positive",
                        "tier": 1,
                        "citation": "external-b",
                        "rationale": "positive control",
                    },
                    {
                        "target": "CD19",
                        "label": "negative",
                        "tier": 1,
                        "citation": "external-c",
                        "rationale": "negative control",
                    },
                    {
                        "target": "CD4",
                        "label": "negative",
                        "tier": 1,
                        "citation": "external-d",
                        "rationale": "negative control",
                    },
                ]
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    config_path = tmp_path / "real.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "version": 1,
                "score": {
                    "definition": "test LCB",
                    "z_value": 1.96,
                    "observability_fraction_cutoff": 0.05,
                },
                "endpoints": [
                    {
                        "study": "BASE_IFN",
                        "stimulus": "IFN_I_18h",
                        "axis": "IFN-I",
                        "time_hours": 18,
                        "role": "primary",
                    },
                    {
                        "study": "BASE_T",
                        "stimulus": "lymphocyte_activation_24h",
                        "axis": "activation",
                        "time_hours": 24,
                        "role": "primary",
                    },
                ],
                "secondary_timecourse": {"stimuli": ["lymphocyte_activation_24h"]},
                "annotation_leakage_control": {
                    "method": "leave-target-out",
                    "marker_targets": [],
                },
                "lineage_mapping_basis": "external test mapping",
                "target_lineages": {
                    "IL2RA": ["T_cell"],
                    "CD274": ["Myeloid"],
                    "CD19": ["B_cell"],
                    "CD4": ["T_cell"],
                },
            },
            sort_keys=False,
        ),
        encoding="utf-8",
    )
    baseline_audit_rows = pd.DataFrame(
        [
            _row("IL2RA", "T_cell", "lymphocyte_activation_24h", 3.0, 0.1, 0.8),
            _row("CD274", "Myeloid", "IFN_I_18h", 2.5, 0.1, 0.8),
            _row("CD19", "B_cell", "IFN_I_18h", 0.0, 0.1, 0.8),
            _row("CD4", "T_cell", "lymphocyte_activation_24h", -0.2, 0.1, 0.8),
        ]
    )
    baseline_cells: list[dict[str, object]] = []
    for target, study, stimulus, cell_type, resting, activated in (
        ("IL2RA", "BASE_T", "lymphocyte_activation_24h", "T_cell", 10, 80),
        ("CD274", "BASE_IFN", "IFN_I_18h", "Myeloid", 5, 50),
        ("CD19", "BASE_IFN", "IFN_I_18h", "B_cell", 50, 50),
        ("CD4", "BASE_T", "lymphocyte_activation_24h", "T_cell", 80, 70),
    ):
        for donor in ("D1", "D2", "D3", "D4"):
            for condition, count in (("resting", resting), ("activated", activated)):
                baseline_cells.append(
                    {
                        "study": study,
                        "donor": donor,
                        "cell_id": f"{study}_{target}_{donor}_{condition}",
                        "cell_type": cell_type,
                        "condition": condition,
                        "stimulus": stimulus,
                        "target": target,
                        "rna_count": count,
                        "library_size": 1000,
                        "adt_count": float("nan"),
                        "adt_library_size": float("nan"),
                    }
                )
    from activation_liability.real_benchmark import run_real_benchmark

    baseline = run_real_benchmark(
        pd.DataFrame(baseline_cells),
        baseline_audit_rows,
        controls_path=controls_path,
        config_path=config_path,
    )
    confirmation_audit_rows = pd.DataFrame(
        [
            _row("IL2RA", "T_cell", "IFN_I_6h_confirmation", 2.2, 0.2, 0.8),
            _row("CD274", "Myeloid", "IFN_I_6h_confirmation", 2.0, 0.2, 0.8),
            _row("CD19", "B_cell", "IFN_I_6h_confirmation", -0.1, 0.1, 0.8),
            _row("CD4", "T_cell", "IFN_I_6h_confirmation", -0.3, 0.1, 0.8),
        ]
    )
    confirmation_cells: list[dict[str, object]] = []
    for target, cell_type, resting, activated in (
        ("IL2RA", "T_cell", 10, 60),
        ("CD274", "Myeloid", 5, 45),
        ("CD19", "B_cell", 50, 45),
        ("CD4", "T_cell", 80, 70),
    ):
        for donor in ("D1", "D2", "D3", "D4"):
            for condition, count in (("resting", resting), ("activated", activated)):
                confirmation_cells.append(
                    {
                        "study": "GSE96583",
                        "donor": donor,
                        "cell_id": f"GSE96583_{target}_{donor}_{condition}",
                        "cell_type": cell_type,
                        "condition": condition,
                        "stimulus": "IFN_I_6h_confirmation",
                        "target": target,
                        "rna_count": count,
                        "library_size": 1000,
                        "adt_count": float("nan"),
                        "adt_library_size": float("nan"),
                    }
                )
    payload = run_gse96583_confirmation(
        pd.DataFrame(confirmation_cells),
        confirmation_audit_rows,
        baseline_benchmark=baseline,
        baseline_audit_rows=baseline_audit_rows,
        controls_path=controls_path,
        config_path=config_path,
    )
    assert payload["status"].startswith("PARTIAL_EXTERNAL_CONFIRMATION")
    assert payload["tier1_ranking_metrics"]["auroc"] == 1.0
    assert payload["frozen_cutoff_classification_tier1"]["specificity"] == 1.0
    assert "No endpoint" in payload["protocol"]["no_tuning_statement"]
    assert payload["target_holdout"]["status"] == "EVALUABLE"
