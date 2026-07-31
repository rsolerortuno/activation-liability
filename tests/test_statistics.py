from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from activation_liability.statistics import (
    benjamini_hochberg,
    meta_analyse,
    paired_effects,
    pseudobulk,
    random_effects_meta,
)


def test_bh_known_values_and_nan() -> None:
    adjusted = benjamini_hochberg([0.01, 0.04, 0.03, np.nan])
    assert adjusted[:3] == pytest.approx([0.03, 0.04, 0.04])
    assert np.isnan(adjusted[3])


def test_random_effects_known_homogeneous_case() -> None:
    result = random_effects_meta(np.array([1.0, 1.0]), np.array([0.1, 0.1]))
    assert result["estimate"] == pytest.approx(1.0)
    assert result["tau2"] == pytest.approx(0.0)
    assert result["i2"] == pytest.approx(0.0)
    assert result["standard_error"] == pytest.approx(np.sqrt(1 / 200))


def test_random_effects_detects_heterogeneity() -> None:
    result = random_effects_meta(np.array([0.0, 2.0]), np.array([0.1, 0.1]))
    assert result["estimate"] == pytest.approx(1.0)
    assert result["tau2"] > 1.0
    assert result["i2"] > 90.0


def test_paired_effect_recovers_exact_signal() -> None:
    rows = []
    for donor in ("D1", "D2", "D3", "D4"):
        for condition, count in (("resting", 100), ("activated", 400)):
            rows.append(
                {
                    "study": "S1",
                    "donor": donor,
                    "cell_id": f"{donor}_{condition}",
                    "cell_type": "T_cell",
                    "condition": condition,
                    "stimulus": "lymphocyte_activation",
                    "target": "X",
                    "rna_count": count,
                    "library_size": 1_000_000,
                }
            )
    cells = pd.DataFrame(rows)
    bulk = pseudobulk(cells)
    assert len(bulk) == 8
    effects = paired_effects(cells)
    assert effects.loc[0, "effect"] == pytest.approx(np.log2(400 + 1) - np.log2(100 + 1))
    assert effects.loc[0, "p_value"] == 0.0
    meta = meta_analyse(effects)
    assert meta.loc[0, "estimate"] == pytest.approx(effects.loc[0, "effect"])


def test_pseudobulk_rejects_missing_columns() -> None:
    with pytest.raises(ValueError, match="missing input columns"):
        pseudobulk(pd.DataFrame({"rna_count": [1]}))


def test_single_study_heterogeneity_is_not_estimable() -> None:
    result = random_effects_meta(np.array([1.0]), np.array([0.2]))
    assert result["estimate"] == pytest.approx(1.0)
    assert np.isnan(result["i2"])
    assert np.isnan(result["tau2"])
