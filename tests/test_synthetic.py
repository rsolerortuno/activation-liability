from __future__ import annotations

from pathlib import Path

import pandas as pd

from activation_liability.statistics import paired_effects
from activation_liability.synthetic import generate_synthetic


def test_synthetic_fixture_is_deterministic(root: Path) -> None:
    first = generate_synthetic(
        root / "data/controls/controls.yaml",
        studies=2,
        donors_per_study=3,
        cells_per_condition=1,
    )
    second = generate_synthetic(
        root / "data/controls/controls.yaml",
        studies=2,
        donors_per_study=3,
        cells_per_condition=1,
    )
    pd.testing.assert_frame_equal(first, second)


def test_planted_signal_and_null_are_numerically_separated(synthetic_cells: pd.DataFrame) -> None:
    subset = synthetic_cells[
        synthetic_cells["target"].isin(["CD69", "CD19"])
        & synthetic_cells["cell_type"].isin(["T_cell", "B_cell"])
    ].copy()
    effects = paired_effects(subset)
    induced = effects[
        (effects["target"] == "CD69")
        & (effects["cell_type"] == "T_cell")
        & (effects["stimulus"] == "lymphocyte_activation")
    ]
    null = effects[
        (effects["target"] == "CD19")
        & (effects["cell_type"] == "B_cell")
        & (effects["stimulus"] == "lymphocyte_activation")
    ]
    assert len(induced) == 3
    assert induced["effect"].mean() > 1.0
    assert abs(null["effect"].mean()) < 0.5
