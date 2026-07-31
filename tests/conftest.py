from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from activation_liability.benchmark import run_benchmark
from activation_liability.synthetic import generate_synthetic

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="session")
def root() -> Path:
    return ROOT


@pytest.fixture(scope="session")
def synthetic_cells(root: Path) -> pd.DataFrame:
    return generate_synthetic(root / "data/controls/controls.yaml")


@pytest.fixture(scope="session")
def benchmark_result(root: Path, synthetic_cells: pd.DataFrame) -> dict[str, object]:
    return run_benchmark(
        synthetic_cells,
        controls_path=root / "data/controls/controls.yaml",
        rules_path=root / "config/evidence_classes.yaml",
    )
