from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from activation_liability.audit import run_audit, validate_input
from activation_liability.evidence import classify_evidence, load_evidence_rules


def test_audit_outputs_claim_contract(root: Path, synthetic_cells: pd.DataFrame) -> None:
    result = run_audit(
        synthetic_cells,
        rules_path=root / "config/evidence_classes.yaml",
    )
    assert result["status"] == "SYNTHETIC"
    assert len(result["target_summary"]) == 38
    assert result["claims"]["coverage"]["targets"] == 38
    unsupported = " ".join(result["claims"]["unsupported_claims"]).lower()
    assert "predicts clinical toxicity" in unsupported
    assert "synthetic benchmark" in unsupported


def test_surface_claim_abstains_without_protein(root: Path) -> None:
    rules = load_evidence_rules(root / "config/evidence_classes.yaml")
    record = {
        "total_pairs": 10,
        "i2": 0.0,
        "positive_fraction_resting": 0.0,
        "positive_fraction_activated": 0.5,
        "relevant_stimulus": True,
        "protein_concordance": "UNAVAILABLE",
        "n_studies": 2,
        "inducibility_lfc": 2.0,
    }
    evidence, reasons = classify_evidence(record, rules, protein_claim_requested=True)
    assert evidence == "INSUFFICIENT"
    assert reasons == ["NO_PROTEIN_CORROBORATION_FOR_PROTEIN_CLAIM"]


def test_undetected_target_is_normal_abstention(root: Path) -> None:
    rules = load_evidence_rules(root / "config/evidence_classes.yaml")
    record = {
        "total_pairs": 4,
        "i2": 0.0,
        "positive_fraction_resting": 0.0,
        "positive_fraction_activated": 0.0,
        "relevant_stimulus": True,
        "protein_concordance": "UNAVAILABLE",
        "n_studies": 1,
        "inducibility_lfc": 0.0,
    }
    evidence, reasons = classify_evidence(record, rules)
    assert evidence == "INSUFFICIENT"
    assert "TARGET_UNDETECTED" in reasons


def test_validate_input_rejects_bad_condition() -> None:
    frame = pd.DataFrame(
        {
            "study": ["S"],
            "donor": ["D"],
            "cell_id": ["C"],
            "cell_type": ["T_cell"],
            "condition": ["case"],
            "stimulus": ["x"],
            "target": ["X"],
            "rna_count": [1],
            "library_size": [100],
        }
    )
    with pytest.raises(ValueError, match="condition"):
        validate_input(frame)
