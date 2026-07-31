"""Versioned evidence-class and abstention rule engine."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

CLASS_ORDER = ("A", "B", "C")


def load_evidence_rules(path: Path) -> dict[str, Any]:
    """Load a versioned YAML evidence rule file."""

    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if (
        not isinstance(payload, dict)
        or "classes" not in payload
        or "hard_abstention" not in payload
    ):
        raise ValueError("invalid evidence rule file")
    return payload


def classify_evidence(
    record: dict[str, Any],
    rules: dict[str, Any],
    *,
    protein_claim_requested: bool = False,
) -> tuple[str, list[str]]:
    """Apply hard abstention then deterministic A/B/C thresholds."""

    hard = rules["hard_abstention"]
    reasons: list[str] = []
    n_studies = int(record["n_studies"])
    i2 = float(record["i2"])
    if int(record["total_pairs"]) < int(hard["min_paired_donors"]):
        reasons.append("FEWER_THAN_MINIMUM_PAIRED_DONORS")
    if n_studies >= 2 and i2 > float(hard["max_i2"]):
        reasons.append("EXCESSIVE_HETEROGENEITY")
    if bool(hard["undetected_both_states"]) and (
        float(record["positive_fraction_resting"]) == 0.0
        and float(record["positive_fraction_activated"]) == 0.0
    ):
        reasons.append("TARGET_UNDETECTED")
    if bool(hard["require_relevant_stimulus"]) and not bool(record["relevant_stimulus"]):
        reasons.append("NO_RELEVANT_STIMULUS_COVERAGE")
    if protein_claim_requested and record["protein_concordance"] != "CONCORDANT":
        reasons.append("NO_PROTEIN_CORROBORATION_FOR_PROTEIN_CLAIM")
    if reasons:
        return "INSUFFICIENT", reasons

    for class_name in CLASS_ORDER:
        threshold = rules["classes"][class_name]
        protein_rule = threshold["protein_concordance"]
        matches_protein = protein_rule == "ANY" or record["protein_concordance"] == protein_rule
        heterogeneity_matches = n_studies < 2 or i2 <= float(threshold["max_i2"])
        if (
            n_studies >= int(threshold["min_studies"])
            and int(record["total_pairs"]) >= int(threshold["min_paired_donors"])
            and heterogeneity_matches
            and float(record["inducibility_lfc"]) >= float(threshold["min_inducibility_lfc"])
            and float(record["positive_fraction_activated"])
            >= float(threshold["min_positive_fraction_activated"])
            and matches_protein
        ):
            return class_name, []
    return "INSUFFICIENT", ["NO_EVIDENCE_CLASS_THRESHOLD_MET"]
