"""Machine-readable claim contract derived from an audit run."""

from __future__ import annotations

from typing import Any

import pandas as pd


def build_claims(audit_rows: pd.DataFrame, target_summary: pd.DataFrame) -> dict[str, Any]:
    """Generate permitted, conditional and unsupported claims from observed coverage."""

    n_targets = int(target_summary["target"].nunique()) if not target_summary.empty else 0
    protein_targets = sorted(
        audit_rows.loc[audit_rows["protein_concordance"] == "CONCORDANT", "target"].unique()
    )
    sufficient_targets = sorted(
        target_summary.loc[target_summary["evidence_class"] != "INSUFFICIENT", "target"].unique()
    )
    abstention_rate = (
        float((target_summary["evidence_class"] == "INSUFFICIENT").mean()) if n_targets else 1.0
    )
    return {
        "schema_version": 1,
        "permitted_claims": [
            {
                "claim": (
                    "Within-study paired RNA inducibility was estimated for covered "
                    "target-cell-stimulus combinations."
                ),
                "condition": (
                    "Valid only for the studies, donors, cell types and stimuli recorded "
                    "in audit.json."
                ),
            },
            {
                "claim": (
                    "Activation expanded or contracted the detected normal-cell footprint "
                    "for evaluated targets."
                ),
                "condition": "Depends on the stated detection and positive-fraction thresholds.",
            },
        ],
        "conditional_claims": [
            {
                "claim": "Matched surface-protein induction corroborates RNA induction.",
                "targets": protein_targets,
                "condition": "Only for target-cell-stimulus rows labelled CONCORDANT.",
            },
            {
                "claim": "There is sufficient evidence for an activation-liability call.",
                "targets": sufficient_targets,
                "condition": "Means evidence class A, B or C; it does not mean clinical toxicity.",
            },
        ],
        "unsupported_claims": [
            "The tool predicts clinical toxicity.",
            "The tool establishes a therapeutic window.",
            "A non-induced target is safe.",
            "RNA induction proves accessible surface protein.",
            "Synthetic benchmark performance is real-data validation.",
        ],
        "coverage": {
            "targets": n_targets,
            "protein_corroborated_targets": len(protein_targets),
            "abstention_rate": abstention_rate,
        },
    }
