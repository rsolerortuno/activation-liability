from __future__ import annotations

from pathlib import Path

from activation_liability.report import render_report


def test_report_renders_status_and_boundary(tmp_path: Path) -> None:
    output = tmp_path / "report.html"
    render_report(
        {
            "status": "SYNTHETIC",
            "target_scores": [
                {
                    "target": "CD69",
                    "score": 3.2,
                    "tier": 1,
                    "protein_concordance": "CONCORDANT",
                }
            ],
        },
        output,
    )
    text = output.read_text()
    assert "SYNTHETIC" in text
    assert "not a toxicity predictor" in text
    assert "CD69" in text
