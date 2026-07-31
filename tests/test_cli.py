from __future__ import annotations

from pathlib import Path

import pandas as pd
from typer.testing import CliRunner

from activation_liability.cli import app

runner = CliRunner()


def test_validate_manifest_cli(root: Path) -> None:
    result = runner.invoke(app, ["validate-manifest", str(root / "data/registry")])
    assert result.exit_code == 0, result.output
    assert "Validated 3 manifests" in result.output


def test_build_and_audit_cli_wiring(root: Path, tmp_path: Path, monkeypatch) -> None:
    tiny = pd.DataFrame(
        {
            "study": ["S"],
            "donor": ["D"],
            "cell_id": ["C"],
            "cell_type": ["T_cell"],
            "condition": ["resting"],
            "stimulus": ["lymphocyte_activation"],
            "target": ["X"],
            "rna_count": [1],
            "library_size": [100],
        }
    )
    monkeypatch.setattr("activation_liability.cli.generate_synthetic", lambda *args, **kwargs: tiny)
    cells = tmp_path / "cells.csv"
    built = runner.invoke(
        app,
        [
            "build",
            "--synthetic",
            "--controls",
            str(root / "data/controls/controls.yaml"),
            "--output",
            str(cells),
        ],
    )
    assert built.exit_code == 0, built.output
    assert cells.exists()

    fake_payload = {
        "target_summary": [{"target": "X", "score": 1.0}],
        "audit_rows": [{"target": "X", "inducibility_lfc": 1.0}],
        "claims": {"unsupported_claims": ["toxicity prediction"]},
    }
    monkeypatch.setattr("activation_liability.cli.run_audit", lambda *args, **kwargs: fake_payload)
    output = tmp_path / "audit"
    audited = runner.invoke(
        app,
        [
            "audit",
            "--input",
            str(cells),
            "--rules",
            str(root / "config/evidence_classes.yaml"),
            "--output",
            str(output),
        ],
    )
    assert audited.exit_code == 0, audited.output
    assert (output / "audit.json").exists()
    assert (output / "claims.json").exists()


def test_tissue_validate_cli_wiring(root: Path, tmp_path: Path, monkeypatch) -> None:
    crohn = tmp_path / "crohn.tar"
    crohn.write_bytes(b"tar")
    psoriasis = tmp_path / "psoriasis"
    psoriasis.mkdir()
    pseudobulk = tmp_path / "pb.csv.gz"
    pd.DataFrame({"x": [1]}).to_csv(pseudobulk, index=False, compression="gzip")
    qc = tmp_path / "qc.json"
    qc.write_text("{}", encoding="utf-8")
    baseline = tmp_path / "baseline.json"
    baseline.write_text("{}", encoding="utf-8")

    captured: dict[str, object] = {}

    def fake_pipeline(**kwargs):
        captured.update(kwargs)
        return {"tissue_benchmark": {"tier1_metrics": {"status": "COMPUTED", "auroc": 0.8}}}

    monkeypatch.setattr("activation_liability.cli.run_tissue_pipeline", fake_pipeline)
    output = tmp_path / "result"
    result = runner.invoke(
        app,
        [
            "tissue-validate",
            "--crohn-tar",
            str(crohn),
            "--psoriasis-directory",
            str(psoriasis),
            "--psoriasis-pseudobulk",
            str(pseudobulk),
            "--psoriasis-qc",
            str(qc),
            "--output",
            str(output),
            "--controls",
            str(root / "data/controls/controls.yaml"),
            "--tissue-config",
            str(root / "config/tissue_extension.yaml"),
            "--baseline-benchmark",
            str(baseline),
        ],
    )
    assert result.exit_code == 0, result.output
    assert captured["crohn_tar"] == crohn
    assert captured["psoriasis_pseudobulk_path"] == pseudobulk
    assert '"auroc": 0.8' in result.output
