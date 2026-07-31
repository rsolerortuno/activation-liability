from __future__ import annotations

from pathlib import Path

import yaml


def test_controls_are_complete_and_not_registry_derived(root: Path) -> None:
    payload = yaml.safe_load((root / "data/controls/controls.yaml").read_text())
    controls = payload["controls"]
    assert len(controls) == 38
    assert {entry["label"] for entry in controls} == {"positive", "negative"}
    assert all(entry["tier"] in {1, 2} for entry in controls)
    registry_accessions = {
        yaml.safe_load(path.read_text())["accession"]
        for path in (root / "data/registry").glob("*.yaml")
    }
    provenance = " ".join(f"{entry['citation']} {entry['rationale']}" for entry in controls)
    assert all(accession not in provenance for accession in registry_accessions)
    assert all("pubmed.ncbi.nlm.nih.gov" in entry["citation"] for entry in controls)
