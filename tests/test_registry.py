from __future__ import annotations

from pathlib import Path

import pytest
import yaml
from pydantic import ValidationError

from activation_liability.models import StudyManifest, load_manifest
from activation_liability.registry import default_manifests, validate_registry


def test_shipped_registry_valid_but_disabled(root: Path) -> None:
    manifests, errors = validate_registry(root / "data/registry")
    assert not errors
    assert len(manifests) == 3
    assert default_manifests(root / "data/registry") == []
    assert all(not manifest.default for manifest in manifests)


def test_unverified_manifest_cannot_be_default(root: Path, tmp_path: Path) -> None:
    payload = yaml.safe_load((root / "data/registry/gse126030.yaml").read_text())
    payload["default"] = True
    path = tmp_path / "bad.yaml"
    path.write_text(yaml.safe_dump(payload))
    with pytest.raises(ValidationError, match="unverified manifest"):
        load_manifest(path)


def test_verified_manifest_rejects_unresolved_fields() -> None:
    with pytest.raises(ValidationError, match="NOT_VERIFIED"):
        StudyManifest.model_validate(
            {
                "accession": "X",
                "source_repository": "TEST",
                "publication_doi": "NOT_VERIFIED",
                "tissue": "blood",
                "disease_context": "healthy",
                "stimulus_axis": "IFN_II",
                "platform": "test",
                "modality": "RNA",
                "donor_count": 3,
                "condition_labels": {"resting": "r", "activated": "a"},
                "licence": "CC0",
                "verified": True,
                "default": False,
                "downloads": [],
            }
        )
