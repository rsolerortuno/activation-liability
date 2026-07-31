from __future__ import annotations

import hashlib
import io
from pathlib import Path

import pytest
import yaml

from activation_liability.fetch import fetch_default_registry


def _manifest(checksum: str) -> dict[str, object]:
    return {
        "accession": "TEST001",
        "source_repository": "TEST",
        "publication_doi": "10.0000/test",
        "tissue": "blood",
        "disease_context": "healthy",
        "stimulus_axis": "IFN_II",
        "platform": "test",
        "modality": "RNA",
        "donor_count": 3,
        "condition_labels": {"resting": "control", "activated": "stimulated"},
        "licence": "CC0-1.0",
        "verified": True,
        "default": True,
        "downloads": [
            {
                "url": "https://example.org/test.csv",
                "sha256": checksum,
                "filename": "test.csv",
            }
        ],
        "notes": "Synthetic fetch adapter test only.",
    }


def test_fetcher_uses_checksum_and_content_address(tmp_path: Path, monkeypatch) -> None:
    content = b"target,count\nX,1\n"
    checksum = hashlib.sha256(content).hexdigest()
    registry = tmp_path / "registry"
    registry.mkdir()
    (registry / "test.yaml").write_text(yaml.safe_dump(_manifest(checksum)))
    monkeypatch.setattr(
        "activation_liability.fetch.urllib.request.urlopen",
        lambda *args, **kwargs: io.BytesIO(content),
    )
    cache = tmp_path / "cache"
    result = fetch_default_registry(registry, cache)
    assert result["status"] == "OK"
    stored = cache / checksum / "test.csv"
    assert stored.read_bytes() == content
    assert (cache / "fetch_index.json").exists()


def test_fetcher_rejects_checksum_mismatch(tmp_path: Path, monkeypatch) -> None:
    registry = tmp_path / "registry"
    registry.mkdir()
    (registry / "test.yaml").write_text(yaml.safe_dump(_manifest("0" * 64)))
    monkeypatch.setattr(
        "activation_liability.fetch.urllib.request.urlopen",
        lambda *args, **kwargs: io.BytesIO(b"wrong"),
    )
    with pytest.raises(ValueError, match="checksum mismatch"):
        fetch_default_registry(registry, tmp_path / "cache")
