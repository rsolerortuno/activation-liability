"""Checksum-pinned direct download support for verified manifests."""

from __future__ import annotations

import shutil
import urllib.request
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from activation_liability.io import sha256_file, write_json
from activation_liability.registry import default_manifests


def fetch_default_registry(registry: Path, cache_dir: Path) -> dict[str, Any]:
    """Fetch only verified/default direct downloads into a content-addressed cache."""

    manifests = default_manifests(registry)
    cache_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, str]] = []
    for manifest in manifests:
        for item in manifest.downloads:
            temporary = cache_dir / f".{item.filename}.partial"
            with (
                urllib.request.urlopen(str(item.url), timeout=120) as response,
                temporary.open("wb") as output,
            ):
                shutil.copyfileobj(response, output)
            actual = sha256_file(temporary)
            if actual.lower() != item.sha256.lower():
                temporary.unlink(missing_ok=True)
                raise ValueError(
                    f"checksum mismatch for {manifest.accession}/{item.filename}: {actual}"
                )
            destination = cache_dir / actual / item.filename
            destination.parent.mkdir(parents=True, exist_ok=True)
            temporary.replace(destination)
            records.append(
                {
                    "accession": manifest.accession,
                    "url": str(item.url),
                    "sha256": actual,
                    "path": str(destination),
                    "retrieved_at": datetime.now(UTC).isoformat(),
                }
            )
    payload: dict[str, Any] = {
        "status": "OK" if records else "NO_VERIFIED_DEFAULT_DOWNLOADS",
        "records": records,
    }
    write_json(cache_dir / "fetch_index.json", payload)
    return payload
