#!/usr/bin/env python3
"""Download public validation files, resume partial transfers, and record checksums.

This bootstrapper intentionally reads a candidate catalog rather than enabled registry manifests.
A successful download does not promote a study to a default: metadata, pairing, licence and
conversion checks remain mandatory.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
import time
import urllib.error
import urllib.request
from collections.abc import Iterable
from pathlib import Path
from typing import Any, cast

import yaml

CHUNK_SIZE = 1024 * 1024
USER_AGENT = "activation-liability-real-data-bootstrap/0.2"


def load_catalog(path: Path) -> list[dict[str, Any]]:
    """Load catalog entries and validate the minimal bootstrap contract."""

    payload = cast(dict[str, Any], yaml.safe_load(path.read_text(encoding="utf-8")))
    datasets = payload.get("datasets")
    if not isinstance(datasets, list):
        raise ValueError("catalog must contain a datasets list")
    normalized: list[dict[str, Any]] = []
    for raw in datasets:
        if not isinstance(raw, dict):
            raise ValueError("every dataset entry must be a mapping")
        accession = raw.get("accession")
        files = raw.get("files")
        if not isinstance(accession, str) or not accession:
            raise ValueError("every dataset requires an accession")
        if not isinstance(files, list):
            raise ValueError(f"{accession}: files must be a list")
        normalized.append(cast(dict[str, Any], raw))
    return normalized


def sha256_file(path: Path) -> str:
    """Return a file SHA-256 digest."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(block)
    return digest.hexdigest()


def download_with_resume(url: str, destination: Path, retries: int = 3) -> None:
    """Download one URL, resuming a `.part` file when the server supports Range."""

    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    for attempt in range(1, retries + 1):
        offset = partial.stat().st_size if partial.exists() else 0
        headers = {"User-Agent": USER_AGENT}
        if offset:
            headers["Range"] = f"bytes={offset}-"
        request = urllib.request.Request(url, headers=headers)
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                status = getattr(response, "status", 200)
                append = offset > 0 and status == 206
                if offset and not append:
                    partial.unlink(missing_ok=True)
                mode = "ab" if append else "wb"
                with partial.open(mode) as handle:
                    shutil.copyfileobj(response, handle, length=CHUNK_SIZE)
            partial.replace(destination)
            return
        except (OSError, urllib.error.URLError) as exc:
            if attempt == retries:
                raise RuntimeError(f"download failed after {retries} attempts: {url}") from exc
            time.sleep(float(attempt))


def choose_datasets(
    catalog: Iterable[dict[str, Any]], accessions: set[str], include_large: bool
) -> list[dict[str, Any]]:
    """Select explicitly named datasets or the recommended small-data set."""

    selected: list[dict[str, Any]] = []
    for dataset in catalog:
        accession = cast(str, dataset["accession"])
        group = dataset.get("download_group")
        explicit = accession in accessions
        recommended_small = not accessions and group == "small"
        if explicit or recommended_small or (not accessions and include_large and group == "large"):
            selected.append(dataset)
    missing = accessions - {cast(str, item["accession"]) for item in selected}
    if missing:
        raise ValueError(f"unknown accessions: {', '.join(sorted(missing))}")
    return selected


def bootstrap_dataset(dataset: dict[str, Any], root: Path) -> dict[str, Any]:
    """Download one dataset and return a provenance record."""

    accession = cast(str, dataset["accession"])
    relative = Path(cast(str, dataset["drive_folder"]))
    output_dir = root / relative
    output_dir.mkdir(parents=True, exist_ok=True)
    files = cast(list[dict[str, Any]], dataset["files"])
    records: list[dict[str, Any]] = []
    if not files:
        return {
            "accession": accession,
            "status": "NO_VERIFIED_FILE_MAPPING",
            "source_record": dataset["source_record"],
            "files": [],
        }
    for spec in files:
        filename = cast(str, spec["filename"])
        url = cast(str, spec["url"])
        destination = output_dir / filename
        if not destination.exists():
            print(f"Downloading {accession}/{filename}", flush=True)
            download_with_resume(url, destination)
        digest = sha256_file(destination)
        records.append(
            {
                "filename": filename,
                "source_url": url,
                "size_bytes": destination.stat().st_size,
                "sha256": digest,
            }
        )
    record = {
        "accession": accession,
        "status": "DOWNLOADED_NOT_PROMOTED",
        "source_record": dataset["source_record"],
        "downloaded_at_unix": int(time.time()),
        "files": records,
        "remaining_blockers": dataset.get("blockers_before_default", []),
    }
    (output_dir / "download_manifest.json").write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return record


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--catalog",
        type=Path,
        default=Path("data/candidates/real_data_catalog.yaml"),
    )
    parser.add_argument("--dest", type=Path, required=True)
    parser.add_argument("--dataset", action="append", default=[])
    parser.add_argument("--include-large", action="store_true")
    parser.add_argument("--list", action="store_true", dest="list_only")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Run the bootstrap command."""

    args = parse_args(argv)
    catalog = load_catalog(args.catalog)
    if args.list_only:
        for dataset in catalog:
            print(
                f"{dataset['accession']}\tpriority={dataset['priority']}\t"
                f"group={dataset['download_group']}\trole={dataset['role']}"
            )
        return 0
    selected = choose_datasets(catalog, set(args.dataset), bool(args.include_large))
    if not selected:
        print("No downloadable datasets selected.", file=sys.stderr)
        return 2
    records = [bootstrap_dataset(dataset, args.dest) for dataset in selected]
    summary = {
        "schema_version": 1,
        "status": "DOWNLOADS_DO_NOT_IMPLY_REGISTRY_PROMOTION",
        "datasets": records,
    }
    args.dest.mkdir(parents=True, exist_ok=True)
    (args.dest / "bootstrap_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
