from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pytest

RESULT_PATTERN = re.compile(r'data-alia-result="(?P<path>[^"]+)">(?P<value>-?[0-9]+(?:\.[0-9]+)?)<')
REAL_RESULT_PATTERN = re.compile(
    r'data-alia-real-result="(?P<path>[^"]+)">(?P<value>-?[0-9]+(?:\.[0-9]+)?)<'
)
CONFIRM_RESULT_PATTERN = re.compile(
    r'data-alia-confirm-result="(?P<path>[^"]+)">(?P<value>-?[0-9]+(?:\.[0-9]+)?)<'
)
TISSUE_RESULT_PATTERN = re.compile(
    r'data-alia-tissue-result="(?P<path>[^"]+)">(?P<value>-?[0-9]+(?:\.[0-9]+)?)<'
)


def _resolve(payload: Any, path: str) -> Any:
    current = payload
    for part in path.split("."):
        length = part.endswith("#len")
        key = part.removesuffix("#len")
        current = current[int(key)] if isinstance(current, list) else current[key]
        if length:
            current = len(current)
    return current


def test_readme_quantitative_claims_match_committed_results(root: Path) -> None:
    readme = (root / "README.md").read_text(encoding="utf-8")
    payload = json.loads((root / "results/synthetic/benchmark.json").read_text())
    matches = list(RESULT_PATTERN.finditer(readme))
    assert len(matches) >= 15
    for match in matches:
        displayed_text = match.group("value")
        displayed = float(displayed_text)
        actual = float(_resolve(payload, match.group("path")))
        decimals = len(displayed_text.split(".", maxsplit=1)[1]) if "." in displayed_text else 0
        assert displayed == pytest.approx(actual, abs=0.5 * 10 ** (-decimals))


def test_readme_real_quantitative_claims_match_committed_results(root: Path) -> None:
    readme = (root / "README.md").read_text(encoding="utf-8")
    payload = json.loads(
        (root / "results/real/public_v0_4_0/benchmark.json").read_text(encoding="utf-8")
    )
    matches = list(REAL_RESULT_PATTERN.finditer(readme))
    assert len(matches) >= 10
    for match in matches:
        displayed_text = match.group("value")
        displayed = float(displayed_text)
        actual = float(_resolve(payload, match.group("path")))
        decimals = len(displayed_text.split(".", maxsplit=1)[1]) if "." in displayed_text else 0
        assert displayed == pytest.approx(actual, abs=0.5 * 10 ** (-decimals))


def test_readme_confirmation_claims_match_committed_results(root: Path) -> None:
    readme = (root / "README.md").read_text(encoding="utf-8")
    payload = json.loads(
        (root / "results/real/public_v0_4_0/gse96583_confirmation_benchmark.json").read_text(
            encoding="utf-8"
        )
    )
    matches = list(CONFIRM_RESULT_PATTERN.finditer(readme))
    assert len(matches) >= 10
    for match in matches:
        displayed_text = match.group("value")
        displayed = float(displayed_text)
        actual = float(_resolve(payload, match.group("path")))
        decimals = len(displayed_text.split(".", maxsplit=1)[1]) if "." in displayed_text else 0
        assert displayed == pytest.approx(actual, abs=0.5 * 10 ** (-decimals))


def test_readme_tissue_claims_match_committed_results(root: Path) -> None:
    readme = (root / "README.md").read_text(encoding="utf-8")
    payload = json.loads(
        (root / "results/real/public_v0_5_0/tissue_benchmark.json").read_text(encoding="utf-8")
    )
    matches = list(TISSUE_RESULT_PATTERN.finditer(readme))
    assert len(matches) >= 15
    for match in matches:
        displayed_text = match.group("value")
        displayed = float(displayed_text)
        actual = float(_resolve(payload, match.group("path")))
        decimals = len(displayed_text.split(".", maxsplit=1)[1]) if "." in displayed_text else 0
        assert displayed == pytest.approx(actual, abs=0.5 * 10 ** (-decimals))


def test_tissue_result_manifest_matches_files(root: Path) -> None:
    import hashlib

    result_root = root / "results/real/public_v0_5_0"
    payload = json.loads((result_root / "tissue_result_manifest.json").read_text(encoding="utf-8"))
    assert len(payload["files"]) >= 10
    for row in payload["files"]:
        path = result_root / row["file_name"]
        assert path.stat().st_size == row["size"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"]
