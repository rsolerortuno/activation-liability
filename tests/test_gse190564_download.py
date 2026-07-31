from __future__ import annotations

import gzip
import importlib.util
import io
import tarfile
from pathlib import Path

import pytest


def _load_validator(root: Path):
    path = root / "scripts" / "validate_gse190564_download.py"
    spec = importlib.util.spec_from_file_location("validate_gse190564_download", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_soft(directory: Path) -> None:
    with gzip.open(directory / "GSE190564_family.soft.gz", "wt") as handle:
        handle.write("^SERIES = GSE190564\n")


def test_validator_detects_required_modalities(root: Path, tmp_path: Path) -> None:
    module = _load_validator(root)
    archive = tmp_path / "GSE190564_processed_data.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        for name in ["Pool1_GEX/matrix.mtx.gz", "Pool1_ADT/counts.tsv.gz"]:
            data = b"test"
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    _write_soft(tmp_path)
    payload = module.validate_gse190564(tmp_path, minimum_archive_size=1)
    assert payload["archive"]["gex_member_count"] == 1
    assert payload["archive"]["adt_member_count"] == 1
    assert payload["status"] == "DOWNLOAD_VALIDATED_NOT_ANALYSED"


def test_validator_rejects_missing_adt(root: Path, tmp_path: Path) -> None:
    module = _load_validator(root)
    archive = tmp_path / "GSE190564_processed_data.tar.gz"
    with tarfile.open(archive, "w:gz") as tar:
        data = b"test"
        info = tarfile.TarInfo("Pool1_GEX/matrix.mtx.gz")
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    _write_soft(tmp_path)
    with pytest.raises(ValueError, match="both GEX and ADT"):
        module.validate_gse190564(tmp_path, minimum_archive_size=1)
