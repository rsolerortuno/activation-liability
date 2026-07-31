from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import tarfile
from pathlib import Path
from typing import Any


def sha256_file(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def validate_gse190564(
    directory: Path, minimum_archive_size: int = 10_000_000_000
) -> dict[str, Any]:
    archive = directory / "GSE190564_processed_data.tar.gz"
    soft = directory / "GSE190564_family.soft.gz"
    if not archive.exists():
        raise FileNotFoundError(archive)
    if not soft.exists():
        raise FileNotFoundError(soft)
    if archive.stat().st_size < minimum_archive_size:
        raise ValueError(f"Archive is unexpectedly small: {archive.stat().st_size:,} bytes")

    with gzip.open(soft, "rt", encoding="utf-8", errors="replace") as handle:
        prefix = handle.read(4096)
    if "GSE190564" not in prefix:
        raise ValueError("SOFT metadata does not identify GSE190564")

    with tarfile.open(archive, mode="r:gz") as tar:
        members = [member.name for member in tar.getmembers() if member.isfile()]
    gex = [name for name in members if "GEX" in name.upper()]
    adt = [name for name in members if "ADT" in name.upper()]
    if not gex or not adt:
        raise ValueError("Archive does not contain both GEX and ADT members")

    payload: dict[str, Any] = {
        "accession": "GSE190564",
        "status": "DOWNLOAD_VALIDATED_NOT_ANALYSED",
        "archive": {
            "file_name": archive.name,
            "size_bytes": archive.stat().st_size,
            "sha256": sha256_file(archive),
            "member_count": len(members),
            "gex_member_count": len(gex),
            "adt_member_count": len(adt),
        },
        "soft": {
            "file_name": soft.name,
            "size_bytes": soft.stat().st_size,
            "sha256": sha256_file(soft),
        },
    }
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate a downloaded GSE190564 bundle")
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = validate_gse190564(args.directory)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
