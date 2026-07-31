#!/usr/bin/env python3
"""Build the paired tissue extension with memory-isolated psoriasis samples."""

from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import tarfile
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any, cast

import pandas as pd

from activation_liability.benchmark import load_controls
from activation_liability.io import write_json
from activation_liability.tissue import (
    PSORIASIS_BASELINE_MAP,
    _aggregate_sample_files,
    load_tissue_config,
)
from activation_liability.tissue_pipeline import run_tissue_pipeline


def _extract_baseline_members(archive_path: Path, output_directory: Path) -> None:
    output_directory.mkdir(parents=True, exist_ok=True)
    wanted = {f"{prefix}.matrix.mtx.gz" for prefix in PSORIASIS_BASELINE_MAP} | {
        f"{prefix}.features.tsv.gz" for prefix in PSORIASIS_BASELINE_MAP
    }
    with tarfile.open(archive_path, "r") as archive:
        for member in archive.getmembers():
            basename = Path(member.name).name
            if not any(basename.endswith(token) for token in wanted):
                continue
            source = archive.extractfile(member)
            if source is None:
                raise ValueError(f"unable to extract {member.name}")
            destination = output_directory / basename
            with destination.open("wb") as stream:
                while chunk := source.read(16 * 1024 * 1024):
                    stream.write(chunk)
    matrices = list(output_directory.glob("*_V1_*.matrix.mtx.gz"))
    features = list(output_directory.glob("*_V1_*.features.tsv.gz"))
    if len(matrices) != 10 or len(features) != 10:
        raise ValueError(
            "expected ten baseline psoriasis matrices and ten feature files; "
            f"found {len(matrices)} and {len(features)}"
        )


def _one_sample(task: dict[str, Any]) -> tuple[pd.DataFrame, dict[str, object]]:
    return _aggregate_sample_files(
        matrix_path=Path(task["matrix_path"]),
        feature_path=Path(task["feature_path"]),
        donor=str(task["donor"]),
        condition=str(task["condition"]),
        study="GSE228421",
        stimulus="tissue_inflammation_psoriasis",
        tissue="skin",
        targets=cast(list[str], task["targets"]),
        marker_map=cast(Mapping[str, Iterable[str]], task["marker_map"]),
        minimum_positive_markers=cast(Mapping[str, int], task["minimum_positive_markers"]),
        qc_parameters=cast(Mapping[str, float], task["qc_parameters"]),
        excluded_lineages=set(cast(list[str], task["excluded_lineages"])),
    )


def _one_sample_to_disk(
    task: dict[str, Any],
    csv_path: Path,
    json_path: Path,
) -> None:
    frame, qc = _one_sample(task)
    frame.to_csv(csv_path, index=False)
    write_json(json_path, qc)


def _build_psoriasis_cache(
    directory: Path,
    *,
    controls_path: Path,
    config_path: Path,
    output_pseudobulk: Path,
    output_qc: Path,
) -> None:
    config = load_tissue_config(config_path)
    controls = load_controls(controls_path)
    targets = sorted(controls["target"].astype(str))
    endpoint = next(item for item in config["endpoints"] if str(item["study"]) == "GSE228421")
    tasks: list[dict[str, Any]] = []
    for prefix, (donor, condition) in PSORIASIS_BASELINE_MAP.items():
        matrices = sorted(directory.glob(f"*_{prefix}.matrix.mtx.gz"))
        features = sorted(directory.glob(f"*_{prefix}.features.tsv.gz"))
        if len(matrices) != 1 or len(features) != 1:
            raise ValueError(f"missing unique matrix/features pair for {prefix}")
        tasks.append(
            {
                "matrix_path": str(matrices[0]),
                "feature_path": str(features[0]),
                "donor": donor,
                "condition": condition,
                "targets": targets,
                "marker_map": config["annotation"]["markers"],
                "minimum_positive_markers": config["annotation"].get(
                    "minimum_positive_markers", {}
                ),
                "qc_parameters": config["qc"],
                "excluded_lineages": endpoint.get("excluded_lineages", []),
            }
        )

    frames: list[pd.DataFrame] = []
    sample_qc: list[dict[str, object]] = []
    sample_cache = output_pseudobulk.parent / "_sample_cache"
    sample_cache.mkdir(parents=True, exist_ok=True)
    context = mp.get_context("spawn")
    for task in tasks:
        prefix = Path(str(task["matrix_path"])).name.removesuffix(".matrix.mtx.gz")
        csv_path = sample_cache / f"{prefix}.csv"
        json_path = sample_cache / f"{prefix}.json"
        if not csv_path.is_file() or not json_path.is_file():
            process = context.Process(
                target=_one_sample_to_disk,
                args=(task, csv_path, json_path),
            )
            process.start()
            process.join()
            if process.exitcode != 0:
                raise RuntimeError(
                    f"isolated psoriasis worker failed for {prefix}: exit {process.exitcode}"
                )
        frames.append(pd.read_csv(csv_path))
        qc_value = json.loads(json_path.read_text(encoding="utf-8"))
        sample_qc.append(cast(dict[str, object], qc_value))

    combined = pd.concat(frames, ignore_index=True).sort_values(
        ["study", "target", "cell_type", "donor", "condition"]
    )
    output_pseudobulk.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(
        output_pseudobulk,
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    write_json(
        output_qc,
        {
            "accession": "GSE228421",
            "analysis_set": "baseline_V1_paired_only",
            "input_mode": "extracted_members_memory_isolated",
            "included_donors": sorted({donor for donor, _ in PSORIASIS_BASELINE_MAP.values()}),
            "excluded_visits": ["V2_day3_lesional", "V3_day14_lesional"],
            "sample_qc": sample_qc,
        },
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--crohn-tar", type=Path, required=True)
    parser.add_argument("--psoriasis-tar", type=Path)
    parser.add_argument("--psoriasis-directory", type=Path, required=True)
    parser.add_argument("--controls", type=Path, default=Path("data/controls/controls.yaml"))
    parser.add_argument("--config", type=Path, default=Path("config/tissue_extension.yaml"))
    parser.add_argument(
        "--baseline-benchmark",
        type=Path,
        default=Path("results/real/public_v0_4_0/benchmark.json"),
    )
    parser.add_argument("--output", type=Path, default=Path("results/real/public_v0_5_0"))
    args = parser.parse_args()

    if args.psoriasis_tar is not None:
        _extract_baseline_members(args.psoriasis_tar, args.psoriasis_directory)

    cache_root = args.output / "_build_cache"
    psoriasis_pseudobulk = cache_root / "gse228421_pseudobulk.csv.gz"
    psoriasis_qc = cache_root / "gse228421_qc.json"
    _build_psoriasis_cache(
        args.psoriasis_directory,
        controls_path=args.controls,
        config_path=args.config,
        output_pseudobulk=psoriasis_pseudobulk,
        output_qc=psoriasis_qc,
    )
    payload = run_tissue_pipeline(
        crohn_tar=args.crohn_tar,
        psoriasis_directory=args.psoriasis_directory,
        psoriasis_pseudobulk_path=psoriasis_pseudobulk,
        psoriasis_qc_path=psoriasis_qc,
        controls_path=args.controls,
        tissue_config_path=args.config,
        baseline_benchmark_path=args.baseline_benchmark,
        output_root=args.output,
    )
    print(json.dumps(payload["tissue_benchmark"]["tier1_metrics"], indent=2))


if __name__ == "__main__":
    main()
