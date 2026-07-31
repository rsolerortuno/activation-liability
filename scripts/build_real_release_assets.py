#!/usr/bin/env python3
# ruff: noqa: E501
"""Build deterministic release-side manifests and an HTML report for real results."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
from pathlib import Path
from typing import Any

SOURCE_PATHS = (
    (
        "GSE157857",
        "01_DISCOVERY_PRIMARY_RNA_ADT/GSE157857_IFNB_CITEseq_PRIMARY/GSE157857_RAW.tar",
        "primary RNA/ADT/HTO source",
    ),
    (
        "GSE178429",
        "02_VALIDATION_IFN_TNF_TLR/GSE178429_MULTI_STIMULUS_scRNAseq/"
        "GSE178429_PBMCs_stim_scRNAseq_cellMeta.txt.gz",
        "cell metadata",
    ),
    (
        "GSE178429",
        "02_VALIDATION_IFN_TNF_TLR/GSE178429_MULTI_STIMULUS_scRNAseq/"
        "GSE178429_PBMCs_stim_scRNAseq_counts.txt.gz",
        "count matrix",
    ),
    (
        "GSE178429",
        "02_VALIDATION_IFN_TNF_TLR/GSE178429_MULTI_STIMULUS_scRNAseq/"
        "GSE178429_PBMCs_stim_scRNAseq_geneNames.txt.gz",
        "gene names",
    ),
    (
        "GSE96583",
        "02_VALIDATION_IFN_TNF_TLR/GSE96583_IFNB_demuxlet_scRNAseq/GSE96583_RAW.tar",
        "control and IFN-beta 6 h matrices",
    ),
    (
        "GSE96583",
        "02_VALIDATION_IFN_TNF_TLR/GSE96583_IFNB_demuxlet_scRNAseq/GSE96583_batch2.genes.tsv.gz",
        "authoritative matrix row order",
    ),
    (
        "GSE96583",
        "02_VALIDATION_IFN_TNF_TLR/GSE96583_IFNB_demuxlet_scRNAseq/"
        "GSE96583_batch2.total.tsne.df.tsv.gz",
        "demuxlet donor and singlet metadata",
    ),
    (
        "GSE140244",
        "04_SORTED_BULK_TIMECOURSES/GSE140244_CD4_ACTIVATION_TIMECOURSE/"
        "GSE140244_rnaseq_gene_counts.txt.gz",
        "bulk RNA-seq counts",
    ),
    (
        "GSE140244",
        "04_SORTED_BULK_TIMECOURSES/GSE140244_CD4_ACTIVATION_TIMECOURSE/"
        "GSE140244_rnaseq_meta_data.txt.gz",
        "bulk RNA-seq metadata",
    ),
)


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    """Read a JSON object."""

    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise TypeError(f"Expected a JSON object in {path}")
    return value


def write_json(path: Path, payload: dict[str, Any]) -> None:
    """Write canonical indented JSON."""

    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def build_source_checksums(data_root: Path, output_root: Path) -> None:
    """Write hashes for every input consumed by the v0.4 pipeline."""

    rows: list[dict[str, Any]] = []
    for accession, relative_path, role in SOURCE_PATHS:
        path = data_root / relative_path
        if not path.is_file():
            raise FileNotFoundError(path)
        rows.append(
            {
                "accession": accession,
                "file_name": path.name,
                "role": role,
                "sha256": sha256(path),
                "size_bytes": path.stat().st_size,
            }
        )
    write_json(
        output_root / "source_checksums.json",
        {
            "code_version": "0.4.0",
            "excluded_inputs": [],
            "inputs": rows,
            "schema_version": 1,
        },
    )


def build_execution_summary(output_root: Path, blockers_path: Path) -> None:
    """Write a compact machine-readable real-run summary."""

    benchmark = read_json(output_root / "benchmark.json")
    confirmation = read_json(output_root / "gse96583_confirmation_benchmark.json")
    qc = read_json(output_root / "ingestion_qc.json")
    confirmation_qc = read_json(output_root / "gse96583_confirmation_ingestion_qc.json")
    blockers = read_json(blockers_path)
    cohorts = [row["accession"] for row in qc.get("cohorts", [])]
    cohorts.append(str(confirmation_qc["accession"]))
    write_json(
        output_root / "execution_summary.json",
        {
            "blockers": blockers["blockers"],
            "cohorts_executed": cohorts,
            "coverage": benchmark["coverage"],
            "expanded_tier1_tier2": benchmark["expanded_tier1_tier2_metrics"],
            "heterogeneity_status": benchmark["heterogeneity"]["status"],
            "diagnostic_holdout_status": benchmark["holdout"]["integrity_status"],
            "external_confirmation": {
                "accession": "GSE96583",
                "status": confirmation["status"],
                "tier1_ranking_metrics": confirmation["tier1_ranking_metrics"],
                "frozen_cutoff_classification_tier1": confirmation[
                    "frozen_cutoff_classification_tier1"
                ],
                "target_holdout_status": confirmation["target_holdout"]["status"],
                "score_correlation_with_v0_3": confirmation["agreement_with_v0_3"][
                    "spearman_score_correlation"
                ],
            },
            "interpretation": (
                "Coverage-limited public-data validation with one independent IFN-I cohort. "
                "The observable control ranking replicated, but the preassigned positive "
                "target holdout was not evaluable at the frozen observability threshold."
            ),
            "primary_endpoints": [
                "GSE157857/IFN_I_18h",
                "GSE178429/IFN_II_6h",
                "GSE178429/TLR_6h",
                "GSE140244/lymphocyte_activation_24h",
            ],
            "confirmation_endpoint": "GSE96583/IFN_I_6h_confirmation",
            "primary_tier1": benchmark["primary_tier1_metrics"],
            "protein_benchmark_status": benchmark["protein_benchmark"]["status"],
            "schema_version": 1,
            "software_version": "0.4.0",
            "status": "COMPUTED_REAL_PUBLIC_DATA_WITH_PARTIAL_EXTERNAL_CONFIRMATION",
        },
    )


def build_report(output_root: Path, blockers_path: Path) -> None:
    """Build a standalone HTML report from committed JSON/CSV artefacts."""

    benchmark = read_json(output_root / "benchmark.json")
    confirmation = read_json(output_root / "gse96583_confirmation_benchmark.json")
    blockers = read_json(blockers_path)
    with (output_root / "target_scores.csv").open(encoding="utf-8", newline="") as stream:
        scores = [row for row in csv.DictReader(stream) if row["observable"] == "True"]
    scores.sort(key=lambda row: float(row["score"]), reverse=True)
    rows = "".join(
        "<tr>"
        + "".join(
            f"<td>{html.escape(str(row[column]))}</td>"
            for column in (
                "target",
                "label",
                "tier",
                "score",
                "driver_stimulus",
                "driver_cell_type",
                "driver_total_pairs",
            )
        )
        + "</tr>"
        for row in scores[:15]
    )
    blocker_rows = "".join(
        "<li><strong>"
        + html.escape(str(row.get("accession", row.get("component", "blocker"))))
        + "</strong>: "
        + html.escape(str(row["reason"]))
        + "</li>"
        for row in blockers["blockers"]
    )
    coverage = benchmark["coverage"]
    primary = benchmark["primary_tier1_metrics"]
    expanded = benchmark["expanded_tier1_tier2_metrics"]
    inflation = benchmark["required_ablations"]["pseudobulk_vs_cell_level"]["inflation_ratio"]
    confirm_metrics = confirmation["tier1_ranking_metrics"]
    confirm_classification = confirmation["frozen_cutoff_classification_tier1"]
    confirm_coverage = confirmation["coverage"]
    holdout_status = confirmation["target_holdout"]["status"]
    correlation = confirmation["agreement_with_v0_3"]["spearman_score_correlation"]
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>activation-liability v0.4.0 — public real-data validation</title>
<style>body{{font-family:Arial,sans-serif;max-width:1050px;margin:40px auto;padding:0 20px;color:#1f2937;line-height:1.5}}h1,h2{{color:#111827}}.warning{{border-left:5px solid #b45309;background:#fffbeb;padding:14px 18px}}.ok{{border-left:5px solid #047857;background:#ecfdf5;padding:14px 18px}}.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px}}.card{{border:1px solid #d1d5db;border-radius:8px;padding:14px}}.metric{{font-size:1.55rem;font-weight:700}}table{{border-collapse:collapse;width:100%;font-size:.9rem}}th,td{{border:1px solid #d1d5db;padding:7px;text-align:left}}th{{background:#f3f4f6}}code{{background:#f3f4f6;padding:2px 4px}}footer{{margin-top:32px;color:#6b7280;font-size:.85rem}}</style></head>
<body><h1>activation-liability v0.4.0</h1><p>Public real-data validation and independent IFN-I replication report</p>
<div class="warning"><strong>Interpretation boundary.</strong> This is a coverage-limited liability benchmark using donor-paired within-study contrasts. It does not predict clinical toxicity. The original target holdout remains diagnostic.</div>
<h2>Primary benchmark</h2><div class="grid">
<div class="card"><div>Observable tier-1 controls</div><div class="metric">{coverage["observable_tier1_controls"]} / {coverage["tier1_controls"]}</div></div>
<div class="card"><div>Tier-1 AUROC</div><div class="metric">{primary["auroc"]:.4f}</div></div>
<div class="card"><div>Tier-1 average precision</div><div class="metric">{primary["average_precision"]:.4f}</div></div>
<div class="card"><div>Expanded AUROC</div><div class="metric">{expanded["auroc"]:.4f}</div></div>
<div class="card"><div>Cell-level inflation</div><div class="metric">{inflation:.2f}×</div></div></div>
<h2>Independent GSE96583 confirmation</h2>
<div class="ok"><strong>What replicated.</strong> The unchanged v0.3 score ranked the {confirm_coverage["observable_tier1_controls"]} observable tier-1 controls with AUROC {confirm_metrics["auroc"]:.4f} and average precision {confirm_metrics["average_precision"]:.4f}. The score correlation with v0.3 was {correlation:.3f}.</div>
<div class="grid"><div class="card"><div>Frozen-cutoff balanced accuracy</div><div class="metric">{confirm_classification["balanced_accuracy"]:.3f}</div></div><div class="card"><div>Sensitivity</div><div class="metric">{confirm_classification["sensitivity"]:.3f}</div></div><div class="card"><div>Specificity</div><div class="metric">{confirm_classification["specificity"]:.3f}</div></div><div class="card"><div>Target holdout</div><div class="metric" style="font-size:1rem">{html.escape(holdout_status)}</div></div></div>
<p>The target holdout was not declared confirmed: IL2RA, PDCD1LG2 and TNFRSF9 all abstained at the frozen 5% observability threshold in this cohort.</p>
<h2>Executed cohorts</h2><ul><li>GSE157857: IFN-beta 18 h, broad myeloid, 3 HTO-resolved donors.</li><li>GSE178429: IFN-gamma 6 h and LPS 6 h, broad PBMC lineages, donor-paired.</li><li>GSE140244: anti-CD3/CD28 24 h, sorted CD4-memory T cells, 24 complete donors.</li><li>GSE96583: IFN-beta 6 h, broad PBMC lineages, 8 demuxlet-resolved SLE donors; frozen external confirmation.</li></ul>
<h2>Top observable primary scores</h2><p>Score = maximum one-sided 95% lower confidence bound across predeclared target-relevant endpoints and lineages.</p><table><thead><tr><th>Target</th><th>Label</th><th>Tier</th><th>Score</th><th>Driver stimulus</th><th>Cell type</th><th>Paired donors</th></tr></thead><tbody>{rows}</tbody></table>
<h2>Guardrails</h2><ul><li>Targets outside covered lineages or below 5% observability abstain.</li><li>CD3E, CD79A and MS4A1 use leave-target-out lineage annotation in both PBMC cohorts.</li><li>The 313 duplicated 10x barcode sequences in GSE96583 are resolved condition-specifically and checked one-to-one.</li><li>No protein-aware control benchmark is claimed because ADT/control overlap is zero.</li></ul>
<h2>Current blockers</h2><ul>{blocker_rows}</ul><h2>Reproducibility</h2><p>Machine-readable inputs and outputs include <code>source_checksums.json</code>, <code>benchmark.json</code>, <code>gse96583_confirmation_benchmark.json</code>, <code>claims.json</code>, and <code>execution_summary.json</code>.</p>
<footer>Generated from committed v0.4.0 result artefacts. Quantitative README claims are checked against result JSON by the test suite.</footer></body></html>"""
    (output_root / "report.html").write_text(document, encoding="utf-8")


def build_result_manifest(output_root: Path) -> None:
    """Hash all release result files except the manifest itself."""

    files = [
        {
            "file_name": path.name,
            "sha256": sha256(path),
            "size_bytes": path.stat().st_size,
        }
        for path in sorted(output_root.iterdir())
        if path.is_file() and path.name != "result_manifest.json"
    ]
    write_json(
        output_root / "result_manifest.json",
        {"files": files, "schema_version": 1, "software_version": "0.4.0"},
    )


def main() -> None:
    """Parse arguments and generate deterministic release assets."""

    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--blockers", type=Path, default=Path("results/real/BLOCKERS.json"))
    args = parser.parse_args()
    build_source_checksums(args.data_root, args.output)
    build_execution_summary(args.output, args.blockers)
    build_report(args.output, args.blockers)
    build_result_manifest(args.output)


if __name__ == "__main__":
    main()
