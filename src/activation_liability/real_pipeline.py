"""End-to-end execution of the public real-data validation layer."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import pandas as pd

from activation_liability.audit import run_audit
from activation_liability.benchmark import run_benchmark
from activation_liability.io import write_json
from activation_liability.real_benchmark import (
    run_gse96583_confirmation,
    run_real_benchmark,
)
from activation_liability.real_data import (
    build_gse96583,
    build_gse140244,
    build_gse157857,
    build_gse178429,
    combine_real_cohorts,
    ensembl_symbol_map_from_gse157857,
    load_control_targets,
)


def _canonicalize_floats(value: Any, *, digits: int = 12) -> Any:
    """Round finite floats before persisted real-result serialization."""

    if isinstance(value, float):
        return round(value, digits) if math.isfinite(value) else value
    if isinstance(value, dict):
        return {key: _canonicalize_floats(item, digits=digits) for key, item in value.items()}
    if isinstance(value, list):
        return [_canonicalize_floats(item, digits=digits) for item in value]
    return value


def run_real_pipeline(
    *,
    data_root: Path,
    output_root: Path,
    controls_path: Path,
    rules_path: Path,
    benchmark_config_path: Path,
) -> dict[str, Any]:
    """Build the frozen cohorts, audit them and run the constrained real benchmark."""

    targets = load_control_targets(controls_path)
    gse157 = build_gse157857(
        data_root / "01_DISCOVERY_PRIMARY_RNA_ADT/GSE157857_IFNB_CITEseq_PRIMARY/GSE157857_RAW.tar",
        targets=targets,
    )
    gse178_dir = data_root / "02_VALIDATION_IFN_TNF_TLR/GSE178429_MULTI_STIMULUS_scRNAseq"
    gse178 = build_gse178429(
        gse178_dir / "GSE178429_PBMCs_stim_scRNAseq_counts.txt.gz",
        gse178_dir / "GSE178429_PBMCs_stim_scRNAseq_geneNames.txt.gz",
        gse178_dir / "GSE178429_PBMCs_stim_scRNAseq_cellMeta.txt.gz",
        targets=targets,
    )
    gse140_dir = data_root / "04_SORTED_BULK_TIMECOURSES/GSE140244_CD4_ACTIVATION_TIMECOURSE"
    gene_map = ensembl_symbol_map_from_gse157857(
        data_root / "01_DISCOVERY_PRIMARY_RNA_ADT/GSE157857_IFNB_CITEseq_PRIMARY/GSE157857_RAW.tar"
    )
    gse140 = build_gse140244(
        gse140_dir / "GSE140244_rnaseq_gene_counts.txt.gz",
        gse140_dir / "GSE140244_rnaseq_meta_data.txt.gz",
        targets=targets,
        ensembl_to_symbol=gene_map,
    )
    baseline = combine_real_cohorts([gse157, gse178, gse140])
    gse96583_dir = data_root / "02_VALIDATION_IFN_TNF_TLR/GSE96583_IFNB_demuxlet_scRNAseq"
    confirmation = build_gse96583(
        gse96583_dir / "GSE96583_RAW.tar",
        gse96583_dir / "GSE96583_batch2.genes.tsv.gz",
        gse96583_dir / "GSE96583_batch2.total.tsne.df.tsv.gz",
        targets=targets,
    )

    output_root.mkdir(parents=True, exist_ok=True)
    baseline.cells.to_csv(
        output_root / "real_cells.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    confirmation.cells.to_csv(
        output_root / "gse96583_confirmation_cells.csv.gz",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    write_json(output_root / "ingestion_qc.json", baseline.qc)
    write_json(output_root / "gse96583_confirmation_ingestion_qc.json", confirmation.qc)

    audit = run_audit(baseline.cells, rules_path=rules_path)
    audit_rows = pd.DataFrame(audit["audit_rows"])
    target_summary = pd.DataFrame(audit["target_summary"])
    write_json(output_root / "audit.json", audit)
    audit_rows.to_csv(output_root / "audit_rows.csv", index=False)
    target_summary.to_csv(output_root / "target_summary.csv", index=False)
    write_json(output_root / "claims.json", audit["claims"])

    confirmation_audit = run_audit(confirmation.cells, rules_path=rules_path)
    confirmation_audit_rows = pd.DataFrame(confirmation_audit["audit_rows"])
    confirmation_target_summary = pd.DataFrame(confirmation_audit["target_summary"])
    write_json(output_root / "gse96583_confirmation_audit.json", confirmation_audit)
    confirmation_audit_rows.to_csv(
        output_root / "gse96583_confirmation_audit_rows.csv", index=False
    )
    confirmation_target_summary.to_csv(
        output_root / "gse96583_confirmation_target_summary.csv", index=False
    )
    write_json(
        output_root / "gse96583_confirmation_claims.json",
        confirmation_audit["claims"],
    )

    unconstrained = run_benchmark(
        baseline.cells,
        controls_path=controls_path,
        rules_path=rules_path,
    )
    write_json(output_root / "unconstrained_diagnostic.json", unconstrained)

    benchmark = run_real_benchmark(
        baseline.cells,
        audit_rows,
        controls_path=controls_path,
        config_path=benchmark_config_path,
    )
    benchmark["required_ablations"] = {
        "within_study_vs_cross_study": unconstrained["ablations"]["within_study_vs_cross_study"],
        "pseudobulk_vs_cell_level": unconstrained["ablations"]["pseudobulk_vs_cell_level"],
        "detection_threshold_sweep": unconstrained["ablations"]["detection_threshold_sweep"],
        "rna_only_vs_protein_corroborated": {
            "status": "UNAVAILABLE_NO_CONTROL_TARGET_OVERLAP",
            "reason": (
                "The GSE157857 ADT panel measures ligands and myeloid markers but none of "
                "the benchmark control targets."
            ),
        },
    }
    confirmation_benchmark = run_gse96583_confirmation(
        confirmation.cells,
        confirmation_audit_rows,
        baseline_benchmark=benchmark,
        baseline_audit_rows=audit_rows,
        controls_path=controls_path,
        config_path=benchmark_config_path,
    )
    benchmark["external_confirmation"] = {
        "status": confirmation_benchmark["status"],
        "accession": "GSE96583",
        "tier1_ranking_metrics": confirmation_benchmark["tier1_ranking_metrics"],
        "frozen_cutoff_classification_tier1": confirmation_benchmark[
            "frozen_cutoff_classification_tier1"
        ],
        "target_holdout_status": confirmation_benchmark["target_holdout"]["status"],
        "coverage": confirmation_benchmark["coverage"],
        "full_result": "gse96583_confirmation_benchmark.json",
    }
    benchmark = _canonicalize_floats(benchmark)
    confirmation_benchmark = _canonicalize_floats(confirmation_benchmark)
    write_json(output_root / "benchmark.json", benchmark)
    write_json(
        output_root / "gse96583_confirmation_benchmark.json",
        confirmation_benchmark,
    )
    pd.DataFrame(benchmark["target_scores"]).to_csv(output_root / "target_scores.csv", index=False)
    pd.DataFrame(benchmark["donor_robustness"]["rows"]).to_csv(
        output_root / "donor_robustness.csv", index=False
    )
    pd.DataFrame(benchmark["timecourse_validation"]["rows"]).to_csv(
        output_root / "timecourse_validation.csv", index=False
    )
    pd.DataFrame(confirmation_benchmark["target_scores"]).to_csv(
        output_root / "gse96583_confirmation_target_scores.csv", index=False
    )
    pd.DataFrame(confirmation_benchmark["donor_robustness"]["rows"]).to_csv(
        output_root / "gse96583_confirmation_donor_robustness.csv", index=False
    )
    return {
        "status": "COMPUTED_REAL_PUBLIC_DATA_WITH_PARTIAL_EXTERNAL_CONFIRMATION",
        "n_baseline_rows": int(len(baseline.cells)),
        "n_confirmation_rows": int(len(confirmation.cells)),
        "n_targets": int(baseline.cells["target"].nunique()),
        "n_baseline_studies": int(baseline.cells["study"].nunique()),
        "primary_tier1_metrics": benchmark["primary_tier1_metrics"],
        "diagnostic_holdout_metrics": benchmark["holdout"]["holdout_metrics"],
        "confirmation_tier1_metrics": confirmation_benchmark["tier1_ranking_metrics"],
        "confirmation_target_holdout_status": confirmation_benchmark["target_holdout"]["status"],
        "coverage": benchmark["coverage"],
        "confirmation_coverage": confirmation_benchmark["coverage"],
        "output_root": str(output_root),
    }


def write_pipeline_summary(payload: dict[str, Any], path: Path) -> None:
    """Write a compact plain JSON summary for orchestration scripts."""

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
