"""Typer command-line interface."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import typer

from activation_liability.audit import run_audit
from activation_liability.benchmark import run_benchmark
from activation_liability.fetch import fetch_default_registry
from activation_liability.io import write_json
from activation_liability.real_pipeline import run_real_pipeline
from activation_liability.registry import validate_registry
from activation_liability.report import render_report
from activation_liability.sensitivity import run_sensitivity
from activation_liability.synthetic import generate_synthetic
from activation_liability.tissue_pipeline import run_tissue_pipeline

app = typer.Typer(
    name="alia",
    help="Audit activation-induced erosion of apparent tumour selectivity.",
    no_args_is_help=True,
)


def _default(path: str) -> Path:
    return Path(path)


@app.command("validate-manifest")
def validate_manifest(
    registry: Path = typer.Argument(_default("data/registry"), exists=True, file_okay=False),
) -> None:
    """Validate all study registry YAML files."""

    manifests, errors = validate_registry(registry)
    if errors:
        for error in errors:
            typer.echo(error, err=True)
        raise typer.Exit(code=1)
    enabled = sum(manifest.default and manifest.verified for manifest in manifests)
    typer.echo(f"Validated {len(manifests)} manifests; {enabled} verified/default entries enabled.")


@app.command()
def build(
    output: Path = typer.Option(_default("data/synthetic_cells.csv"), "--output", "-o"),
    synthetic: bool = typer.Option(False, "--synthetic", help="Generate deterministic fixtures."),
    controls: Path = typer.Option(_default("data/controls/controls.yaml"), exists=True),
    seed: int = typer.Option(20260730),
    studies: int = typer.Option(3, min=1, hidden=True),
    donors_per_study: int = typer.Option(5, min=1, hidden=True),
    cells_per_condition: int = typer.Option(12, min=1, hidden=True),
) -> None:
    """Build a standardised analysis table; currently supports synthetic fixtures."""

    if not synthetic:
        typer.echo("Only --synthetic is available in the sandbox build.", err=True)
        raise typer.Exit(code=2)
    frame = generate_synthetic(
        controls,
        seed=seed,
        studies=studies,
        donors_per_study=donors_per_study,
        cells_per_condition=cells_per_condition,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, index=False)
    typer.echo(f"Wrote {len(frame):,} rows to {output}.")


@app.command()
def audit(
    input_path: Path = typer.Option(..., "--input", "-i", exists=True, dir_okay=False),
    output: Path = typer.Option(_default("results/audit"), "--output", "-o"),
    rules: Path = typer.Option(_default("config/evidence_classes.yaml"), exists=True),
    detection_threshold: int = typer.Option(0, min=0),
    positive_fraction_cutoff: float = typer.Option(0.10, min=0.0, max=1.0),
) -> None:
    """Run the primary within-study donor-pseudobulk audit."""

    cells = pd.read_csv(input_path)
    payload = run_audit(
        cells,
        rules_path=rules,
        detection_count_threshold=detection_threshold,
        positive_fraction_cutoff=positive_fraction_cutoff,
    )
    output.mkdir(parents=True, exist_ok=True)
    write_json(output / "audit.json", payload)
    write_json(output / "claims.json", payload["claims"])
    pd.DataFrame(payload["audit_rows"]).to_csv(output / "audit_rows.csv", index=False)
    pd.DataFrame(payload["target_summary"]).to_csv(output / "target_summary.csv", index=False)
    typer.echo(f"Audit completed for {len(payload['target_summary'])} targets in {output}.")


@app.command()
def benchmark(
    input_path: Path = typer.Option(..., "--input", "-i", exists=True, dir_okay=False),
    output: Path = typer.Option(_default("results/synthetic"), "--output", "-o"),
    controls: Path = typer.Option(_default("data/controls/controls.yaml"), exists=True),
    rules: Path = typer.Option(_default("config/evidence_classes.yaml"), exists=True),
) -> None:
    """Run tiered controls, holdout metrics and required ablations."""

    cells = pd.read_csv(input_path)
    payload = run_benchmark(cells, controls_path=controls, rules_path=rules)
    output.mkdir(parents=True, exist_ok=True)
    write_json(output / "benchmark.json", payload)
    pd.DataFrame(payload["target_scores"]).to_csv(output / "target_scores.csv", index=False)
    write_json(output / "claims.json", payload["claims"])
    typer.echo(json.dumps(payload["holdout"]["holdout_metrics"], indent=2, sort_keys=True))


@app.command("real-validate")
def real_validate(
    data_root: Path = typer.Option(..., "--data-root", exists=True, file_okay=False),
    output: Path = typer.Option(_default("results/real/public_v0_4_0"), "--output", "-o"),
    controls: Path = typer.Option(_default("data/controls/controls.yaml"), exists=True),
    rules: Path = typer.Option(_default("config/evidence_classes.yaml"), exists=True),
    benchmark_config: Path = typer.Option(_default("config/real_benchmark.yaml"), exists=True),
) -> None:
    """Build and benchmark the frozen public real-data cohorts."""

    payload = run_real_pipeline(
        data_root=data_root,
        output_root=output,
        controls_path=controls,
        rules_path=rules,
        benchmark_config_path=benchmark_config,
    )
    typer.echo(json.dumps(payload, indent=2, sort_keys=True))


@app.command("tissue-validate")
def tissue_validate(
    crohn_tar: Path = typer.Option(..., "--crohn-tar", exists=True, dir_okay=False),
    psoriasis_directory: Path = typer.Option(
        ..., "--psoriasis-directory", exists=True, file_okay=False
    ),
    psoriasis_pseudobulk: Path | None = typer.Option(
        None, "--psoriasis-pseudobulk", exists=True, dir_okay=False
    ),
    psoriasis_qc: Path | None = typer.Option(None, "--psoriasis-qc", exists=True, dir_okay=False),
    output: Path = typer.Option(_default("results/real/public_v0_5_0"), "--output", "-o"),
    controls: Path = typer.Option(_default("data/controls/controls.yaml"), exists=True),
    tissue_config: Path = typer.Option(_default("config/tissue_extension.yaml"), exists=True),
    baseline_benchmark: Path = typer.Option(
        _default("results/real/public_v0_4_0/benchmark.json"), exists=True
    ),
) -> None:
    """Run paired Crohn and psoriasis tissue-inflammation validation."""

    payload = run_tissue_pipeline(
        crohn_tar=crohn_tar,
        psoriasis_directory=psoriasis_directory,
        psoriasis_pseudobulk_path=psoriasis_pseudobulk,
        psoriasis_qc_path=psoriasis_qc,
        controls_path=controls,
        tissue_config_path=tissue_config,
        baseline_benchmark_path=baseline_benchmark,
        output_root=output,
    )
    typer.echo(json.dumps(payload["tissue_benchmark"]["tier1_metrics"], indent=2))


@app.command()
def sensitivity(
    input_path: Path = typer.Option(..., "--input", "-i", exists=True, dir_okay=False),
    output: Path = typer.Option(_default("results/sensitivity.json"), "--output", "-o"),
    rules: Path = typer.Option(_default("config/evidence_classes.yaml"), exists=True),
) -> None:
    """Run detection, fraction and LFC threshold sensitivity."""

    cells = pd.read_csv(input_path)
    payload = run_sensitivity(cells, rules_path=rules)
    write_json(output, payload)
    typer.echo(f"Sensitivity: {payload['ranking_stability']} ({payload['n_settings']} settings).")


@app.command()
def report(
    results: Path = typer.Option(..., "--results", exists=True, dir_okay=False),
    output: Path = typer.Option(_default("results/report.html"), "--output", "-o"),
) -> None:
    """Render an honest, dependency-free HTML report."""

    payload = json.loads(results.read_text(encoding="utf-8"))
    render_report(payload, output)
    typer.echo(f"Wrote {output}.")


@app.command()
def fetch(
    registry: Path = typer.Option(_default("data/registry"), exists=True, file_okay=False),
    cache_dir: Path = typer.Option(_default(".cache/alia")),
) -> None:
    """Fetch checksum-pinned files for verified/default manifests only."""

    payload = fetch_default_registry(registry, cache_dir)
    typer.echo(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    app()
