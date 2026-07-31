from __future__ import annotations


def test_primary_benchmark_recovers_controls(benchmark_result: dict[str, object]) -> None:
    result = benchmark_result
    metrics = result["all_tier1_metrics"]
    assert metrics["auroc"] >= 0.90
    assert metrics["average_precision"] >= 0.90


def test_holdout_was_not_used_for_cutoff(benchmark_result: dict[str, object]) -> None:
    holdout = benchmark_result["holdout"]
    assert set(holdout["train_positive_targets"]).isdisjoint(holdout["holdout_positive_targets"])
    assert len(holdout["holdout_positive_targets"]) == 3
    assert holdout["holdout_metrics"]["auroc"] >= 0.90


def test_invalid_cross_study_contrast_degrades(benchmark_result: dict[str, object]) -> None:
    ablation = benchmark_result["ablations"]["within_study_vs_cross_study"]
    assert (
        ablation["primary_within_study"]["auroc"]
        > ablation["INVALID_FOR_INFERENCE_cross_study"]["auroc"]
    )
    assert ablation["INVALID_FOR_INFERENCE_cross_study"]["auroc"] < 0.85


def test_cell_level_testing_inflates_calls(benchmark_result: dict[str, object]) -> None:
    ablation = benchmark_result["ablations"]["pseudobulk_vs_cell_level"]
    assert ablation["pseudobulk_significant_calls"] > 0
    assert (
        ablation["INVALID_FOR_INFERENCE_cell_level_significant_calls"]
        > ablation["pseudobulk_significant_calls"]
    )
    assert ablation["inflation_ratio"] > 2.0


def test_sensitivity_is_explicit(benchmark_result: dict[str, object]) -> None:
    sensitivity = benchmark_result["ablations"]["detection_threshold_sweep"]
    assert sensitivity["n_settings"] == 27
    assert sensitivity["ranking_stability"] in {"STABLE", "UNSTABLE"}
    assert 0.0 <= sensitivity["median_spearman"] <= 1.0


def test_protein_corroboration_improves_planted_rna_artifact(
    benchmark_result: dict[str, object],
) -> None:
    ablation = benchmark_result["ablations"]["rna_only_vs_protein_corroborated"]
    assert ablation["protein_corroborated_score"]["auroc"] > ablation["rna_only"]["auroc"]
    assert (
        ablation["protein_corroborated_score"]["average_precision"]
        > ablation["rna_only"]["average_precision"]
    )


def test_abstention_improves_non_abstained_accuracy(benchmark_result: dict[str, object]) -> None:
    curve = benchmark_result["holdout"]["coverage_vs_accuracy"]
    assert curve[1]["coverage"] < curve[0]["coverage"]
    assert curve[1]["accuracy"] > curve[0]["accuracy"]
