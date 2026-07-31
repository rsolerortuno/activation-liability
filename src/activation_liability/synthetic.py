"""Synthetic fixtures with planted activation, batch, donor and dropout effects."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml

CELL_TYPES = (
    "T_cell",
    "B_cell",
    "Plasma_cell",
    "Myeloid",
    "Endothelial",
    "Fibroblast",
    "Epithelial",
)
STIMULUS_BY_CELL = {
    "T_cell": "lymphocyte_activation",
    "B_cell": "lymphocyte_activation",
    "Plasma_cell": "lymphocyte_activation",
    "Myeloid": "IFN_II",
    "Endothelial": "IFN_II",
    "Fibroblast": "tissue_damage",
    "Epithelial": "IFN_II",
}
LINEAGE_TARGETS = {
    "B_cell": {"CD19", "MS4A1", "CD79A"},
    "Plasma_cell": {"TNFRSF17", "GPRC5D"},
    "T_cell": {"CD3E", "CD4", "CD8A"},
    "Epithelial": {"CLDN18", "CLDN6", "DLL3", "ERBB2", "FOLR1"},
}
T_CELL_INDUCED = {
    "IL2RA",
    "TNFRSF9",
    "TNFRSF4",
    "ICOS",
    "PDCD1",
    "CTLA4",
    "LAG3",
    "HAVCR2",
    "TIGIT",
    "CD38",
    "CD69",
    "TNFRSF18",
    "ENTPD1",
    "CD70",
}
IFN_INDUCED = {
    "CD274",
    "PDCD1LG2",
    "MICA",
    "MICB",
    "HLA-DRA",
    "ICAM1",
    "VCAM1",
    "SIGLEC1",
    "IDO1",
}
DAMAGE_INDUCED = {"FAP", "TNFRSF12A"}


def control_targets(path: Path) -> tuple[list[str], dict[str, str]]:
    """Read target names and labels without using them for statistical fitting."""

    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    labels = {entry["target"]: entry["label"] for entry in payload["controls"]}
    return sorted(labels), labels


def _is_induced(target: str, cell_type: str) -> bool:
    if cell_type == "T_cell" and target in T_CELL_INDUCED:
        return True
    if cell_type in {"Myeloid", "Endothelial", "Epithelial"} and target in IFN_INDUCED:
        return True
    return cell_type == "Fibroblast" and target in DAMAGE_INDUCED


def generate_synthetic(
    controls_path: Path,
    *,
    seed: int = 20260730,
    studies: int = 3,
    donors_per_study: int = 5,
    cells_per_condition: int = 12,
) -> pd.DataFrame:
    """Generate a long-form count table with analytically planted qualitative truth.

    The study/target batch multiplier is deliberately much larger than donor noise.
    It cancels in the paired within-study contrast but corrupts invalid cross-study
    absolute contrasts.
    """

    targets, labels = control_targets(controls_path)
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []

    batch_log2 = {
        (study, target): rng.normal(0.0, 2.2) for study in range(studies) for target in targets
    }
    donor_log2 = {
        (study, donor, target): rng.normal(0.0, 0.30)
        for study in range(studies)
        for donor in range(donors_per_study)
        for target in targets
    }

    for study in range(studies):
        study_name = f"SYNTH_{study + 1}"
        for donor in range(donors_per_study):
            donor_name = f"D{donor + 1}"
            for cell_type in CELL_TYPES:
                stimulus = STIMULUS_BY_CELL[cell_type]
                for condition in ("resting", "activated"):
                    for cell_index in range(cells_per_condition):
                        cell_id = (
                            f"{study_name}_{donor_name}_{cell_type}_{condition}_{cell_index:03d}"
                        )
                        library_size = int(rng.lognormal(np.log(9000), 0.18))
                        adt_library_size = int(rng.lognormal(np.log(1800), 0.15))
                        for target in targets:
                            baseline = 0.06
                            if target in LINEAGE_TARGETS.get(cell_type, set()):
                                baseline = 8.0
                            elif labels[target] == "positive":
                                baseline = 0.12

                            activation_log2 = 0.0
                            if _is_induced(target, cell_type) and condition == "activated":
                                activation_log2 = (
                                    2.0
                                    if target
                                    in {
                                        "IL2RA",
                                        "TNFRSF9",
                                        "TNFRSF4",
                                        "CD69",
                                        "CD274",
                                        "PDCD1LG2",
                                        "HLA-DRA",
                                        "ICAM1",
                                        "IDO1",
                                    }
                                    else 1.35
                                )
                            # Planted RNA-only technical artefact: activated T cells
                            # show a strong CD3E compositional/count shift without a
                            # matching surface-protein change.
                            if (
                                target == "CD3E"
                                and cell_type == "T_cell"
                                and condition == "activated"
                            ):
                                activation_log2 = 4.05

                            mean_count = (
                                baseline
                                * 2.0 ** batch_log2[(study, target)]
                                * 2.0 ** donor_log2[(study, donor, target)]
                                * 2.0**activation_log2
                                * library_size
                                / 9000.0
                            )
                            rna_count = int(rng.poisson(max(mean_count, 0.0001)))

                            protein_measured = target in {
                                "IL2RA",
                                "TNFRSF9",
                                "TNFRSF4",
                                "CD69",
                                "CD274",
                                "PDCD1LG2",
                                "ICAM1",
                                "CD19",
                                "MS4A1",
                                "TNFRSF17",
                                "CD3E",
                                "ERBB2",
                            }
                            adt_count: float
                            if protein_measured:
                                adt_base = 0.12
                                if target in LINEAGE_TARGETS.get(cell_type, set()):
                                    adt_base = 5.0
                                protein_lfc = activation_log2 * 0.90
                                # One planted RNA-only discordance demonstrates the flag.
                                if target == "ICAM1":
                                    protein_lfc = 0.15 if condition == "activated" else 0.0
                                if target == "CD3E":
                                    protein_lfc = 0.0
                                adt_mean = (
                                    adt_base
                                    * 2.0 ** (0.35 * batch_log2[(study, target)])
                                    * 2.0**protein_lfc
                                    * adt_library_size
                                    / 1800.0
                                )
                                adt_count = float(rng.poisson(max(adt_mean, 0.0001)))
                            else:
                                adt_count = float("nan")

                            tumour_reference_cpm = 125.0 if labels[target] == "negative" else 80.0
                            rows.append(
                                {
                                    "study": study_name,
                                    "donor": donor_name,
                                    "cell_id": cell_id,
                                    "cell_type": cell_type,
                                    "condition": condition,
                                    "stimulus": stimulus,
                                    "target": target,
                                    "rna_count": rna_count,
                                    "library_size": library_size,
                                    "adt_count": adt_count,
                                    "adt_library_size": adt_library_size,
                                    "tumour_reference_cpm": tumour_reference_cpm,
                                }
                            )
    frame = pd.DataFrame.from_records(rows)
    return frame.sort_values(
        ["study", "donor", "cell_type", "condition", "cell_id", "target"]
    ).reset_index(drop=True)
