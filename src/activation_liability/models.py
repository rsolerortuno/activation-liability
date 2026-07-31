"""Validated configuration and registry models."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, HttpUrl, model_validator

StimulusAxis = Literal[
    "IFN_I",
    "IFN_II",
    "TLR",
    "TNF",
    "lymphocyte_activation",
    "tissue_damage",
    "other",
]
Modality = Literal["RNA", "RNA+ADT", "sorted_bulk"]


class DownloadSpec(BaseModel):
    """A checksum-pinned direct download."""

    model_config = ConfigDict(extra="forbid")

    url: HttpUrl
    sha256: str = Field(pattern=r"^[a-fA-F0-9]{64}$")
    filename: str = Field(min_length=1)


class StudyManifest(BaseModel):
    """One real-data study registry entry."""

    model_config = ConfigDict(extra="forbid")

    accession: str = Field(min_length=1)
    source_repository: str = Field(min_length=1)
    publication_doi: str = Field(min_length=1)
    tissue: str = Field(min_length=1)
    disease_context: str = Field(min_length=1)
    stimulus_axis: StimulusAxis
    platform: str = Field(min_length=1)
    modality: Modality
    donor_count: int = Field(ge=0)
    condition_labels: dict[str, str]
    licence: str = Field(min_length=1)
    verified: bool
    default: bool = False
    downloads: list[DownloadSpec] = Field(default_factory=list)
    notes: str = ""

    @model_validator(mode="after")
    def enforce_truthful_defaults(self) -> StudyManifest:
        required = {"resting", "activated"}
        if set(self.condition_labels) != required:
            raise ValueError("condition_labels must contain exactly resting and activated")
        if self.default and not self.verified:
            raise ValueError("an unverified manifest cannot be enabled by default")
        if self.verified:
            unresolved = {
                self.publication_doi,
                self.licence,
                *self.condition_labels.values(),
            }
            if any(value.startswith("NOT_VERIFIED") for value in unresolved):
                raise ValueError("verified manifests cannot contain NOT_VERIFIED fields")
        return self


def load_manifest(path: Path) -> StudyManifest:
    """Load and validate a manifest YAML file."""

    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    return StudyManifest.model_validate(payload)
