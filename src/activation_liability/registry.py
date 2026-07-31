"""Study registry validation."""

from __future__ import annotations

from pathlib import Path

from pydantic import ValidationError

from activation_liability.models import StudyManifest, load_manifest


def validate_registry(directory: Path) -> tuple[list[StudyManifest], list[str]]:
    """Validate every YAML manifest, returning valid entries and readable errors."""

    manifests: list[StudyManifest] = []
    errors: list[str] = []
    for path in sorted(directory.glob("*.yaml")):
        try:
            manifests.append(load_manifest(path))
        except (ValidationError, ValueError, OSError) as exc:
            errors.append(f"{path}: {exc}")
    if not list(directory.glob("*.yaml")):
        errors.append(f"{directory}: no YAML manifests found")
    return manifests, errors


def default_manifests(directory: Path) -> list[StudyManifest]:
    """Return only verified, explicitly enabled entries."""

    manifests, errors = validate_registry(directory)
    if errors:
        raise ValueError("\n".join(errors))
    return [manifest for manifest in manifests if manifest.verified and manifest.default]
