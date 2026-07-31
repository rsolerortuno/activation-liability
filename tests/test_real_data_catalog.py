from pathlib import Path

import yaml

CATALOG = Path("data/candidates/real_data_catalog.yaml")


def test_primary_catalog_has_paired_rna_adt_study() -> None:
    payload = yaml.safe_load(CATALOG.read_text(encoding="utf-8"))
    datasets = {item["accession"]: item for item in payload["datasets"]}
    primary = datasets["GSE157857"]
    assert primary["priority"] == 1
    assert primary["modality"] == "RNA+ADT"
    assert primary["donors_reported"] >= 3
    assert "paired" in primary["pairing"]


def test_downloadable_catalog_entries_are_https_and_not_registry_defaults() -> None:
    payload = yaml.safe_load(CATALOG.read_text(encoding="utf-8"))
    for dataset in payload["datasets"]:
        assert "default" not in dataset
        for file_spec in dataset["files"]:
            assert file_spec["url"].startswith("https://")
            assert file_spec["filename"]
