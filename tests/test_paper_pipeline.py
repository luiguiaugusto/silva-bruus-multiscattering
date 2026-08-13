from __future__ import annotations

from copy import deepcopy
from pathlib import Path

import pytest

from acoustic_ms.paper_pipeline import (
    ManifestValidationError,
    load_json_yaml,
    validate_manifest,
    validate_manifest_file,
)


ROOT = Path(__file__).resolve().parents[1]
SCHEMAS = ROOT / "campaigns" / "schemas"
TEMPLATES = ROOT / "campaigns" / "templates"


@pytest.mark.parametrize(
    ("name", "kind"),
    [
        ("campaign_manifest.example.yaml", "campaign"),
        ("figure_manifest.example.yaml", "figure"),
    ],
)
def test_manifest_examples_satisfy_their_schemas(name: str, kind: str) -> None:
    document = validate_manifest_file(TEMPLATES / name, kind=kind)
    assert document["schema_version"] == "1.0.0"


def test_campaign_schema_rejects_missing_required_provenance() -> None:
    document = load_json_yaml(TEMPLATES / "campaign_manifest.example.yaml")
    schema = load_json_yaml(SCHEMAS / "campaign_manifest.schema.json")
    invalid = deepcopy(document)
    del invalid["provenance"]["git_commit"]
    with pytest.raises(ManifestValidationError, match="git_commit"):
        validate_manifest(invalid, schema)


def test_campaign_schema_rejects_unstable_case_id() -> None:
    document = load_json_yaml(TEMPLATES / "campaign_manifest.example.yaml")
    schema = load_json_yaml(SCHEMAS / "campaign_manifest.schema.json")
    invalid = deepcopy(document)
    invalid["cases"][0]["case_id"] = "P1 contains spaces"
    with pytest.raises(ManifestValidationError, match="case_id"):
        validate_manifest(invalid, schema)


@pytest.mark.parametrize(
    ("mutator", "message"),
    [
        (
            lambda value: value["physical"].__setitem__("ka", 0.2),
            "radius_m \\* k_rad_m",
        ),
        (
            lambda value: value["numerical"].update(
                {"lmax_min": 14, "lmax_max": 13}
            ),
            "lmax_min",
        ),
        (
            lambda value: value["cases"].append(deepcopy(value["cases"][0])),
            "case_id values must be unique",
        ),
        (
            lambda value: value.update({"status": "preregistered"}),
            "TBD is allowed only",
        ),
    ],
)
def test_campaign_semantic_mutations_fail(mutator, message: str) -> None:
    document = load_json_yaml(TEMPLATES / "campaign_manifest.example.yaml")
    schema = load_json_yaml(SCHEMAS / "campaign_manifest.schema.json")
    mutator(document)
    with pytest.raises(ManifestValidationError, match=message):
        validate_manifest(document, schema)


def test_figure_schema_rejects_duplicate_or_unknown_formats() -> None:
    document = load_json_yaml(TEMPLATES / "figure_manifest.example.yaml")
    schema = load_json_yaml(SCHEMAS / "figure_manifest.schema.json")
    for formats in (["png", "png"], ["eps"]):
        invalid = deepcopy(document)
        invalid["formats"] = formats
        with pytest.raises(ManifestValidationError, match="formats"):
            validate_manifest(invalid, schema)


def test_figure_schema_rejects_panel_with_unknown_source_table() -> None:
    document = load_json_yaml(TEMPLATES / "figure_manifest.example.yaml")
    schema = load_json_yaml(SCHEMAS / "figure_manifest.schema.json")
    document["panels"][0]["source_table"] = "missing_table"
    with pytest.raises(ManifestValidationError, match="unknown source table"):
        validate_manifest(document, schema)


def test_validator_rejects_unknown_manifest_kind() -> None:
    with pytest.raises(ValueError, match="kind"):
        validate_manifest_file(
            TEMPLATES / "campaign_manifest.example.yaml", kind="unknown"
        )
