from __future__ import annotations

from copy import deepcopy
import json
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
MULTICASE_EXAMPLE = TEMPLATES / "campaign_manifest.multicase.example.yaml"


@pytest.mark.parametrize(
    ("name", "kind", "version"),
    [
        ("campaign_manifest.example.yaml", "campaign", "1.0.0"),
        ("campaign_manifest.multicase.example.yaml", "campaign", "1.1.0"),
        ("figure_manifest.example.yaml", "figure", "1.0.0"),
    ],
)
def test_manifest_examples_satisfy_their_schemas(
    name: str, kind: str, version: str
) -> None:
    document = validate_manifest_file(TEMPLATES / name, kind=kind)
    assert document["schema_version"] == version


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


@pytest.mark.parametrize(
    ("mutator", "message"),
    [
        (
            lambda value: value["cases"][0]["parameters"].update({"ka": 0.2}),
            r"physical\.radius_m \* k_rad_m",
        ),
        (
            lambda value: value["cases"][1].update({"case_order": 1}),
            "case_order",
        ),
        (
            lambda value: value["cases"][0].update({"enabled": True}),
            "planned campaigns",
        ),
        (
            lambda value: value["cases"][1]["parameters"].update(
                {"twin_case_id": "missing_case"}
            ),
            "unknown case",
        ),
        (
            lambda value: value["cases"][0]["parameters"].update(
                {
                    "material_model": "rigid",
                    "f1": 1.0,
                    "f0_applicable": True,
                }
            ),
            "rigid material requires false",
        ),
        (
            lambda value: value["cases"][0]["parameters"].update(
                {
                    "material_model": "rigid",
                    "f0": 0.2,
                    "f1": 1.0,
                    "f0_applicable": False,
                }
            ),
            "rigid API sentinel must be zero",
        ),
        (
            lambda value: value["numerical"].update({"minimum_stop_lmax": 22}),
            "minimum_stop_lmax",
        ),
    ],
)
def test_campaign_v1_1_semantic_mutations_fail(mutator, message: str) -> None:
    document = load_json_yaml(MULTICASE_EXAMPLE)
    schema = load_json_yaml(SCHEMAS / "campaign_manifest.v1.1.schema.json")
    mutator(document)
    with pytest.raises(ManifestValidationError, match=message):
        validate_manifest(document, schema)


def test_validator_rejects_unknown_campaign_schema_version(tmp_path: Path) -> None:
    document = load_json_yaml(MULTICASE_EXAMPLE)
    document["schema_version"] = "1.2.0"
    manifest = tmp_path / "unknown.yaml"
    manifest.write_text(json.dumps(document), encoding="utf-8")
    with pytest.raises(ManifestValidationError, match="unsupported"):
        validate_manifest_file(manifest, kind="campaign")


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
