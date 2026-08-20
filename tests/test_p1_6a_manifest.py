"""Response-blind tests for the P1.6A confirmatory lock."""

from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path

import pytest

from acoustic_ms.p1_campaign import load_p1_6_configuration
from acoustic_ms.paper_pipeline import (
    ManifestValidationError,
    P1_FROZEN_MANIFEST_SHA256,
    P1_HISTORICAL_MANIFEST_SHA256,
    load_json_yaml,
    manifest_file_sha256,
    manifest_sha256,
    validate_manifest_file,
    verify_manifest_sha256,
)


ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = ROOT / "campaigns" / "p1" / "campaign_manifest.yaml"
PILOT = ROOT / "campaigns" / "p1" / "pilot_manifest.yaml"
PILOT_OUTPUT = ROOT / "campaigns" / "p1" / "pilot"
NEW_SHA256 = "3a63fd66501f8a7ec967ba26fbb8a46f8219fcd65ef1aca4c3ae999803ace6fe"
P1_4_SHA256 = "9d360de6e61d901cff3f84c477f367773251103db12386dbb8156bd1ec2addca"
IMMUTABLE_PROJECTION_SHA256 = (
    "48a3d19b6f3412d6b96e4a028e77e0fda4fca023111073f48e393c460d0a8137"
)
PILOT_BYTES_SHA256 = "37db8cb4d0342fec4a72d18f3dce87e56f21970dc419abb3dd242c976d97361e"
PILOT_ARTIFACT_SHA256 = {
    "data_raw.csv": "a4416cae58654371ddcf680ce1a8470ab227c58760b8e1d507893e91883574da",
    "data_derived.csv": "ccd1f7a1aac92a25c51dfb822530fc55f2c27d77c0d85bbcc4215397a3bf2026",
    "data_plot.csv": "9a94fb1203ae122f89a3eb3f49074bea89f1fd7879224b58c6af7e9cafc38424",
    "failures.csv": "cd60af766e7340aa04e1b3a1fb2f4b7948f7901163ea75fb5ac42ef4e93e3e8f",
    "performance.csv": "9bb573c524a31a183856289610dc91478d14c8cacb254c9f3a85a2ad00048222",
}


def _encode(document: dict[str, object]) -> bytes:
    return (json.dumps(document, indent=2, ensure_ascii=False) + "\n").encode(
        "utf-8"
    )


def _rehash(document: dict[str, object]) -> bytes:
    value = deepcopy(document)
    value["provenance"]["manifest_sha256"] = "0" * 64
    digest = manifest_sha256(_encode(value))
    value["provenance"]["manifest_sha256"] = digest
    return _encode(value)


def _immutable_projection(document: dict[str, object]) -> bytes:
    value = deepcopy(document)
    value.pop("provenance")
    value.pop("resources")
    for case in value["cases"]:
        case.pop("enabled")
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")


def test_new_lock_resources_enablement_and_budget_are_exact() -> None:
    campaign = validate_manifest_file(CAMPAIGN, kind="campaign")
    cases = campaign["cases"]

    assert campaign["provenance"]["manifest_sha256"] == NEW_SHA256
    assert P1_FROZEN_MANIFEST_SHA256["p1_dimer_confirmatory"] == NEW_SHA256
    assert manifest_file_sha256(CAMPAIGN) == NEW_SHA256
    assert verify_manifest_sha256(
        CAMPAIGN.read_bytes(), expected_sha256=NEW_SHA256
    ) == NEW_SHA256
    assert P1_HISTORICAL_MANIFEST_SHA256[
        "p1_dimer_confirmatory_p1_4"
    ] == P1_4_SHA256
    assert len(cases) == 102
    assert [case["case_order"] for case in cases] == list(range(1, 103))
    assert all(case["enabled"] for case in cases)
    assert sum(case["parameters"]["evidence_role"] == "primary" for case in cases) == 96
    assert sum(
        case["parameters"]["evidence_role"] == "rotational_audit"
        for case in cases
    ) == 6
    assert sum(
        case["parameters"]["include_in_scientific_tables"] for case in cases
    ) == 96
    assert campaign["resources"] == {
        "worker_count": 1,
        "blas_threads": 1,
        "peak_rss_bytes_per_case": 4 * 1024**3,
        "wall_seconds_per_case": 1800,
        "wall_seconds_campaign": 64800,
        "limits_status": "frozen",
    }
    assert 102 * 1800 == 183600
    assert 64800 == 18 * 3600
    assert 102 * (21 - 2 + 1) == 2040


def test_only_authorized_manifest_fields_changed_since_p1_4() -> None:
    campaign = load_json_yaml(CAMPAIGN)
    assert sha256(_immutable_projection(campaign)).hexdigest() == (
        IMMUTABLE_PROJECTION_SHA256
    )


def test_exact_twelve_strata_and_six_audit_twins_are_frozen() -> None:
    campaign = load_json_yaml(CAMPAIGN)
    cases = campaign["cases"]
    by_id = {case["case_id"]: case for case in cases}
    primaries = [
        case for case in cases if case["parameters"]["evidence_role"] == "primary"
    ]
    audits = [
        case
        for case in cases
        if case["parameters"]["evidence_role"] == "rotational_audit"
    ]
    strata = {
        (case["parameters"]["ka"], case["parameters"]["material_id"])
        for case in primaries
    }

    assert len(strata) == 12
    assert len(audits) == 6
    assert all(audit["parameters"]["twin_case_id"] in by_id for audit in audits)
    assert all(
        by_id[audit["parameters"]["twin_case_id"]]["parameters"][
            "evidence_role"
        ]
        == "primary"
        for audit in audits
    )


@pytest.mark.parametrize("mutation", ["bytes", "order", "enabled", "resources", "twin"])
def test_public_lock_rejects_every_frozen_mutation(
    tmp_path: Path, mutation: str
) -> None:
    document = load_json_yaml(CAMPAIGN)
    if mutation == "bytes":
        path = tmp_path / "mutated.yaml"
        path.write_bytes(CAMPAIGN.read_bytes().replace(b'"title": ', b'"title":  ', 1))
        message = "stored hash"
    else:
        if mutation == "order":
            document["cases"][0], document["cases"][1] = (
                document["cases"][1],
                document["cases"][0],
            )
        elif mutation == "enabled":
            document["cases"][0]["enabled"] = False
        elif mutation == "resources":
            document["resources"]["wall_seconds_campaign"] += 1
        else:
            document["cases"][96]["parameters"]["twin_case_id"] = document[
                "cases"
            ][1]["case_id"]
        path = tmp_path / f"{mutation}.yaml"
        path.write_bytes(_rehash(document))
        message = "frozen public lock|case_order|twin physical"
    with pytest.raises(ManifestValidationError, match=message):
        validate_manifest_file(path, kind="campaign")


def test_pilot_manifest_and_five_artifacts_are_byte_immutable() -> None:
    assert sha256(PILOT.read_bytes()).hexdigest() == PILOT_BYTES_SHA256
    assert {
        path.name: sha256(path.read_bytes()).hexdigest()
        for path in sorted(PILOT_OUTPUT.iterdir())
    } == PILOT_ARTIFACT_SHA256


def test_configuration_validation_performs_no_campaign_solve_or_output() -> None:
    configuration = load_p1_6_configuration(ROOT)
    assert len(configuration.cases) == 102
    assert configuration.manifest["resources"]["limits_status"] == "frozen"
    assert not any((ROOT / relative).exists() for relative in (
        "campaigns/p1/data_raw.csv",
        "campaigns/p1/data_derived.csv",
        "campaigns/p1/data_plot.csv",
        "campaigns/p1/failures.csv",
        "campaigns/p1/performance.csv",
    ))
