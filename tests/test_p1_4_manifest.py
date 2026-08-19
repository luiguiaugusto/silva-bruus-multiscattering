from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path

import pytest

from acoustic_ms.paper_pipeline import (
    ManifestValidationError,
    P1_FROZEN_MANIFEST_SHA256,
    load_json_yaml,
    manifest_file_sha256,
    manifest_sha256,
    validate_executable_manifest_file,
    validate_manifest_file,
    verify_manifest_sha256,
)


ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = ROOT / "campaigns" / "p1" / "campaign_manifest.yaml"
PILOT = ROOT / "campaigns" / "p1" / "pilot_manifest.yaml"
MULTICASE_EXAMPLE = (
    ROOT / "campaigns" / "templates" / "campaign_manifest.multicase.example.yaml"
)
ZERO_SHA256 = "0" * 64


def _encode(document: dict[str, object]) -> bytes:
    return (
        json.dumps(document, indent=2, ensure_ascii=False) + "\n"
    ).encode("utf-8")


def _rehash(document: dict[str, object]) -> bytes:
    updated = deepcopy(document)
    updated["provenance"]["manifest_sha256"] = ZERO_SHA256
    normalized = _encode(updated)
    updated["provenance"]["manifest_sha256"] = manifest_sha256(normalized)
    return _encode(updated)


@pytest.mark.parametrize(
    ("path", "campaign_id"),
    [
        (CAMPAIGN, "p1_dimer_confirmatory"),
        (PILOT, "p1_dimer_resource_pilot"),
    ],
)
def test_frozen_manifest_hashes_recompute_from_exact_bytes(
    path: Path, campaign_id: str
) -> None:
    exact_bytes = path.read_bytes()
    document = validate_manifest_file(path, kind="campaign")
    expected = P1_FROZEN_MANIFEST_SHA256[campaign_id]

    assert document["provenance"]["manifest_sha256"] == expected
    assert manifest_sha256(exact_bytes) == expected
    assert manifest_file_sha256(path) == expected
    assert verify_manifest_sha256(
        exact_bytes, expected_sha256=expected
    ) == expected


def test_hash_normalizes_only_its_own_value() -> None:
    exact_bytes = PILOT.read_bytes()
    expected = P1_FROZEN_MANIFEST_SHA256["p1_dimer_resource_pilot"]
    zeroed_field = exact_bytes.replace(
        expected.encode("ascii"), ZERO_SHA256.encode("ascii"), 1
    )

    assert manifest_sha256(zeroed_field) == expected
    with pytest.raises(ManifestValidationError, match="stored hash"):
        verify_manifest_sha256(zeroed_field)


def test_any_non_hash_byte_mutation_invalidates_manifest() -> None:
    exact_bytes = CAMPAIGN.read_bytes()
    mutated = exact_bytes.replace(b'"title": ', b'"title":  ', 1)

    assert manifest_sha256(mutated) != manifest_sha256(exact_bytes)
    with pytest.raises(ManifestValidationError, match="stored hash"):
        verify_manifest_sha256(mutated)


def test_rehashed_physical_mutation_is_rejected_by_public_lock(
    tmp_path: Path,
) -> None:
    document = load_json_yaml(CAMPAIGN)
    parameters = document["cases"][2]["parameters"]
    parameters["ka"] = 0.051
    parameters["k_rad_m"] = 0.051
    mutated = _rehash(document)

    assert verify_manifest_sha256(mutated) == manifest_sha256(mutated)
    with pytest.raises(ManifestValidationError, match="frozen public lock"):
        verify_manifest_sha256(
            mutated,
            expected_sha256=P1_FROZEN_MANIFEST_SHA256[
                "p1_dimer_confirmatory"
            ],
        )

    path = tmp_path / "mutated-physical.yaml"
    path.write_bytes(mutated)
    with pytest.raises(ManifestValidationError, match="frozen public lock"):
        validate_manifest_file(path, kind="campaign")


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("duplicate", "case_id values must be unique"),
        ("reorder", "case_order"),
        ("twin", "twin physical parameters differ"),
    ],
)
def test_structural_mutations_are_rejected_before_execution(
    tmp_path: Path, mutation: str, message: str
) -> None:
    document = load_json_yaml(CAMPAIGN)
    if mutation == "duplicate":
        document["cases"][1]["case_id"] = document["cases"][0]["case_id"]
    elif mutation == "reorder":
        document["cases"][0], document["cases"][1] = (
            document["cases"][1],
            document["cases"][0],
        )
    else:
        audit = document["cases"][96]
        audit["parameters"]["twin_case_id"] = document["cases"][1]["case_id"]

    path = tmp_path / f"{mutation}.yaml"
    path.write_bytes(_rehash(document))
    with pytest.raises(ManifestValidationError, match=message):
        validate_manifest_file(path, kind="campaign")


@pytest.mark.parametrize("field", ["manifest_sha256", "task"])
def test_nonplanned_manifest_rejects_tbd_everywhere(
    tmp_path: Path, field: str
) -> None:
    document = load_json_yaml(PILOT)
    document["provenance"][field] = "TBD"
    path = tmp_path / f"tbd-{field}.yaml"
    path.write_bytes(_encode(document))

    with pytest.raises(ManifestValidationError, match="TBD is allowed only"):
        validate_manifest_file(path, kind="campaign")


def test_ids_twins_enablement_and_scientific_exclusion_are_frozen() -> None:
    campaign = validate_manifest_file(CAMPAIGN, kind="campaign")
    pilot = validate_manifest_file(PILOT, kind="campaign")
    cases = campaign["cases"]
    ids = tuple(case["case_id"] for case in cases)
    audits = [
        case
        for case in cases
        if case["parameters"]["evidence_role"] == "rotational_audit"
    ]
    by_id = {case["case_id"]: case for case in cases}

    assert len(ids) == len(set(ids)) == 102
    assert tuple(case["case_order"] for case in cases) == tuple(range(1, 103))
    assert len(audits) == 6
    assert all(
        audit["parameters"]["twin_case_id"] in by_id
        and by_id[audit["parameters"]["twin_case_id"]]["parameters"][
            "evidence_role"
        ] == "primary"
        for audit in audits
    )
    assert not any(case["enabled"] for case in cases)
    assert [case["case_id"] for case in pilot["cases"] if case["enabled"]] == [
        "p1_pilot_rigid_ka010_d0210_t000"
    ]
    assert pilot["cases"][0]["parameters"][
        "include_in_scientific_tables"
    ] is False
    assert pilot["cases"][0]["case_id"] not in ids
    assert campaign["resources"]["limits_status"] == "provisional"
    assert pilot["resources"]["limits_status"] == "provisional"


def test_execution_guard_authorizes_only_the_p1_5_pilot() -> None:
    pilot = validate_executable_manifest_file(
        PILOT,
        expected_campaign_id="p1_dimer_resource_pilot",
    )
    assert [case["case_id"] for case in pilot["cases"] if case["enabled"]] == [
        "p1_pilot_rigid_ka010_d0210_t000"
    ]

    with pytest.raises(ManifestValidationError, match="no enabled cases"):
        validate_executable_manifest_file(
            CAMPAIGN,
            expected_campaign_id="p1_dimer_confirmatory",
        )
    with pytest.raises(ManifestValidationError, match="status"):
        validate_executable_manifest_file(
            MULTICASE_EXAMPLE,
            expected_campaign_id="p1_multicase_example",
        )


def test_execution_guard_rejects_invalid_hash(tmp_path: Path) -> None:
    mutated = PILOT.read_bytes().replace(b'"title": ', b'"title":  ', 1)
    path = tmp_path / "invalid-pilot.yaml"
    path.write_bytes(mutated)

    with pytest.raises(ManifestValidationError, match="stored hash"):
        validate_executable_manifest_file(
            path,
            expected_campaign_id="p1_dimer_resource_pilot",
        )
