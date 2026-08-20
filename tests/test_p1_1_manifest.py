from __future__ import annotations

from math import pi
from pathlib import Path

import pytest

from acoustic_ms.paper_pipeline import validate_manifest_file


ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = ROOT / "campaigns" / "p1" / "campaign_manifest.yaml"
PILOT = ROOT / "campaigns" / "p1" / "pilot_manifest.yaml"

KA_LEVELS = ((0.05, "ka005"), (0.1, "ka010"))
MATERIALS = (
    ("fluid_f000_f1010", "fluid", 0.0, True, 0.1),
    ("fluid_f000_f1040", "fluid", 0.0, True, 0.4),
    ("fluid_f000_f1060", "fluid", 0.0, True, 0.6),
    ("fluid_f020_f1060", "fluid", 0.2, True, 0.6),
    ("fluid_f000_f1080", "fluid", 0.0, True, 0.8),
    ("rigid_boundary", "rigid", 0.0, False, 1.0),
)
DISTANCES = (
    (2.1, "d0210"),
    (2.25, "d0225"),
    (2.5, "d0250"),
    (3.0, "d0300"),
    (4.0, "d0400"),
    (6.0, "d0600"),
    (8.0, "d0800"),
    (10.0, "d1000"),
)
AUDIT_CASES = (
    ("fluid_f000_f1010", 0.05, "ka005", 2.1, "d0210"),
    ("fluid_f000_f1040", 0.05, "ka005", 10.0, "d1000"),
    ("fluid_f000_f1060", 0.1, "ka010", 2.1, "d0210"),
    ("fluid_f020_f1060", 0.1, "ka010", 10.0, "d1000"),
    ("fluid_f000_f1080", 0.05, "ka005", 2.1, "d0210"),
    ("rigid_boundary", 0.1, "ka010", 10.0, "d1000"),
)


def _primary_id(ka_tag: str, material_id: str, distance_tag: str) -> str:
    return f"p1_dimer_{ka_tag}_{material_id}_{distance_tag}_t000"


def _expected_case_ids() -> tuple[str, ...]:
    primary = tuple(
        _primary_id(ka_tag, material_id, distance_tag)
        for _, ka_tag in KA_LEVELS
        for material_id, *_ in MATERIALS
        for _, distance_tag in DISTANCES
    )
    audits = tuple(
        f"p1_dimer_{ka_tag}_{material_id}_{distance_tag}_t045_audit"
        for material_id, _, ka_tag, _, distance_tag in AUDIT_CASES
    )
    return primary + audits


def test_confirmatory_manifest_freezes_exact_grid_ids_and_order() -> None:
    document = validate_manifest_file(CAMPAIGN, kind="campaign")
    cases = document["cases"]
    expected_ids = _expected_case_ids()

    assert document["schema_version"] == "1.1.0"
    assert document["classification"] == "confirmatory_new"
    assert document["status"] == "preregistered"
    assert document["provenance"]["manifest_sha256"] == "3a63fd66501f8a7ec967ba26fbb8a46f8219fcd65ef1aca4c3ae999803ace6fe"
    assert len(cases) == 102
    assert tuple(case["case_id"] for case in cases) == expected_ids
    assert tuple(case["case_order"] for case in cases) == tuple(range(1, 103))
    assert all(case["enabled"] is True for case in cases)

    primary = [
        case
        for case in cases
        if case["parameters"]["evidence_role"] == "primary"
    ]
    audits = [
        case
        for case in cases
        if case["parameters"]["evidence_role"] == "rotational_audit"
    ]
    assert len(primary) == 96
    assert len(audits) == 6
    assert all(
        case["parameters"]["include_in_scientific_tables"] is True
        for case in primary
    )
    assert all(
        case["parameters"]["include_in_scientific_tables"] is False
        for case in audits
    )


def test_confirmatory_manifest_has_exact_per_case_physical_grid() -> None:
    document = validate_manifest_file(CAMPAIGN, kind="campaign")
    primary_parameters = [
        case["parameters"]
        for case in document["cases"]
        if case["parameters"]["evidence_role"] == "primary"
    ]

    expected = {
        (ka, material_id, model, f0, f0_applicable, f1, distance)
        for ka, _ in KA_LEVELS
        for material_id, model, f0, f0_applicable, f1 in MATERIALS
        for distance, _ in DISTANCES
    }
    observed = {
        (
            parameters["ka"],
            parameters["material_id"],
            parameters["material_model"],
            parameters["f0"],
            parameters["f0_applicable"],
            parameters["f1"],
            parameters["distance_ratio"],
        )
        for parameters in primary_parameters
    }
    assert observed == expected
    assert all(
        parameters["k_rad_m"] == parameters["ka"]
        and parameters["theta_rad"] == 0.0
        for parameters in primary_parameters
    )


def test_rotational_audits_have_exact_zero_angle_twins() -> None:
    document = validate_manifest_file(CAMPAIGN, kind="campaign")
    by_id = {case["case_id"]: case for case in document["cases"]}
    audits = [
        case
        for case in document["cases"]
        if case["parameters"]["evidence_role"] == "rotational_audit"
    ]

    assert {
        (
            case["parameters"]["material_id"],
            case["parameters"]["ka"],
            case["parameters"]["distance_ratio"],
        )
        for case in audits
    } == {
        (material_id, ka, distance)
        for material_id, ka, _, distance, _ in AUDIT_CASES
    }
    for audit in audits:
        parameters = audit["parameters"]
        twin = by_id[parameters["twin_case_id"]]["parameters"]
        assert parameters["theta_rad"] == pytest.approx(pi / 4.0)
        assert twin["theta_rad"] == 0.0
        assert parameters["ka"] == twin["ka"]
        assert parameters["distance_ratio"] == twin["distance_ratio"]
        for field in (
            "k_rad_m",
            "material_id",
            "material_model",
            "f0",
            "f0_applicable",
            "f1",
        ):
            assert parameters[field] == twin[field]


def test_rigid_cases_mark_f0_as_api_sentinel_not_physical_contrast() -> None:
    document = validate_manifest_file(CAMPAIGN, kind="campaign")
    rigid = [
        case["parameters"]
        for case in document["cases"]
        if case["parameters"]["material_model"] == "rigid"
    ]

    assert len(rigid) == 17
    assert all(parameters["material_id"] == "rigid_boundary" for parameters in rigid)
    assert all(parameters["f1"] == 1.0 for parameters in rigid)
    assert all(parameters["f0"] == 0.0 for parameters in rigid)
    assert all(parameters["f0_applicable"] is False for parameters in rigid)


def test_numerical_and_p1_6a_resource_policy_is_frozen_and_enabled() -> None:
    document = validate_manifest_file(CAMPAIGN, kind="campaign")
    numerical = document["numerical"]
    resources = document["resources"]

    assert numerical["lmax_min"] == 2
    assert numerical["lmax_max"] == 21
    assert numerical["minimum_stop_lmax"] == 5
    assert numerical["convergence_tolerance"] == 1.0e-5
    assert numerical["consecutive_passes"] == 2
    assert numerical["require_all_applicable_channels"] is True
    assert numerical["required_channels"] == [
        "total",
        "interaction",
        "external_scattered",
        "scattered_scattered",
    ]
    assert resources == {
        "worker_count": 1,
        "blas_threads": 1,
        "peak_rss_bytes_per_case": 4 * 1024**3,
        "wall_seconds_per_case": 30 * 60,
        "wall_seconds_campaign": 18 * 60 * 60,
        "limits_status": "frozen",
    }


def test_resource_pilot_is_separate_enabled_p1_5_development_evidence() -> None:
    campaign = validate_manifest_file(CAMPAIGN, kind="campaign")
    pilot = validate_manifest_file(PILOT, kind="campaign")
    case = pilot["cases"][0]
    parameters = case["parameters"]

    assert pilot["classification"] == "development"
    assert pilot["status"] == "preregistered"
    assert pilot["provenance"]["manifest_sha256"] == "d8f56ce20f6f0821d84fd6f36e1f76c855f63f55d809ba9a7201ba52097a43bf"
    assert len(pilot["cases"]) == 1
    assert case["case_id"] == "p1_pilot_rigid_ka010_d0210_t000"
    assert case["enabled"] is True
    assert parameters["material_model"] == "rigid"
    assert parameters["f1"] == 1.0
    assert parameters["f0_applicable"] is False
    assert parameters["ka"] == parameters["k_rad_m"] == 0.1
    assert parameters["distance_ratio"] == 2.1
    assert parameters["theta_rad"] == 0.0
    assert parameters["evidence_role"] == "resource_pilot"
    assert parameters["include_in_scientific_tables"] is False
    assert case["case_id"] not in {
        item["case_id"] for item in campaign["cases"]
    }


def test_p1_3a_rigid_pair_is_exactly_the_translated_p1_5_pilot() -> None:
    pilot = validate_manifest_file(PILOT, kind="campaign")
    case = pilot["cases"][0]
    parameters = case["parameters"]

    p1_3a = {
        "ka": 0.1,
        "f0": 0.0,
        "f1": 1.0,
        "distance_ratio": 2.1,
        "theta_rad": 0.0,
    }
    assert {field: parameters[field] for field in p1_3a} == p1_3a
    assert parameters["material_model"] == "rigid"
    assert parameters["f0_applicable"] is False

    p1_3a_positions = ((0.0, 0.0, 0.0), (2.1, 0.0, 0.0))
    pilot_positions = ((-1.05, 0.0, 0.0), (1.05, 0.0, 0.0))

    def relative(points: tuple[tuple[float, ...], ...]) -> tuple[float, ...]:
        return tuple(
            points[1][axis] - points[0][axis]
            for axis in range(3)
        )

    assert relative(p1_3a_positions) == relative(pilot_positions)
    assert pilot["geometry"]["positions_source"] == (
        "r1=(-1.05,0,0), r2=(1.05,0,0); fixed development pilot only"
    )
