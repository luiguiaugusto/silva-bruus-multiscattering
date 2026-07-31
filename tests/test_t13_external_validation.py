"""Artifact and scientific-protocol tests for the revealed T13 holdout."""

from __future__ import annotations

import csv
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import subprocess

import numpy as np
import pytest

from acoustic_ms import (
    EXPECTED_CASE_IDS,
    cluster_family,
    external_eligibility_mask,
    maximum_geometric_coupling,
    nodal_pair_force_on_probe,
    normalized_rms_error_xyz,
    solve_multipolar_nodal_interaction_forces,
    spectral_radius_l1,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
PHASE_A_COMMIT = "af29faf89cb8f8c6883cc8bea0d44073e7caf020"
PHASE_A_HASHES = {
    "t13_holdout_manifest.csv": "25d79db59d9dd6d52c5674d0a64fe2fea351cf213a0cdcd92b45845a9ecc2b38",
    "t13_frozen_predictions.csv": "581a748dca2e5d161890284fca673ed20f2a4fbcbc0ff356d5d31db6ec8ac9c2",
    "t13_frozen_protocol.csv": "eb1878e3425ede7a2b599fd20f63550d2fdb23d177264a043d443694907dc650",
}


def _read(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _vectors(rows: list[dict[str, str]], prefix: str) -> np.ndarray:
    return np.asarray([
        [float(row[f"{prefix}_{axis}"]) for axis in ("x", "y", "z")]
        for row in rows
    ])


def test_phase_a_artifacts_are_byte_identical_to_published_commit():
    for name, expected_hash in PHASE_A_HASHES.items():
        path = DATA / name
        relative = str(path.relative_to(ROOT))
        published = subprocess.run(
            ["git", "show", f"{PHASE_A_COMMIT}:{relative}"], cwd=ROOT,
            check=True, capture_output=True,
        ).stdout
        assert sha256(path.read_bytes()).hexdigest() == expected_hash
        assert path.read_bytes() == published


def test_manifest_raw_and_summaries_have_exact_frozen_identity():
    manifest = _read("t13_holdout_manifest.csv")
    raw = _read("t13_model_e_convergence.csv")
    summary = _read("t13_case_summary.csv")
    assert tuple(row["case_id"] for row in manifest) == EXPECTED_CASE_IDS
    assert tuple(row["case_id"] for row in summary) == EXPECTED_CASE_IDS
    assert set(row["case_id"] for row in raw) == set(EXPECTED_CASE_IDS)
    assert len(manifest) == len(summary) == 24
    assert len(raw) == 205
    assert all(row["eligible"] == "true" for row in summary)


def test_all_six_strata_reconstruct_and_predictors_are_geometric():
    manifest = _read("t13_holdout_manifest.csv")
    assert {row["stratum"] for row in manifest} == {
        "n6_linear", "n6_compact", "n6_irregular",
        "n10_linear", "n10_compact", "n10_irregular",
    }
    for row in manifest:
        positions = cluster_family(
            int(row["particle_count"]), row["family"], float(row["distance_ratio"])
        )
        stored = np.asarray([
            [float(value) for value in particle.split(":")]
            for particle in row["coordinates_xyz"].split(";")
        ])
        assert np.array_equal(positions, stored)
        assert maximum_geometric_coupling(
            positions, float(row["radius"]), float(row["f1"])
        ) == pytest.approx(float(row["lambda_max"]), rel=3e-14, abs=3e-15)


def test_lambda_and_rho_are_invariant_under_rigid_motion_and_permutation():
    row = _read("t13_holdout_manifest.csv")[5]
    positions = cluster_family(6, row["family"], float(row["distance_ratio"]))
    angle = 0.731
    rotation = np.asarray([
        [np.cos(angle), -np.sin(angle), 0.0],
        [np.sin(angle), np.cos(angle), 0.0],
        [0.0, 0.0, 1.0],
    ])
    transformed = (positions @ rotation.T) + np.asarray([1.7, -2.3, 0.0])
    transformed = transformed[np.asarray([4, 0, 5, 2, 1, 3])]
    original_lambda = maximum_geometric_coupling(positions, 1.0, float(row["f1"]))
    transformed_lambda = maximum_geometric_coupling(transformed, 1.0, float(row["f1"]))
    original_rho = spectral_radius_l1(solve_multipolar_nodal_interaction_forces(
        positions, 0.1, 1.0, 1.0, 0.0, float(row["f1"]), 1
    ))
    transformed_rho = spectral_radius_l1(solve_multipolar_nodal_interaction_forces(
        transformed, 0.1, 1.0, 1.0, 0.0, float(row["f1"]), 1
    ))
    assert transformed_lambda == pytest.approx(original_lambda, rel=3e-14, abs=3e-15)
    assert transformed_rho == pytest.approx(original_rho, rel=3e-12, abs=3e-13)


def test_frozen_models_a_and_d_are_audited_against_public_apis():
    manifest = _read("t13_holdout_manifest.csv")[4]
    case_id = manifest["case_id"]
    positions = cluster_family(6, manifest["family"], float(manifest["distance_ratio"]))
    f1 = float(manifest["f1"])
    model_a = np.zeros((6, 2))
    for first, second in combinations(range(6), 2):
        model_a[first] += nodal_pair_force_on_probe(
            positions[first, :2], positions[second, :2], 0.1, 1.0, 1.0, f1
        )
        model_a[second] += nodal_pair_force_on_probe(
            positions[second, :2], positions[first, :2], 0.1, 1.0, 1.0, f1
        )
    model_d = solve_multipolar_nodal_interaction_forces(
        positions, 0.1, 1.0, 1.0, 0.0, f1, int(manifest["reference_lmax_d"])
    ).forces_xy
    frozen = [row for row in _read("t13_forces.csv") if row["case_id"] == case_id]
    assert np.allclose(model_a, _vectors(frozen, "model_a")[:, :2], rtol=5e-12, atol=5e-14)
    assert np.allclose(model_d, _vectors(frozen, "model_d")[:, :2], rtol=5e-12, atol=5e-14)


def test_interaction_force_is_the_exclusive_external_reference():
    summary = _read("t13_case_summary.csv")[0]
    rows = [row for row in _read("t13_forces.csv") if row["case_id"] == summary["case_id"]]
    model_a = _vectors(rows, "model_a")
    interaction = _vectors(rows, "model_e_interaction")
    interaction_error, applicable = normalized_rms_error_xyz(interaction, model_a)
    assert applicable
    assert interaction_error == pytest.approx(float(summary["epsilon_a_e"]), rel=2e-15)
    assert {row["response_source"] for row in _read("t13_external_predictions.csv")} == {"Model_E_interaction_force_only"}


def test_convergence_sequence_obeys_two_step_and_extension_protocol():
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in _read("t13_model_e_convergence.csv"):
        grouped.setdefault(row["case_id"], []).append(row)
    for rows in grouped.values():
        rows.sort(key=lambda row: int(row["lmax"]))
        orders = [int(row["lmax"]) for row in rows]
        assert orders == list(range(2, orders[-1] + 1))
        assert orders[-1] <= 21
        final = rows[-1]
        assert final["interaction_confirmed"] == "true"
        if orders[-1] < 13:
            assert orders[-1] >= 5
            assert all(final[f"{channel}_confirmed"] == "true" for channel in (
                "total", "interaction", "external_scattered", "scattered_scattered"
            ))
        if orders[-1] > 13:
            assert final["extended_beyond_13"] == "true"


def test_nonconfirmed_or_failed_cases_are_explicitly_excluded():
    mask = external_eligibility_mask(
        [True, False, True, True],
        [True, True, False, True],
        [True, True, True, False],
    )
    assert mask.tolist() == [True, False, False, False]
    with pytest.raises(ValueError):
        external_eligibility_mask([True], [True, False], [True])


def test_metrics_and_gate_use_all_and_only_eligible_cases():
    summaries = _read("t13_case_summary.csv")
    eligible = sum(row["eligible"] == "true" for row in summaries)
    metrics = _read("t13_metrics.csv")
    global_rows = [row for row in metrics if row["scope_type"] == "global"]
    assert len(global_rows) == 2
    assert {int(row["point_count"]) for row in global_rows} == {eligible}
    gate = _read("t13_gate.csv")
    decision = [row for row in gate if row["stage"] == "decision"]
    assert len(decision) == 1
    assert decision[0]["decision"] == "PASS_T13_EXTERNAL_VALIDATION_LAMBDA_MAX"
    assert all(row["passed"] == "true" for row in gate)


def test_all_t13_tables_are_finite_and_deterministically_ordered():
    names = (
        "t13_model_e_convergence.csv", "t13_forces.csv", "t13_case_summary.csv",
        "t13_external_predictions.csv", "t13_metrics.csv",
        "t13_threshold_audit.csv", "t13_gate.csv",
    )
    for name in names:
        text = (DATA / name).read_text(encoding="utf-8")
        assert "nan" not in text.lower()
        assert "inf" not in text.lower()
        assert text.endswith("\n")
    summary = _read("t13_case_summary.csv")
    assert [int(row["holdout_order"]) for row in summary] == list(range(1, 25))
    forces = _read("t13_forces.csv")
    keys = [(int(row["holdout_order"]), int(row["particle_index"])) for row in forces]
    assert keys == sorted(keys)
