"""Unit and scientific checks for the preregistered T12 sentinel analysis."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import sys

import numpy as np
import pytest

from acoustic_ms import (
    compare_model_e_forces,
    normalized_rms_error_xyz,
    rms_vector_magnitude_xyz,
    solve_model_e_nodal,
    symmetric_rms_error_xyz,
)


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import analyze_t12_model_e_sentinels as t12  # noqa: E402


EXPECTED_IDS = tuple(item.case_id for item in t12.SENTINELS)


def test_preregistered_manifest_has_exact_order_and_strata() -> None:
    assert len(EXPECTED_IDS) == len(set(EXPECTED_IDS)) == 28
    assert all(item.particle_count <= 4 for item in t12.SENTINELS)
    strata = {(item.particle_count, item.family) for item in t12.SENTINELS}
    assert strata == {
        (2, "pair"),
        (3, "compact"), (3, "irregular"), (3, "linear"),
        (4, "compact"), (4, "irregular"), (4, "linear"),
    }
    for stratum in strata:
        selected = [
            item for item in t12.SENTINELS
            if (item.particle_count, item.family) == stratum
        ]
        assert [item.rho_band for item in selected] == [1, 2, 3, 4]


def test_frozen_fit_and_thresholds_are_literal_t08_values() -> None:
    assert t12.FROZEN_PREFACTOR == 2.6353684041458636
    assert t12.FROZEN_EXPONENT == 1.1088518115798773
    assert t12.FROZEN_THRESHOLDS == (
        (0.01, 0.0053990295322641655),
        (0.05, 0.02000077753569526),
        (0.10, 0.03914887870730305),
    )
    t12._frozen_fit_and_thresholds()


def test_frozen_positions_and_models_a_d_are_reproduced() -> None:
    sentinel = t12.SENTINELS[0]
    cases, forces = t12._load_t08()
    rows = forces[sentinel.case_id]
    positions = np.asarray(
        [[float(row[name]) for name in ("x", "y", "z")] for row in rows]
    )
    np.testing.assert_allclose(
        positions,
        t12.cluster_family(
            sentinel.particle_count, sentinel.family, sentinel.distance_ratio
        ),
        rtol=0.0,
        atol=5e-15,
    )
    frozen_a = np.asarray([[float(row["a_x"]), float(row["a_y"])] for row in rows])
    frozen_d = np.asarray([[float(row["d_x"]), float(row["d_y"])] for row in rows])
    np.testing.assert_allclose(
        t12._model_a(positions, sentinel.f1), frozen_a, rtol=5e-12, atol=5e-14
    )
    recalculated_d = t12.solve_multipolar_nodal_interaction_forces(
        positions, t12.KA, t12.RADIUS, t12.ENERGY_DENSITY, t12.F0,
        sentinel.f1, int(cases[sentinel.case_id]["reference_lmax"]),
    ).forces_xy
    np.testing.assert_allclose(
        recalculated_d, frozen_d, rtol=5e-12, atol=5e-14
    )


@pytest.mark.parametrize(
    "function",
    (rms_vector_magnitude_xyz, symmetric_rms_error_xyz, normalized_rms_error_xyz),
)
def test_xyz_metrics_reject_invalid_vectors(function) -> None:
    with pytest.raises(ValueError):
        function(np.ones((2, 2))) if function is rms_vector_magnitude_xyz else function(np.ones((2, 2)), np.ones((2, 2)))


def test_xyz_rms_and_errors_use_the_z_component() -> None:
    reference = np.array([[3.0, 4.0, 12.0], [0.0, 0.0, 0.0]])
    model = np.array([[3.0, 4.0, 0.0], [0.0, 0.0, 0.0]])
    assert rms_vector_magnitude_xyz(reference) == pytest.approx(np.sqrt(169.0 / 2.0))
    error, applicable = normalized_rms_error_xyz(reference, model)
    assert applicable
    assert error == pytest.approx(np.sqrt(144.0 / 169.0))
    expected_symmetric = 2.0 * np.sqrt(144.0 / 2.0) / (
        np.sqrt(169.0 / 2.0) + np.sqrt(25.0 / 2.0)
    )
    assert symmetric_rms_error_xyz(reference, model) == pytest.approx(expected_symmetric)


def test_unresolved_xyz_reference_is_flagged_without_nan() -> None:
    zeros = np.zeros((2, 3))
    error, applicable = normalized_rms_error_xyz(zeros, zeros)
    assert error == 0.0
    assert not applicable
    assert symmetric_rms_error_xyz(zeros, zeros) == 0.0


def test_signed_a_d_e_decomposition_closes_particle_by_particle() -> None:
    a = np.array([[1.0, -2.0], [-1.0, 2.0]])
    d = np.array([[1.2, -1.9], [-1.2, 1.9]])
    external_scattered = np.array([[1.3, -1.8, 0.1], [-1.3, 1.8, -0.1]])
    scattered_scattered = np.array([[0.05, -0.02, 0.03], [-0.05, 0.02, -0.03]])
    interaction = external_scattered + scattered_scattered
    result = compare_model_e_forces(a, d, interaction, external_scattered, scattered_scattered)
    assert result.decomposition_max_abs_error < 1e-15
    assert result.decomposition_relative_error < 1e-15
    assert result.max_abs_interaction_fz == pytest.approx(0.13)
    correction = rms_vector_magnitude_xyz(
        interaction - np.column_stack((a, np.zeros(len(a))))
    )
    components = (
        result.rms_d_minus_a
        + result.rms_mie_external_correction
        + result.rms_scattered_scattered
    )
    assert result.cancellation_ratio == pytest.approx(components / correction)


@dataclass(frozen=True)
class _DummyResult:
    lmax: int
    total_forces_xyz: np.ndarray
    interaction_forces_xyz: np.ndarray
    external_scattered_forces_xyz: np.ndarray
    scattered_scattered_forces_xyz: np.ndarray


def test_confirmation_requires_two_successive_applicable_changes() -> None:
    results = []
    for order, value in ((2, 1.0), (3, 1.1), (4, 1.100001), (5, 1.1000015)):
        vectors = np.array([[value, 0.0, 0.0]])
        results.append(_DummyResult(order, vectors, vectors, vectors, vectors))
    assert t12._minimum_confirmation(results[:3], "interaction_forces_xyz") == 0
    assert t12._minimum_confirmation(results, "interaction_forces_xyz") == 5


def test_unconfirmed_status_is_not_described_as_divergent() -> None:
    source = (SCRIPTS / "analyze_t12_model_e_sentinels.py").read_text(encoding="utf-8")
    assert '"unconfirmed"' in source
    assert "divergent" not in source


def test_t12_path_does_not_use_forbidden_linear_algebra_or_regression() -> None:
    sources = [
        (SCRIPTS / "analyze_t12_model_e_sentinels.py").read_text(encoding="utf-8"),
        (Path(__file__).resolve().parents[1] / "src/acoustic_ms/model_e_comparison.py").read_text(encoding="utf-8"),
    ]
    combined = "\n".join(sources)
    for forbidden in ("np.linalg.inv", "np.linalg.pinv", "np.linalg.lstsq", "polyfit(", "curve_fit("):
        assert forbidden not in combined


def test_model_e_production_solver_and_force_closure() -> None:
    positions = np.array([[-1.5, 0.0, 0.0], [1.5, 0.0, 0.0]])
    result = solve_model_e_nodal(positions, 0.1, 1.0, 1.0, 0.0, 0.4, 3)
    assert result.solution.production_solver == "balanced_sqrt"
    np.testing.assert_allclose(
        result.interaction_forces_xyz,
        result.external_scattered_forces_xyz + result.scattered_scattered_forces_xyz,
        rtol=0.0,
        atol=2e-15,
    )


def test_model_e_rotation_and_permutation_covariance() -> None:
    positions = np.array([[-1.5, 0.0, 0.0], [1.5, 0.0, 0.0]])
    angle = 0.37
    rotation = np.array([
        [np.cos(angle), -np.sin(angle), 0.0],
        [np.sin(angle), np.cos(angle), 0.0],
        [0.0, 0.0, 1.0],
    ])
    baseline = solve_model_e_nodal(positions, 0.1, 1.0, 1.0, 0.0, 0.4, 3)
    rotated = solve_model_e_nodal(positions @ rotation.T, 0.1, 1.0, 1.0, 0.0, 0.4, 3)
    np.testing.assert_allclose(
        rotated.interaction_forces_xyz,
        baseline.interaction_forces_xyz @ rotation.T,
        rtol=3e-12,
        atol=3e-14,
    )
    permuted = solve_model_e_nodal(positions[::-1], 0.1, 1.0, 1.0, 0.0, 0.4, 3)
    np.testing.assert_allclose(
        permuted.interaction_forces_xyz,
        baseline.interaction_forces_xyz[::-1],
        rtol=3e-12,
        atol=3e-14,
    )
