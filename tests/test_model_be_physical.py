"""Response-blind P1.3 physical validation of the Model-B_E composition."""

from __future__ import annotations

from functools import lru_cache
import json

import numpy as np

from acoustic_ms import (
    corrected_nodal_pair_forces,
    evaluate_model_e_numerical_diagnostics,
    nodal_pair_forces,
    rms_vector_magnitude_xyz,
    solve_model_be_nodal,
    solve_model_e_nodal,
)


RADIUS = 1.0
ENERGY_DENSITY = 1.0
CONVERGENCE_TOLERANCE = 1.0e-5
NUMERICAL_RELATIVE_BUDGET = 5.0e-10
ROUND_OFF_FACTOR = 512.0
FLOAT_EPSILON = np.finfo(float).eps


def _dimer_positions(distance: float, theta: float) -> np.ndarray:
    axis = np.array([np.cos(theta), np.sin(theta), 0.0])
    return np.vstack((-0.5 * distance * axis, 0.5 * distance * axis))


def _rms_error(actual: np.ndarray, expected: np.ndarray) -> tuple[float, float, float]:
    absolute = rms_vector_magnitude_xyz(actual - expected)
    scale = max(
        rms_vector_magnitude_xyz(actual),
        rms_vector_magnitude_xyz(expected),
    )
    relative = absolute / scale if scale > 0.0 else absolute
    return absolute, relative, scale


def _numerical_tolerance(scale: float) -> float:
    return (NUMERICAL_RELATIVE_BUDGET + ROUND_OFF_FACTOR * FLOAT_EPSILON) * scale


def _solve_be(
    positions: np.ndarray,
    *,
    ka: float,
    f0: float,
    f1: float,
):
    return solve_model_be_nodal(
        positions,
        ka / RADIUS,
        RADIUS,
        ENERGY_DENSITY,
        f0,
        f1,
        convergence_tolerance=CONVERGENCE_TOLERANCE,
    )


def _required_forces(result) -> np.ndarray:
    assert result.eligible, result.failure_reason
    assert result.forces_xyz is not None
    return np.asarray(result.forces_xyz)


def _pair_summary(record) -> dict[str, object]:
    assert record.diagnostics is not None
    diagnostics = record.diagnostics
    return {
        "pair": list(record.particle_indices),
        "final_lmax": record.final_lmax,
        "eligible": record.eligible,
        "condition": diagnostics.balanced_condition_number,
        "backward_error": diagnostics.balanced_backward_error,
        "incident_closure": diagnostics.effective_incident_closure_error,
        "scattering_closure": diagnostics.scattering_closure_error,
        "decomposition_residual": diagnostics.force_decomposition_residual,
        "max_abs_fz": diagnostics.max_abs_fz,
        "fz_tolerance": diagnostics.fz_tolerance,
        "diagnostics_passed": diagnostics.passed,
        "channels": {
            channel.channel: {
                "applicable": channel.applicable,
                "confirmed": channel.confirmed,
                "confirmation_lmax": channel.confirmation_lmax,
                "final_successive_change": channel.final_successive_change,
                "final_absolute_change": channel.final_absolute_change,
                "final_window": [
                    {
                        "lmax": step.lmax,
                        "successive_change": step.successive_change,
                        "applicable": step.applicable,
                    }
                    for step in channel.history[-2:]
                ],
            }
            for channel in record.convergence
        },
    }


@lru_cache(maxsize=1)
def _physical_evidence() -> dict[str, object]:
    geo_positions = _dimer_positions(2.7, 0.41)
    geo = _solve_be(geo_positions, ka=0.04, f0=0.0, f1=0.35)
    geo_forces = _required_forces(geo)

    rotation_angle = 0.73
    rotation = np.array(
        [
            [np.cos(rotation_angle), -np.sin(rotation_angle), 0.0],
            [np.sin(rotation_angle), np.cos(rotation_angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    reflection = np.diag([-1.0, 1.0, 1.0])
    geo_rotated = _solve_be(
        geo_positions @ rotation.T,
        ka=0.04,
        f0=0.0,
        f1=0.35,
    )
    geo_reflected = _solve_be(
        geo_positions @ reflection.T,
        ka=0.04,
        f0=0.0,
        f1=0.35,
    )

    low_positions = _dimer_positions(2.7, 0.41)
    analytic_low = _solve_be(low_positions, ka=0.02, f0=0.0, f1=0.35)

    null_positions = _dimer_positions(3.2, 0.29)
    null = _solve_be(null_positions, ka=0.07, f0=0.0, f1=0.0)

    weak_positions = _dimer_positions(3.3, 0.27)
    weak = tuple(
        _solve_be(weak_positions, ka=0.06, f0=0.0, f1=f1)
        for f1 in (0.32, 0.08, 0.02)
    )

    n3_positions = np.array(
        [
            [-1.55, -0.45, 0.0],
            [1.25, -0.55, 0.0],
            [0.20, 2.05, 0.0],
        ]
    )
    permutation = np.array([2, 0, 1])
    n3 = _solve_be(n3_positions, ka=0.075, f0=0.12, f1=0.37)
    n3_permuted = _solve_be(
        n3_positions[permutation],
        ka=0.075,
        f0=0.12,
        f1=0.37,
    )

    geo_record = geo.pair_ledger[0]
    assert geo_record.final_lmax is not None
    direct_geo = solve_model_e_nodal(
        geo_positions,
        0.04,
        RADIUS,
        ENERGY_DENSITY,
        0.0,
        0.35,
        geo_record.final_lmax,
    )
    direct_geo_diagnostics = evaluate_model_e_numerical_diagnostics(direct_geo)

    analytic: dict[str, dict[str, float]] = {}
    for label, ka, result, positions in (
        ("DEV-GEO", 0.04, geo, geo_positions),
        ("DEV-EQ30-LOW", 0.02, analytic_low, low_positions),
    ):
        first, second = corrected_nodal_pair_forces(
            positions[0, :2],
            positions[1, :2],
            ka / RADIUS,
            RADIUS,
            ENERGY_DENSITY,
            0.35,
        )
        oracle = np.column_stack((np.vstack((first, second)), np.zeros(2)))
        absolute, relative, scale = _rms_error(_required_forces(result), oracle)
        analytic[label] = {
            "ka": ka,
            "absolute_error": absolute,
            "relative_error": relative,
            "scale": scale,
            "ka_squared_envelope": ka**2,
        }

    weak_metrics: list[dict[str, float]] = []
    for f1, result in zip((0.32, 0.08, 0.02), weak):
        first, second = nodal_pair_forces(
            weak_positions[0, :2],
            weak_positions[1, :2],
            0.06,
            RADIUS,
            ENERGY_DENSITY,
            f1,
        )
        model_a = np.column_stack((np.vstack((first, second)), np.zeros(2)))
        absolute, relative, scale = _rms_error(_required_forces(result), model_a)
        weak_metrics.append(
            {
                "f1": f1,
                "absolute_correction": absolute,
                "relative_correction": relative,
                "scale": scale,
            }
        )

    n3_forces = _required_forces(n3)
    common_lmax = max(
        record.final_lmax
        for record in n3.pair_ledger
        if record.final_lmax is not None
    )
    common_forces = np.zeros_like(n3_forces)
    common_diagnostics = []
    for record in n3.pair_ledger:
        common_result = solve_model_e_nodal(
            record.positions_xyz,
            0.075,
            RADIUS,
            ENERGY_DENSITY,
            0.12,
            0.37,
            common_lmax,
        )
        common_diagnostics.append(
            evaluate_model_e_numerical_diagnostics(common_result)
        )
        first, second = record.particle_indices
        common_forces[first] += common_result.interaction_forces_xyz[0]
        common_forces[second] += common_result.interaction_forces_xyz[1]
    common_absolute, common_relative, common_scale = _rms_error(
        n3_forces,
        common_forces,
    )

    identity = _rms_error(geo_forces, direct_geo.interaction_forces_xyz)
    rotation_error = _rms_error(
        _required_forces(geo_rotated),
        geo_forces @ rotation.T,
    )
    reflection_error = _rms_error(
        _required_forces(geo_reflected),
        geo_forces @ reflection.T,
    )
    permutation_error = _rms_error(
        _required_forces(n3_permuted),
        n3_forces[permutation],
    )

    dimer_axis = (
        geo_positions[1] - geo_positions[0]
    ) / np.linalg.norm(geo_positions[1] - geo_positions[0])
    radial_projection = np.outer(geo_forces @ dimer_axis, dimer_axis)
    dimer_scale = rms_vector_magnitude_xyz(geo_forces)
    dimer_physics = {
        "scale": dimer_scale,
        "planarity_residual": rms_vector_magnitude_xyz(
            np.column_stack((np.zeros((2, 2)), geo_forces[:, 2]))
        ),
        "radiality_residual": rms_vector_magnitude_xyz(
            geo_forces - radial_projection
        ),
        "action_reaction_residual": float(np.linalg.norm(np.sum(geo_forces, axis=0))),
    }

    cases = {
        "DEV-GEO": geo,
        "DEV-GEO-R": geo_rotated,
        "DEV-GEO-M": geo_reflected,
        "DEV-EQ30-LOW": analytic_low,
        "DEV-NULL": null,
        "DEV-WEAK-1": weak[0],
        "DEV-WEAK-2": weak[1],
        "DEV-WEAK-3": weak[2],
        "DEV-N3": n3,
        "DEV-N3-P": n3_permuted,
    }
    return {
        "cases": cases,
        "identity": identity,
        "direct_geo_diagnostics": direct_geo_diagnostics,
        "rotation": rotation_error,
        "reflection": reflection_error,
        "permutation": permutation_error,
        "dimer_physics": dimer_physics,
        "analytic": analytic,
        "weak": weak_metrics,
        "common_lmax": common_lmax,
        "common_forces": common_forces,
        "common_diagnostics": tuple(common_diagnostics),
        "common_error": (common_absolute, common_relative, common_scale),
    }


def test_real_pair_eligibility_diagnostics_and_final_convergence_window() -> None:
    evidence = _physical_evidence()
    for result in evidence["cases"].values():
        assert result.eligible
        assert result.failure_stage is None
        for record in result.pair_ledger:
            assert record.eligible
            assert record.converged
            assert record.final_lmax is not None
            assert 5 <= record.final_lmax <= 21
            assert record.diagnostics is not None
            assert record.diagnostics.passed

    geo_record = evidence["cases"]["DEV-GEO"].pair_ledger[0]
    applicable_count = 0
    for channel in geo_record.convergence:
        if not channel.applicable:
            continue
        applicable_count += 1
        assert channel.confirmed
        final_window = channel.history[-2:]
        assert len(final_window) == 2
        assert all(step.applicable for step in final_window)
        assert all(
            step.successive_change <= CONVERGENCE_TOLERANCE
            for step in final_window
        )
    assert applicable_count > 0


def test_n2_identity_with_real_model_e_at_converged_order() -> None:
    evidence = _physical_evidence()
    absolute, _, scale = evidence["identity"]
    assert evidence["direct_geo_diagnostics"].passed
    assert absolute <= _numerical_tolerance(scale)


def test_planar_rotation_reflection_and_particle_permutation_covariance() -> None:
    evidence = _physical_evidence()
    for key in ("rotation", "reflection", "permutation"):
        absolute, _, scale = evidence[key]
        assert absolute <= _numerical_tolerance(scale)


def test_applicable_dimer_planarity_radiality_and_action_reaction() -> None:
    physics = _physical_evidence()["dimer_physics"]
    tolerance = _numerical_tolerance(physics["scale"])
    assert physics["planarity_residual"] <= tolerance
    assert physics["radiality_residual"] <= tolerance
    assert physics["action_reaction_residual"] <= tolerance


def test_null_channels_are_inapplicable_and_exactly_zero() -> None:
    result = _physical_evidence()["cases"]["DEV-NULL"]
    forces = _required_forces(result)
    assert not np.any(forces)
    record = result.pair_ledger[0]
    assert record.final_lmax == 5
    assert all(not channel.applicable for channel in record.convergence)
    assert all(not channel.confirmed for channel in record.convergence)


def test_corrected_fifth_order_formula_in_frozen_asymptotic_domain() -> None:
    analytic = _physical_evidence()["analytic"]
    high = analytic["DEV-GEO"]
    low = analytic["DEV-EQ30-LOW"]
    assert high["relative_error"] <= high["ka_squared_envelope"]
    assert low["relative_error"] <= low["ka_squared_envelope"]
    assert low["relative_error"] < high["relative_error"]
    ratio = low["relative_error"] / high["relative_error"]
    assert 0.15 <= ratio <= 0.35


def test_model_be_minus_model_a_tends_to_zero_with_weak_dipole_contrast() -> None:
    metrics = _physical_evidence()["weak"]
    absolute = [item["absolute_correction"] for item in metrics]
    relative = [item["relative_correction"] for item in metrics]
    assert absolute[0] > absolute[1] > absolute[2]
    assert relative[0] > relative[1] > relative[2]
    for previous, current, item in zip(absolute, absolute[1:], metrics):
        assert previous - current > _numerical_tolerance(item["scale"])


def test_n3_common_order_audit_with_real_model_e() -> None:
    evidence = _physical_evidence()
    assert all(item.passed for item in evidence["common_diagnostics"])
    absolute, _, scale = evidence["common_error"]
    truncation_budget = (
        10.0 * CONVERGENCE_TOLERANCE + ROUND_OFF_FACTOR * FLOAT_EPSILON
    ) * scale
    assert absolute <= truncation_budget


def test_p1_3_evidence_snapshot() -> None:
    evidence = _physical_evidence()
    summary = {
        "cases": {
            label: [_pair_summary(record) for record in result.pair_ledger]
            for label, result in evidence["cases"].items()
        },
        "identity": {
            "absolute_error": evidence["identity"][0],
            "relative_error": evidence["identity"][1],
        },
        "rotation": {
            "absolute_error": evidence["rotation"][0],
            "relative_error": evidence["rotation"][1],
        },
        "reflection": {
            "absolute_error": evidence["reflection"][0],
            "relative_error": evidence["reflection"][1],
        },
        "permutation": {
            "absolute_error": evidence["permutation"][0],
            "relative_error": evidence["permutation"][1],
        },
        "dimer_physics": evidence["dimer_physics"],
        "analytic": evidence["analytic"],
        "weak": evidence["weak"],
        "common_order": {
            "lmax": evidence["common_lmax"],
            "absolute_error": evidence["common_error"][0],
            "relative_error": evidence["common_error"][1],
            "scale": evidence["common_error"][2],
            "truncation_budget": (
                10.0 * CONVERGENCE_TOLERANCE
                + ROUND_OFF_FACTOR * FLOAT_EPSILON
            ) * evidence["common_error"][2],
            "diagnostics": [
                {
                    "condition": item.balanced_condition_number,
                    "backward_error": item.balanced_backward_error,
                    "incident_closure": item.effective_incident_closure_error,
                    "scattering_closure": item.scattering_closure_error,
                    "decomposition_residual": item.force_decomposition_residual,
                    "max_abs_fz": item.max_abs_fz,
                    "passed": item.passed,
                }
                for item in evidence["common_diagnostics"]
            ],
        },
    }
    encoded = json.dumps(summary, sort_keys=True)
    print("P1_3_EVIDENCE=" + encoded)
    assert '\"DEV-GEO\"' in encoded
