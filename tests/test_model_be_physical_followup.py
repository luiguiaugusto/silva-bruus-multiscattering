"""P1.3a development follow-up for the failed fifth-order comparison."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import json

import numpy as np
import pytest

from acoustic_ms import (
    complete_radiation_force_from_bsc,
    corrected_nodal_pair_force_magnitude,
    evaluate_model_e_numerical_diagnostics,
    mie_scattering_coefficients_from_contrasts,
    mode_count,
    mode_index,
    modes,
    nodal_standing_wave_coefficients,
    rayleigh_multipolar_scattering_coefficients,
    rms_vector_magnitude_xyz,
    separation_coefficient,
    solve_model_be_nodal,
    solve_model_e_nodal,
    solve_multipolar_nodal_interaction_forces,
)
from acoustic_ms.special import spherical_harmonic


RADIUS = 1.0
ENERGY_DENSITY = 1.0
DISTANCE = 2.7
F0 = 0.0
F1 = 0.35
LMAX = 5
KA_VALUES = np.array([0.08, 0.04, 0.02, 0.01])
CONVERGENCE_TOLERANCE = 1.0e-5
NUMERICAL_RELATIVE_BUDGET = 5.0e-10
ROUND_OFF_FACTOR = 512.0
FLOAT_EPSILON = np.finfo(float).eps


@dataclass(frozen=True)
class ProjectedDimerAudit:
    """One test-only modal projection of the dimer scattering system."""

    space: str
    coefficient_model: str
    selected_modes: tuple[tuple[int, int], ...]
    effective_incident_coefficients: np.ndarray
    scattered_coefficients: np.ndarray
    external_incident_coefficients: np.ndarray
    scattering_by_ell: np.ndarray
    total_forces_xyz: np.ndarray
    external_forces_xyz: np.ndarray
    interaction_forces_xyz: np.ndarray
    external_scattered_forces_xyz: np.ndarray
    scattered_scattered_forces_xyz: np.ndarray
    appendix_b_forces_xy: np.ndarray
    balanced_condition_number: float
    balanced_backward_error: float
    effective_incident_closure_error: float
    scattering_closure_error: float
    force_decomposition_residual: float


@dataclass(frozen=True)
class StrictArticleAudit:
    """Independent one-particle reduction used in the article's Eq. (30)."""

    force: float
    coefficients: np.ndarray
    condition_number: float
    residual_relative: float


def _positions(distance: float = DISTANCE) -> np.ndarray:
    return np.array([[-distance / 2.0, 0.0, 0.0], [distance / 2.0, 0.0, 0.0]])


def _numerical_tolerance(scale: float) -> float:
    return (NUMERICAL_RELATIVE_BUDGET + ROUND_OFF_FACTOR * FLOAT_EPSILON) * scale


def _harmonic_identity_metrics() -> dict[str, float]:
    maximum_error = 0.0
    maximum_tolerance = 0.0
    for theta, phi in ((0.37, 0.23), (0.91, 1.17), (1.42, 2.03)):
        for ell in range(LMAX + 1):
            for m in range(-ell, ell + 1):
                left = spherical_harmonic(ell, m, np.pi - theta, phi)
                right = (-1) ** (ell + m) * spherical_harmonic(
                    ell, m, theta, phi
                )
                tolerance = (
                    ROUND_OFF_FACTOR
                    * FLOAT_EPSILON
                    * max(1.0, abs(left), abs(right))
                )
                maximum_error = max(maximum_error, float(abs(left - right)))
                maximum_tolerance = max(maximum_tolerance, float(tolerance))
    return {"maximum_error": maximum_error, "maximum_tolerance": maximum_tolerance}


def _backward_error(matrix: np.ndarray, solution: np.ndarray, rhs: np.ndarray) -> float:
    numerator = float(np.linalg.norm(matrix @ solution - rhs))
    denominator = (
        float(np.linalg.norm(matrix)) * float(np.linalg.norm(solution))
        + float(np.linalg.norm(rhs))
    )
    return 0.0 if denominator == 0.0 else numerator / denominator


def _closure_error(residual: np.ndarray, *terms: np.ndarray) -> float:
    denominator = sum(float(np.linalg.norm(term)) for term in terms)
    return 0.0 if denominator == 0.0 else float(np.linalg.norm(residual)) / denominator


def _selected_modes(space: str) -> tuple[tuple[int, int], ...]:
    if space == "planar_complete":
        return tuple((ell, m) for ell, m in modes(LMAX) if (ell + m) % 2 == 1)
    if space == "article_reduced":
        return tuple(
            (ell, m)
            for ell, m in modes(LMAX)
            if ell % 2 == 1 and m % 2 == 0
        )
    raise ValueError(f"unknown audit space: {space}")


def _projected_dimer_audit(
    ka: float,
    *,
    space: str,
    coefficient_model: str,
) -> ProjectedDimerAudit:
    """Solve a selected modal subspace using only established primitives."""

    selected = _selected_modes(space)
    full_modes = modes(LMAX)
    selected_indices = tuple(full_modes.index(mode) for mode in selected)
    if coefficient_model == "mie_exact":
        scattering = mie_scattering_coefficients_from_contrasts(
            ka, F0, F1, LMAX
        )
    elif coefficient_model == "appendix_a":
        scattering = rayleigh_multipolar_scattering_coefficients(
            ka, F0, F1, LMAX
        )
    else:
        raise ValueError(f"unknown coefficient model: {coefficient_model}")

    particle_positions = _positions()
    active_count = len(selected)
    particle_count = len(particle_positions)
    local_scattering = np.asarray(
        [scattering[ell] for ell, _ in selected], dtype=complex
    )
    global_scattering = np.tile(local_scattering, particle_count)
    translation = np.zeros(
        (particle_count * active_count, particle_count * active_count),
        dtype=complex,
    )
    for target in range(particle_count):
        for source in range(particle_count):
            if target == source:
                continue
            row_offset = target * active_count
            column_offset = source * active_count
            for row, (target_ell, target_m) in enumerate(selected):
                for column, (source_ell, source_m) in enumerate(selected):
                    translation[row_offset + row, column_offset + column] = (
                        separation_coefficient(
                            target_ell,
                            target_m,
                            source_ell,
                            source_m,
                            ka,
                            particle_positions[target],
                            particle_positions[source],
                        )
                    )

    external_local = nodal_standing_wave_coefficients(LMAX)
    external_active = external_local[np.asarray(selected_indices)]
    external_global = np.tile(external_active, particle_count)
    square_root = np.sqrt(global_scattering)
    identity = np.eye(particle_count * active_count, dtype=complex)
    balanced_system = (
        identity
        - square_root[:, None] * translation * square_root[None, :]
    )
    balanced_rhs = square_root * external_global
    balanced_solution = np.linalg.solve(balanced_system, balanced_rhs)
    scattered_active = square_root * balanced_solution
    effective_active = external_global + translation @ scattered_active

    full_count = mode_count(LMAX)
    external = np.tile(external_local, (particle_count, 1))
    local_full = np.zeros((particle_count, full_count), dtype=complex)
    scattered = np.zeros((particle_count, full_count), dtype=complex)
    active_array = np.asarray(selected_indices)
    scattered[:, active_array] = scattered_active.reshape(
        particle_count, active_count
    )
    for target in range(particle_count):
        for source in range(particle_count):
            if target == source:
                continue
            for target_index, (target_ell, target_m) in enumerate(full_modes):
                for source_index, (source_ell, source_m) in enumerate(selected):
                    local_full[target, target_index] += separation_coefficient(
                        target_ell,
                        target_m,
                        source_ell,
                        source_m,
                        ka,
                        particle_positions[target],
                        particle_positions[source],
                    ) * scattered[source, active_array[source_index]]
    effective = external + local_full
    force_scattering = np.array(scattering, copy=True)
    if space == "article_reduced":
        force_scattering[0::2] = 0.0

    total = np.vstack(
        [
            complete_radiation_force_from_bsc(
                row, force_scattering, ka, ENERGY_DENSITY
            )
            for row in effective
        ]
    )
    external_force = np.vstack(
        [
            complete_radiation_force_from_bsc(
                row, force_scattering, ka, ENERGY_DENSITY
            )
            for row in external
        ]
    )
    scattered_incident = effective - external
    scattered_scattered = np.vstack(
        [
            complete_radiation_force_from_bsc(
                row, force_scattering, ka, ENERGY_DENSITY
            )
            for row in scattered_incident
        ]
    )
    interaction = total - external_force
    external_scattered = interaction - scattered_scattered
    reconstruction = external_force + external_scattered + scattered_scattered
    decomposition_scale = max(float(np.linalg.norm(total)), np.finfo(float).eps)
    decomposition_residual = float(
        np.linalg.norm(total - reconstruction) / decomposition_scale
    )

    local_scattered = effective - external
    b2m1 = local_scattered[:, mode_index(2, -1)]
    b21 = local_scattered[:, mode_index(2, 1)]
    prefactor = (
        np.sqrt(30.0 * np.pi)
        * ka
        * RADIUS**3
        * ENERGY_DENSITY
        / 15.0
    )
    appendix_b = np.empty((particle_count, 2), dtype=float)
    appendix_b[:, 0] = prefactor * np.real(F1 * (b2m1 - b21))
    appendix_b[:, 1] = prefactor * np.real(-1j * F1 * (b21 + b2m1))

    return ProjectedDimerAudit(
        space=space,
        coefficient_model=coefficient_model,
        selected_modes=selected,
        effective_incident_coefficients=effective,
        scattered_coefficients=scattered,
        external_incident_coefficients=external,
        scattering_by_ell=np.asarray(scattering),
        total_forces_xyz=total,
        external_forces_xyz=external_force,
        interaction_forces_xyz=interaction,
        external_scattered_forces_xyz=external_scattered,
        scattered_scattered_forces_xyz=scattered_scattered,
        appendix_b_forces_xy=appendix_b,
        balanced_condition_number=float(np.linalg.cond(balanced_system)),
        balanced_backward_error=_backward_error(
            balanced_system, balanced_solution, balanced_rhs
        ),
        effective_incident_closure_error=_closure_error(
            effective_active - external_global - translation @ scattered_active,
            effective_active,
            external_global,
            translation @ scattered_active,
        ),
        scattering_closure_error=_closure_error(
            scattered_active - global_scattering * effective_active,
            scattered_active,
            global_scattering * effective_active,
        ),
        force_decomposition_residual=decomposition_residual,
    )


def _strict_article_audit(ka: float) -> StrictArticleAudit:
    """Evaluate the established one-particle Appendix-B reduction."""

    particle_positions = _positions()
    selected = tuple(
        (ell, m) for ell in (1, 3, 5) for m in range(0, ell + 1, 2)
    )
    scattering = rayleigh_multipolar_scattering_coefficients(
        ka, F0, F1, LMAX
    )
    incident = nodal_standing_wave_coefficients(LMAX)
    translation = np.empty((len(selected), len(selected)), dtype=complex)
    external_active = np.empty(len(selected), dtype=complex)
    for row, (ell, m) in enumerate(selected):
        external_active[row] = incident[mode_index(ell, m)]
        for column, (source_ell, source_m) in enumerate(selected):
            translated = separation_coefficient(
                ell,
                m,
                source_ell,
                source_m,
                ka,
                particle_positions[0],
                particle_positions[1],
            )
            if source_m:
                translated += separation_coefficient(
                    ell,
                    m,
                    source_ell,
                    -source_m,
                    ka,
                    particle_positions[0],
                    particle_positions[1],
                )
            translation[row, column] = translated
    local_scattering = np.asarray(
        [scattering[ell] for ell, _ in selected], dtype=complex
    )
    square_root = np.sqrt(local_scattering)
    matrix = (
        np.eye(len(selected), dtype=complex)
        - square_root[:, None] * translation * square_root[None, :]
    )
    rhs = square_root * external_active
    coefficients = square_root * np.linalg.solve(matrix, rhs)
    b21 = 0.0j
    for (ell, m), coefficient in zip(selected, coefficients):
        b21 += separation_coefficient(
            2,
            1,
            ell,
            m,
            ka,
            particle_positions[0],
            particle_positions[1],
        ) * coefficient
        if m:
            b21 += separation_coefficient(
                2,
                1,
                ell,
                -m,
                ka,
                particle_positions[0],
                particle_positions[1],
            ) * coefficient
    force = float(-2.0 * np.sqrt(30.0 * np.pi) * ka * F1 * np.real(b21) / 15.0)
    physical_matrix = (
        np.eye(len(selected), dtype=complex)
        - local_scattering[:, None] * translation
    )
    physical_rhs = local_scattering * external_active
    residual = float(
        np.linalg.norm(physical_matrix @ coefficients - physical_rhs)
        / max(float(np.linalg.norm(physical_rhs)), np.finfo(float).eps)
    )
    return StrictArticleAudit(
        force=force,
        coefficients=coefficients,
        condition_number=float(np.linalg.cond(matrix)),
        residual_relative=residual,
    )


def _signed_radial(forces: np.ndarray) -> float:
    return float(np.asarray(forces)[0, 0])


def _scalar_errors(left: float, right: float) -> tuple[float, float, float]:
    absolute = abs(left - right)
    scale = max(abs(left), abs(right))
    relative = absolute / scale if scale > 0.0 else absolute
    return absolute, relative, scale


def _series_metrics(left: np.ndarray, right: np.ndarray) -> dict[str, object]:
    absolute = np.abs(left - right)
    scale = np.maximum(np.abs(left), np.abs(right))
    relative = np.divide(
        absolute,
        scale,
        out=absolute.copy(),
        where=scale > 0.0,
    )

    def local_orders(values: np.ndarray) -> list[float | None]:
        orders: list[float | None] = []
        for index in range(len(values) - 1):
            threshold = max(
                _numerical_tolerance(scale[index]),
                _numerical_tolerance(scale[index + 1]),
            )
            if values[index] <= threshold or values[index + 1] <= threshold:
                orders.append(None)
            else:
                orders.append(
                    float(np.log(values[index] / values[index + 1]) / np.log(2.0))
                )
        return orders

    applicable = absolute > np.asarray([_numerical_tolerance(item) for item in scale])
    global_absolute_order = (
        float(np.polyfit(np.log(KA_VALUES), np.log(absolute), 1)[0])
        if np.all(applicable)
        else None
    )
    global_relative_order = (
        float(np.polyfit(np.log(KA_VALUES), np.log(relative), 1)[0])
        if np.all(applicable) and np.all(relative > 0.0)
        else None
    )
    return {
        "absolute_error": absolute,
        "relative_error": relative,
        "scale": scale,
        "absolute_order": local_orders(absolute),
        "relative_order": local_orders(relative),
        "global_absolute_order": global_absolute_order,
        "global_relative_order": global_relative_order,
    }


def _required_be_forces(result) -> np.ndarray:
    assert result.eligible, result.failure_reason
    assert result.forces_xyz is not None
    return np.asarray(result.forces_xyz)


def _common_order_case(
    positions: np.ndarray,
    ka: float,
    f0: float,
    f1: float,
) -> dict[str, object]:
    independent = solve_model_be_nodal(
        positions,
        ka,
        RADIUS,
        ENERGY_DENSITY,
        f0,
        f1,
        convergence_tolerance=CONVERGENCE_TOLERANCE,
    )
    final_orders = tuple(record.final_lmax for record in independent.pair_ledger)
    assert all(order is not None for order in final_orders)
    common_lmax = max(int(order) for order in final_orders if order is not None)
    if not independent.eligible or independent.forces_xyz is None:
        return {
            "result": independent,
            "final_orders": final_orders,
            "common_lmax": common_lmax,
            "common_forces": None,
            "diagnostics": tuple(),
            "absolute_error": None,
            "relative_error": None,
            "scale": None,
            "budget": None,
            "failure_reason": independent.failure_reason,
        }
    independent_forces = np.asarray(independent.forces_xyz)
    common_forces = np.zeros_like(independent_forces)
    diagnostics = []
    for record in independent.pair_ledger:
        direct = solve_model_e_nodal(
            record.positions_xyz,
            ka,
            RADIUS,
            ENERGY_DENSITY,
            f0,
            f1,
            common_lmax,
        )
        diagnostics.append(evaluate_model_e_numerical_diagnostics(direct))
        first, second = record.particle_indices
        common_forces[first] += direct.interaction_forces_xyz[0]
        common_forces[second] += direct.interaction_forces_xyz[1]
    absolute = rms_vector_magnitude_xyz(common_forces - independent_forces)
    scale = max(
        rms_vector_magnitude_xyz(common_forces),
        rms_vector_magnitude_xyz(independent_forces),
    )
    relative = absolute / scale if scale > 0.0 else absolute
    budget = (
        10.0 * CONVERGENCE_TOLERANCE
        + ROUND_OFF_FACTOR * FLOAT_EPSILON
    ) * scale
    return {
        "result": independent,
        "final_orders": final_orders,
        "common_lmax": common_lmax,
        "common_forces": common_forces,
        "diagnostics": tuple(diagnostics),
        "absolute_error": absolute,
        "relative_error": relative,
        "scale": scale,
        "budget": budget,
        "failure_reason": None,
    }


@lru_cache(maxsize=1)
def _common_order_evidence() -> dict[str, object]:
    return _common_order_case(
        np.array(
            [[0.0, 0.0, 0.0], [2.1, 0.0, 0.0], [0.0, 8.0, 0.0]]
        ),
        0.1,
        0.0,
        1.0,
    )


@lru_cache(maxsize=1)
def _fallback_common_order_evidence() -> dict[str, object]:
    return _common_order_case(
        np.array(
            [[0.0, 0.0, 0.0], [2.7, 0.0, 0.0], [0.0, 8.0, 0.0]]
        ),
        0.04,
        0.0,
        0.35,
    )


@lru_cache(maxsize=1)
def _followup_evidence() -> dict[str, object]:
    cases = []
    for ka in KA_VALUES:
        positions = _positions()
        model_be = solve_model_be_nodal(
            positions,
            ka,
            RADIUS,
            ENERGY_DENSITY,
            F0,
            F1,
            convergence_tolerance=CONVERGENCE_TOLERANCE,
        )
        _required_be_forces(model_be)
        record = model_be.pair_ledger[0]
        assert record.final_lmax is not None
        converged = solve_model_e_nodal(
            positions,
            ka,
            RADIUS,
            ENERGY_DENSITY,
            F0,
            F1,
            record.final_lmax,
        )
        fixed = solve_model_e_nodal(
            positions,
            ka,
            RADIUS,
            ENERGY_DENSITY,
            F0,
            F1,
            LMAX,
        )
        exact_planar = _projected_dimer_audit(
            ka, space="planar_complete", coefficient_model="mie_exact"
        )
        exact_reduced = _projected_dimer_audit(
            ka, space="article_reduced", coefficient_model="mie_exact"
        )
        rayleigh_planar = _projected_dimer_audit(
            ka, space="planar_complete", coefficient_model="appendix_a"
        )
        rayleigh_reduced = _projected_dimer_audit(
            ka, space="article_reduced", coefficient_model="appendix_a"
        )
        model_d = solve_multipolar_nodal_interaction_forces(
            positions,
            ka,
            RADIUS,
            ENERGY_DENSITY,
            F0,
            F1,
            LMAX,
        )
        strict = _strict_article_audit(ka)
        equation_30 = -corrected_nodal_pair_force_magnitude(
            ka, RADIUS, DISTANCE, ENERGY_DENSITY, F1
        )
        q_indices = np.asarray(
            [
                mode_index(ell, m)
                for ell, m in modes(LMAX)
                if ell % 2 == 0 and m % 2 != 0
            ]
        )
        exact_q_norm = float(
            np.linalg.norm(exact_planar.scattered_coefficients[:, q_indices])
        )
        exact_total_norm = float(
            np.linalg.norm(exact_planar.scattered_coefficients)
        )
        rayleigh_q_norm = float(
            np.linalg.norm(rayleigh_planar.scattered_coefficients[:, q_indices])
        )
        rayleigh_total_norm = float(
            np.linalg.norm(rayleigh_planar.scattered_coefficients)
        )
        values = {
            "e_conv_interaction": _signed_radial(converged.interaction_forces_xyz),
            "e_conv_external_scattered": _signed_radial(
                converged.external_scattered_forces_xyz
            ),
            "e_conv_scattered_scattered": _signed_radial(
                converged.scattered_scattered_forces_xyz
            ),
            "e_l5_mie_p_external_scattered": _signed_radial(
                fixed.external_scattered_forces_xyz
            ),
            "e_l5_mie_r_external_scattered": _signed_radial(
                exact_reduced.external_scattered_forces_xyz
            ),
            "e_l5_rayleigh_p_external_scattered": _signed_radial(
                rayleigh_planar.external_scattered_forces_xyz
            ),
            "e_l5_rayleigh_r_external_scattered": _signed_radial(
                rayleigh_reduced.external_scattered_forces_xyz
            ),
            "e_l5_rayleigh_p_appendix_b": _signed_radial(
                rayleigh_planar.appendix_b_forces_xy
            ),
            "e_l5_rayleigh_r_appendix_b": _signed_radial(
                rayleigh_reduced.appendix_b_forces_xy
            ),
            "model_d_planar_appendix_b": _signed_radial(model_d.forces_xy),
            "strict_article": strict.force,
            "equation_30": equation_30,
        }
        cases.append(
            {
                "ka": float(ka),
                "model_be": model_be,
                "converged": converged,
                "fixed": fixed,
                "exact_planar": exact_planar,
                "exact_reduced": exact_reduced,
                "rayleigh_planar": rayleigh_planar,
                "rayleigh_reduced": rayleigh_reduced,
                "model_d": model_d,
                "strict": strict,
                "exact_q_norm": exact_q_norm,
                "exact_q_ratio": exact_q_norm / exact_total_norm,
                "rayleigh_q_norm": rayleigh_q_norm,
                "rayleigh_q_ratio": rayleigh_q_norm / rayleigh_total_norm,
                "values": values,
            }
        )

    bridge_pairs = {
        "scattered_scattered_channel": (
            "e_conv_interaction",
            "e_conv_external_scattered",
        ),
        "converged_vs_l5": (
            "e_conv_external_scattered",
            "e_l5_mie_p_external_scattered",
        ),
        "mie_modal_p_vs_r": (
            "e_l5_mie_p_external_scattered",
            "e_l5_mie_r_external_scattered",
        ),
        "reduced_mie_vs_appendix_a": (
            "e_l5_mie_r_external_scattered",
            "e_l5_rayleigh_r_external_scattered",
        ),
        "reduced_complete_vs_appendix_b": (
            "e_l5_rayleigh_r_external_scattered",
            "e_l5_rayleigh_r_appendix_b",
        ),
        "appendix_b_vs_equation_30": (
            "e_l5_rayleigh_r_appendix_b",
            "equation_30",
        ),
        "planar_mie_vs_appendix_a": (
            "e_l5_mie_p_external_scattered",
            "e_l5_rayleigh_p_external_scattered",
        ),
        "rayleigh_modal_p_vs_r_complete": (
            "e_l5_rayleigh_p_external_scattered",
            "e_l5_rayleigh_r_external_scattered",
        ),
        "rayleigh_modal_p_vs_r_appendix_b": (
            "e_l5_rayleigh_p_appendix_b",
            "e_l5_rayleigh_r_appendix_b",
        ),
        "global_reduced_vs_strict": (
            "e_l5_rayleigh_r_appendix_b",
            "strict_article",
        ),
        "strict_vs_equation_30": ("strict_article", "equation_30"),
    }
    bridges = {}
    for name, (left_key, right_key) in bridge_pairs.items():
        left = np.asarray([case["values"][left_key] for case in cases])
        right = np.asarray([case["values"][right_key] for case in cases])
        bridges[name] = {
            "left_key": left_key,
            "right_key": right_key,
            "left": left,
            "right": right,
            **_series_metrics(left, right),
        }

    coefficient_errors = {}
    for ell in range(1, LMAX + 1):
        exact = np.asarray(
            [case["exact_planar"].scattering_by_ell[ell] for case in cases]
        )
        rayleigh = np.asarray(
            [case["rayleigh_planar"].scattering_by_ell[ell] for case in cases]
        )
        absolute = np.abs(exact - rayleigh)
        relative = np.abs(exact / rayleigh - 1.0)
        coefficient_errors[ell] = {
            "absolute": absolute,
            "relative": relative,
            "global_order": float(
                np.polyfit(np.log(KA_VALUES), np.log(relative), 1)[0]
            ),
        }

    hierarchical_closure = []
    for case in cases:
        values = case["values"]
        sequence = np.asarray(
            [
                values["e_conv_interaction"],
                values["e_conv_external_scattered"],
                values["e_l5_mie_p_external_scattered"],
                values["e_l5_mie_r_external_scattered"],
                values["e_l5_rayleigh_r_external_scattered"],
                values["e_l5_rayleigh_r_appendix_b"],
                values["equation_30"],
            ]
        )
        signed_bridges = sequence[:-1] - sequence[1:]
        closure = abs(float(np.sum(signed_bridges) - (sequence[0] - sequence[-1])))
        scale = float(np.max(np.abs(sequence)))
        hierarchical_closure.append(
            {
                "signed_bridges": signed_bridges,
                "closure_error": closure,
                "scale": scale,
            }
        )

    return {
        "cases": tuple(cases),
        "bridges": bridges,
        "coefficient_errors": coefficient_errors,
        "hierarchical_closure": tuple(hierarchical_closure),
    }


def test_spherical_harmonic_planar_reflection_identity() -> None:
    metrics = _harmonic_identity_metrics()
    assert metrics["maximum_error"] <= metrics["maximum_tolerance"]


def test_projected_full_space_reproduces_model_e_and_channel_closure() -> None:
    for case in _followup_evidence()["cases"]:
        fixed = case["fixed"]
        projected = case["exact_planar"]
        diagnostics = evaluate_model_e_numerical_diagnostics(fixed)
        assert diagnostics.passed
        for attribute in (
            "total_forces_xyz",
            "interaction_forces_xyz",
            "external_scattered_forces_xyz",
            "scattered_scattered_forces_xyz",
        ):
            actual = np.asarray(getattr(projected, attribute))
            expected = np.asarray(getattr(fixed, attribute))
            absolute = rms_vector_magnitude_xyz(actual - expected)
            scale = max(
                rms_vector_magnitude_xyz(actual),
                rms_vector_magnitude_xyz(expected),
            )
            assert absolute <= _numerical_tolerance(scale)
        closure = (
            projected.interaction_forces_xyz
            - projected.external_scattered_forces_xyz
            - projected.scattered_scattered_forces_xyz
        )
        scale = rms_vector_magnitude_xyz(projected.interaction_forces_xyz)
        assert rms_vector_magnitude_xyz(closure) <= _numerical_tolerance(scale)


def test_appendix_a_exact_mie_relative_corrections_have_order_two() -> None:
    errors = _followup_evidence()["coefficient_errors"]
    for ell in range(1, LMAX + 1):
        relative = errors[ell]["relative"]
        assert np.all(relative[:-1] > relative[1:])
        assert errors[ell]["global_order"] == pytest.approx(2.0, abs=0.05)


def test_article_reduced_global_solution_matches_independent_branch() -> None:
    evidence = _followup_evidence()
    bridge = evidence["bridges"]["global_reduced_vs_strict"]
    for absolute, scale in zip(bridge["absolute_error"], bridge["scale"]):
        assert absolute <= _numerical_tolerance(scale)
    for case in evidence["cases"]:
        model_d_force = case["values"]["model_d_planar_appendix_b"]
        projected_force = case["values"]["e_l5_rayleigh_p_appendix_b"]
        absolute, _, scale = _scalar_errors(model_d_force, projected_force)
        assert absolute <= _numerical_tolerance(scale)


def test_equivalent_article_branch_to_equation_30_has_order_two() -> None:
    bridge = _followup_evidence()["bridges"]["strict_vs_equation_30"]
    relative = bridge["relative_error"]
    assert np.all(relative[:-1] > relative[1:])
    assert all(
        order is not None and 1.7 <= order <= 2.3
        for order in bridge["relative_order"]
    )
    assert 1.7 <= bridge["global_relative_order"] <= 2.3


def test_even_ell_odd_m_sector_is_resolved_and_changes_the_force() -> None:
    evidence = _followup_evidence()
    for index in (1, 2):
        case = evidence["cases"][index]
        assert case["exact_q_ratio"] > ROUND_OFF_FACTOR * FLOAT_EPSILON
        full_force = _signed_radial(
            case["exact_planar"].interaction_forces_xyz
        )
        reduced_force = _signed_radial(
            case["exact_reduced"].interaction_forces_xyz
        )
        absolute, _, scale = _scalar_errors(full_force, reduced_force)
        assert absolute > _numerical_tolerance(scale)


def test_hierarchical_signed_bridges_close_the_original_difference() -> None:
    evidence = _followup_evidence()
    for item in evidence["hierarchical_closure"]:
        assert item["closure_error"] <= _numerical_tolerance(item["scale"])
        assert np.all(np.isfinite(item["signed_bridges"]))


def test_all_followup_solves_pass_frozen_numerical_diagnostics() -> None:
    evidence = _followup_evidence()
    for case in evidence["cases"]:
        assert case["model_be"].eligible
        assert evaluate_model_e_numerical_diagnostics(case["converged"]).passed
        assert evaluate_model_e_numerical_diagnostics(case["fixed"]).passed
        for projected in (
            case["exact_planar"],
            case["exact_reduced"],
            case["rayleigh_planar"],
            case["rayleigh_reduced"],
        ):
            assert np.isfinite(projected.balanced_condition_number)
            assert projected.balanced_condition_number < 10.0
            assert projected.balanced_backward_error <= 1.0e-12
            assert projected.effective_incident_closure_error <= 1.0e-12
            assert projected.scattering_closure_error <= 1.0e-12
            assert projected.force_decomposition_residual <= 1.0e-12
        assert case["strict"].condition_number < 10.0
        assert case["strict"].residual_relative <= 1.0e-12


def test_initial_common_order_candidate_records_nonconvergence() -> None:
    common = _common_order_evidence()
    assert not common["result"].eligible
    assert common["failure_reason"] is not None
    failed = [
        record
        for record in common["result"].pair_ledger
        if not record.eligible
    ]
    assert [record.particle_indices for record in failed] == [(0, 1)]
    assert failed[0].failure_stage == "convergence"
    assert failed[0].final_lmax == 21


def test_common_order_fallback_uses_distinct_independent_orders() -> None:
    common = _fallback_common_order_evidence()
    assert common["result"].eligible, common["failure_reason"]
    assert len(set(common["final_orders"])) >= 2
    assert all(record.eligible for record in common["result"].pair_ledger)
    assert all(item.passed for item in common["diagnostics"])
    assert common["absolute_error"] <= common["budget"]


def test_p1_3a_evidence_snapshot() -> None:
    evidence = _followup_evidence()
    summary = {
        "harmonic_identity": _harmonic_identity_metrics(),
        "cases": [
            {
                "ka": case["ka"],
                "final_lmax": case["model_be"].pair_ledger[0].final_lmax,
                "exact_q_norm": case["exact_q_norm"],
                "exact_q_ratio": case["exact_q_ratio"],
                "rayleigh_q_norm": case["rayleigh_q_norm"],
                "rayleigh_q_ratio": case["rayleigh_q_ratio"],
                "values": case["values"],
                "projected_diagnostics": {
                    name: {
                        "condition": projected.balanced_condition_number,
                        "backward_error": projected.balanced_backward_error,
                        "incident_closure": projected.effective_incident_closure_error,
                        "scattering_closure": projected.scattering_closure_error,
                        "force_decomposition": projected.force_decomposition_residual,
                    }
                    for name, projected in (
                        ("mie_p", case["exact_planar"]),
                        ("mie_r", case["exact_reduced"]),
                        ("rayleigh_p", case["rayleigh_planar"]),
                        ("rayleigh_r", case["rayleigh_reduced"]),
                    )
                },
                "strict_diagnostics": {
                    "condition": case["strict"].condition_number,
                    "residual": case["strict"].residual_relative,
                },
            }
            for case in evidence["cases"]
        ],
        "bridges": {
            name: {
                key: (
                    value.tolist()
                    if isinstance(value, np.ndarray)
                    else value
                )
                for key, value in bridge.items()
                if key not in {"left", "right"}
            }
            for name, bridge in evidence["bridges"].items()
        },
        "coefficient_errors": {
            str(ell): {
                "absolute": item["absolute"].tolist(),
                "relative": item["relative"].tolist(),
                "global_order": item["global_order"],
            }
            for ell, item in evidence["coefficient_errors"].items()
        },
        "hierarchical_closure": [
            {
                "signed_bridges": item["signed_bridges"].tolist(),
                "closure_error": item["closure_error"],
                "scale": item["scale"],
            }
            for item in evidence["hierarchical_closure"]
        ],
    }
    encoded = json.dumps(summary, sort_keys=True)
    print("P1_3A_EVIDENCE=" + encoded)
    assert '"bridges"' in encoded


def test_p1_3a_common_order_snapshot() -> None:
    def summarize(common_order: dict[str, object]) -> dict[str, object]:
        return {
            "eligible": common_order["result"].eligible,
            "failure_reason": common_order["failure_reason"],
            "final_orders": list(common_order["final_orders"]),
            "common_lmax": common_order["common_lmax"],
            "absolute_error": common_order["absolute_error"],
            "relative_error": common_order["relative_error"],
            "scale": common_order["scale"],
            "budget": common_order["budget"],
            "conditions": [
                item.balanced_condition_number
                for item in common_order["diagnostics"]
            ],
            "pairs": [
                {
                    "pair": list(record.particle_indices),
                    "final_lmax": record.final_lmax,
                    "eligible": record.eligible,
                    "failure_stage": record.failure_stage,
                    "final_changes": {
                        channel.channel: channel.final_successive_change
                        for channel in record.convergence
                    },
                }
                for record in common_order["result"].pair_ledger
            ],
        }

    summary = {
        "initial": summarize(_common_order_evidence()),
        "fallback": summarize(_fallback_common_order_evidence()),
    }
    encoded = json.dumps(summary, sort_keys=True)
    print("P1_3A_COMMON_ORDER=" + encoded)
    assert '"fallback"' in encoded
