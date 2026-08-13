"""Exact-Mie global multiple scattering with the complete radiation force."""

from __future__ import annotations

from dataclasses import dataclass
import numbers

import numpy as np

from .complete_force import complete_radiation_force_from_bsc
from .mie_multiparticle import MieMultiparticleSolution, solve_mie_multiparticle_nodal
from .model_e_comparison import rms_vector_magnitude_xyz


@dataclass(frozen=True)
class ModelENodalResult:
    """Model-E coefficients, force channels, and decomposition diagnostic."""

    solution: MieMultiparticleSolution
    total_forces_xyz: np.ndarray
    external_forces_xyz: np.ndarray
    interaction_forces_xyz: np.ndarray
    external_scattered_forces_xyz: np.ndarray
    scattered_scattered_forces_xyz: np.ndarray
    scattered_incident_coefficients: np.ndarray
    decomposition_residual: float
    lmax: int


@dataclass(frozen=True)
class ModelENumericalDiagnostics:
    """Reusable established numerical gates for one Model-E result."""

    production_solver: str
    balanced_condition_number: float
    balanced_backward_error: float
    effective_incident_closure_error: float
    scattering_closure_error: float
    force_decomposition_residual: float
    max_abs_fz: float
    fz_tolerance: float
    finite: bool
    mode_dimension_consistent: bool
    planar_symmetry_pass: bool
    passed: bool


def _energy(value: object) -> float:
    if isinstance(value, (bool, np.bool_, complex, np.complexfloating)) or not isinstance(value, numbers.Real):
        raise TypeError("energy_density must be a real scalar")
    result = float(value)
    if not np.isfinite(result) or result < 0.0:
        raise ValueError("energy_density must be finite and non-negative")
    return result


def evaluate_model_e_numerical_diagnostics(
    result: ModelENodalResult,
    *,
    diagnostic_tolerance: float = 1.0e-12,
    maximum_balanced_condition: float = 10.0,
    planar_tolerance_factor: float = 128.0,
) -> ModelENumericalDiagnostics:
    """Apply the established T13--T14.1 finite, closure, and planar gates."""

    thresholds = np.asarray(
        [
            diagnostic_tolerance,
            maximum_balanced_condition,
            planar_tolerance_factor,
        ],
        dtype=float,
    )
    if not np.all(np.isfinite(thresholds)) or np.any(thresholds <= 0.0):
        raise ValueError("diagnostic thresholds must be finite and positive")

    solution = result.solution
    force_arrays = (
        np.asarray(result.total_forces_xyz, dtype=float),
        np.asarray(result.interaction_forces_xyz, dtype=float),
        np.asarray(result.external_scattered_forces_xyz, dtype=float),
        np.asarray(result.scattered_scattered_forces_xyz, dtype=float),
    )
    reference_shape = force_arrays[0].shape
    force_shapes_valid = (
        len(reference_shape) == 2
        and reference_shape[1:] == (3,)
        and reference_shape[0] > 0
        and all(values.shape == reference_shape for values in force_arrays)
    )
    force_finite = force_shapes_valid and all(
        np.all(np.isfinite(values)) for values in force_arrays
    )
    diagnostic_values = (
        float(solution.balanced_condition_number),
        float(solution.balanced_backward_error),
        float(solution.effective_incident_closure_error),
        float(solution.scattering_closure_error),
        float(result.decomposition_residual),
    )
    finite = bool(
        force_finite
        and np.all(np.isfinite(solution.scattered_coefficients))
        and all(np.isfinite(value) for value in diagnostic_values)
    )
    if force_finite:
        maximum_force_rms = max(
            rms_vector_magnitude_xyz(values) for values in force_arrays
        )
        max_abs_fz = max(
            float(np.max(np.abs(values[:, 2]))) for values in force_arrays
        )
        fz_tolerance = (
            planar_tolerance_factor * np.finfo(float).eps * maximum_force_rms
        )
        planar = max_abs_fz <= fz_tolerance
    else:
        max_abs_fz = float("inf")
        fz_tolerance = 0.0
        planar = False

    particle_count = len(force_arrays[0]) if force_shapes_valid else 0
    expected_dimension = particle_count * len(solution.active_modes)
    mode_consistent = bool(
        force_shapes_valid
        and solution.balanced_system_matrix.shape
        == (expected_dimension, expected_dimension)
        and len(solution.modes) == (int(result.lmax) + 1) ** 2
    )
    passed = bool(
        solution.production_solver == "balanced_sqrt"
        and finite
        and solution.balanced_condition_number < maximum_balanced_condition
        and solution.balanced_backward_error < diagnostic_tolerance
        and solution.effective_incident_closure_error < diagnostic_tolerance
        and solution.scattering_closure_error < diagnostic_tolerance
        and result.decomposition_residual < diagnostic_tolerance
        and mode_consistent
        and planar
    )
    return ModelENumericalDiagnostics(
        production_solver=str(solution.production_solver),
        balanced_condition_number=diagnostic_values[0],
        balanced_backward_error=diagnostic_values[1],
        effective_incident_closure_error=diagnostic_values[2],
        scattering_closure_error=diagnostic_values[3],
        force_decomposition_residual=diagnostic_values[4],
        max_abs_fz=max_abs_fz,
        fz_tolerance=fz_tolerance,
        finite=finite,
        mode_dimension_consistent=mode_consistent,
        planar_symmetry_pass=planar,
        passed=passed,
    )


def solve_model_e_nodal(
    positions_xyz: object,
    k: object,
    radius: object,
    energy_density: object,
    f0: object,
    f1: object,
    lmax: object,
    *,
    use_planar_symmetry: bool = True,
) -> ModelENodalResult:
    """Solve exact-Mie Model E and separate complete force contributions."""

    if isinstance(lmax, (bool, np.bool_)) or not isinstance(lmax, numbers.Integral):
        raise TypeError("lmax must be an integer")
    if int(lmax) < 2:
        raise ValueError("Model E requires lmax >= 2 for the complete force")
    energy = _energy(energy_density)
    solution = solve_mie_multiparticle_nodal(
        positions_xyz, k, radius, f0, f1, int(lmax),
        use_planar_symmetry=use_planar_symmetry,
    )
    particle_count = len(solution.effective_incident_coefficients)
    total = np.empty((particle_count, 3), dtype=float)
    external = np.empty_like(total)
    scattered_scattered = np.empty_like(total)
    scattered_incident = solution.effective_incident_coefficients - solution.external_incident_coefficients
    for particle in range(particle_count):
        total[particle] = complete_radiation_force_from_bsc(
            solution.effective_incident_coefficients[particle], solution.scattering_coefficients, k, energy
        )
        external[particle] = complete_radiation_force_from_bsc(
            solution.external_incident_coefficients[particle], solution.scattering_coefficients, k, energy
        )
        scattered_scattered[particle] = complete_radiation_force_from_bsc(
            scattered_incident[particle], solution.scattering_coefficients, k, energy
        )
    interaction = total - external
    external_scattered = interaction - scattered_scattered
    reconstruction = external + external_scattered + scattered_scattered
    residual = float(
        np.linalg.norm(total - reconstruction)
        / max(float(np.linalg.norm(total)), np.finfo(float).eps)
    )
    return ModelENodalResult(
        solution=solution,
        total_forces_xyz=total,
        external_forces_xyz=external,
        interaction_forces_xyz=interaction,
        external_scattered_forces_xyz=external_scattered,
        scattered_scattered_forces_xyz=scattered_scattered,
        scattered_incident_coefficients=scattered_incident,
        decomposition_residual=residual,
        lmax=int(lmax),
    )
