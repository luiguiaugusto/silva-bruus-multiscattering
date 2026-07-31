"""Three-dimensional comparison metrics for the T12 Model E sentinels."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]


def _finite_vectors(vectors: ArrayLike, *, name: str) -> FloatArray:
    values = np.asarray(vectors, dtype=float)
    if values.ndim != 2 or values.shape[1] != 3:
        raise ValueError(f"{name} must have shape (N, 3)")
    if values.shape[0] == 0:
        raise ValueError(f"{name} must contain at least one vector")
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{name} must contain only finite values")
    return values


def rms_vector_magnitude_xyz(vectors: ArrayLike) -> float:
    """Return the RMS of three-dimensional vector magnitudes."""

    values = _finite_vectors(vectors, name="vectors")
    return float(np.sqrt(np.mean(np.sum(values * values, axis=1))))


def _resolved_scale(*vectors: FloatArray) -> tuple[float, float]:
    scale = max(rms_vector_magnitude_xyz(value) for value in vectors)
    tolerance = 128.0 * np.finfo(float).eps * scale
    return scale, tolerance


def normalized_rms_error_xyz(reference: ArrayLike, model: ArrayLike) -> tuple[float, bool]:
    """Return ``RMS(model-reference) / RMS(reference)`` and applicability."""

    reference_values = _finite_vectors(reference, name="reference")
    model_values = _finite_vectors(model, name="model")
    if reference_values.shape != model_values.shape:
        raise ValueError("reference and model must have matching shapes")
    _, tolerance = _resolved_scale(reference_values, model_values)
    denominator = rms_vector_magnitude_xyz(reference_values)
    absolute = rms_vector_magnitude_xyz(model_values - reference_values)
    if denominator <= tolerance:
        return absolute, False
    return absolute / denominator, True


def symmetric_rms_error_xyz(reference: ArrayLike, model: ArrayLike) -> float:
    """Return the bounded symmetric RMS difference of two vector fields."""

    reference_values = _finite_vectors(reference, name="reference")
    model_values = _finite_vectors(model, name="model")
    if reference_values.shape != model_values.shape:
        raise ValueError("reference and model must have matching shapes")
    denominator = (
        rms_vector_magnitude_xyz(reference_values)
        + rms_vector_magnitude_xyz(model_values)
    )
    _, tolerance = _resolved_scale(reference_values, model_values)
    if denominator <= tolerance:
        return 0.0
    return 2.0 * rms_vector_magnitude_xyz(model_values - reference_values) / denominator


@dataclass(frozen=True)
class ModelEForceComparison:
    """Signed-vector decomposition and RMS diagnostics for one sentinel."""

    epsilon_a_e: float
    epsilon_a_e_applicable: bool
    epsilon_d_e: float
    epsilon_d_e_applicable: bool
    epsilon_a_external_scattered: float
    epsilon_a_external_scattered_applicable: bool
    epsilon_d_external_scattered: float
    epsilon_d_external_scattered_applicable: bool
    symmetric_a_e: float
    symmetric_d_e: float
    rms_model_e_interaction: float
    rms_model_e_external_scattered: float
    rms_model_e_scattered_scattered: float
    rms_d_minus_a: float
    rms_mie_external_correction: float
    rms_scattered_scattered: float
    x_d_minus_a: float
    x_mie_external: float
    x_scattered_scattered: float
    cancellation_ratio: float
    cancellation_ratio_applicable: bool
    decomposition_max_abs_error: float
    decomposition_relative_error: float
    max_abs_interaction_fz: float


def compare_model_e_forces(
    model_a_forces_xy: ArrayLike,
    model_d_forces_xy: ArrayLike,
    model_e_interaction_forces_xyz: ArrayLike,
    model_e_external_scattered_forces_xyz: ArrayLike,
    model_e_scattered_scattered_forces_xyz: ArrayLike,
) -> ModelEForceComparison:
    """Compare Models A/D with the full three-dimensional Model E interaction."""

    a_xy = np.asarray(model_a_forces_xy, dtype=float)
    d_xy = np.asarray(model_d_forces_xy, dtype=float)
    if a_xy.ndim != 2 or a_xy.shape[1] != 2 or a_xy.shape[0] == 0:
        raise ValueError("model_a_forces_xy must have shape (N, 2)")
    if d_xy.shape != a_xy.shape:
        raise ValueError("model_d_forces_xy must match model_a_forces_xy")
    if not np.all(np.isfinite(a_xy)) or not np.all(np.isfinite(d_xy)):
        raise ValueError("Models A and D must contain only finite values")

    interaction = _finite_vectors(
        model_e_interaction_forces_xyz,
        name="model_e_interaction_forces_xyz",
    )
    external_scattered = _finite_vectors(
        model_e_external_scattered_forces_xyz,
        name="model_e_external_scattered_forces_xyz",
    )
    scattered_scattered = _finite_vectors(
        model_e_scattered_scattered_forces_xyz,
        name="model_e_scattered_scattered_forces_xyz",
    )
    if interaction.shape[0] != a_xy.shape[0]:
        raise ValueError("all models must contain the same number of particles")
    if external_scattered.shape != interaction.shape or scattered_scattered.shape != interaction.shape:
        raise ValueError("all Model E force fields must have matching shapes")

    a = np.column_stack((a_xy, np.zeros(a_xy.shape[0], dtype=float)))
    d = np.column_stack((d_xy, np.zeros(d_xy.shape[0], dtype=float)))
    epsilon_a_e, epsilon_a_e_applicable = normalized_rms_error_xyz(interaction, a)
    epsilon_d_e, epsilon_d_e_applicable = normalized_rms_error_xyz(interaction, d)
    epsilon_a_ext, epsilon_a_ext_applicable = normalized_rms_error_xyz(
        external_scattered,
        a,
    )
    epsilon_d_ext, epsilon_d_ext_applicable = normalized_rms_error_xyz(
        external_scattered,
        d,
    )

    rms_interaction = rms_vector_magnitude_xyz(interaction)
    rms_external_scattered = rms_vector_magnitude_xyz(external_scattered)
    rms_scattered_scattered = rms_vector_magnitude_xyz(scattered_scattered)
    d_minus_a = d - a
    mie_external = external_scattered - d
    identity_lhs = interaction - a
    identity_rhs = d_minus_a + mie_external + scattered_scattered
    identity_error = identity_lhs - identity_rhs
    identity_scale, identity_tolerance = _resolved_scale(identity_lhs, identity_rhs)
    identity_absolute = rms_vector_magnitude_xyz(identity_error)
    identity_relative = (
        identity_absolute / identity_scale
        if identity_scale > identity_tolerance
        else identity_absolute
    )

    rms_d_minus_a = rms_vector_magnitude_xyz(d_minus_a)
    rms_mie_external = rms_vector_magnitude_xyz(mie_external)
    denominator = rms_interaction
    _, denominator_tolerance = _resolved_scale(interaction, a, d)
    if denominator > denominator_tolerance:
        x_d_minus_a = rms_d_minus_a / denominator
        x_mie_external = rms_mie_external / denominator
        x_scattered_scattered = rms_scattered_scattered / denominator
    else:
        x_d_minus_a = 0.0
        x_mie_external = 0.0
        x_scattered_scattered = 0.0

    component_sum = rms_d_minus_a + rms_mie_external + rms_scattered_scattered
    correction_rms = rms_vector_magnitude_xyz(identity_lhs)
    _, cancellation_tolerance = _resolved_scale(
        d_minus_a,
        mie_external,
        scattered_scattered,
    )
    cancellation_applicable = (
        component_sum > cancellation_tolerance
        and correction_rms > cancellation_tolerance
    )
    cancellation = (
        component_sum / correction_rms
        if cancellation_applicable
        else 0.0
    )

    return ModelEForceComparison(
        epsilon_a_e=epsilon_a_e,
        epsilon_a_e_applicable=epsilon_a_e_applicable,
        epsilon_d_e=epsilon_d_e,
        epsilon_d_e_applicable=epsilon_d_e_applicable,
        epsilon_a_external_scattered=epsilon_a_ext,
        epsilon_a_external_scattered_applicable=epsilon_a_ext_applicable,
        epsilon_d_external_scattered=epsilon_d_ext,
        epsilon_d_external_scattered_applicable=epsilon_d_ext_applicable,
        symmetric_a_e=symmetric_rms_error_xyz(interaction, a),
        symmetric_d_e=symmetric_rms_error_xyz(interaction, d),
        rms_model_e_interaction=rms_interaction,
        rms_model_e_external_scattered=rms_external_scattered,
        rms_model_e_scattered_scattered=rms_scattered_scattered,
        rms_d_minus_a=rms_d_minus_a,
        rms_mie_external_correction=rms_mie_external,
        rms_scattered_scattered=rms_scattered_scattered,
        x_d_minus_a=x_d_minus_a,
        x_mie_external=x_mie_external,
        x_scattered_scattered=x_scattered_scattered,
        cancellation_ratio=cancellation,
        cancellation_ratio_applicable=cancellation_applicable,
        decomposition_max_abs_error=float(np.max(np.abs(identity_error))),
        decomposition_relative_error=identity_relative,
        max_abs_interaction_fz=float(np.max(np.abs(interaction[:, 2]))),
    )
