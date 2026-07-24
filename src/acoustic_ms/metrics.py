"""Vector-error metrics for Model A/B/C force comparisons."""

import numpy as np


def _arrays(reference: object, model: object) -> tuple[np.ndarray, np.ndarray]:
    reference, model = np.asarray(reference, dtype=float), np.asarray(model, dtype=float)
    if reference.shape != model.shape or reference.ndim != 2 or reference.shape[1] != 2:
        raise ValueError("reference and model must share shape (N, 2)")
    if not np.all(np.isfinite(reference)) or not np.all(np.isfinite(model)):
        raise ValueError("force arrays must be finite")
    return reference, model


def _norms_and_zero_tolerance(
    reference: object, model: object
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    """Validate force arrays and return their norms with a shared zero tolerance."""
    reference, model = _arrays(reference, model)
    ref_norm = np.linalg.norm(reference, axis=1)
    model_norm = np.linalg.norm(model, axis=1)
    scale = max(float(np.max(ref_norm, initial=0.0)), float(np.max(model_norm, initial=0.0)))
    tolerance = 128 * np.finfo(float).eps * scale
    return reference, model, ref_norm, model_norm, tolerance


def symmetric_particle_errors(reference: object, model: object) -> np.ndarray:
    """Return bounded symmetric vector error per particle, in [0, 2]."""
    reference, model, ref_norm, model_norm, tolerance = _norms_and_zero_tolerance(reference, model)
    denominator = ref_norm + model_norm
    both_zero = (ref_norm <= tolerance) & (model_norm <= tolerance)
    errors = np.divide(
        2 * np.linalg.norm(reference - model, axis=1),
        denominator,
        out=np.zeros_like(denominator),
        where=~both_zero,
    )
    return np.clip(errors, 0.0, 2.0)


def rms_relative_error(reference: object, model: object) -> float:
    """Return global RMS relative error, or infinity for nonzero model vs zero ref."""
    reference, model = _arrays(reference, model)
    reference_squared = float(np.sum(reference**2))
    model_squared = float(np.sum(model**2))
    if reference_squared == 0.0:
        return 0.0 if model_squared == 0.0 else float("inf")
    return float(np.sqrt(np.sum((reference - model) ** 2) / reference_squared))


def rms_vector_magnitude(vectors: object) -> float:
    """Return the RMS magnitude of a finite ``(N, 2)`` force array."""
    vectors = np.asarray(vectors, dtype=float)
    if vectors.ndim != 2 or vectors.shape[1] != 2 or len(vectors) == 0:
        raise ValueError("vectors must have nonempty shape (N, 2)")
    if not np.all(np.isfinite(vectors)):
        raise ValueError("vectors must be finite")
    return float(np.sqrt(np.mean(np.sum(vectors**2, axis=1))))


def angular_errors_degrees(reference: object, model: object) -> np.ndarray:
    """Return per-particle angle in degrees; undefined zero-force directions are NaN."""
    reference, model, ref_norm, model_norm, tolerance = _norms_and_zero_tolerance(reference, model)
    result = np.full(len(reference), np.nan)
    valid = (ref_norm > tolerance) & (model_norm > tolerance)
    cosine = np.sum(reference[valid] * model[valid], axis=1) / (ref_norm[valid] * model_norm[valid])
    result[valid] = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
    return result
