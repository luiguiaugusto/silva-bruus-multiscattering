"""Vector-error metrics for Model A/B/C force comparisons."""

import numpy as np


def _arrays(reference: object, model: object) -> tuple[np.ndarray, np.ndarray]:
    reference, model = np.asarray(reference, dtype=float), np.asarray(model, dtype=float)
    if reference.shape != model.shape or reference.ndim != 2 or reference.shape[1] != 2:
        raise ValueError("reference and model must share shape (N, 2)")
    if not np.all(np.isfinite(reference)) or not np.all(np.isfinite(model)):
        raise ValueError("force arrays must be finite")
    return reference, model


def symmetric_particle_errors(reference: object, model: object) -> np.ndarray:
    """Return bounded symmetric vector error per particle, in [0, 2]."""
    reference, model = _arrays(reference, model)
    scale = max(float(np.max(np.linalg.norm(reference, axis=1))), float(np.max(np.linalg.norm(model, axis=1))), 1.0)
    threshold = 128 * np.finfo(float).eps * scale
    ref_norm, model_norm = np.linalg.norm(reference, axis=1), np.linalg.norm(model, axis=1)
    denominator = ref_norm + model_norm
    errors = np.divide(2 * np.linalg.norm(reference - model, axis=1), denominator, out=np.zeros_like(denominator), where=denominator > threshold)
    return np.clip(errors, 0.0, 2.0)


def rms_relative_error(reference: object, model: object) -> float:
    """Return global RMS relative error, or infinity for nonzero model vs zero ref."""
    reference, model = _arrays(reference, model)
    reference_squared = float(np.sum(reference**2))
    model_squared = float(np.sum(model**2))
    if reference_squared == 0.0:
        return 0.0 if model_squared == 0.0 else float("inf")
    return float(np.sqrt(np.sum((reference - model) ** 2) / reference_squared))


def angular_errors_degrees(reference: object, model: object) -> np.ndarray:
    """Return per-particle angle in degrees; undefined zero-force directions are NaN."""
    reference, model = _arrays(reference, model)
    ref_norm, model_norm = np.linalg.norm(reference, axis=1), np.linalg.norm(model, axis=1)
    result = np.full(len(reference), np.nan)
    valid = (ref_norm > 0.0) & (model_norm > 0.0)
    cosine = np.sum(reference[valid] * model[valid], axis=1) / (ref_norm[valid] * model_norm[valid])
    result[valid] = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
    return result
