"""Dimensionless material contrast factors for a small sphere."""

import numpy as np


def _finite_scalar(value: float, name: str) -> float:
    """Return a finite scalar float or raise ``ValueError``."""
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a real finite scalar") from exc
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def monopole_contrast(compressibility_ratio: float) -> float:
    """Return ``f0 = 1 - kappa_p / kappa_0`` (dimensionless)."""
    ratio = _finite_scalar(compressibility_ratio, "compressibility_ratio")
    return 1.0 - ratio


def dipole_contrast(density_ratio: float) -> float:
    """Return ``f1`` for positive ``rho_p / rho_0`` (dimensionless)."""
    ratio = _finite_scalar(density_ratio, "density_ratio")
    if ratio <= 0.0:
        raise ValueError("density_ratio must be positive")
    return 2.0 * (ratio - 1.0) / (2.0 * ratio + 1.0)
