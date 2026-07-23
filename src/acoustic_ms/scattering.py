"""Rayleigh scattering coefficients for a compressible sphere."""

import numpy as np


def rayleigh_scattering_coefficients(ka: float, f0: float, f1: float) -> np.ndarray:
    """Return dominant complex ``[s0, s1]`` Rayleigh coefficients."""
    ka, f0, f1 = (float(value) for value in (ka, f0, f1))
    if not all(np.isfinite((ka, f0, f1))) or not 0.0 < ka <= 0.1 or not -2.0 <= f1 <= 1.0:
        raise ValueError("require finite 0 < ka <= 0.1, finite f0, and -2 <= f1 <= 1")
    return np.array([-1j * f0 * ka**3 / 3.0, 1j * f1 * ka**3 / 6.0], dtype=complex)
