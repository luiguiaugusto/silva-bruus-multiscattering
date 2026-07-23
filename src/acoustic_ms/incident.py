"""Beam-shape coefficients for the nodal standing wave sin(k z)."""

import numpy as np

from .multipoles import mode_count, mode_index


def nodal_standing_wave_coefficients(lmax: int) -> np.ndarray:
    """Return external coefficients in the standard multipole ordering."""
    if not isinstance(lmax, int) or lmax < 0:
        raise ValueError("lmax must be a non-negative integer")
    coefficients = np.zeros(mode_count(lmax), dtype=complex)
    for ell in range(1, lmax + 1, 2):
        coefficients[mode_index(ell, 0)] = (-1) ** ((ell - 1) // 2) * np.sqrt(4.0 * np.pi * (2 * ell + 1))
    return coefficients
