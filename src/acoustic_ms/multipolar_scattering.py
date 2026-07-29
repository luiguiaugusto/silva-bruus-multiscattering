"""Leading Rayleigh scattering coefficients through arbitrary multipole order."""

import numpy as np


def _odd_double_factorial(value: int) -> int:
    """Return an odd double factorial as an exact Python integer."""
    result = 1
    for factor in range(1, value + 1, 2):
        result *= factor
    return result


def rayleigh_multipolar_scattering_coefficients(
    ka: object,
    f0: object,
    f1: object,
    lmax: int,
) -> np.ndarray:
    """Return leading Rayleigh coefficients ``s_ell`` for ``0 <= ell <= lmax``.

    The monopole retains the established ``-i f0 (ka)^3 / 3`` expression.
    Every positive order uses the leading small-``ka`` term derived from the
    exact single-sphere coefficient.  These are not exact finite-frequency
    T-matrix coefficients.
    """
    if not isinstance(lmax, int) or isinstance(lmax, bool) or lmax < 1:
        raise ValueError("lmax must be a positive integer")
    try:
        ka, f0, f1 = (float(value) for value in (ka, f0, f1))
    except (TypeError, ValueError) as exc:
        raise ValueError("ka, f0, and f1 must be real scalars") from exc
    if not np.all(np.isfinite((ka, f0, f1))):
        raise ValueError("ka, f0, and f1 must be finite")
    if not 0.0 < ka <= 0.1:
        raise ValueError("the Rayleigh regime requires 0 < ka <= 0.1")
    if not -2.0 <= f1 <= 1.0:
        raise ValueError("f1 must lie in the physical interval [-2, 1]")

    coefficients = np.empty(lmax + 1, dtype=complex)
    coefficients[0] = -1j * f0 * ka**3 / 3.0
    for ell in range(1, lmax + 1):
        odd_lower = _odd_double_factorial(2 * ell - 1)
        odd_upper = _odd_double_factorial(2 * ell + 1)
        material = 2 * (2 * ell + 1) - (ell - 1) * f1
        coefficient = 3 * ell * f1 / (odd_lower * odd_upper * material)
        coefficients[ell] = 1j * coefficient * ka ** (2 * ell + 1)
    return coefficients
