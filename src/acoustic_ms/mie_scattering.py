r"""Exact partial-wave coefficients for a lossless fluid sphere.

The coefficients use the project's :math:`e^{-i\omega t}` convention and an
outgoing spherical Hankel function of the first kind.  This module describes
an isolated sphere only; it is intentionally not connected to Model D.
"""

from __future__ import annotations

import numbers

import numpy as np
from scipy.special import spherical_jn, spherical_yn


def _real_scalar(name: str, value: object) -> float:
    if isinstance(value, (bool, np.bool_, complex, np.complexfloating)):
        raise TypeError(f"{name} must be a real scalar")
    if not isinstance(value, numbers.Real):
        raise TypeError(f"{name} must be a real scalar")
    result = float(value)
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _multipole_order(lmax: object) -> int:
    if isinstance(lmax, (bool, np.bool_)) or not isinstance(
        lmax, numbers.Integral
    ):
        raise TypeError("lmax must be an integer")
    result = int(lmax)
    if result < 0:
        raise ValueError("lmax must be non-negative")
    return result


def material_ratios_from_contrasts(
    f0: object, f1: object
) -> tuple[float, float, float]:
    """Return ``(rho_p/rho_0, kappa_p/kappa_0, c_p/c_0)``.

    This conversion applies to a finite-density lossless fluid sphere.  The
    exact rigid value ``f1 == 1`` is handled by
    :func:`mie_scattering_coefficients_from_contrasts` and is rejected here
    because its density ratio is infinite.
    """

    monopole = _real_scalar("f0", f0)
    dipole = _real_scalar("f1", f1)
    if monopole >= 1.0:
        raise ValueError("f0 must be smaller than 1")
    if not -2.0 < dipole < 1.0:
        raise ValueError("finite-density conversion requires -2 < f1 < 1")
    compressibility_ratio = 1.0 - monopole
    density_ratio = (2.0 + dipole) / (2.0 * (1.0 - dipole))
    sound_speed_ratio = 1.0 / np.sqrt(
        density_ratio * compressibility_ratio
    )
    return density_ratio, compressibility_ratio, float(sound_speed_ratio)


def _hankel_values(
    ell: np.ndarray, argument: float, *, derivative: bool = False
) -> np.ndarray:
    return spherical_jn(ell, argument, derivative=derivative) + 1j * spherical_yn(
        ell, argument, derivative=derivative
    )


def fluid_sphere_mie_scattering_coefficients(
    ka: object,
    density_ratio: object,
    compressibility_ratio: object,
    lmax: object,
) -> np.ndarray:
    """Return exact lossless-fluid-sphere coefficients through ``lmax``.

    The result is indexed by multipole order ``ell=0, ..., lmax``.  The
    routine accepts any positive size parameter that SciPy can represent;
    choosing a sufficiently large truncation remains the caller's
    responsibility outside the project's Rayleigh domain.
    """

    size = _real_scalar("ka", ka)
    density = _real_scalar("density_ratio", density_ratio)
    compressibility = _real_scalar(
        "compressibility_ratio", compressibility_ratio
    )
    order = _multipole_order(lmax)
    if size <= 0.0:
        raise ValueError("ka must be positive")
    if density <= 0.0 or compressibility <= 0.0:
        raise ValueError("material ratios must be positive")
    if density == 1.0 and compressibility == 1.0:
        return np.zeros(order + 1, dtype=complex)

    ell = np.arange(order + 1)
    internal_size = size * np.sqrt(density * compressibility)
    beta = np.sqrt(compressibility / density)
    exterior = spherical_jn(ell, size)
    exterior_derivative = spherical_jn(ell, size, derivative=True)
    interior = spherical_jn(ell, internal_size)
    interior_derivative = spherical_jn(
        ell, internal_size, derivative=True
    )
    hankel = _hankel_values(ell, size)
    hankel_derivative = _hankel_values(ell, size, derivative=True)
    numerator = (
        beta * exterior * interior_derivative
        - interior * exterior_derivative
    )
    denominator = (
        beta * hankel * interior_derivative
        - interior * hankel_derivative
    )
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        coefficients = -numerator / denominator
    if not np.all(np.isfinite(coefficients)):
        raise FloatingPointError(
            "Mie coefficients are not finite at the requested parameters"
        )
    return np.asarray(coefficients, dtype=complex)


def rigid_sphere_scattering_coefficients(
    ka: object, lmax: object
) -> np.ndarray:
    """Return the exact sound-hard (infinite-density) sphere coefficients."""

    size = _real_scalar("ka", ka)
    order = _multipole_order(lmax)
    if size <= 0.0:
        raise ValueError("ka must be positive")
    ell = np.arange(order + 1)
    numerator = spherical_jn(ell, size, derivative=True)
    denominator = _hankel_values(ell, size, derivative=True)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        coefficients = -numerator / denominator
    if not np.all(np.isfinite(coefficients)):
        raise FloatingPointError(
            "rigid-sphere coefficients are not finite at the requested parameters"
        )
    return np.asarray(coefficients, dtype=complex)


def mie_scattering_coefficients_from_contrasts(
    ka: object, f0: object, f1: object, lmax: object
) -> np.ndarray:
    """Return exact fluid or rigid coefficients from Silva--Bruus contrasts."""

    monopole = _real_scalar("f0", f0)
    dipole = _real_scalar("f1", f1)
    if monopole >= 1.0:
        raise ValueError("f0 must be smaller than 1")
    if not -2.0 < dipole <= 1.0:
        raise ValueError("f1 must satisfy -2 < f1 <= 1")
    if dipole == 1.0:
        return rigid_sphere_scattering_coefficients(ka, lmax)
    density, compressibility, _ = material_ratios_from_contrasts(
        monopole, dipole
    )
    return fluid_sphere_mie_scattering_coefficients(
        ka, density, compressibility, lmax
    )
