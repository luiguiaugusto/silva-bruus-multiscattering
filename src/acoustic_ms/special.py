"""Special functions using the Condon--Shortley spherical-harmonic convention."""

import numpy as np
from scipy.special import spherical_jn, spherical_yn, sph_harm_y


def spherical_hankel1(ell: int, x: object, derivative: bool = False) -> object:
    """Return outgoing ``h_ell^(1)(x)`` or its derivative with respect to ``x``."""
    if not isinstance(ell, int) or ell < 0:
        raise ValueError("ell must be a non-negative integer")
    values = np.asarray(x)
    if not np.all(np.isfinite(values)):
        raise ValueError("x must be finite")
    return spherical_jn(ell, x, derivative=derivative) + 1j * spherical_yn(ell, x, derivative=derivative)


def spherical_harmonic(ell: int, m: int, theta: object, phi: object) -> object:
    """Return complex orthonormal ``Y_ell^m(theta, phi)``."""
    if not isinstance(ell, int) or not isinstance(m, int) or ell < 0 or abs(m) > ell:
        raise ValueError("invalid spherical-harmonic indices")
    if not np.all(np.isfinite(theta)) or not np.all(np.isfinite(phi)):
        raise ValueError("theta and phi must be finite")
    return sph_harm_y(ell, m, theta, phi)


def cartesian_to_spherical(vector_xyz: object) -> tuple[float, float, float]:
    """Return ``(r, theta, phi)`` for a finite three-component Cartesian vector."""
    vector = np.asarray(vector_xyz, dtype=float)
    if vector.shape != (3,) or not np.all(np.isfinite(vector)):
        raise ValueError("vector_xyz must be a finite vector with shape (3,)")
    radius = float(np.linalg.norm(vector))
    if radius == 0.0:
        return 0.0, 0.0, 0.0
    theta = float(np.arccos(np.clip(vector[2] / radius, -1.0, 1.0)))
    phi = float(np.mod(np.arctan2(vector[1], vector[0]), 2.0 * np.pi))
    return radius, theta, phi
