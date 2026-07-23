"""Cached Gaunt coefficients for complex spherical harmonics."""

from functools import lru_cache
from math import pi, sqrt

from sympy.physics.wigner import wigner_3j


@lru_cache(maxsize=None)
def gaunt_coefficient(n_source: int, m_source: int, q: int, mu: int, n_target: int, m_target: int) -> float:
    """Return the cached complex-spherical-harmonic Gaunt integral."""
    integers = (n_source, m_source, q, mu, n_target, m_target)
    if any(not isinstance(value, int) for value in integers):
        raise ValueError("Gaunt indices must be integers")
    if n_source < 0 or q < 0 or n_target < 0 or abs(m_source) > n_source or abs(mu) > q or abs(m_target) > n_target:
        return 0.0
    if m_source + mu - m_target != 0 or not abs(n_source - q) <= n_target <= n_source + q or (n_source + q + n_target) % 2:
        return 0.0
    factor = sqrt((2 * n_source + 1) * (2 * q + 1) * (2 * n_target + 1) / (4 * pi))
    value = factor * wigner_3j(n_source, q, n_target, 0, 0, 0) * wigner_3j(n_source, q, n_target, m_source, mu, -m_target)
    return float(value)
