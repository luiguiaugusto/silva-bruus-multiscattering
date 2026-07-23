"""Outgoing-to-regular translation, oriented target <- source."""

import numpy as np

from .gaunt import gaunt_coefficient
from .multipoles import mode_count, modes
from .special import cartesian_to_spherical, spherical_hankel1, spherical_harmonic


def _position(value: object, name: str) -> np.ndarray:
    position = np.asarray(value, dtype=float)
    if position.shape != (3,) or not np.all(np.isfinite(position)):
        raise ValueError(f"{name} must be a finite vector with shape (3,)")
    return position


def separation_coefficient(n_target: int, m_target: int, n_source: int, m_source: int, k: float, target_position_xyz: object, source_position_xyz: object) -> complex:
    """Return one translation coefficient with rows target and columns source.

    The displacement is explicitly ``source_position - target_position``.
    """
    if not np.isfinite(k) or k <= 0.0:
        raise ValueError("k must be finite and positive")
    target = _position(target_position_xyz, "target_position_xyz")
    source = _position(source_position_xyz, "source_position_xyz")
    radius, theta, phi = cartesian_to_spherical(source - target)
    if radius == 0.0:
        raise ValueError("source and target positions must not coincide")
    if min(n_target, n_source) < 0 or abs(m_target) > n_target or abs(m_source) > n_source:
        raise ValueError("invalid multipole indices")
    q0 = max(abs(n_target - n_source), abs(m_target - m_source))
    if q0 % 2 != (n_target + n_source) % 2:
        q0 += 1
    total = 0.0j
    mu = m_target - m_source
    for q in range(q0, n_target + n_source + 1, 2):
        gaunt = gaunt_coefficient(n_source, m_source, q, mu, n_target, m_target)
        total += (1j**q) * spherical_hankel1(q, k * radius) * np.conj(spherical_harmonic(q, mu, theta, phi)) * gaunt
    return complex(4.0 * np.pi * (-1) ** (n_target + n_source + m_target) * 1j ** (n_target - n_source) * total)


def translation_matrix(k: float, target_position_xyz: object, source_position_xyz: object, target_lmax: int, source_lmax: int | None = None) -> np.ndarray:
    """Return outgoing-source to regular-target matrix, oriented target <- source."""
    if source_lmax is None:
        source_lmax = target_lmax
    target_modes = modes(target_lmax)
    source_modes = modes(source_lmax)
    matrix = np.empty((mode_count(target_lmax), mode_count(source_lmax)), dtype=complex)
    for row, (n_target, m_target) in enumerate(target_modes):
        for column, (n_source, m_source) in enumerate(source_modes):
            matrix[row, column] = separation_coefficient(n_target, m_target, n_source, m_source, k, target_position_xyz, source_position_xyz)
    return matrix
