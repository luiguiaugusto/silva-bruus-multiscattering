"""Complete partial-wave radiation force for one spherical particle."""

from __future__ import annotations

import numbers

import numpy as np

from .multipoles import mode_count, mode_index


def _real_scalar(name: str, value: object, *, positive: bool) -> float:
    if isinstance(value, (bool, np.bool_, complex, np.complexfloating)) or not isinstance(value, numbers.Real):
        raise TypeError(f"{name} must be a real scalar")
    result = float(value)
    if not np.isfinite(result) or (result <= 0.0 if positive else result < 0.0):
        qualifier = "positive" if positive else "non-negative"
        raise ValueError(f"{name} must be finite and {qualifier}")
    return result


def complete_radiation_force_from_bsc(
    effective_incident_coefficients: object,
    scattering_coefficients: object,
    k: object,
    energy_density: object,
) -> np.ndarray:
    r"""Return the complete real Cartesian radiation force.

    The BSCs are the effective incident coefficients ``b_nm`` through
    ``lmax`` and scattering coefficients are indexed by multipole order.
    With the project convention ``E_LAS = 2 E0``, the implemented coupling
    is ``Gamma_n = s_n + s_{n+1}^* + 2 s_n s_{n+1}^*``.
    """

    wave_number = _real_scalar("k", k, positive=True)
    energy = _real_scalar("energy_density", energy_density, positive=False)
    bsc = np.asarray(effective_incident_coefficients, dtype=complex)
    scattering = np.asarray(scattering_coefficients, dtype=complex)
    if bsc.ndim != 1 or scattering.ndim != 1:
        raise ValueError("BSCs and scattering coefficients must be one-dimensional")
    if len(scattering) < 2 or len(bsc) != mode_count(len(scattering) - 1):
        raise ValueError("BSC length must equal (lmax + 1)**2")
    if not np.all(np.isfinite(bsc)) or not np.all(np.isfinite(scattering)):
        raise ValueError("BSCs and scattering coefficients must be finite")
    if energy == 0.0:
        return np.zeros(3)

    transverse = 0.0j
    longitudinal_sum = 0.0j
    for ell in range(len(scattering) - 1):
        gamma = (
            scattering[ell]
            + np.conj(scattering[ell + 1])
            + 2.0 * scattering[ell] * np.conj(scattering[ell + 1])
        )
        for m in range(-ell, ell + 1):
            current = bsc[mode_index(ell, m)]
            transverse_weight = np.sqrt(
                (ell + m + 1) * (ell + m + 2) / ((2 * ell + 1) * (2 * ell + 3))
            )
            transverse += transverse_weight * (
                gamma * current * np.conj(bsc[mode_index(ell + 1, m + 1)])
                + np.conj(gamma)
                * np.conj(bsc[mode_index(ell, -m)])
                * bsc[mode_index(ell + 1, -m - 1)]
            )
            longitudinal_weight = np.sqrt(
                (ell - m + 1) * (ell + m + 1) / ((2 * ell + 1) * (2 * ell + 3))
            )
            longitudinal_sum += (
                longitudinal_weight
                * gamma
                * current
                * np.conj(bsc[mode_index(ell + 1, m)])
            )
    transverse *= 1j * energy / wave_number**2
    longitudinal = 2.0 * energy / wave_number**2 * np.imag(longitudinal_sum)
    return np.array([transverse.real, transverse.imag, longitudinal], dtype=float)
