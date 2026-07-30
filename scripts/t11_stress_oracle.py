"""Independent stress-tensor quadrature used only by T11 validation."""

from __future__ import annotations

import numpy as np
from scipy.special import sph_harm_y, spherical_jn, spherical_yn


def _mode_index(ell: int, m: int) -> int:
    return ell * ell + ell + m


def stress_tensor_force(
    effective_incident_coefficients: object,
    scattering_coefficients: object,
    k: float,
    energy_density: float,
    integration_radius: float,
    theta_order: int,
    phi_count: int,
) -> np.ndarray:
    """Integrate the dimensionless inviscid radiation-stress traction."""

    bsc = np.asarray(effective_incident_coefficients, dtype=complex)
    scattering = np.asarray(scattering_coefficients, dtype=complex)
    lmax = len(scattering) - 1
    if bsc.shape != ((lmax + 1) ** 2,):
        raise ValueError("incompatible BSC and scattering-coefficient lengths")
    if not np.all(np.isfinite(bsc)) or not np.all(np.isfinite(scattering)):
        raise ValueError("coefficients must be finite")
    if not np.isfinite(k) or k <= 0.0 or not np.isfinite(energy_density) or energy_density < 0.0:
        raise ValueError("require finite k > 0 and energy density >= 0")
    if not np.isfinite(integration_radius) or integration_radius <= 0.0:
        raise ValueError("integration radius must be positive")
    if theta_order < 2 or phi_count < 4:
        raise ValueError("quadrature orders are too small")

    mu, mu_weights = np.polynomial.legendre.leggauss(theta_order)
    theta = np.arccos(mu)[:, None]
    phi = (2.0 * np.pi * np.arange(phi_count) / phi_count)[None, :]
    theta_grid = np.broadcast_to(theta, (theta_order, phi_count))
    phi_grid = np.broadcast_to(phi, (theta_order, phi_count))
    sin_theta = np.sqrt(1.0 - mu * mu)[:, None]
    cos_theta = mu[:, None]
    cos_phi = np.cos(phi_grid)
    sin_phi = np.sin(phi_grid)
    er = np.stack(
        (
            sin_theta * cos_phi,
            sin_theta * sin_phi,
            np.broadcast_to(cos_theta, theta_grid.shape),
        ),
        axis=-1,
    )
    etheta = np.stack(
        (
            cos_theta * cos_phi,
            cos_theta * sin_phi,
            -np.broadcast_to(sin_theta, theta_grid.shape),
        ),
        axis=-1,
    )
    ephi = np.stack((-sin_phi, cos_phi, np.zeros_like(phi_grid)), axis=-1)

    x = k * integration_radius
    psi = np.zeros(theta_grid.shape, dtype=complex)
    radial_gradient = np.zeros_like(psi)
    theta_gradient = np.zeros_like(psi)
    phi_gradient = np.zeros_like(psi)
    for ell in range(lmax + 1):
        radial = spherical_jn(ell, x) + scattering[ell] * (
            spherical_jn(ell, x) + 1j * spherical_yn(ell, x)
        )
        radial_derivative = spherical_jn(ell, x, derivative=True) + scattering[ell] * (
            spherical_jn(ell, x, derivative=True)
            + 1j * spherical_yn(ell, x, derivative=True)
        )
        for m in range(-ell, ell + 1):
            coefficient = bsc[_mode_index(ell, m)]
            if coefficient == 0.0:
                continue
            harmonic = sph_harm_y(ell, m, theta_grid, phi_grid)
            if m < ell:
                theta_derivative = (
                    m * cos_theta / sin_theta * harmonic
                    + np.sqrt((ell - m) * (ell + m + 1))
                    * np.exp(-1j * phi_grid)
                    * sph_harm_y(ell, m + 1, theta_grid, phi_grid)
                )
            else:
                theta_derivative = m * cos_theta / sin_theta * harmonic
            psi += coefficient * radial * harmonic
            radial_gradient += coefficient * radial_derivative * harmonic
            theta_gradient += coefficient * radial * theta_derivative / x
            phi_gradient += coefficient * radial * (1j * m) * harmonic / (x * sin_theta)

    gradient = (
        radial_gradient[..., None] * er
        + theta_gradient[..., None] * etheta
        + phi_gradient[..., None] * ephi
    )
    gradient_squared = np.sum(np.abs(gradient) ** 2, axis=-1)
    traction = -(
        gradient_squared - np.abs(psi) ** 2
    )[..., None] * er + 2.0 * np.real(gradient * np.conj(radial_gradient)[..., None])
    angular_weights = mu_weights[:, None] * (2.0 * np.pi / phi_count)
    integral = np.sum(traction * angular_weights[..., None], axis=(0, 1))
    return np.asarray(-energy_density * integration_radius**2 * integral, dtype=float)
