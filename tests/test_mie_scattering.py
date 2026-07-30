import numpy as np
import pytest
from scipy.special import spherical_jn, spherical_yn

from acoustic_ms.mie_scattering import (
    fluid_sphere_mie_scattering_coefficients,
    material_ratios_from_contrasts,
    mie_scattering_coefficients_from_contrasts,
    rigid_sphere_scattering_coefficients,
)
from acoustic_ms.multipolar_scattering import (
    rayleigh_multipolar_scattering_coefficients,
)


def _hankel(ell, value, derivative=False):
    return spherical_jn(ell, value, derivative=derivative) + 1j * spherical_yn(
        ell, value, derivative=derivative
    )


def _boundary_oracle(ka, density_ratio, compressibility_ratio, ell):
    internal_ka = ka * np.sqrt(density_ratio * compressibility_ratio)
    beta = np.sqrt(compressibility_ratio / density_ratio)
    matrix = np.array(
        [
            [_hankel(ell, ka), -spherical_jn(ell, internal_ka)],
            [
                _hankel(ell, ka, derivative=True),
                -beta * spherical_jn(ell, internal_ka, derivative=True),
            ],
        ],
        dtype=complex,
    )
    rhs = -np.array(
        [
            spherical_jn(ell, ka),
            spherical_jn(ell, ka, derivative=True),
        ],
        dtype=complex,
    )
    return np.linalg.solve(matrix, rhs)


def _boundary_residuals(
    ka, density_ratio, compressibility_ratio, ell, scattered, internal
):
    internal_ka = ka * np.sqrt(density_ratio * compressibility_ratio)
    beta = np.sqrt(compressibility_ratio / density_ratio)
    pressure = (
        spherical_jn(ell, ka)
        + scattered * _hankel(ell, ka)
        - internal * spherical_jn(ell, internal_ka)
    )
    velocity = (
        spherical_jn(ell, ka, derivative=True)
        + scattered * _hankel(ell, ka, derivative=True)
        - beta
        * internal
        * spherical_jn(ell, internal_ka, derivative=True)
    )
    return pressure, velocity


def test_material_ratios_and_exact_matching():
    density, compressibility, sound_speed = material_ratios_from_contrasts(
        0.25, 0.4
    )
    assert density == pytest.approx(2.0)
    assert compressibility == pytest.approx(0.75)
    assert sound_speed == pytest.approx(1 / np.sqrt(1.5))
    coefficients = fluid_sphere_mie_scattering_coefficients(0.1, 1.0, 1.0, 5)
    assert np.array_equal(coefficients, np.zeros(6, dtype=complex))


@pytest.mark.parametrize("ka", [0.01, 0.05, 0.1])
@pytest.mark.parametrize("contrasts", [(0.25, 0.4), (0.0, 0.8)])
def test_fluid_coefficients_match_independent_boundary_oracle(ka, contrasts):
    density, compressibility, _ = material_ratios_from_contrasts(*contrasts)
    actual = fluid_sphere_mie_scattering_coefficients(
        ka, density, compressibility, 5
    )
    expected = np.array(
        [
            _boundary_oracle(ka, density, compressibility, ell)[0]
            for ell in range(6)
        ]
    )
    assert np.allclose(actual, expected, rtol=4e-12, atol=2e-25)


@pytest.mark.parametrize("contrasts", [(0.25, 0.4), (-0.2, 0.8)])
def test_pressure_and_radial_velocity_boundary_residuals(contrasts):
    density, compressibility, _ = material_ratios_from_contrasts(*contrasts)
    for ka in (0.01, 0.05, 0.1):
        coefficients = fluid_sphere_mie_scattering_coefficients(
            ka, density, compressibility, 5
        )
        for ell, coefficient in enumerate(coefficients):
            _, internal = _boundary_oracle(
                ka, density, compressibility, ell
            )
            pressure, velocity = _boundary_residuals(
                ka,
                density,
                compressibility,
                ell,
                coefficient,
                internal,
            )
            assert abs(pressure) < 8e-16
            assert abs(velocity) < 8e-16


def test_rayleigh_limit_and_first_relative_correction_order():
    f0, f1 = 0.2, 0.6
    sizes = np.array([0.02, 0.01, 0.005])
    errors = {ell: [] for ell in range(6)}
    for ka in sizes:
        exact = mie_scattering_coefficients_from_contrasts(ka, f0, f1, 5)
        rayleigh = rayleigh_multipolar_scattering_coefficients(
            ka, f0, f1, 5
        )
        for ell in range(6):
            assert np.allclose(exact[ell] / rayleigh[ell], 1.0, rtol=8e-4)
            errors[ell].append(abs(exact[ell] / rayleigh[ell] - 1.0))
    for ell in range(6):
        slope = np.polyfit(np.log(sizes), np.log(errors[ell]), 1)[0]
        assert slope == pytest.approx(2.0, abs=0.03)


def test_rayleigh_powers_and_published_odd_order_coefficients():
    f1 = 0.8
    high = rayleigh_multipolar_scattering_coefficients(0.02, 0.0, f1, 5)
    low = rayleigh_multipolar_scattering_coefficients(0.01, 0.0, f1, 5)
    for ell in range(1, 6):
        assert low[ell] / high[ell] == pytest.approx(2.0 ** (-(2 * ell + 1)))
    ka = 0.07
    expected = {
        1: 1j * f1 * ka**3 / 6,
        3: 1j * f1 * ka**7 / (350 * (7 - f1)),
        5: 1j * f1 * ka**11 / (1309770 * (11 - 2 * f1)),
    }
    coefficients = rayleigh_multipolar_scattering_coefficients(
        ka, 0.0, f1, 5
    )
    for ell, value in expected.items():
        assert coefficients[ell] == pytest.approx(value, rel=2e-15)


def test_rigid_formula_and_finite_density_limit():
    ka = 0.073
    rigid = rigid_sphere_scattering_coefficients(ka, 5)
    direct = np.array(
        [
            -spherical_jn(ell, ka, derivative=True)
            / _hankel(ell, ka, derivative=True)
            for ell in range(6)
        ]
    )
    assert np.array_equal(rigid, direct)
    errors = []
    for density in (1e8, 1e10, 1e12):
        fluid = fluid_sphere_mie_scattering_coefficients(
            ka, density, 1.0, 5
        )
        errors.append(np.linalg.norm(fluid - rigid) / np.linalg.norm(rigid))
    assert errors[2] < errors[1] < errors[0]
    assert errors[-1] < 2e-5
    assert np.array_equal(
        mie_scattering_coefficients_from_contrasts(ka, 0.0, 1.0, 5),
        rigid,
    )


def test_positive_order_rigid_rayleigh_limit():
    for ka in (0.01, 0.005):
        exact = rigid_sphere_scattering_coefficients(ka, 5)
        rayleigh = rayleigh_multipolar_scattering_coefficients(
            ka, 0.0, 1.0, 5
        )
        assert np.allclose(
            exact[1:] / rayleigh[1:], 1.0, rtol=2e-4, atol=0.0
        )


@pytest.mark.parametrize(
    "density,compressibility", [(0.7, 1.3), (2.0, 0.75), (8.0, 0.2)]
)
def test_lossless_unitarity(density, compressibility):
    coefficients = fluid_sphere_mie_scattering_coefficients(
        0.37, density, compressibility, 8
    )
    defect = coefficients.real + np.abs(coefficients) ** 2
    assert np.max(np.abs(defect)) < 3e-17


@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf, 1 + 0j, True])
def test_physical_scalars_reject_non_real_or_non_finite_values(value):
    with pytest.raises((TypeError, ValueError)):
        fluid_sphere_mie_scattering_coefficients(value, 1.0, 1.0, 1)
    with pytest.raises((TypeError, ValueError)):
        material_ratios_from_contrasts(value, 0.2)


@pytest.mark.parametrize("ka", [0.0, -0.1])
def test_nonpositive_size_is_rejected(ka):
    with pytest.raises(ValueError):
        rigid_sphere_scattering_coefficients(ka, 1)


@pytest.mark.parametrize(
    "density,compressibility", [(0.0, 1.0), (-1.0, 1.0), (1.0, 0.0), (1.0, -1.0)]
)
def test_nonpositive_material_ratios_are_rejected(density, compressibility):
    with pytest.raises(ValueError):
        fluid_sphere_mie_scattering_coefficients(
            0.1, density, compressibility, 1
        )


@pytest.mark.parametrize("f0,f1", [(1.0, 0.0), (1.1, 0.0), (0.0, -2.0), (0.0, -3.0), (0.0, 1.1)])
def test_invalid_contrast_domains_are_rejected(f0, f1):
    with pytest.raises(ValueError):
        mie_scattering_coefficients_from_contrasts(0.1, f0, f1, 1)


@pytest.mark.parametrize("lmax", [-1, 1.2, True, 1 + 0j])
def test_invalid_lmax_is_rejected(lmax):
    with pytest.raises((TypeError, ValueError)):
        rigid_sphere_scattering_coefficients(0.1, lmax)


def test_lmax_zero_is_supported():
    assert fluid_sphere_mie_scattering_coefficients(0.1, 2.0, 0.75, 0).shape == (1,)
