import numpy as np
from scipy.special import spherical_jn, spherical_yn

from acoustic_ms.special import cartesian_to_spherical, spherical_hankel1, spherical_harmonic


def test_hankel_derivative_wronskian_and_harmonic_conjugacy():
    x = 1.7
    ell = 2
    assert np.allclose(spherical_hankel1(ell, x), spherical_jn(ell, x) + 1j * spherical_yn(ell, x))
    assert np.isclose(spherical_jn(ell, x) * spherical_yn(ell, x, derivative=True) - spherical_jn(ell, x, derivative=True) * spherical_yn(ell, x), 1 / x**2)
    assert np.allclose(spherical_harmonic(3, -2, 0.8, 1.2), spherical_harmonic(3, 2, 0.8, 1.2).conjugate())


def test_cartesian_angles():
    assert cartesian_to_spherical([0, 0, 2]) == (2.0, 0.0, 0.0)
    radius, theta, phi = cartesian_to_spherical([1, 1, 1])
    assert np.isclose(radius, np.sqrt(3)) and np.isclose(theta, np.arccos(1 / np.sqrt(3))) and np.isclose(phi, np.pi / 4)


def test_hankel_derivative_public_api_matches_recurrence():
    for ell, x in ((1, 0.9), (2, 1.7), (4, 3.2)):
        expected = spherical_hankel1(ell - 1, x) - (ell + 1) * spherical_hankel1(ell, x) / x
        assert np.allclose(spherical_hankel1(ell, x, derivative=True), expected, rtol=2e-13, atol=2e-13)


def test_harmonic_conjugacy_includes_odd_m_sign():
    theta, phi = 0.8, 1.2
    assert np.allclose(spherical_harmonic(3, -1, theta, phi), -spherical_harmonic(3, 1, theta, phi).conjugate())
    assert np.allclose(spherical_harmonic(3, -2, theta, phi), spherical_harmonic(3, 2, theta, phi).conjugate())


def test_spherical_harmonic_orthonormality_by_independent_quadrature():
    from scipy.special import roots_legendre, sph_harm_y
    u, weights = roots_legendre(80)
    theta = np.arccos(u)[:, None]
    phi = (2 * np.pi * np.arange(160) / 160)[None, :]
    measure = weights[:, None] * (2 * np.pi / 160)
    def inner(left, right):
        return np.sum(sph_harm_y(*left, theta, phi) * sph_harm_y(*right, theta, phi).conjugate() * measure)
    assert np.allclose(inner((3, 1), (3, 1)), 1.0, rtol=2e-12, atol=2e-12)
    assert np.allclose(inner((2, -1), (2, -1)), 1.0, rtol=2e-12, atol=2e-12)
    assert abs(inner((3, 1), (2, 1))) < 2e-12
