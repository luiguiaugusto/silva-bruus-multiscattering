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
