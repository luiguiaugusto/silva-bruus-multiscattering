import numpy as np
import pytest
from scipy.special import spherical_jn

from acoustic_ms.special import cartesian_to_spherical, spherical_hankel1, spherical_harmonic
from acoustic_ms.translation import separation_coefficient, translation_matrix


def _outgoing(n, m, k, point, source):
    radius, theta, phi = cartesian_to_spherical(np.asarray(point) - np.asarray(source))
    return spherical_hankel1(n, k * radius) * spherical_harmonic(n, m, theta, phi)


def _regular_sum(n_source, m_source, k, point, target, source, lmax):
    matrix = translation_matrix(k, target, source, lmax, n_source)
    radius, theta, phi = cartesian_to_spherical(np.asarray(point) - np.asarray(target))
    return sum(matrix[ell * ell + ell + m, n_source * n_source + n_source + m_source] * spherical_jn(ell, k * radius) * spherical_harmonic(ell, m, theta, phi) for ell in range(lmax + 1) for m in range(-ell, ell + 1))


def test_analytic_planar_dipole_translation():
    k, distance = .4, 2.3
    actual = separation_coefficient(1, 0, 1, 0, k, [0, 0, 0], [distance, 0, 0])
    assert np.allclose(actual, spherical_hankel1(0, k * distance) + spherical_hankel1(2, k * distance))


@pytest.mark.parametrize("source_mode", [(0, 0), (1, 1)])
def test_direct_reexpansion_theorem(source_mode):
    target = np.array([.2, -.1, .3]); source = np.array([2.1, .7, 1.4]); point = np.array([.31, -.04, .38]); k = .8
    n, m = source_mode
    expected = _outgoing(n, m, k, point, source)
    actual = _regular_sum(n, m, k, point, target, source, 10)
    assert abs(actual - expected) / abs(expected) < 1e-9


def test_translation_shape_and_coincident_rejection():
    assert translation_matrix(.4, [0, 0, 0], [1, 0, 0], 1, 2).shape == (4, 9)
    with pytest.raises(ValueError): translation_matrix(.4, [0, 0, 0], [0, 0, 0], 1)


def test_reexpansion_error_converges_with_test_truncation():
    target = np.array([.2, -.1, .3]); source = np.array([2.1, .7, 1.4]); point = np.array([.31, -.04, .38]); k = .8
    expected = _outgoing(1, 1, k, point, source)
    errors = [abs(_regular_sum(1, 1, k, point, target, source, lmax) - expected) / abs(expected) for lmax in (2, 4, 6, 8, 10)]
    assert all(later < earlier for earlier, later in zip(errors, errors[1:]))
    assert errors[-1] < 1e-9
