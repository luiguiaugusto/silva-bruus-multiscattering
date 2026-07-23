import numpy as np

from acoustic_ms.gaunt import gaunt_coefficient


def test_known_gaunt_values_and_selection_zeroes():
    assert np.isclose(gaunt_coefficient(1, 0, 0, 0, 1, 0), 1 / np.sqrt(4 * np.pi))
    assert np.isclose(gaunt_coefficient(1, 0, 2, 0, 1, 0), np.sqrt(5) / (5 * np.sqrt(np.pi)))
    assert gaunt_coefficient(1, 1, 1, 1, 1, 0) == 0.0
    assert gaunt_coefficient(1, 0, 1, 0, 1, 0) == 0.0


def _gaunt_quadrature(n_source, m_source, q, mu, n_target, m_target):
    from scipy.special import roots_legendre, sph_harm_y
    u, weights = roots_legendre(100)
    theta = np.arccos(u)[:, None]
    phi = (2 * np.pi * np.arange(192) / 192)[None, :]
    measure = weights[:, None] * (2 * np.pi / 192)
    return np.sum(sph_harm_y(n_source, m_source, theta, phi) * sph_harm_y(q, mu, theta, phi) * sph_harm_y(n_target, -m_target, theta, phi) * measure)


def test_nontrivial_gaunt_coefficients_match_independent_quadrature():
    cases = ((1, 1, 1, -1, 0, 0), (1, 0, 2, 0, 1, 0), (2, 1, 2, -1, 2, 0))
    for indices in cases:
        actual = gaunt_coefficient(*indices)
        expected = _gaunt_quadrature(*indices)
        assert abs(actual) > 1e-12
        assert np.allclose(actual, expected, rtol=3e-12, atol=3e-12)
