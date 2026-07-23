import numpy as np
import pytest

from acoustic_ms.scattering import rayleigh_scattering_coefficients


def test_rayleigh_coefficients_and_cubic_scaling():
    s0, s1 = rayleigh_scattering_coefficients(0.05, 0.6, -0.4)
    assert np.isclose(s0, -1j * .6 / 3 * .05**3)
    assert np.isclose(s1, 1j * -.4 / 6 * .05**3)
    assert np.allclose(rayleigh_scattering_coefficients(.05, 0, 0), 0)
    assert np.allclose(rayleigh_scattering_coefficients(.1, .6, -.4), 8 * np.array([s0, s1]))


def test_invalid_rayleigh_parameters():
    with pytest.raises(ValueError): rayleigh_scattering_coefficients(.11, 0, 0)
