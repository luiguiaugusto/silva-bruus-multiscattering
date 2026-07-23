import numpy as np

from acoustic_ms.incident import nodal_standing_wave_coefficients


def test_nodal_coefficients():
    coefficients = nodal_standing_wave_coefficients(3)
    assert np.isclose(coefficients[2], np.sqrt(12 * np.pi))
    assert np.isclose(coefficients[12], -np.sqrt(28 * np.pi))
    assert np.count_nonzero(coefficients) == 2
