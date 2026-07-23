import numpy as np

from acoustic_ms.gaunt import gaunt_coefficient


def test_known_gaunt_values_and_selection_zeroes():
    assert np.isclose(gaunt_coefficient(1, 0, 0, 0, 1, 0), 1 / np.sqrt(4 * np.pi))
    assert np.isclose(gaunt_coefficient(1, 0, 2, 0, 1, 0), np.sqrt(5) / (5 * np.sqrt(np.pi)))
    assert gaunt_coefficient(1, 1, 1, 1, 1, 0) == 0.0
    assert gaunt_coefficient(1, 0, 1, 0, 1, 0) == 0.0
