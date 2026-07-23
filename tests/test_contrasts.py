import numpy as np
import pytest

from acoustic_ms.contrasts import dipole_contrast, monopole_contrast


def test_monopole_contrast_is_zero_for_matched_compressibility():
    assert monopole_contrast(1.0) == 0.0


def test_dipole_contrast_is_zero_for_matched_density():
    assert dipole_contrast(1.0) == 0.0


def test_dipole_contrast_has_correct_physical_limits():
    assert np.isclose(dipole_contrast(1.0e15), 1.0, rtol=0.0, atol=1.0e-14)
    assert np.isclose(dipole_contrast(1.0e-15), -2.0, rtol=0.0, atol=1.0e-14)


@pytest.mark.parametrize("ratio", [0.0, -1.0, np.nan, np.inf])
def test_invalid_density_ratio_is_rejected(ratio):
    with pytest.raises(ValueError):
        dipole_contrast(ratio)
