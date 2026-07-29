"""Independent strict dimer reduction used by the fifth-order Eq. (30)."""

import numpy as np

from acoustic_ms import (
    corrected_nodal_pair_force_magnitude,
    mode_index,
    nodal_standing_wave_coefficients,
    rayleigh_multipolar_scattering_coefficients,
    separation_coefficient,
)


def _strict_reduced_force(ka, distance, f1):
    positions = np.array([[-distance / 2, 0.0, 0.0], [distance / 2, 0.0, 0.0]])
    selected = [(ell, m) for ell in (1, 3, 5) for m in range(0, ell + 1, 2)]
    scattering = rayleigh_multipolar_scattering_coefficients(ka, 0.0, f1, 5)
    incident = nodal_standing_wave_coefficients(5)
    matrix = np.eye(len(selected), dtype=complex)
    rhs = np.zeros(len(selected), dtype=complex)
    for row, (ell, m) in enumerate(selected):
        rhs[row] = scattering[ell] * incident[mode_index(ell, m)]
        for column, (source_ell, source_m) in enumerate(selected):
            translated = separation_coefficient(
                ell, m, source_ell, source_m, ka, positions[0], positions[1]
            )
            if source_m:
                translated += separation_coefficient(
                    ell, m, source_ell, -source_m, ka, positions[0], positions[1]
                )
            matrix[row, column] -= scattering[ell] * translated
    coefficients = np.linalg.solve(matrix, rhs)
    b21 = 0j
    for (ell, m), coefficient in zip(selected, coefficients):
        b21 += separation_coefficient(
            2, 1, ell, m, ka, positions[0], positions[1]
        ) * coefficient
        if m:
            b21 += separation_coefficient(
                2, 1, ell, -m, ka, positions[0], positions[1]
            ) * coefficient
    return float(-2 * np.sqrt(30 * np.pi) * ka * f1 * np.real(b21) / 15)


def test_equation_30_has_the_expected_asymptotic_trend():
    errors = []
    for ka in (0.1, 0.05, 0.025):
        reduced = _strict_reduced_force(ka, 2.1, 0.8)
        equation_30 = -corrected_nodal_pair_force_magnitude(ka, 1.0, 2.1, 1.0, 0.8)
        errors.append(abs(reduced - equation_30) / abs(reduced))
    np.testing.assert_allclose(
        errors,
        [0.0017061326608159644, 0.00043074064107634445, 0.0001079493043921584],
        rtol=3e-11,
        atol=3e-13,
    )
    assert errors[0] > errors[1] > errors[2]
