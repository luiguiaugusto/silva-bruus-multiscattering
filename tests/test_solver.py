import numpy as np
import pytest

from acoustic_ms.scattering import rayleigh_scattering_coefficients
from acoustic_ms.special import spherical_hankel1
from acoustic_ms.solver import solve_rayleigh_nodal


def test_single_particle_analytical_benchmark():
    ka, f1 = .08, .7
    solution = solve_rayleigh_nodal([[0, 0, 0]], ka, 1, -.3, f1)
    expected = 1j * f1 * np.sqrt(np.pi / 3) * ka**3
    assert abs(solution.coefficients[0, 2] - expected) / abs(expected) < 1e-13
    assert np.count_nonzero(solution.coefficients) == 1
    assert solution.residual_relative < 1e-13


@pytest.mark.parametrize("ka, separation", [(.03, 20.), (.08, 4.), (.1, 2.)])
def test_two_particle_analytical_benchmark(ka, separation):
    f1 = .6; k = ka
    solution = solve_rayleigh_nodal([[-separation / 2, 0, 0], [separation / 2, 0, 0]], k, 1, .2, f1)
    s1 = rayleigh_scattering_coefficients(ka, .2, f1)[1]
    c = spherical_hankel1(0, k * separation) + spherical_hankel1(2, k * separation)
    expected = s1 * np.sqrt(12 * np.pi) / (1 - s1 * c)
    assert np.max(abs(solution.coefficients[:, 2] - expected)) / abs(expected) < 1e-12
    assert solution.residual_relative < 1e-12


def test_solver_validation_and_three_particle_symmetry_structure():
    solution = solve_rayleigh_nodal([[0, 0, 0], [3, .5, 0], [-1, 2.5, 0]], .08, 1, -.5, .8)
    assert solution.residual_relative < 1e-12
    assert np.max(abs(solution.coefficients[:, [0, 1, 3]])) < 1e-12
    with pytest.raises(ValueError): solve_rayleigh_nodal([[0, 0, .1]], .08, 1, 0, 0)
    with pytest.raises(ValueError): solve_rayleigh_nodal([[0, 0, 0]], .2, 1, 0, 0)


def test_neumann_partial_sums_converge_to_dense_solution():
    positions = [[-2., 0., 0.], [2., 0., 0.]]
    solution = solve_rayleigh_nodal(positions, .08, 1., -.2, .8)
    rhs = solution.right_hand_side
    operator = np.eye(len(rhs), dtype=complex) - solution.system_matrix
    partial = rhs.copy(); term = rhs.copy(); errors = []
    exact = solution.coefficients.reshape(-1)
    for _ in range(6):
        errors.append(np.linalg.norm(partial - exact))
        term = operator @ term
        partial += term
    assert all(later < earlier for earlier, later in zip(errors, errors[1:]))


def test_three_particle_rotation_and_permutation_preserve_s10_up_to_relabeling():
    positions = np.array([[0., 0., 0.], [3., .4, 0.], [-.8, 2.6, 0.]])
    base = solve_rayleigh_nodal(positions, .08, 1., -.5, .8).coefficients[:, 2]
    angle = .63; rotation = np.array([[np.cos(angle), -np.sin(angle), 0.], [np.sin(angle), np.cos(angle), 0.], [0., 0., 1.]])
    rotated = solve_rayleigh_nodal(positions @ rotation.T, .08, 1., -.5, .8).coefficients[:, 2]
    order = np.array([2, 0, 1])
    permuted = solve_rayleigh_nodal(positions[order], .08, 1., -.5, .8).coefficients[:, 2]
    assert np.allclose(rotated, base, rtol=1e-12, atol=1e-14)
    assert np.allclose(permuted, base[order], rtol=1e-12, atol=1e-14)


def test_small_kd_rescattering_asymptotics():
    radius, separation, ka, f1 = 1.0, 2.0, 1e-3, .8
    k = ka / radius
    s1 = rayleigh_scattering_coefficients(ka, 0.0, f1)[1]
    coupling = s1 * (spherical_hankel1(0, k * separation) + spherical_hankel1(2, k * separation))
    expected_coupling = f1 / 2 * (radius / separation) ** 3
    assert np.isclose(coupling, expected_coupling, rtol=3e-6)
    single = solve_rayleigh_nodal([[0., 0., 0.]], k, radius, 0.0, f1).coefficients[0, 2]
    pair = solve_rayleigh_nodal([[-separation / 2, 0., 0.], [separation / 2, 0., 0.]], k, radius, 0.0, f1).coefficients[0, 2]
    expected_ratio = 1 / (1 - expected_coupling)
    assert np.isclose(pair / single, expected_ratio, rtol=2e-7)
