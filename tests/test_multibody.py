"""T05 independent scalar-oracle and trimer Model A/B/C validation."""

import pytest

import numpy as np
from scipy.special import spherical_jn, spherical_yn

from acoustic_ms import (
    angular_errors_degrees, compare_nodal_force_models, equilateral_trimer,
    linear_trimer, rms_relative_error, rms_vector_magnitude, scalene_trimer, symmetric_particle_errors,
)


def _hankel(ell, x, derivative=False):
    return spherical_jn(ell, x, derivative=derivative) + 1j * spherical_yn(ell, x, derivative=derivative)


def _scalar_oracle(positions, k, radius, energy, f1):
    positions = np.asarray(positions, dtype=float); ka = k * radius
    s1 = 1j * f1 * ka**3 / 6; a10 = np.sqrt(12 * np.pi); count = len(positions)
    matrix = np.eye(count, dtype=complex)
    for target in range(count):
        for source in range(count):
            if source != target:
                x = k * np.linalg.norm(positions[target] - positions[source])
                matrix[target, source] = -s1 * (_hankel(0, x) + _hankel(2, x))
    sigma = np.linalg.solve(matrix, np.full(count, s1 * a10, dtype=complex))
    forces = np.zeros((count, 2)); factor = 4 * np.pi * k * radius**3 * energy * np.sqrt(3 / (4 * np.pi))
    for target in range(count):
        for source in range(count):
            if target != source:
                displacement = positions[target] - positions[source]; distance = np.linalg.norm(displacement); x = k * distance
                q = _hankel(1, x, True) / x - _hankel(1, x) / x**2
                forces[target] += factor * np.real(np.conj(f1) * sigma[source] * q) * displacement[:2] / distance
    return sigma, forces


def _check_oracle(positions):
    result = compare_nodal_force_models(positions, .1, 1, 1, 0, .8)
    sigma, forces = _scalar_oracle(positions, .1, 1, 1, .8)
    assert np.allclose(result.global_result.solution.coefficients[:, 2], sigma, rtol=3e-12, atol=3e-13)
    assert np.max(abs(result.global_result.solution.coefficients[:, [0, 1, 3]])) < 3e-13
    assert np.allclose(result.model_c_forces_xy, forces, rtol=3e-12, atol=3e-13)
    assert np.allclose(result.model_c_forces_xy - result.model_a_forces_xy, result.two_body_correction_xy + result.irreducible_multibody_xy)
    return result


def test_n1_n2_reductions_and_metrics_basics():
    one = compare_nodal_force_models([[0, 0, 0]], .1, 1, 1, 0, .8)
    assert np.allclose(one.model_a_forces_xy, 0) and np.allclose(one.model_b_forces_xy, 0) and np.allclose(one.model_c_forces_xy, 0)
    two = compare_nodal_force_models([[-1.1, 0, 0], [1.1, 0, 0]], .1, 1, 1, 0, .8)
    assert np.allclose(two.model_b_forces_xy, two.model_c_forces_xy, rtol=2e-12, atol=2e-13)
    assert np.allclose(two.irreducible_multibody_xy, 0, atol=2e-13)
    vectors = np.array([[1., 0.], [0., 0.]])
    assert np.allclose(symmetric_particle_errors(vectors, vectors), 0)
    assert np.isclose(rms_relative_error(vectors, -vectors), 2)
    assert np.allclose(angular_errors_degrees(np.array([[1.,0.],[-1.,0.],[1.,0.]]), np.array([[1.,0.],[1.,0.],[0.,1.]])), [0.,180.,90.])
    assert np.isnan(angular_errors_degrees(vectors, vectors)[1])


def test_scalar_oracle_and_canonical_regressions():
    chain = _check_oracle(linear_trimer(2.1)); equilateral = _check_oracle(equilateral_trimer(2.1)); scalene = _check_oracle(scalene_trimer(2.2))
    assert np.allclose(chain.model_c_forces_xy[:, 0], [0.725204866335105, 0., -0.725204866335105], rtol=5e-12, atol=3e-13)
    assert np.allclose(equilateral.model_c_forces_xy, [[1.028117630785575,.593583990892652],[-1.028117630785576,.593583990892652],[0.,-1.187167981785305]], rtol=5e-12, atol=3e-13)
    assert np.allclose(scalene.model_c_forces_xy, [[.632645020072720,.330698096562575],[-.657831717956925,.154039050571838],[.019458368752579,-.493350707384327]], rtol=5e-12, atol=3e-13)
    assert np.isclose(rms_relative_error(chain.model_c_forces_xy, chain.model_a_forces_xy), .0831973786371517, rtol=5e-12)
    assert np.isclose(rms_relative_error(scalene.model_c_forces_xy, scalene.model_b_forces_xy), .0270655878714466, rtol=5e-12)


def test_permutation_translation_rotation_scaling_and_geometric_symmetries():
    positions = scalene_trimer(2.2); base = _check_oracle(positions)
    order = np.array([2, 0, 1]); permuted = compare_nodal_force_models(positions[order], .1, 1, 1, 0, .8)
    for name in ("model_a_forces_xy", "model_b_forces_xy", "model_c_forces_xy", "two_body_correction_xy", "irreducible_multibody_xy"):
        assert np.allclose(getattr(permuted, name), getattr(base, name)[order], rtol=2e-12, atol=2e-13)
    assert np.allclose(permuted.global_result.solution.coefficients, base.global_result.solution.coefficients[order], rtol=2e-12, atol=2e-13)
    shifted = compare_nodal_force_models(positions + [3, -2, 0], .1, 1, 1, 0, .8)
    angle=.43; rotation=np.array([[np.cos(angle),-np.sin(angle),0],[np.sin(angle),np.cos(angle),0],[0,0,1]])
    rotated=compare_nodal_force_models(positions @ rotation.T,.1,1,1,0,.8)
    for name in ("model_a_forces_xy", "model_b_forces_xy", "model_c_forces_xy", "two_body_correction_xy", "irreducible_multibody_xy"):
        assert np.allclose(getattr(shifted,name), getattr(base,name), rtol=2e-12, atol=2e-13)
        assert np.allclose(getattr(rotated,name), getattr(base,name) @ rotation[:2,:2].T, rtol=2e-12, atol=2e-13)
    scale=2.3; scaled=compare_nodal_force_models(positions*scale,.1/scale,scale,1,0,.8)
    assert np.allclose(scaled.model_c_forces_xy, scale**2*base.model_c_forces_xy, rtol=2e-12, atol=2e-13)
    chain=compare_nodal_force_models(linear_trimer(2.1),.1,1,1,0,.8)
    assert np.linalg.norm(chain.model_c_forces_xy[1]) < 2e-13 and np.allclose(chain.model_c_forces_xy[0],-chain.model_c_forces_xy[2],atol=2e-13)
    equilateral=compare_nodal_force_models(equilateral_trimer(2.1),.1,1,1,0,.8)
    assert np.allclose(np.linalg.norm(equilateral.model_c_forces_xy,axis=1),np.linalg.norm(equilateral.model_c_forces_xy[0]),rtol=2e-12)


def test_energy_contrast_f0_and_weak_coupling_limits():
    positions=scalene_trimer(2.2); base=compare_nodal_force_models(positions,.1,1,1,0,.8)
    assert np.allclose(compare_nodal_force_models(positions,.1,1,2,0,.8).model_c_forces_xy,2*base.model_c_forces_xy)
    assert np.allclose(compare_nodal_force_models(positions,.1,1,1,5,.8).model_c_forces_xy,base.model_c_forces_xy)
    assert np.allclose(compare_nodal_force_models(positions,.1,1,1,0,0).model_c_forces_xy,0)
    weak=compare_nodal_force_models(scalene_trimer(8),.05,1,1,0,.1)
    assert rms_relative_error(weak.model_b_forces_xy,weak.model_a_forces_xy) < .01
    assert rms_relative_error(weak.model_c_forces_xy,weak.model_b_forces_xy) < .01


FIELDS = (
    "model_a_forces_xy",
    "model_b_forces_xy",
    "model_c_forces_xy",
    "two_body_correction_xy",
    "irreducible_multibody_xy",
)


def test_metrics_share_a_global_scale_without_an_absolute_floor():
    reference = np.array([[1.0, 0.0], [-9.2e-17, 8.1e-17], [-1.0, 0.0]])
    model = np.array([[0.9, 0.0], [0.0, 8.0e-17], [-0.9, 0.0]])
    symmetric = symmetric_particle_errors(reference, model)
    angles = angular_errors_degrees(reference, model)
    assert np.isclose(symmetric[0], 2 / 19)
    assert symmetric[1] == 0.0
    assert np.isnan(angles[1])
    assert np.isclose(angles[0], 0.0)

    zeros = np.zeros((1, 2))
    assert symmetric_particle_errors(zeros, zeros)[0] == 0.0
    assert np.isnan(angular_errors_degrees(zeros, zeros)[0])
    assert np.isclose(symmetric_particle_errors([[1.0, 0.0]], zeros)[0], 2.0)
    assert np.isnan(angular_errors_degrees([[1.0, 0.0]], zeros)[0])
    assert np.isclose(angular_errors_degrees([[1.0, 0.0]], [[0.0, 1.0]])[0], 90.0)

    tiny_reference = np.array([[1e-20, 0.0]])
    tiny_model = np.array([[-1e-20, 0.0]])
    assert symmetric_particle_errors(tiny_reference, tiny_model)[0] == 2.0
    assert np.isclose(angular_errors_degrees(tiny_reference, tiny_model)[0], 180.0)
    assert np.isclose(rms_vector_magnitude([[3.0, 4.0], [0.0, 0.0]]), np.sqrt(25 / 2))


def test_all_fields_obey_permutation_scaling_energy_and_contrast_invariances():
    positions = scalene_trimer(2.2)
    base = compare_nodal_force_models(positions, .1, 1, 1, 0, .8)
    for order in (np.array([2, 0, 1]), np.array([1, 2, 0])):
        permuted = compare_nodal_force_models(positions[order], .1, 1, 1, 0, .8)
        for name in FIELDS:
            assert np.allclose(getattr(permuted, name), getattr(base, name)[order], rtol=2e-12, atol=2e-13)
        assert np.allclose(permuted.global_result.solution.coefficients, base.global_result.solution.coefficients[order], rtol=2e-12, atol=2e-13)

    scale = 2.3
    scaled = compare_nodal_force_models(positions * scale, .1 / scale, scale, 1, 0, .8)
    doubled_energy = compare_nodal_force_models(positions, .1, 1, 2, 0, .8)
    f0_changed = compare_nodal_force_models(positions, .1, 1, 1, 5, .8)
    zero_energy = compare_nodal_force_models(positions, .1, 1, 0, 0, .8)
    zero_f1 = compare_nodal_force_models(positions, .1, 1, 1, 0, 0)
    for name in FIELDS:
        field = getattr(base, name)
        assert np.allclose(getattr(scaled, name), scale**2 * field, rtol=2e-12, atol=2e-13)
        assert np.allclose(getattr(doubled_energy, name), 2 * field, rtol=2e-12, atol=2e-13)
        assert np.allclose(getattr(f0_changed, name), field, rtol=2e-12, atol=2e-13)
        assert np.allclose(getattr(zero_energy, name), 0, atol=2e-13)
        assert np.allclose(getattr(zero_f1, name), 0, atol=2e-13)
    with pytest.raises(ValueError):
        compare_nodal_force_models(positions, .1, 1, 1, 0, .8, lmax=2)


def test_chain_and_equilateral_symmetries_hold_for_every_force_field():
    chain = compare_nodal_force_models(linear_trimer(2.1), .1, 1, 1, 0, .8)
    for name in FIELDS:
        field = getattr(chain, name)
        scale = max(np.linalg.norm(field, axis=1).max(), 1e-300)
        tolerance = 3e-12 * scale
        assert np.linalg.norm(field[1]) <= tolerance
        assert np.allclose(field[0], -field[2], atol=tolerance, rtol=0)
        assert np.max(np.abs(field[:, 1])) <= tolerance
        assert np.linalg.norm(field.sum(axis=0)) <= tolerance

    positions = equilateral_trimer(2.1)
    equilateral = compare_nodal_force_models(positions, .1, 1, 1, 0, .8)
    angle = 2 * np.pi / 3
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0.0], [np.sin(angle), np.cos(angle), 0.0], [0.0, 0.0, 1.0]])
    rotated = compare_nodal_force_models(positions @ rotation.T, .1, 1, 1, 0, .8)
    order = np.array([1, 2, 0])
    for name in FIELDS:
        field = getattr(equilateral, name)
        scale = max(np.linalg.norm(field, axis=1).max(), 1e-300)
        tolerance = 3e-12 * scale
        assert np.allclose(np.linalg.norm(field, axis=1), np.linalg.norm(field[0]), rtol=3e-12, atol=tolerance)
        assert np.max(np.abs(positions[:, 0] * field[:, 1] - positions[:, 1] * field[:, 0])) <= tolerance
        assert np.linalg.norm(field.sum(axis=0)) <= tolerance
        assert np.allclose(getattr(rotated, name), field @ rotation[:2, :2].T, rtol=3e-12, atol=tolerance)
        assert np.allclose(getattr(rotated, name), field[order], rtol=3e-12, atol=tolerance)
