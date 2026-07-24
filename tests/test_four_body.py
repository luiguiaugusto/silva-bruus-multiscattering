"""T06 independent-oracle and connected quartet-expansion validation."""

from itertools import combinations

import numpy as np
import pytest
from scipy.special import spherical_jn, spherical_yn

from acoustic_ms import (
    decompose_nodal_quartet,
    irregular_quartet,
    linear_quartet,
    rms_relative_error,
    rms_vector_magnitude,
    square_quartet,
)


FIELDS = (
    "model_a_forces_xy",
    "model_b_forces_xy",
    "model_c_forces_xy",
    "two_body_correction_xy",
    "collective_correction_xy",
    "irreducible_three_body_sum_xy",
    "up_to_three_body_forces_xy",
    "irreducible_four_body_xy",
)
TRIPLETS = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))


def _hankel(ell, x, derivative=False):
    return spherical_jn(ell, x, derivative=derivative) + 1j * spherical_yn(
        ell, x, derivative=derivative
    )


def _scalar_oracle(positions, k, radius, energy_density, f1):
    positions = np.asarray(positions, dtype=float)
    count = len(positions)
    s1 = 1j * f1 * (k * radius) ** 3 / 6
    matrix = np.eye(count, dtype=complex)
    for target in range(count):
        for source in range(count):
            if target != source:
                x = k * np.linalg.norm(positions[target] - positions[source])
                matrix[target, source] = -s1 * (_hankel(0, x) + _hankel(2, x))
    sigma = np.linalg.solve(
        matrix, np.full(count, s1 * np.sqrt(12 * np.pi), dtype=complex)
    )
    forces = np.zeros((count, 2))
    prefactor = 4 * np.pi * k * radius**3 * energy_density * np.sqrt(3 / (4 * np.pi))
    for target in range(count):
        for source in range(count):
            if target != source:
                displacement = positions[target] - positions[source]
                distance = np.linalg.norm(displacement)
                x = k * distance
                q = _hankel(1, x, True) / x - _hankel(1, x) / x**2
                forces[target] += (
                    prefactor
                    * np.real(np.conj(f1) * sigma[source] * q)
                    * displacement[:2]
                    / distance
                )
    return sigma, forces


def _oracle_expansion(positions, k=.1, radius=1, energy_density=1, f1=.8):
    positions = np.asarray(positions)
    full_sigma, full_force = _scalar_oracle(
        positions, k, radius, energy_density, f1
    )
    pair_forces = {}
    model_b = np.zeros((4, 2))
    for pair in combinations(range(4), 2):
        _, force = _scalar_oracle(
            positions[list(pair)], k, radius, energy_density, f1
        )
        pair_forces[pair] = force
        model_b[list(pair)] += force
    by_triplet = np.zeros((4, 4, 2))
    triplet_sigmas = []
    triplet_forces = []
    for row, triplet in enumerate(TRIPLETS):
        sigma, force = _scalar_oracle(
            positions[list(triplet)], k, radius, energy_density, f1
        )
        triplet_sigmas.append(sigma)
        triplet_forces.append(force)
        pair_sum = np.zeros((3, 2))
        for local_first, local_second in combinations(range(3), 2):
            global_pair = tuple(sorted((triplet[local_first], triplet[local_second])))
            pair_force = pair_forces[global_pair]
            if global_pair != (triplet[local_first], triplet[local_second]):
                pair_force = pair_force[::-1]
            pair_sum[local_first] += pair_force[0]
            pair_sum[local_second] += pair_force[1]
        by_triplet[row, list(triplet)] = force - pair_sum
    three_sum = by_triplet.sum(axis=0)
    four_body = full_force - model_b - three_sum
    return full_sigma, full_force, model_b, tuple(triplet_sigmas), tuple(triplet_forces), by_triplet, three_sum, four_body


def _decompose(positions, **kwargs):
    parameters = dict(k=.1, radius=1, energy_density=1, f0=0, f1=.8)
    parameters.update(kwargs)
    return decompose_nodal_quartet(positions, **parameters)


def test_quartet_geometry_definitions_and_invalid_distances():
    for geometry in (linear_quartet, square_quartet, irregular_quartet):
        positions = geometry(2.1)
        assert positions.shape == (4, 3)
        assert np.allclose(positions.mean(axis=0), 0, atol=2e-16)
        assert np.all(positions[:, 2] == 0)
        distances = [np.linalg.norm(positions[i] - positions[j]) for i, j in combinations(range(4), 2)]
        assert np.isclose(min(distances), 2.1)
        for invalid in (0, -1, np.nan, np.inf):
            with pytest.raises(ValueError):
                geometry(invalid)


def test_quartet_input_validation():
    good = irregular_quartet(2.1)
    for positions in (good[:3], np.zeros((4, 2)), np.vstack([good, good[0]])):
        with pytest.raises(ValueError):
            _decompose(positions)
    invalid = good.copy(); invalid[0, 0] = np.nan
    with pytest.raises(ValueError):
        _decompose(invalid)
    off_plane = good.copy(); off_plane[0, 2] = 1e-3
    with pytest.raises(ValueError):
        _decompose(off_plane)
    overlap = good.copy(); overlap[1] = overlap[0]
    with pytest.raises(ValueError):
        _decompose(overlap)
    with pytest.raises(ValueError):
        _decompose(linear_quartet(1.9))
    for lmax in (0, 2, 3, 5):
        with pytest.raises(ValueError):
            _decompose(good, lmax=lmax)


def test_triplet_order_embedding_and_vector_identities():
    result = _decompose(irregular_quartet(2.1))
    assert result.triplet_indices == TRIPLETS
    assert result.irreducible_three_body_by_triplet_xy.shape == (4, 4, 2)
    for row, (indices, comparison) in enumerate(
        zip(result.triplet_indices, result.triplet_comparisons, strict=True)
    ):
        missing = next(index for index in range(4) if index not in indices)
        assert np.array_equal(result.irreducible_three_body_by_triplet_xy[row, missing], np.zeros(2))
        assert np.allclose(
            result.irreducible_three_body_by_triplet_xy[row, list(indices)],
            comparison.irreducible_multibody_xy,
        )
    assert np.allclose(
        result.model_c_forces_xy,
        result.model_b_forces_xy
        + result.irreducible_three_body_sum_xy
        + result.irreducible_four_body_xy,
        atol=2e-16,
    )
    assert np.allclose(
        result.collective_correction_xy,
        result.irreducible_three_body_sum_xy + result.irreducible_four_body_xy,
        atol=2e-16,
    )
    assert np.allclose(
        result.model_c_forces_xy - result.model_a_forces_xy,
        result.two_body_correction_xy
        + result.irreducible_three_body_sum_xy
        + result.irreducible_four_body_xy,
        atol=2e-16,
    )
    closed = result.model_c_forces_xy.copy()
    for indices, comparison in zip(
        result.triplet_indices, result.triplet_comparisons, strict=True
    ):
        closed[list(indices)] -= comparison.model_c_forces_xy
    closed += result.model_b_forces_xy
    assert np.allclose(closed, result.irreducible_four_body_xy, atol=4e-16)


@pytest.mark.parametrize(
    "geometry", (linear_quartet, square_quartet, irregular_quartet)
)
def test_independent_scalar_oracle_for_quartet_triplets_and_connected_terms(geometry):
    positions = geometry(2.1)
    result = _decompose(positions)
    oracle = _oracle_expansion(positions)
    sigma, full_force, model_b, triplet_sigmas, triplet_forces, by_triplet, three_sum, four_body = oracle
    coefficients = result.full_comparison.global_result.solution.coefficients
    assert coefficients.shape == (4, 4)
    assert np.allclose(coefficients[:, 2], sigma, rtol=3e-12, atol=3e-13)
    assert np.max(np.abs(coefficients[:, [0, 1, 3]])) < 3e-13
    assert np.allclose(result.model_c_forces_xy, full_force, rtol=3e-12, atol=3e-13)
    assert np.allclose(result.model_b_forces_xy, model_b, rtol=3e-12, atol=3e-13)
    for comparison, expected_sigma, expected_force in zip(
        result.triplet_comparisons, triplet_sigmas, triplet_forces, strict=True
    ):
        actual = comparison.global_result.solution.coefficients
        assert np.allclose(actual[:, 2], expected_sigma, rtol=3e-12, atol=3e-13)
        assert np.max(np.abs(actual[:, [0, 1, 3]])) < 3e-13
        assert np.allclose(comparison.model_c_forces_xy, expected_force, rtol=3e-12, atol=3e-13)
    assert np.allclose(result.irreducible_three_body_by_triplet_xy, by_triplet, rtol=3e-12, atol=3e-13)
    assert np.allclose(result.irreducible_three_body_sum_xy, three_sum, rtol=3e-12, atol=3e-13)
    assert np.allclose(result.irreducible_four_body_xy, four_body, rtol=3e-12, atol=3e-13)


def test_canonical_regressions():
    expected = {
        linear_quartet: (
            [[.741001670038274,0],[.069982023490876,0],[-.069982023490876,0],[-.741001670038274,0]],
            [[.037276624590095,0],[.029982785116269,0],[-.029982785116269,0],[-.037276624590095,0]],
            [[.001539740231837,0],[-.000193609810481,0],[.000193609810481,0],[-.001539740231837,0]],
            (.09977255014752409,.06573916615217347,.002085002107029354,.03382686908479156,.001097334210769815),
        ),
        square_quartet: (
            [[.821848100064453,.821848100064453],[-.821848100064453,.821848100064453],[-.821848100064453,-.821848100064453],[.821848100064453,-.821848100064453]],
            [[.052751430705498,.052751430705498],[-.052751430705498,.052751430705498],[-.052751430705498,-.052751430705498],[.052751430705498,-.052751430705498]],
            [[.002227631154520,.002227631154520],[-.002227631154520,.002227631154520],[-.002227631154520,-.002227631154520],[.002227631154520,-.002227631154520]],
            (.10419668455946075,.06689686555910586,.002710514454368684,.07460178873829999,.003150346190687219),
        ),
        irregular_quartet: (
            [[.808294694345220,.322600068278836],[-.645766967125507,.494718667472942],[.473022318883894,-.482298140039304],[-.630849957580654,-.344985947033440]],
            [[.036631993918007,.016467994763673],[-.022090087132087,.023128523747475],[.019230383145640,-.027601910171466],[-.029203334725065,-.021760503592003]],
            [[.001207441363755,.000463073332837],[-.000687638003254,.000451324799227],[.000720308361924,-.000617437338051],[-.001108978405968,-.000496416862659]],
            (.07877328473634558,.04751654171343999,.001405346638645700,.03568622194513429,.001086879416669834),
        ),
    }
    for geometry, (force, three, four, metrics) in expected.items():
        result = _decompose(geometry(2.1))
        assert np.allclose(result.model_c_forces_xy, force, rtol=5e-12, atol=3e-13)
        assert np.allclose(result.irreducible_three_body_sum_xy, three, rtol=5e-12, atol=3e-13)
        assert np.allclose(result.irreducible_four_body_xy, four, rtol=5e-12, atol=3e-13)
        actual = (
            rms_relative_error(result.model_c_forces_xy, result.model_a_forces_xy),
            rms_relative_error(result.model_c_forces_xy, result.model_b_forces_xy),
            rms_relative_error(result.model_c_forces_xy, result.up_to_three_body_forces_xy),
            rms_vector_magnitude(result.irreducible_three_body_sum_xy),
            rms_vector_magnitude(result.irreducible_four_body_xy),
        )
        assert np.allclose(actual, metrics, rtol=5e-12, atol=3e-13)


def test_permutation_maps_all_fields_coefficients_and_identified_triplets():
    positions = irregular_quartet(2.1)
    base = _decompose(positions)
    order = np.array([3, 1, 0, 2])
    permuted = _decompose(positions[order])
    for name in FIELDS:
        assert np.allclose(getattr(permuted, name), getattr(base, name)[order], rtol=3e-12, atol=3e-13)
    assert np.allclose(
        permuted.full_comparison.global_result.solution.coefficients,
        base.full_comparison.global_result.solution.coefficients[order],
        rtol=3e-12,
        atol=3e-13,
    )
    for row, indices in enumerate(permuted.triplet_indices):
        original_set = tuple(sorted(order[list(indices)]))
        original_row = base.triplet_indices.index(original_set)
        assert np.allclose(
            permuted.irreducible_three_body_by_triplet_xy[row],
            base.irreducible_three_body_by_triplet_xy[original_row, order],
            rtol=3e-12,
            atol=3e-13,
        )


def test_translation_rotation_scaling_energy_and_contrasts():
    positions = irregular_quartet(2.1)
    base = _decompose(positions)
    shifted = _decompose(positions + [2.3, -1.7, 0])
    angle = .37
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    rotated = _decompose(positions @ rotation.T)
    scale = 2.4
    scaled = _decompose(positions * scale, k=.1 / scale, radius=scale)
    doubled = _decompose(positions, energy_density=2)
    changed_f0 = _decompose(positions, f0=5)
    zero_energy = _decompose(positions, energy_density=0)
    zero_f1 = _decompose(positions, f1=0)
    for name in FIELDS:
        field = getattr(base, name)
        assert np.allclose(getattr(shifted, name), field, rtol=3e-12, atol=3e-13)
        assert np.allclose(getattr(rotated, name), field @ rotation[:2, :2].T, rtol=3e-12, atol=3e-13)
        assert np.allclose(getattr(scaled, name), scale**2 * field, rtol=3e-12, atol=3e-13)
        assert np.allclose(getattr(doubled, name), 2 * field, rtol=3e-12, atol=3e-13)
        assert np.allclose(getattr(changed_f0, name), field, rtol=3e-12, atol=3e-13)
        assert np.allclose(getattr(zero_energy, name), 0, atol=3e-13)
        assert np.allclose(getattr(zero_f1, name), 0, atol=3e-13)


def test_chain_and_square_symmetries_for_all_vector_fields():
    chain = _decompose(linear_quartet(2.1))
    for name in FIELDS:
        field = getattr(chain, name)
        tolerance = 4e-12 * max(np.linalg.norm(field, axis=1).max(), 1e-300)
        assert np.allclose(field[0], -field[3], atol=tolerance, rtol=0)
        assert np.allclose(field[1], -field[2], atol=tolerance, rtol=0)
        assert np.max(np.abs(field[:, 1])) <= tolerance
        assert np.linalg.norm(field.sum(axis=0)) <= tolerance

    positions = square_quartet(2.1)
    square = _decompose(positions)
    rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    rotated = _decompose(positions @ rotation.T)
    order = np.array([1, 2, 3, 0])
    for name in FIELDS:
        field = getattr(square, name)
        tolerance = 4e-12 * max(np.linalg.norm(field, axis=1).max(), 1e-300)
        assert np.allclose(np.linalg.norm(field, axis=1), np.linalg.norm(field[0]), rtol=4e-12, atol=tolerance)
        cross = positions[:, 0] * field[:, 1] - positions[:, 1] * field[:, 0]
        assert np.max(np.abs(cross)) <= tolerance
        assert np.linalg.norm(field.sum(axis=0)) <= tolerance
        assert np.allclose(getattr(rotated, name), field @ rotation[:2, :2].T, rtol=4e-12, atol=tolerance)
        assert np.allclose(getattr(rotated, name), field[order], rtol=4e-12, atol=tolerance)


def test_weak_coupling_limit():
    result = _decompose(irregular_quartet(8), k=.05, f1=.1)
    assert np.isclose(
        rms_relative_error(result.model_c_forces_xy, result.model_b_forces_xy),
        1.13e-4,
        rtol=.02,
    )
    assert np.isclose(
        rms_relative_error(result.model_c_forces_xy, result.up_to_three_body_forces_xy),
        7.75e-9,
        rtol=.03,
    )
