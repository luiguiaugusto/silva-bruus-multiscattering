"""Geometry validation for the thirteen deterministic T08 families."""

import numpy as np
import pytest

from acoustic_ms.cluster_families import (
    cluster_family,
    compact_cluster,
    enumerate_transferability_configurations,
    irregular_cluster,
    linear_cluster,
)


def _minimum_distance(positions):
    distances = np.linalg.norm(positions[:, None] - positions[None, :], axis=2)
    return np.min(distances[np.triu_indices(len(positions), 1)])


@pytest.mark.parametrize(
    "particle_count,family",
    [(2, "pair")] + [(n, family) for n in (3, 4, 6, 10) for family in ("linear", "compact", "irregular")],
)
def test_family_geometry_contract(particle_count, family):
    positions = cluster_family(particle_count, family, 2.1)
    assert positions.shape == (particle_count, 3)
    np.testing.assert_allclose(positions.mean(axis=0), 0, atol=3e-16)
    assert np.count_nonzero(positions[:, 2]) == 0
    assert _minimum_distance(positions) == pytest.approx(2.1, rel=3e-15)
    assert _minimum_distance(positions) >= 2.0
    assert np.array_equal(positions, cluster_family(particle_count, family, 2.1))


def test_exact_configuration_count_and_order():
    configurations = enumerate_transferability_configurations()
    assert len(configurations) == 312
    assert len({item.case_id for item in configurations}) == 312
    assert configurations[0].case_id == "n2_pair_f0.1_d2.1"
    assert configurations[-1].case_id == "n10_irregular_f1.0_d10.0"
    assert sum(item.split == "calibration" for item in configurations) == 168
    assert sum(item.split == "holdout" for item in configurations) == 144


def test_configuration_order_is_explicit_and_deterministic():
    geometries = [(2, "pair")]
    geometries.extend(
        (n, family)
        for n in (3, 4, 6, 10)
        for family in ("linear", "compact", "irregular")
    )
    expected = [
        f"n{n}_{family}_f{f1:.1f}_d{distance:.1f}"
        for n, family in geometries
        for f1 in (0.1, 0.4, 0.8, 1.0)
        for distance in (2.1, 2.5, 3.0, 4.0, 6.0, 10.0)
    ]
    first = [item.case_id for item in enumerate_transferability_configurations()]
    second = [item.case_id for item in enumerate_transferability_configurations()]
    assert first == expected == second


@pytest.mark.parametrize("builder", [linear_cluster, compact_cluster, irregular_cluster])
def test_geometry_distance_validation(builder):
    for distance in (0, -1, np.nan, np.inf):
        with pytest.raises(ValueError):
            builder(6, distance)


def test_family_validation():
    with pytest.raises(ValueError): cluster_family(5, "linear", 2.1)
    with pytest.raises(ValueError): cluster_family(4, "unknown", 2.1)
    with pytest.raises(ValueError): compact_cluster(5, 2.1)
    with pytest.raises(ValueError): irregular_cluster(5, 2.1)
