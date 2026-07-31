"""Independent unit tests for the reusable T12.1 diagnostics."""

import numpy as np
import pytest

from acoustic_ms import (
    convergence_tail_diagnostics,
    fit_log_linear,
    leave_group_out_folds,
    mechanism_diagnostics,
    out_of_fold_metrics,
    vector_field_amplitude_ratio,
    vector_field_cosine,
    vector_field_inner_product,
    vector_field_projection,
)


def test_inner_product_and_analytic_cosines():
    x = np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    y = np.array([[0.0, 2.0, 0.0], [0.0, 2.0, 0.0]])
    assert vector_field_inner_product(x, x) == 1.0
    assert vector_field_cosine(x, x).value == pytest.approx(1.0)
    assert vector_field_cosine(x, -x).value == pytest.approx(-1.0)
    assert vector_field_cosine(x, y).value == pytest.approx(0.0)


def test_null_cosine_and_projection_are_explicitly_inapplicable():
    zero = np.zeros((2, 3))
    x = np.ones((2, 3))
    assert not vector_field_cosine(zero, x).applicable
    assert not vector_field_projection(x, zero).applicable


def test_amplitude_ratio_has_no_absolute_floor():
    numerator = np.array([[2.0e-20, 0.0, 0.0]])
    denominator = np.array([[1.0e-20, 0.0, 0.0]])
    diagnostic = vector_field_amplitude_ratio(numerator, denominator)
    assert diagnostic.applicable
    assert diagnostic.value == pytest.approx(2.0)


def test_mechanism_identity_and_signed_projection_sum():
    a = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    c_d = np.array([[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]])
    c_m = np.array([[-0.02, 0.01, 0.0], [0.02, -0.01, 0.0]])
    c_s = np.array([[0.01, 0.02, 0.0], [-0.01, -0.02, 0.0]])
    d = a + c_d
    external = d + c_m
    interaction = external + c_s
    result = mechanism_diagnostics(a, d, interaction, external, c_s)
    assert result.closure_rms < 1.0e-16
    assert result.projection_sum.applicable
    assert result.projection_sum.value == pytest.approx(1.0)


@pytest.mark.parametrize(
    "bad",
    [np.ones(3), np.ones((2, 2)), np.array([[np.nan, 0.0, 0.0]])],
)
def test_vector_diagnostics_validate_inputs(bad):
    with pytest.raises(ValueError):
        vector_field_inner_product(bad, bad)


def test_log_linear_fit_recovers_synthetic_power_law():
    x = np.geomspace(0.01, 2.0, 30)
    y = 3.25 * x**1.7
    result = fit_log_linear(x, y)
    assert result.prefactor == pytest.approx(3.25, rel=2.0e-15)
    assert result.coefficient == pytest.approx(1.7, rel=2.0e-15)


@pytest.mark.parametrize(
    ("x", "y"),
    [
        ([1.0], [2.0]),
        ([1.0, 2.0], [1.0]),
        ([0.0, 1.0], [1.0, 2.0]),
        ([1.0, np.inf], [1.0, 2.0]),
        ([1.0, 1.0], [1.0, 2.0]),
    ],
)
def test_log_linear_fit_rejects_invalid_inputs(x, y):
    with pytest.raises(ValueError):
        fit_log_linear(x, y)


def test_leave_group_out_folds_are_deterministic_exhaustive_and_leak_free():
    groups = ["n2_pair", "n3_compact", "n2_pair", "n4_linear", "n3_compact"]
    first = leave_group_out_folds(groups)
    second = leave_group_out_folds(groups)
    assert [fold[0] for fold in first] == ["n2_pair", "n3_compact", "n4_linear"]
    assert [fold[0] for fold in second] == [fold[0] for fold in first]
    covered = []
    for group, train, test in first:
        assert not set(train).intersection(test)
        assert all(groups[index] == group for index in test)
        assert all(groups[index] != group for index in train)
        covered.extend(test.tolist())
    assert sorted(covered) == list(range(len(groups)))


def test_out_of_fold_metrics_match_direct_calculation():
    observed = np.array([1.0, 2.0, 4.0, 8.0])
    predicted = np.array([2.0, 1.0, 8.0, 4.0])
    result = out_of_fold_metrics(observed, predicted)
    assert result.rmse_log == pytest.approx(np.log(2.0))
    assert result.median_factor == pytest.approx(2.0)
    assert result.fraction_within_factor_2 == 1.0


def test_convergence_requires_two_successive_applicable_changes():
    result = convergence_tail_diagnostics(
        [14, 15, 16, 17],
        [2.0e-5, 9.0e-6, 8.0e-6, 7.0e-6],
        [True, True, True, True],
    )
    assert result.confirmation_order == 16
    assert result.classification == "directly_confirmed"


def test_convergence_unconfirmed_at_maximum_order_and_oscillation():
    result = convergence_tail_diagnostics(
        [17, 18, 19, 20, 21],
        [2.0e-5, 3.0e-5, 2.5e-5, 3.5e-5, 3.0e-5],
        [True] * 5,
    )
    assert result.confirmation_order is None
    assert result.classification == "unconfirmed_at_21"
    assert result.oscillatory
