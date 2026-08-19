"""Unit tests for the P1.2 independently converged Model-E pair baseline."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from acoustic_ms import (
    evaluate_model_e_numerical_diagnostics,
    solve_model_be_nodal,
)


POSITIONS = np.array(
    [
        [0.0, 0.0, 0.0],
        [3.0, 1.0, 0.0],
        [7.0, -2.0, 0.0],
    ]
)


def _key(positions: object) -> tuple[tuple[float, ...], ...]:
    return tuple(tuple(float(value) for value in row) for row in positions)


def _oriented_forces(positions: np.ndarray, order: int, *, varying: bool) -> np.ndarray:
    displacement = positions[1] - positions[0]
    first = displacement.copy()
    second = np.array(
        [
            positions[0, 0] + positions[1, 0],
            positions[0, 1] - positions[1, 1],
            0.0,
        ]
    )
    scale = float(order) if varying else 1.0
    return scale * np.vstack((first, second))


def _fake_result(
    positions: np.ndarray,
    order: int,
    *,
    varying: bool = False,
    condition: float = 1.25,
) -> SimpleNamespace:
    interaction = _oriented_forces(positions, order, varying=varying)
    zero = np.zeros((2, 3), dtype=float)
    active_modes = ((1, -1), (1, 0))
    dimension = 2 * len(active_modes)
    solution = SimpleNamespace(
        balanced_condition_number=condition,
        balanced_backward_error=1.0e-15,
        effective_incident_closure_error=2.0e-15,
        scattering_closure_error=3.0e-15,
        scattered_coefficients=np.zeros((2, (order + 1) ** 2), dtype=complex),
        active_modes=active_modes,
        modes=tuple(range((order + 1) ** 2)),
        balanced_system_matrix=np.eye(dimension),
        production_solver="balanced_sqrt",
    )
    return SimpleNamespace(
        solution=solution,
        total_forces_xyz=interaction.copy(),
        external_forces_xyz=zero.copy(),
        interaction_forces_xyz=interaction.copy(),
        external_scattered_forces_xyz=interaction.copy(),
        scattered_scattered_forces_xyz=zero.copy(),
        decomposition_residual=4.0e-15,
        lmax=order,
    )


class FakePairSolver:
    def __init__(
        self,
        *,
        fail_key: tuple[tuple[float, ...], ...] | None = None,
        fail_lmax: int = 3,
        bad_diagnostics_key: tuple[tuple[float, ...], ...] | None = None,
        varying_key: tuple[tuple[float, ...], ...] | None = None,
    ) -> None:
        self.fail_key = fail_key
        self.fail_lmax = fail_lmax
        self.bad_diagnostics_key = bad_diagnostics_key
        self.varying_key = varying_key
        self.calls: list[tuple[tuple[tuple[float, ...], ...], int]] = []

    def __call__(
        self,
        positions: object,
        k: float,
        radius: float,
        energy: float,
        f0: float,
        f1: float,
        order: int,
    ) -> SimpleNamespace:
        del k, radius, energy, f0, f1
        values = np.asarray(positions, dtype=float)
        pair_key = _key(values)
        self.calls.append((pair_key, order))
        if pair_key == self.fail_key and order == self.fail_lmax:
            raise RuntimeError("injected pair failure")
        condition = 11.0 if pair_key == self.bad_diagnostics_key else 1.25
        return _fake_result(
            values,
            order,
            varying=pair_key == self.varying_key,
            condition=condition,
        )


class SettlingPairSolver(FakePairSolver):
    def __init__(self, settling_key: tuple[tuple[float, ...], ...]) -> None:
        super().__init__()
        self.settling_key = settling_key

    def __call__(
        self,
        positions: object,
        k: float,
        radius: float,
        energy: float,
        f0: float,
        f1: float,
        order: int,
    ) -> SimpleNamespace:
        result = super().__call__(positions, k, radius, energy, f0, f1, order)
        if _key(positions) == self.settling_key:
            scale = float(min(order, 4))
            for attribute in (
                "total_forces_xyz",
                "interaction_forces_xyz",
                "external_scattered_forces_xyz",
            ):
                setattr(result, attribute, scale * getattr(result, attribute))
        return result


class NonMonotonicPairSolver(FakePairSolver):
    def __call__(
        self,
        positions: object,
        k: float,
        radius: float,
        energy: float,
        f0: float,
        f1: float,
        order: int,
    ) -> SimpleNamespace:
        result = super().__call__(positions, k, radius, energy, f0, f1, order)
        scale = 1.0 if order <= 4 else 2.0
        for attribute in (
            "total_forces_xyz",
            "interaction_forces_xyz",
            "external_scattered_forces_xyz",
        ):
            setattr(result, attribute, scale * getattr(result, attribute))
        return result


def _interaction_convergence(record):
    return next(
        channel
        for channel in record.convergence
        if channel.channel == "interaction"
    )


def _solve(solver: FakePairSolver, **kwargs):
    return solve_model_be_nodal(
        POSITIONS,
        0.1,
        1.0,
        1.0,
        0.0,
        0.8,
        lmax_max=kwargs.pop("lmax_max", 6),
        solver=solver,
        **kwargs,
    )


def test_vector_sum_orientation_pair_order_and_applicability() -> None:
    solver = FakePairSolver()
    result = _solve(solver)

    assert result.eligible
    assert result.failure_stage is None
    assert result.pair_count == 3
    assert tuple(record.particle_indices for record in result.pair_ledger) == (
        (0, 1),
        (0, 2),
        (1, 2),
    )
    np.testing.assert_array_equal(
        result.forces_xyz,
        np.array(
            [
                [10.0, -1.0, 0.0],
                [7.0, -4.0, 0.0],
                [17.0, 5.0, 0.0],
            ]
        ),
    )
    assert result.forces_xyz is not None
    assert not result.forces_xyz.flags.writeable

    expected_pair_positions = (
        POSITIONS[[0, 1]],
        POSITIONS[[0, 2]],
        POSITIONS[[1, 2]],
    )
    for record, expected in zip(result.pair_ledger, expected_pair_positions):
        np.testing.assert_array_equal(record.positions_xyz, expected)
        np.testing.assert_array_equal(
            record.interaction_forces_xyz,
            _oriented_forces(expected, 5, varying=False),
        )
        assert record.attempted_lmax == (2, 3, 4, 5)
        assert record.evaluated_lmax == (2, 3, 4, 5)
        assert record.final_lmax == 5
        by_channel = {item.channel: item for item in record.convergence}
        assert by_channel["interaction"].applicable
        assert by_channel["interaction"].confirmed
        assert by_channel["interaction"].confirmation_lmax == 4
        assert not by_channel["scattered_scattered"].applicable
        assert not by_channel["scattered_scattered"].confirmed

    expected_calls = tuple(
        (_key(pair), order)
        for pair in expected_pair_positions
        for order in (2, 3, 4, 5)
    )
    assert tuple(solver.calls) == expected_calls


def test_each_pair_converges_to_its_own_final_order() -> None:
    settling_pair = _key(POSITIONS[[0, 2]])
    result = _solve(SettlingPairSolver(settling_pair), lmax_max=7)

    assert result.eligible
    assert tuple(record.final_lmax for record in result.pair_ledger) == (5, 6, 5)
    assert tuple(
        next(
            channel.confirmation_lmax
            for channel in record.convergence
            if channel.channel == "interaction"
        )
        for record in result.pair_ledger
    ) == (4, 6, 4)


def test_reopened_final_window_does_not_authorize_stop() -> None:
    solver = NonMonotonicPairSolver()
    result = solve_model_be_nodal(
        POSITIONS[:2],
        0.1,
        1.0,
        1.0,
        0.0,
        0.8,
        lmax_max=6,
        solver=solver,
    )

    record = result.pair_ledger[0]
    interaction = _interaction_convergence(record)
    assert record.final_lmax == 6
    assert record.attempted_lmax == (2, 3, 4, 5, 6)
    assert interaction.confirmation_lmax == 4
    assert all(step.applicable for step in interaction.history[-2:])
    assert interaction.history[-2].successive_change > 1.0e-5
    assert interaction.history[-1].successive_change <= 1.0e-5
    assert not interaction.confirmed
    assert any(
        channel.applicable and not channel.confirmed
        for channel in record.convergence
    )
    assert not result.eligible
    assert result.forces_xyz is None


def test_reopened_channel_stops_only_after_final_window_reconfirms() -> None:
    solver = NonMonotonicPairSolver()
    result = solve_model_be_nodal(
        POSITIONS[:2],
        0.1,
        1.0,
        1.0,
        0.0,
        0.8,
        lmax_max=8,
        solver=solver,
    )

    record = result.pair_ledger[0]
    interaction = _interaction_convergence(record)
    assert record.final_lmax == 7
    assert record.attempted_lmax == (2, 3, 4, 5, 6, 7)
    assert interaction.confirmation_lmax == 4
    assert all(step.applicable for step in interaction.history[-2:])
    assert all(
        step.successive_change <= 1.0e-5
        for step in interaction.history[-2:]
    )
    assert interaction.confirmed
    assert result.eligible
    assert all(
        not channel.applicable or channel.confirmed
        for channel in record.convergence
    )


def test_solver_failure_is_explicit_and_later_pairs_are_still_audited() -> None:
    failed_pair = _key(POSITIONS[[0, 1]])
    solver = FakePairSolver(fail_key=failed_pair)
    result = _solve(solver)

    assert not result.eligible
    assert result.forces_xyz is None
    assert result.failure_stage == "pair_eligibility"
    assert "1:(0, 1)=pair_solver" in result.failure_reason
    failed = result.pair_ledger[0]
    assert failed.final_lmax == 2
    assert failed.failed_lmax == 3
    assert failed.attempted_lmax == (2, 3)
    assert failed.evaluated_lmax == (2,)
    assert failed.interaction_forces_xyz is None
    assert failed.failure_stage == "pair_solver"
    assert failed.failure_reason == "RuntimeError: injected pair failure"
    assert all(record.eligible for record in result.pair_ledger[1:])
    assert any(call[0] == _key(POSITIONS[[1, 2]]) for call in solver.calls)


def test_nonconvergence_and_numerical_gate_failure_remain_distinct() -> None:
    varying = _key(POSITIONS[[0, 2]])
    bad_diagnostics = _key(POSITIONS[[1, 2]])
    solver = FakePairSolver(
        varying_key=varying,
        bad_diagnostics_key=bad_diagnostics,
    )
    result = _solve(solver, lmax_max=5)

    assert not result.eligible
    assert result.forces_xyz is None
    assert result.pair_ledger[0].eligible
    unconverged = result.pair_ledger[1]
    assert not unconverged.converged
    assert unconverged.failure_stage == "convergence"
    assert unconverged.interaction_forces_xyz is not None
    diagnostic_failure = result.pair_ledger[2]
    assert diagnostic_failure.converged
    assert diagnostic_failure.failure_stage == "numerical_diagnostics"
    assert diagnostic_failure.diagnostics is not None
    assert not diagnostic_failure.diagnostics.passed


def _snapshot(result) -> tuple:
    return (
        None if result.forces_xyz is None else result.forces_xyz.tobytes(),
        result.eligible,
        result.failure_stage,
        result.failure_reason,
        tuple(
            (
                record.pair_order,
                record.particle_indices,
                record.positions_xyz.tobytes(),
                record.interaction_forces_xyz.tobytes(),
                record.attempted_lmax,
                record.evaluated_lmax,
                record.final_lmax,
                tuple(
                    (
                        channel.channel,
                        channel.applicable,
                        channel.confirmed,
                        channel.confirmation_lmax,
                        tuple(channel.history),
                    )
                    for channel in record.convergence
                ),
            )
            for record in result.pair_ledger
        ),
    )


def test_repeated_fake_solver_runs_are_deterministic() -> None:
    first = _solve(FakePairSolver())
    second = _solve(FakePairSolver())
    assert _snapshot(first) == _snapshot(second)


@pytest.mark.parametrize(
    ("owner", "attribute", "value"),
    [
        ("solution", "production_solver", "legacy"),
        ("solution", "balanced_condition_number", 10.0),
        ("solution", "balanced_backward_error", 1.0e-12),
        ("solution", "effective_incident_closure_error", 1.0e-12),
        ("solution", "scattering_closure_error", 1.0e-12),
        ("result", "decomposition_residual", 1.0e-12),
        (
            "solution",
            "scattered_coefficients",
            np.array([[complex(float("nan"), 0.0)]]),
        ),
        ("solution", "balanced_system_matrix", np.eye(3)),
        (
            "result",
            "interaction_forces_xyz",
            np.array([[1.0, 0.0, 1.0e-8], [-1.0, 0.0, 0.0]]),
        ),
        ("result", "total_forces_xyz", np.zeros((1, 3))),
    ],
)
def test_each_established_model_e_gate_is_enforced(
    owner, attribute, value
) -> None:
    result = _fake_result(POSITIONS[:2], 5)
    target = result.solution if owner == "solution" else result
    setattr(target, attribute, value)

    diagnostics = evaluate_model_e_numerical_diagnostics(result)
    assert not diagnostics.passed


@pytest.mark.parametrize(
    ("options", "message"),
    [
        ({"diagnostic_tolerance": 0.0}, "finite and positive"),
        ({"maximum_balanced_condition": float("inf")}, "finite and positive"),
        ({"planar_tolerance_factor": -1.0}, "finite and positive"),
    ],
)
def test_numerical_gate_thresholds_are_validated(options, message) -> None:
    with pytest.raises(ValueError, match=message):
        evaluate_model_e_numerical_diagnostics(
            _fake_result(POSITIONS[:2], 5),
            **options,
        )


@pytest.mark.parametrize(
    ("positions", "kwargs", "message"),
    [
        (np.zeros((1, 3)), {}, "at least two"),
        (np.array([[0.0, 0.0, 0.0], [1.9, 0.0, 0.0]]), {}, "separation"),
        (np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 1.0e-4]]), {}, "nodal plane"),
        (POSITIONS, {"lmax_min": 1}, "at least 2"),
        (POSITIONS, {"lmax_min": 6, "lmax_max": 5}, "must not exceed"),
        (POSITIONS, {"minimum_stop_lmax": 22}, "within the lmax range"),
    ],
)
def test_input_validation_precedes_solver_calls(positions, kwargs, message) -> None:
    solver = FakePairSolver()
    with pytest.raises((TypeError, ValueError), match=message):
        solve_model_be_nodal(
            positions,
            0.1,
            1.0,
            1.0,
            0.0,
            0.8,
            solver=solver,
            **kwargs,
        )
    assert solver.calls == []


@pytest.mark.parametrize(
    ("options", "exception", "message"),
    [
        ({"solver": None}, TypeError, "callable"),
        ({"lmax_min": 2.0}, TypeError, "integer"),
        ({"convergence_tolerance": 0.0}, ValueError, "positive"),
        (
            {"lmax_min": 4, "minimum_stop_lmax": 3},
            ValueError,
            "within the lmax range",
        ),
    ],
)
def test_configuration_validation_precedes_solver_calls(
    options, exception, message
) -> None:
    supplied = dict(options)
    fake = FakePairSolver()
    solver = supplied.pop("solver", fake)
    with pytest.raises(exception, match=message):
        solve_model_be_nodal(
            POSITIONS,
            0.1,
            1.0,
            1.0,
            0.0,
            0.8,
            solver=solver,
            **supplied,
        )
    assert fake.calls == []


def test_malformed_injected_result_is_an_explicit_pair_failure() -> None:
    calls: list[tuple[tuple[tuple[float, ...], ...], int]] = []

    def malformed_solver(
        positions, k, radius, energy, f0, f1, order
    ) -> SimpleNamespace:
        del k, radius, energy, f0, f1
        values = np.asarray(positions, dtype=float)
        calls.append((_key(values), order))
        result = _fake_result(values, order)
        result.lmax = order + 1
        return result

    result = _solve(malformed_solver)
    assert not result.eligible
    assert result.forces_xyz is None

    assert len(calls) == result.pair_count == 3
    assert all(record.failure_stage == "pair_solver" for record in result.pair_ledger)
    assert all("requested order" in record.failure_reason for record in result.pair_ledger)


def test_injected_result_lmax_must_be_an_integer() -> None:
    def float_order_solver(
        positions, k, radius, energy, f0, f1, order
    ) -> SimpleNamespace:
        del k, radius, energy, f0, f1
        result = _fake_result(np.asarray(positions, dtype=float), order)
        result.lmax = float(order)
        return result

    result = _solve(float_order_solver)
    assert not result.eligible
    assert result.forces_xyz is None
    assert all(record.failure_stage == "pair_solver" for record in result.pair_ledger)
    assert all("must be an integer" in record.failure_reason for record in result.pair_ledger)



def test_rigid_boundary_requires_the_documented_f0_sentinel() -> None:
    solver = FakePairSolver()
    with pytest.raises(ValueError, match="sentinel f0=0"):
        solve_model_be_nodal(
            POSITIONS,
            0.1,
            1.0,
            1.0,
            0.2,
            1.0,
            solver=solver,
        )
    assert solver.calls == []
