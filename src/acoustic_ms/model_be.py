"""Independently converged complete-dimer Model-E baseline ``B_E``."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
import numbers
from typing import Callable

import numpy as np
from numpy.typing import NDArray

from .external_validation import minimum_two_step_confirmation, successive_change
from .mie_multiparticle import _positive_real
from .mie_scattering import mie_scattering_coefficients_from_contrasts
from .model_e import (
    ModelENodalResult,
    ModelENumericalDiagnostics,
    _energy,
    evaluate_model_e_numerical_diagnostics,
    solve_model_e_nodal,
)
from .solver import _validate_positions


FloatArray = NDArray[np.float64]
ModelEPairSolver = Callable[..., ModelENodalResult]

_FORCE_CHANNELS = (
    ("total", "total_forces_xyz"),
    ("interaction", "interaction_forces_xyz"),
    ("external_scattered", "external_scattered_forces_xyz"),
    ("scattered_scattered", "scattered_scattered_forces_xyz"),
)


@dataclass(frozen=True)
class ModelBEChannelStep:
    """One order's successive-change evidence for one force channel."""

    lmax: int
    successive_change: float
    absolute_change: float
    applicable: bool


@dataclass(frozen=True)
class ModelBEChannelConvergence:
    """Independent convergence summary for one dimer force channel.

    ``confirmed`` describes only the two most recent changes at the current
    order.  ``confirmation_lmax`` retains the first historical confirmation.
    """

    channel: str
    applicable: bool
    confirmed: bool
    confirmation_lmax: int | None
    final_successive_change: float
    final_absolute_change: float
    history: tuple[ModelBEChannelStep, ...]


@dataclass(frozen=True)
class ModelBEPairRecord:
    """Deterministic ledger entry for one isolated unordered pair."""

    pair_order: int
    particle_indices: tuple[int, int]
    positions_xyz: FloatArray
    interaction_forces_xyz: FloatArray | None
    attempted_lmax: tuple[int, ...]
    evaluated_lmax: tuple[int, ...]
    final_lmax: int | None
    failed_lmax: int | None
    convergence: tuple[ModelBEChannelConvergence, ...]
    diagnostics: ModelENumericalDiagnostics | None
    converged: bool
    eligible: bool
    failure_stage: str | None
    failure_reason: str | None


@dataclass(frozen=True)
class ModelBEResult:
    r"""Complete pairwise Model-E baseline.

    ``forces_xyz`` is available only when every ledger entry is eligible.  It
    is never a partial sum: an ineligible pair makes the global field ``None``.
    """

    forces_xyz: FloatArray | None
    pair_ledger: tuple[ModelBEPairRecord, ...]
    particle_count: int
    pair_count: int
    eligible: bool
    failure_stage: str | None
    failure_reason: str | None


def _readonly_float_array(values: object) -> FloatArray:
    result = np.array(values, dtype=float, copy=True)
    result.setflags(write=False)
    return result


def _integer(name: str, value: object, *, minimum: int) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(
        value, numbers.Integral
    ):
        raise TypeError(f"{name} must be an integer")
    result = int(value)
    if result < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    return result


def _positive_tolerance(name: str, value: object) -> float:
    if isinstance(value, (bool, np.bool_, complex, np.complexfloating)) or not isinstance(
        value, numbers.Real
    ):
        raise TypeError(f"{name} must be a real scalar")
    result = float(value)
    if not np.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return result


def _validate_order_result(result: object, order: int) -> None:
    if not hasattr(result, "lmax"):
        raise ValueError("pair solver result lmax must equal the requested order")
    result_order = _integer("pair solver result lmax", result.lmax, minimum=2)
    if result_order != order:
        raise ValueError("pair solver result lmax must equal the requested order")
    for _, attribute in _FORCE_CHANNELS:
        values = np.asarray(getattr(result, attribute), dtype=float)
        if values.shape != (2, 3):
            raise ValueError(
                f"pair solver {attribute} must have shape (2, 3)"
            )


def _channel_convergence(
    results: list[ModelENodalResult],
    tolerance: float,
) -> tuple[ModelBEChannelConvergence, ...]:
    summaries: list[ModelBEChannelConvergence] = []
    orders = [int(result.lmax) for result in results]
    for channel, attribute in _FORCE_CHANNELS:
        steps: list[ModelBEChannelStep] = []
        changes: list[float] = []
        applicable: list[bool] = []
        for index, result in enumerate(results):
            if index == 0:
                change, flag, absolute = 0.0, False, 0.0
            else:
                change, flag, absolute = successive_change(
                    getattr(result, attribute),
                    getattr(results[index - 1], attribute),
                )
            changes.append(change)
            applicable.append(flag)
            steps.append(
                ModelBEChannelStep(
                    lmax=int(result.lmax),
                    successive_change=float(change),
                    absolute_change=float(absolute),
                    applicable=bool(flag),
                )
            )
        confirmation = minimum_two_step_confirmation(
            changes,
            applicable,
            orders,
            tolerance=tolerance,
        )
        final_confirmed = bool(
            len(changes) >= 2
            and applicable[-2]
            and applicable[-1]
            and changes[-2] <= tolerance
            and changes[-1] <= tolerance
        )
        final = steps[-1]
        summaries.append(
            ModelBEChannelConvergence(
                channel=channel,
                applicable=any(applicable),
                confirmed=final_confirmed,
                confirmation_lmax=confirmation or None,
                final_successive_change=final.successive_change,
                final_absolute_change=final.absolute_change,
                history=tuple(steps),
            )
        )
    return tuple(summaries)


def _all_applicable_channels_confirmed(
    convergence: tuple[ModelBEChannelConvergence, ...],
) -> bool:
    return all(
        not channel.applicable or channel.confirmed for channel in convergence
    )


def _failure_text(error: Exception) -> str:
    return f"{type(error).__name__}: {error}"


def solve_model_be_nodal(
    positions_xyz: object,
    k: object,
    radius: object,
    energy_density: object,
    f0: object,
    f1: object,
    *,
    lmax_min: object = 2,
    lmax_max: object = 21,
    minimum_stop_lmax: object = 5,
    convergence_tolerance: object = 1.0e-5,
    solver: ModelEPairSolver = solve_model_e_nodal,
) -> ModelBEResult:
    r"""Return independently converged isolated-dimer Model-E forces.

    For every unordered pair ``i < j``, ``solver`` receives the original two
    nodal-plane positions in that order.  The pair's Model-E *interaction*
    forces are accumulated only after every pair passes convergence and the
    established numerical gates.  Solver failures are retained and do not
    prevent later pairs from being audited.
    """

    if not callable(solver):
        raise TypeError("solver must be callable")
    minimum_order = _integer("lmax_min", lmax_min, minimum=2)
    maximum_order = _integer("lmax_max", lmax_max, minimum=2)
    stop_order = _integer("minimum_stop_lmax", minimum_stop_lmax, minimum=2)
    if minimum_order > maximum_order:
        raise ValueError("lmax_min must not exceed lmax_max")
    if not minimum_order <= stop_order <= maximum_order:
        raise ValueError("minimum_stop_lmax must lie within the lmax range")
    tolerance = _positive_tolerance(
        "convergence_tolerance", convergence_tolerance
    )

    wave_number = _positive_real("k", k)
    sphere_radius = _positive_real("radius", radius)
    energy = _energy(energy_density)
    # Reuse the exact-Mie contrast validator before an injected solver runs.
    mie_scattering_coefficients_from_contrasts(
        wave_number * sphere_radius,
        f0,
        f1,
        minimum_order,
    )
    monopole = float(f0)
    dipole = float(f1)
    if dipole == 1.0 and monopole != 0.0:
        raise ValueError("rigid Model E requires the API sentinel f0=0")
    positions = _validate_positions(positions_xyz, sphere_radius)
    if len(positions) < 2:
        raise ValueError("Model B_E requires at least two particles")

    ledger: list[ModelBEPairRecord] = []
    for pair_order, (first, second) in enumerate(
        combinations(range(len(positions)), 2),
        start=1,
    ):
        indices = (first, second)
        pair_positions = _readonly_float_array(positions[[first, second]])
        attempted: list[int] = []
        results: list[ModelENodalResult] = []
        diagnostics: list[ModelENumericalDiagnostics] = []
        failure_stage: str | None = None
        failure_reason: str | None = None
        failed_lmax: int | None = None

        for order in range(minimum_order, maximum_order + 1):
            attempted.append(order)
            try:
                result = solver(
                    pair_positions,
                    wave_number,
                    sphere_radius,
                    energy,
                    monopole,
                    dipole,
                    order,
                )
                _validate_order_result(result, order)
                order_diagnostics = evaluate_model_e_numerical_diagnostics(result)
            except Exception as error:
                failure_stage = "pair_solver"
                failure_reason = _failure_text(error)
                failed_lmax = order
                break
            if not order_diagnostics.finite:
                failure_stage = "numerical_diagnostics"
                failure_reason = "non-finite Model-E pair result"
                failed_lmax = order
                break
            results.append(result)
            diagnostics.append(order_diagnostics)
            convergence = _channel_convergence(results, tolerance)
            if (
                order >= stop_order
                and _all_applicable_channels_confirmed(convergence)
            ):
                break

        convergence = (
            _channel_convergence(results, tolerance) if results else tuple()
        )
        converged = bool(
            results and _all_applicable_channels_confirmed(convergence)
        )
        final_diagnostics = diagnostics[-1] if diagnostics else None
        pair_forces: FloatArray | None = None
        if failure_stage is None and results:
            pair_forces = _readonly_float_array(
                results[-1].interaction_forces_xyz
            )
            reasons: list[str] = []
            stages: list[str] = []
            if not converged:
                stages.append("convergence")
                reasons.append(
                    "not all applicable channels confirmed by lmax_max"
                )
            if final_diagnostics is None or not final_diagnostics.passed:
                stages.append("numerical_diagnostics")
                reasons.append("final Model-E numerical gates did not pass")
            if stages:
                failure_stage = "+".join(stages)
                failure_reason = "; ".join(reasons)

        eligible = bool(
            failure_stage is None
            and pair_forces is not None
            and converged
            and final_diagnostics is not None
            and final_diagnostics.passed
        )
        ledger.append(
            ModelBEPairRecord(
                pair_order=pair_order,
                particle_indices=indices,
                positions_xyz=pair_positions,
                interaction_forces_xyz=pair_forces,
                attempted_lmax=tuple(attempted),
                evaluated_lmax=tuple(int(result.lmax) for result in results),
                final_lmax=int(results[-1].lmax) if results else None,
                failed_lmax=failed_lmax,
                convergence=convergence,
                diagnostics=final_diagnostics,
                converged=converged,
                eligible=eligible,
                failure_stage=failure_stage,
                failure_reason=failure_reason,
            )
        )

    globally_eligible = all(record.eligible for record in ledger)
    forces: FloatArray | None = None
    failure_stage = None
    failure_reason = None
    if globally_eligible:
        accumulated = np.zeros((len(positions), 3), dtype=float)
        for record in ledger:
            assert record.interaction_forces_xyz is not None
            first, second = record.particle_indices
            accumulated[first] += record.interaction_forces_xyz[0]
            accumulated[second] += record.interaction_forces_xyz[1]
        forces = _readonly_float_array(accumulated)
    else:
        failure_stage = "pair_eligibility"
        failure_reason = "; ".join(
            f"{record.pair_order}:{record.particle_indices}={record.failure_stage}"
            for record in ledger
            if not record.eligible
        )

    return ModelBEResult(
        forces_xyz=forces,
        pair_ledger=tuple(ledger),
        particle_count=len(positions),
        pair_count=len(ledger),
        eligible=globally_eligible,
        failure_stage=failure_stage,
        failure_reason=failure_reason,
    )
