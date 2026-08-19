"""Deterministic, single-use runner for the preregistered P1.5 pilot."""

from __future__ import annotations

from contextlib import contextmanager
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
from hashlib import sha256
import io
import json
import os
from pathlib import Path
import platform
import resource
import shutil
import signal
import sys
import tempfile
import time
from types import TracebackType
from typing import Any, Callable, Iterator, Mapping

import numpy as np
import scipy

from . import mie_multiparticle as _mie_module
from . import model_e as _model_e_module
from .external_validation import canonical_coordinate_hash
from .model_be import ModelBEResult, solve_model_be_nodal
from .model_e import (
    ModelENodalResult,
    evaluate_model_e_numerical_diagnostics,
    solve_model_e_nodal,
)
from .paper_pipeline import (
    P1_FROZEN_MANIFEST_SHA256,
    validate_executable_manifest_file,
    validate_manifest_file,
)


PILOT_CAMPAIGN_ID = "p1_dimer_resource_pilot"
PILOT_CASE_ID = "p1_pilot_rigid_ka010_d0210_t000"
CONFIRMATORY_CAMPAIGN_ID = "p1_dimer_confirmatory"
PILOT_MANIFEST_RELATIVE = Path("campaigns/p1/pilot_manifest.yaml")
CONFIRMATORY_MANIFEST_RELATIVE = Path("campaigns/p1/campaign_manifest.yaml")
PILOT_OUTPUT_RELATIVE = Path("campaigns/p1/pilot")
PERFORMANCE_RELATIVE = PILOT_OUTPUT_RELATIVE / "performance.csv"
_ARTIFACT_NAMES = (
    "data_raw.csv",
    "data_derived.csv",
    "data_plot.csv",
    "failures.csv",
    "performance.csv",
)
_CHANNELS = (
    ("total", "total_forces_xyz"),
    ("interaction", "interaction_forces_xyz"),
    ("external_scattered", "external_scattered_forces_xyz"),
    ("scattered_scattered", "scattered_scattered_forces_xyz"),
)


class PilotExecutionError(RuntimeError):
    """Raised when the frozen pilot cannot be executed or published safely."""


class _PilotWallTimeout(Exception):
    """Internal exception used to serialize the frozen wall limit."""


@dataclass(frozen=True)
class PilotConfiguration:
    """Validated response-blind configuration for the only P1.5 case."""

    root: Path
    pilot_manifest: Mapping[str, Any]
    confirmatory_manifest: Mapping[str, Any]
    case: Mapping[str, Any]
    positions_xyz: np.ndarray
    output_directory: Path


@dataclass(frozen=True)
class PilotExecutionSummary:
    """Published P1.5 resource outcome."""

    decision: str
    stop_reason: str
    eligible: bool
    converged: bool
    attempted_lmax: tuple[int, ...]
    evaluated_lmax: tuple[int, ...]
    final_lmax: int | None
    total_wall_seconds: float
    peak_rss_bytes: int
    artifact_sha256: Mapping[str, str]


@dataclass
class _ResourceState:
    wall_limit_exceeded: bool = False
    memory_limit_exceeded: bool = False
    effective_address_space_limit_bytes: int | None = None


@dataclass
class _Attempt:
    lmax: int
    result: ModelENodalResult | None
    order_wall_seconds: float
    assembly_diagnostics_seconds: float
    linear_solve_seconds: float
    force_postprocess_seconds: float
    peak_rss_bytes: int
    error: str | None


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


def _utc_text(value: datetime) -> str:
    if value.tzinfo is None or value.utcoffset() is None:
        raise ValueError("UTC clock must return a timezone-aware datetime")
    return value.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")


def _peak_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if sys.platform == "darwin" else value * 1024


def _format(value: object) -> object:
    if value is None:
        return ""
    if isinstance(value, (bool, np.bool_)):
        return str(bool(value)).lower()
    if isinstance(value, (float, np.floating)):
        return format(float(value), ".17g")
    return value


def _csv_bytes(fields: tuple[str, ...], rows: list[Mapping[str, object]]) -> bytes:
    if not rows:
        raise PilotExecutionError("pilot artifacts cannot contain an empty table")
    unexpected = {
        key for row in rows for key in row if key not in fields
    }
    missing = {
        key for row in rows for key in fields if key not in row
    }
    if unexpected or missing:
        raise PilotExecutionError(
            f"inconsistent CSV fields: unexpected={sorted(unexpected)}, "
            f"missing={sorted(missing)}"
        )
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow({key: _format(row[key]) for key in fields})
    return stream.getvalue().encode("utf-8")


def _read_csv_bytes(payload: bytes) -> list[dict[str, str]]:
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise PilotExecutionError("pilot CSV is not UTF-8") from exc
    rows = list(csv.DictReader(io.StringIO(text, newline="")))
    if not rows:
        raise PilotExecutionError("pilot CSV contains no data rows")
    return rows


def _default_environment() -> dict[str, str]:
    thread_keys = (
        "OPENBLAS_NUM_THREADS",
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
    )
    return {
        "python_version": platform.python_version(),
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "pythonhashseed": os.environ.get("PYTHONHASHSEED", ""),
        **{key.lower(): os.environ.get(key, "") for key in thread_keys},
    }


def _require_exact_outputs(manifest: Mapping[str, Any]) -> None:
    expected = {
        "raw": "campaigns/p1/pilot/data_raw.csv",
        "derived": "campaigns/p1/pilot/data_derived.csv",
        "plot_ready": "campaigns/p1/pilot/data_plot.csv",
        "failure_log": "campaigns/p1/pilot/failures.csv",
    }
    if manifest["outputs"] != expected:
        raise PilotExecutionError("pilot output paths differ from the frozen paths")


def load_p1_5_configuration(root: str | Path) -> PilotConfiguration:
    """Validate both frozen manifests and return the only executable case."""

    repository = Path(root).resolve()
    pilot_path = repository / PILOT_MANIFEST_RELATIVE
    confirmatory_path = repository / CONFIRMATORY_MANIFEST_RELATIVE
    pilot = validate_executable_manifest_file(
        pilot_path,
        expected_campaign_id=PILOT_CAMPAIGN_ID,
    )
    confirmatory = validate_manifest_file(confirmatory_path, kind="campaign")
    if confirmatory["campaign_id"] != CONFIRMATORY_CAMPAIGN_ID:
        raise PilotExecutionError("unexpected confirmatory campaign identity")
    if any(case["enabled"] for case in confirmatory["cases"]):
        raise PilotExecutionError("the 102-case confirmatory campaign is blocked")
    if len(confirmatory["cases"]) != 102:
        raise PilotExecutionError("the confirmatory campaign must retain 102 cases")
    enabled = [case for case in pilot["cases"] if case["enabled"]]
    if len(enabled) != 1 or enabled[0]["case_id"] != PILOT_CASE_ID:
        raise PilotExecutionError("only the frozen P1.5 pilot case may execute")
    if pilot["provenance"]["manifest_sha256"] != P1_FROZEN_MANIFEST_SHA256[
        PILOT_CAMPAIGN_ID
    ]:
        raise PilotExecutionError("pilot manifest hash differs from the public lock")
    if confirmatory["provenance"]["manifest_sha256"] != (
        P1_FROZEN_MANIFEST_SHA256[CONFIRMATORY_CAMPAIGN_ID]
    ):
        raise PilotExecutionError(
            "confirmatory manifest hash differs from the public lock"
        )
    _require_exact_outputs(pilot)

    case = enabled[0]
    parameters = case["parameters"]
    exact_parameters = {
        "ka": 0.1,
        "k_rad_m": 0.1,
        "material_model": "rigid",
        "f0": 0,
        "f0_applicable": False,
        "f1": 1,
        "distance_ratio": 2.1,
        "theta_rad": 0,
        "evidence_role": "resource_pilot",
        "include_in_scientific_tables": False,
    }
    if any(parameters[key] != value for key, value in exact_parameters.items()):
        raise PilotExecutionError("pilot physical parameters differ from P1.4")
    numerical = pilot["numerical"]
    if numerical != {
        "model": "E",
        "lmax_min": 2,
        "lmax_max": 21,
        "minimum_stop_lmax": 5,
        "convergence_tolerance": 1.0e-5,
        "consecutive_passes": 2,
        "required_channels": [channel for channel, _ in _CHANNELS],
        "require_all_applicable_channels": True,
        "failure_policy": (
            "retain every attempt; no imputation or physics/tolerance/solver "
            "retry; record controlled stage and reason"
        ),
    }:
        raise PilotExecutionError("pilot numerical policy differs from P1.4")
    resources = pilot["resources"]
    if resources["worker_count"] != 1 or resources["blas_threads"] != 1:
        raise PilotExecutionError("P1.5 requires one worker and one BLAS thread")
    if resources["wall_seconds_per_case"] != 1800:
        raise PilotExecutionError("P1.5 wall limit must remain 1800 seconds")
    if resources["peak_rss_bytes_per_case"] != 4 * 1024**3:
        raise PilotExecutionError("P1.5 memory limit must remain 4 GiB")
    if resources["limits_status"] != "provisional":
        raise PilotExecutionError("P1.5 limits must remain provisional")

    radius = float(pilot["physical"]["radius_m"])
    distance = radius * float(parameters["distance_ratio"])
    angle = float(parameters["theta_rad"])
    displacement = np.array(
        [distance * np.cos(angle), distance * np.sin(angle), 0.0],
        dtype=float,
    )
    positions = np.vstack((-0.5 * displacement, 0.5 * displacement))
    positions.setflags(write=False)
    return PilotConfiguration(
        root=repository,
        pilot_manifest=pilot,
        confirmatory_manifest=confirmatory,
        case=case,
        positions_xyz=positions,
        output_directory=repository / PILOT_OUTPUT_RELATIVE,
    )


@contextmanager
def _instrument_model_e(
    clock: Callable[[], float],
) -> Iterator[dict[str, float]]:
    phases = {"linear": 0.0, "force": 0.0}
    original_linear_solve = _mie_module.np.linalg.solve
    original_force = _model_e_module.complete_radiation_force_from_bsc

    def measured_linear_solve(*args: object, **kwargs: object) -> object:
        started = clock()
        try:
            return original_linear_solve(*args, **kwargs)
        finally:
            phases["linear"] += max(0.0, clock() - started)

    def measured_force(*args: object, **kwargs: object) -> object:
        started = clock()
        try:
            return original_force(*args, **kwargs)
        finally:
            phases["force"] += max(0.0, clock() - started)

    _mie_module.np.linalg.solve = measured_linear_solve
    _model_e_module.complete_radiation_force_from_bsc = measured_force
    try:
        yield phases
    finally:
        _mie_module.np.linalg.solve = original_linear_solve
        _model_e_module.complete_radiation_force_from_bsc = original_force


class _TimedSolver:
    def __init__(
        self,
        solver: Callable[..., ModelENodalResult],
        *,
        clock: Callable[[], float],
        rss_reader: Callable[[], int],
        state: _ResourceState,
        memory_limit_bytes: int,
    ) -> None:
        self.solver = solver
        self.clock = clock
        self.rss_reader = rss_reader
        self.state = state
        self.memory_limit_bytes = memory_limit_bytes
        self.attempts: list[_Attempt] = []

    def __call__(self, *args: object, **kwargs: object) -> ModelENodalResult:
        lmax = int(args[6] if len(args) > 6 else kwargs["lmax"])
        started = self.clock()
        phases = {"linear": 0.0, "force": 0.0}
        result: ModelENodalResult | None = None
        error_text: str | None = None
        try:
            if self.solver is solve_model_e_nodal:
                with _instrument_model_e(self.clock) as phases:
                    result = self.solver(*args, **kwargs)
            else:
                result = self.solver(*args, **kwargs)
            return result
        except MemoryError as exc:
            self.state.memory_limit_exceeded = True
            error_text = f"{type(exc).__name__}: {exc}"
            raise
        except _PilotWallTimeout as exc:
            self.state.wall_limit_exceeded = True
            error_text = f"{type(exc).__name__}: {exc}"
            raise
        except Exception as exc:
            error_text = f"{type(exc).__name__}: {exc}"
            raise
        finally:
            elapsed = max(0.0, self.clock() - started)
            peak = int(self.rss_reader())
            if peak > self.memory_limit_bytes:
                self.state.memory_limit_exceeded = True
            linear = float(phases["linear"])
            force = float(phases["force"])
            if self.solver is solve_model_e_nodal:
                assembly = max(0.0, elapsed - linear - force)
            else:
                assembly = 0.0
                linear = elapsed
                force = 0.0
            self.attempts.append(
                _Attempt(
                    lmax=lmax,
                    result=result,
                    order_wall_seconds=elapsed,
                    assembly_diagnostics_seconds=assembly,
                    linear_solve_seconds=linear,
                    force_postprocess_seconds=force,
                    peak_rss_bytes=peak,
                    error=error_text,
                )
            )


@contextmanager
def _resource_limits(
    wall_seconds: int,
    memory_bytes: int,
    state: _ResourceState,
) -> Iterator[None]:
    old_limit = resource.getrlimit(resource.RLIMIT_AS)
    old_handler = signal.getsignal(signal.SIGALRM)
    old_timer = signal.getitimer(signal.ITIMER_REAL)

    def wall_handler(_signum: int, _frame: object) -> None:
        state.wall_limit_exceeded = True
        raise _PilotWallTimeout(f"wall limit of {wall_seconds} seconds exceeded")

    infinity = resource.RLIM_INFINITY
    old_soft, old_hard = old_limit
    if old_soft == infinity or old_soft > memory_bytes:
        effective_soft = memory_bytes
    else:
        effective_soft = int(old_soft)
    state.effective_address_space_limit_bytes = effective_soft
    try:
        resource.setrlimit(resource.RLIMIT_AS, (effective_soft, old_hard))
        signal.signal(signal.SIGALRM, wall_handler)
        signal.setitimer(signal.ITIMER_REAL, float(wall_seconds))
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, *old_timer)
        signal.signal(signal.SIGALRM, old_handler)
        resource.setrlimit(resource.RLIMIT_AS, old_limit)


def _confirmation_at(history: list[Mapping[str, object]], index: int) -> bool:
    if index < 2:
        return False
    latest = history[index]
    previous = history[index - 1]
    return bool(
        latest["applicable"]
        and previous["applicable"]
        and float(latest["successive_change"]) <= 1.0e-5
        and float(previous["successive_change"]) <= 1.0e-5
    )


def _outcome(
    result: ModelBEResult | None,
    timed_solver: _TimedSolver,
    state: _ResourceState,
    total_wall_seconds: float,
    wall_limit_seconds: int,
) -> tuple[str, str, str | None]:
    if total_wall_seconds > wall_limit_seconds:
        state.wall_limit_exceeded = True
    if state.wall_limit_exceeded:
        return (
            "NO_GO_P1.6_RESOURCE_LIMIT",
            "wall_limit_exceeded",
            f"wall limit of {wall_limit_seconds} seconds exceeded",
        )
    if state.memory_limit_exceeded:
        return (
            "NO_GO_P1.6_RESOURCE_LIMIT",
            "memory_limit_exceeded",
            "4 GiB memory limit exceeded",
        )
    if result is None:
        return (
            "INCONCLUSIVE_P1.5",
            "infrastructure_failure",
            "pilot did not return a Model B_E ledger",
        )
    if any(attempt.error is not None for attempt in timed_solver.attempts):
        return (
            "INCONCLUSIVE_P1.5",
            "infrastructure_failure",
            next(
                attempt.error
                for attempt in timed_solver.attempts
                if attempt.error is not None
            ),
        )
    pair = result.pair_ledger[0]
    if result.eligible:
        return "GO_P1.6A_BLIND_FREEZE", "all_channels_confirmed", None
    if pair.final_lmax == 21 and not pair.converged:
        return (
            "GO_P1.6A_BLIND_FREEZE",
            "unconfirmed_at_21",
            "unconfirmed_at_21",
        )
    return (
        "GO_P1.6A_BLIND_FREEZE",
        pair.failure_stage or "controlled_ineligibility",
        pair.failure_reason or "controlled ineligible resource result",
    )


_RAW_FIELDS = (
    "schema_version", "campaign_id", "case_id", "solve_id",
    "classification", "created_utc", "completed_utc", "git_commit",
    "campaign_manifest_path", "campaign_manifest_sha256",
    "confirmatory_manifest_sha256", "command", "environment_json",
    "radius_m", "k_rad_m", "ka", "energy_density_j_m3", "f0",
    "f0_applicable", "f1", "material_model", "temporal_convention",
    "particle_count", "particle_index", "family", "position_x_m",
    "position_y_m", "position_z_m", "minimum_distance_m",
    "distance_ratio", "coordinate_sha256", "model", "force_channel",
    "lmax", "full_modes_per_particle", "active_modes_per_particle",
    "system_dimension", "solver_name", "force_x_n", "force_y_n",
    "force_z_n", "force_x_over_a2e0", "force_y_over_a2e0",
    "force_z_over_a2e0", "successive_change", "absolute_change",
    "change_applicable", "confirmation_lmax", "confirmed", "finite",
    "diagnostics_pass", "eligible", "failure_stage",
    "ineligibility_reason", "stop_reason", "evidence_role",
    "include_in_scientific_tables",
)

_PERFORMANCE_FIELDS = (
    "schema_version", "campaign_id", "case_id", "solve_id",
    "classification", "started_utc", "completed_utc", "git_commit",
    "campaign_manifest_path", "campaign_manifest_sha256",
    "confirmatory_manifest_sha256", "command", "environment_json",
    "worker_count", "blas_threads", "wall_limit_seconds",
    "memory_limit_bytes", "effective_address_space_limit_bytes", "lmax",
    "attempt_success", "assembly_diagnostics_seconds",
    "linear_solve_seconds", "force_postprocess_seconds",
    "order_wall_seconds", "accumulated_order_wall_seconds",
    "case_wall_seconds", "peak_rss_bytes", "full_modes_per_particle",
    "active_modes_per_particle", "system_dimension", "production_solver",
    "balanced_condition_number", "balanced_backward_error",
    "effective_incident_closure_error", "scattering_closure_error",
    "force_decomposition_residual", "max_abs_fz", "fz_tolerance",
    "finite", "mode_dimension_consistent", "planar_symmetry_pass",
    "diagnostics_pass", "attempt_error", "converged", "eligible",
    "stop_reason", "final_reason", "decision",
)

_FAILURE_FIELDS = (
    "schema_version", "campaign_id", "case_id", "solve_id",
    "classification", "created_utc", "completed_utc", "git_commit",
    "campaign_manifest_sha256", "attempted_lmax", "evaluated_lmax",
    "final_lmax", "converged", "eligible", "failure_stage",
    "failure_reason", "stop_reason", "decision", "evidence_role",
    "include_in_scientific_tables",
)

_DERIVED_FIELDS = (
    "schema_version", "analysis_id", "created_utc", "git_commit",
    "source_path", "source_sha256", "campaign_id", "case_id",
    "classification", "eligible", "inclusion_rule", "model",
    "reference_model", "metric", "value", "unit", "applicable", "reason",
)

_PLOT_FIELDS = (
    "schema_version", "figure_id", "panel_id", "series_id", "point_order",
    "case_id", "x_name", "x_value", "x_unit", "y_name", "y_value",
    "y_unit", "xerr_low", "xerr_high", "yerr_low", "yerr_high",
    "marker", "color", "linestyle", "label", "eligible", "annotation",
)


def _build_primary_tables(
    configuration: PilotConfiguration,
    timed_solver: _TimedSolver,
    result: ModelBEResult | None,
    *,
    solve_id: str,
    started_utc: str,
    completed_utc: str,
    source_commit: str,
    command: str,
    environment_json: str,
    total_wall_seconds: float,
    state: _ResourceState,
    decision: str,
    stop_reason: str,
    final_reason: str | None,
) -> dict[str, bytes]:
    pilot = configuration.pilot_manifest
    case = configuration.case
    parameters = case["parameters"]
    physical = pilot["physical"]
    numerical = pilot["numerical"]
    resources = pilot["resources"]
    pair = result.pair_ledger[0] if result is not None else None
    convergence_by_channel = {
        item.channel: item for item in pair.convergence
    } if pair is not None else {}
    coordinate_hash = canonical_coordinate_hash(configuration.positions_xyz)
    minimum_distance = float(
        np.linalg.norm(configuration.positions_xyz[1] - configuration.positions_xyz[0])
    )
    raw_rows: list[dict[str, object]] = []
    performance_rows: list[dict[str, object]] = []
    accumulated = 0.0

    for attempt in timed_solver.attempts:
        attempt_result = attempt.result
        accumulated += attempt.order_wall_seconds
        diagnostics = (
            evaluate_model_e_numerical_diagnostics(attempt_result)
            if attempt_result is not None
            else None
        )
        if attempt_result is None:
            full_modes = active_modes = system_dimension = None
            solver_name = None
        else:
            solution = attempt_result.solution
            full_modes = len(solution.modes)
            active_modes = len(solution.active_modes)
            system_dimension = int(solution.balanced_system_matrix.shape[0])
            solver_name = str(solution.production_solver)

        for channel, attribute in _CHANNELS:
            summary = convergence_by_channel.get(channel)
            history = [] if summary is None else [
                {
                    "lmax": step.lmax,
                    "successive_change": step.successive_change,
                    "absolute_change": step.absolute_change,
                    "applicable": step.applicable,
                }
                for step in summary.history
            ]
            history_index = next(
                (
                    index
                    for index, item in enumerate(history)
                    if int(item["lmax"]) == attempt.lmax
                ),
                None,
            )
            step = history[history_index] if history_index is not None else None
            confirmation_lmax = (
                summary.confirmation_lmax
                if summary is not None
                and summary.confirmation_lmax is not None
                and summary.confirmation_lmax <= attempt.lmax
                else None
            )
            confirmed = bool(
                history_index is not None
                and _confirmation_at(history, history_index)
            )
            force_values = (
                np.asarray(getattr(attempt_result, attribute), dtype=float)
                if attempt_result is not None
                else None
            )
            for particle_index, position in enumerate(configuration.positions_xyz):
                force = (
                    force_values[particle_index]
                    if force_values is not None
                    else (None, None, None)
                )
                scale = float(physical["radius_m"]) ** 2 * float(
                    physical["energy_density_j_m3"]
                )
                normalized = (
                    tuple(float(value) / scale for value in force)
                    if force_values is not None
                    else (None, None, None)
                )
                raw_rows.append(
                    {
                        "schema_version": "1.0.0",
                        "campaign_id": pilot["campaign_id"],
                        "case_id": case["case_id"],
                        "solve_id": solve_id,
                        "classification": pilot["classification"],
                        "created_utc": started_utc,
                        "completed_utc": completed_utc,
                        "git_commit": source_commit,
                        "campaign_manifest_path": str(PILOT_MANIFEST_RELATIVE),
                        "campaign_manifest_sha256": pilot["provenance"]["manifest_sha256"],
                        "confirmatory_manifest_sha256": configuration.confirmatory_manifest["provenance"]["manifest_sha256"],
                        "command": command,
                        "environment_json": environment_json,
                        "radius_m": physical["radius_m"],
                        "k_rad_m": parameters["k_rad_m"],
                        "ka": parameters["ka"],
                        "energy_density_j_m3": physical["energy_density_j_m3"],
                        "f0": parameters["f0"],
                        "f0_applicable": parameters["f0_applicable"],
                        "f1": parameters["f1"],
                        "material_model": parameters["material_model"],
                        "temporal_convention": physical["temporal_convention"],
                        "particle_count": case["particle_count"],
                        "particle_index": particle_index,
                        "family": case["family"],
                        "position_x_m": position[0],
                        "position_y_m": position[1],
                        "position_z_m": position[2],
                        "minimum_distance_m": minimum_distance,
                        "distance_ratio": parameters["distance_ratio"],
                        "coordinate_sha256": coordinate_hash,
                        "model": "B_E",
                        "force_channel": channel,
                        "lmax": attempt.lmax,
                        "full_modes_per_particle": full_modes,
                        "active_modes_per_particle": active_modes,
                        "system_dimension": system_dimension,
                        "solver_name": solver_name,
                        "force_x_n": force[0],
                        "force_y_n": force[1],
                        "force_z_n": force[2],
                        "force_x_over_a2e0": normalized[0],
                        "force_y_over_a2e0": normalized[1],
                        "force_z_over_a2e0": normalized[2],
                        "successive_change": None if step is None else step["successive_change"],
                        "absolute_change": None if step is None else step["absolute_change"],
                        "change_applicable": False if step is None else step["applicable"],
                        "confirmation_lmax": confirmation_lmax,
                        "confirmed": confirmed,
                        "finite": False if diagnostics is None else diagnostics.finite,
                        "diagnostics_pass": False if diagnostics is None else diagnostics.passed,
                        "eligible": False if result is None else result.eligible,
                        "failure_stage": None if pair is None else pair.failure_stage,
                        "ineligibility_reason": final_reason,
                        "stop_reason": stop_reason,
                        "evidence_role": parameters["evidence_role"],
                        "include_in_scientific_tables": parameters["include_in_scientific_tables"],
                    }
                )

        performance_rows.append(
            {
                "schema_version": "1.0.0",
                "campaign_id": pilot["campaign_id"],
                "case_id": case["case_id"],
                "solve_id": solve_id,
                "classification": pilot["classification"],
                "started_utc": started_utc,
                "completed_utc": completed_utc,
                "git_commit": source_commit,
                "campaign_manifest_path": str(PILOT_MANIFEST_RELATIVE),
                "campaign_manifest_sha256": pilot["provenance"]["manifest_sha256"],
                "confirmatory_manifest_sha256": configuration.confirmatory_manifest["provenance"]["manifest_sha256"],
                "command": command,
                "environment_json": environment_json,
                "worker_count": resources["worker_count"],
                "blas_threads": resources["blas_threads"],
                "wall_limit_seconds": resources["wall_seconds_per_case"],
                "memory_limit_bytes": resources["peak_rss_bytes_per_case"],
                "effective_address_space_limit_bytes": state.effective_address_space_limit_bytes,
                "lmax": attempt.lmax,
                "attempt_success": attempt.result is not None and attempt.error is None,
                "assembly_diagnostics_seconds": attempt.assembly_diagnostics_seconds,
                "linear_solve_seconds": attempt.linear_solve_seconds,
                "force_postprocess_seconds": attempt.force_postprocess_seconds,
                "order_wall_seconds": attempt.order_wall_seconds,
                "accumulated_order_wall_seconds": accumulated,
                "case_wall_seconds": total_wall_seconds,
                "peak_rss_bytes": attempt.peak_rss_bytes,
                "full_modes_per_particle": full_modes,
                "active_modes_per_particle": active_modes,
                "system_dimension": system_dimension,
                "production_solver": solver_name,
                "balanced_condition_number": None if diagnostics is None else diagnostics.balanced_condition_number,
                "balanced_backward_error": None if diagnostics is None else diagnostics.balanced_backward_error,
                "effective_incident_closure_error": None if diagnostics is None else diagnostics.effective_incident_closure_error,
                "scattering_closure_error": None if diagnostics is None else diagnostics.scattering_closure_error,
                "force_decomposition_residual": None if diagnostics is None else diagnostics.force_decomposition_residual,
                "max_abs_fz": None if diagnostics is None else diagnostics.max_abs_fz,
                "fz_tolerance": None if diagnostics is None else diagnostics.fz_tolerance,
                "finite": False if diagnostics is None else diagnostics.finite,
                "mode_dimension_consistent": False if diagnostics is None else diagnostics.mode_dimension_consistent,
                "planar_symmetry_pass": False if diagnostics is None else diagnostics.planar_symmetry_pass,
                "diagnostics_pass": False if diagnostics is None else diagnostics.passed,
                "attempt_error": attempt.error,
                "converged": False if pair is None else pair.converged,
                "eligible": False if result is None else result.eligible,
                "stop_reason": stop_reason,
                "final_reason": final_reason,
                "decision": decision,
            }
        )

    if not timed_solver.attempts:
        raise PilotExecutionError("pilot produced no attempted orders")
    failure_rows = [
        {
            "schema_version": "1.0.0",
            "campaign_id": pilot["campaign_id"],
            "case_id": case["case_id"],
            "solve_id": solve_id,
            "classification": pilot["classification"],
            "created_utc": started_utc,
            "completed_utc": completed_utc,
            "git_commit": source_commit,
            "campaign_manifest_sha256": pilot["provenance"]["manifest_sha256"],
            "attempted_lmax": ";".join(str(item.lmax) for item in timed_solver.attempts),
            "evaluated_lmax": "" if pair is None else ";".join(str(value) for value in pair.evaluated_lmax),
            "final_lmax": None if pair is None else pair.final_lmax,
            "converged": False if pair is None else pair.converged,
            "eligible": False if result is None else result.eligible,
            "failure_stage": None if pair is None else pair.failure_stage,
            "failure_reason": final_reason,
            "stop_reason": stop_reason,
            "decision": decision,
            "evidence_role": parameters["evidence_role"],
            "include_in_scientific_tables": parameters["include_in_scientific_tables"],
        }
    ]
    return {
        "data_raw.csv": _csv_bytes(_RAW_FIELDS, raw_rows),
        "performance.csv": _csv_bytes(_PERFORMANCE_FIELDS, performance_rows),
        "failures.csv": _csv_bytes(_FAILURE_FIELDS, failure_rows),
    }


def derive_p1_5_artifacts(
    raw_bytes: bytes,
    performance_bytes: bytes,
    failure_bytes: bytes,
) -> dict[str, bytes]:
    """Generate resource-only derived and plot tables without any solver."""

    raw = _read_csv_bytes(raw_bytes)
    performance = _read_csv_bytes(performance_bytes)
    failures = _read_csv_bytes(failure_bytes)
    first = raw[0]
    final = performance[-1]
    failure = failures[0]
    source_hash = sha256(raw_bytes).hexdigest()
    eligible = final["eligible"] == "true"
    reason = failure["failure_reason"]
    metrics = (
        ("pilot_case_wall_seconds", final["case_wall_seconds"], "s"),
        ("pilot_peak_rss_bytes", max(int(row["peak_rss_bytes"]) for row in performance), "byte"),
        ("pilot_attempt_count", len(performance), "count"),
        ("pilot_final_lmax", failure["final_lmax"], "1"),
        (
            "pilot_wall_limit_fraction",
            float(final["case_wall_seconds"]) / float(final["wall_limit_seconds"]),
            "1",
        ),
        (
            "pilot_memory_limit_fraction",
            max(int(row["peak_rss_bytes"]) for row in performance)
            / float(final["memory_limit_bytes"]),
            "1",
        ),
    )
    derived_rows = [
        {
            "schema_version": "1.0.0",
            "analysis_id": "p1_5_resource_derivation_v1",
            "created_utc": first["completed_utc"],
            "git_commit": first["git_commit"],
            "source_path": str(PILOT_OUTPUT_RELATIVE / "data_raw.csv"),
            "source_sha256": source_hash,
            "campaign_id": first["campaign_id"],
            "case_id": first["case_id"],
            "classification": first["classification"],
            "eligible": eligible,
            "inclusion_rule": "resource_evidence_only;exclude_from_scientific_tables",
            "model": "B_E",
            "reference_model": "",
            "metric": metric,
            "value": value,
            "unit": unit,
            "applicable": True,
            "reason": reason,
        }
        for metric, value, unit in metrics
    ]
    plot_rows: list[dict[str, object]] = []
    series = (
        ("order_wall_seconds", "order_wall_seconds", "s", "circle", "#0072B2"),
        ("peak_rss_bytes", "peak_rss_bytes", "byte", "square", "#D55E00"),
    )
    for series_id, column, unit, marker, color in series:
        for point_order, row in enumerate(performance, start=1):
            plot_rows.append(
                {
                    "schema_version": "1.0.0",
                    "figure_id": "p1_5_resource_pilot",
                    "panel_id": "resource_profile",
                    "series_id": series_id,
                    "point_order": point_order,
                    "case_id": row["case_id"],
                    "x_name": "lmax",
                    "x_value": row["lmax"],
                    "x_unit": "1",
                    "y_name": column,
                    "y_value": row[column],
                    "y_unit": unit,
                    "xerr_low": "",
                    "xerr_high": "",
                    "yerr_low": "",
                    "yerr_high": "",
                    "marker": marker,
                    "color": color,
                    "linestyle": "solid",
                    "label": series_id,
                    "eligible": eligible,
                    "annotation": reason,
                }
            )
    return {
        "data_derived.csv": _csv_bytes(_DERIVED_FIELDS, derived_rows),
        "data_plot.csv": _csv_bytes(_PLOT_FIELDS, plot_rows),
    }


def _write_bytes_fsynced(path: Path, payload: bytes) -> None:
    with path.open("xb") as stream:
        stream.write(payload)
        stream.flush()
        os.fsync(stream.fileno())


def _publish_directory_atomic(
    output_directory: Path,
    artifacts: Mapping[str, bytes],
) -> None:
    if output_directory.exists():
        raise FileExistsError(
            f"pilot output already exists; second execution refused: {output_directory}"
        )
    if set(artifacts) != set(_ARTIFACT_NAMES):
        raise PilotExecutionError("pilot publication requires exactly five artifacts")
    parent = output_directory.parent
    parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=".p1-5-pilot-", dir=parent))
    published = False
    try:
        for name in _ARTIFACT_NAMES:
            _write_bytes_fsynced(temporary / name, artifacts[name])
        os.rename(temporary, output_directory)
        published = True
        directory_fd = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        if not published and temporary.exists():
            shutil.rmtree(temporary)


def execute_p1_5_pilot(
    root: str | Path,
    *,
    source_commit: str,
    command: str,
    solver: Callable[..., ModelENodalResult] = solve_model_e_nodal,
    clock: Callable[[], float] = time.perf_counter,
    utc_now: Callable[[], datetime] = _utc_now,
    rss_reader: Callable[[], int] = _peak_rss_bytes,
    environment: Mapping[str, str] | None = None,
    enforce_resource_limits: bool = True,
) -> PilotExecutionSummary:
    """Execute and atomically publish the single frozen P1.5 attempt."""

    if not callable(solver):
        raise TypeError("solver must be callable")
    if len(source_commit) != 40 or any(
        character not in "0123456789abcdef" for character in source_commit
    ):
        raise ValueError("source_commit must be a lowercase 40-character SHA")
    if not command.strip():
        raise ValueError("command must be non-empty")
    configuration = load_p1_5_configuration(root)
    if configuration.output_directory.exists():
        raise FileExistsError(
            "pilot output exists; overwrite and second execution are forbidden"
        )
    pilot = configuration.pilot_manifest
    resources = pilot["resources"]
    parameters = configuration.case["parameters"]
    environment_payload = dict(
        _default_environment() if environment is None else environment
    )
    environment_json = json.dumps(
        environment_payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    )
    solve_id = sha256(
        "|".join(
            (
                PILOT_CAMPAIGN_ID,
                PILOT_CASE_ID,
                source_commit,
                pilot["provenance"]["manifest_sha256"],
            )
        ).encode("ascii")
    ).hexdigest()
    started_utc = _utc_text(utc_now())
    state = _ResourceState()
    timed_solver = _TimedSolver(
        solver,
        clock=clock,
        rss_reader=rss_reader,
        state=state,
        memory_limit_bytes=int(resources["peak_rss_bytes_per_case"]),
    )
    result: ModelBEResult | None = None
    runner_error: str | None = None
    run_started = clock()
    try:
        context = (
            _resource_limits(
                int(resources["wall_seconds_per_case"]),
                int(resources["peak_rss_bytes_per_case"]),
                state,
            )
            if enforce_resource_limits
            else _null_context()
        )
        with context:
            result = solve_model_be_nodal(
                configuration.positions_xyz,
                parameters["k_rad_m"],
                pilot["physical"]["radius_m"],
                pilot["physical"]["energy_density_j_m3"],
                parameters["f0"],
                parameters["f1"],
                lmax_min=pilot["numerical"]["lmax_min"],
                lmax_max=pilot["numerical"]["lmax_max"],
                minimum_stop_lmax=pilot["numerical"]["minimum_stop_lmax"],
                convergence_tolerance=pilot["numerical"]["convergence_tolerance"],
                solver=timed_solver,
            )
    except (_PilotWallTimeout, MemoryError, OSError, RuntimeError) as exc:
        runner_error = f"{type(exc).__name__}: {exc}"
    total_wall_seconds = max(0.0, clock() - run_started)
    completed_utc = _utc_text(utc_now())
    peak_rss_bytes = max(
        [int(rss_reader()), *(attempt.peak_rss_bytes for attempt in timed_solver.attempts)]
    )
    if peak_rss_bytes > int(resources["peak_rss_bytes_per_case"]):
        state.memory_limit_exceeded = True
    decision, stop_reason, final_reason = _outcome(
        result,
        timed_solver,
        state,
        total_wall_seconds,
        int(resources["wall_seconds_per_case"]),
    )
    if runner_error is not None and final_reason is None:
        final_reason = runner_error
    primary = _build_primary_tables(
        configuration,
        timed_solver,
        result,
        solve_id=solve_id,
        started_utc=started_utc,
        completed_utc=completed_utc,
        source_commit=source_commit,
        command=command,
        environment_json=environment_json,
        total_wall_seconds=total_wall_seconds,
        state=state,
        decision=decision,
        stop_reason=stop_reason,
        final_reason=final_reason,
    )
    derived = derive_p1_5_artifacts(
        primary["data_raw.csv"],
        primary["performance.csv"],
        primary["failures.csv"],
    )
    artifacts = {**primary, **derived}
    _publish_directory_atomic(configuration.output_directory, artifacts)
    hashes = {
        str(PILOT_OUTPUT_RELATIVE / name): sha256(payload).hexdigest()
        for name, payload in artifacts.items()
    }
    pair = result.pair_ledger[0] if result is not None else None
    return PilotExecutionSummary(
        decision=decision,
        stop_reason=stop_reason,
        eligible=False if result is None else result.eligible,
        converged=False if pair is None else pair.converged,
        attempted_lmax=tuple(attempt.lmax for attempt in timed_solver.attempts),
        evaluated_lmax=tuple() if pair is None else pair.evaluated_lmax,
        final_lmax=None if pair is None else pair.final_lmax,
        total_wall_seconds=total_wall_seconds,
        peak_rss_bytes=peak_rss_bytes,
        artifact_sha256=hashes,
    )


class _null_context:
    def __enter__(self) -> None:
        return None

    def __exit__(
        self,
        _exc_type: type[BaseException] | None,
        _exc_value: BaseException | None,
        _traceback: TracebackType | None,
    ) -> None:
        return None


def verify_p1_5_derivations(root: str | Path) -> Mapping[str, str]:
    """Regenerate derived tables twice and verify every published byte."""

    configuration = load_p1_5_configuration(root)
    output = configuration.output_directory
    observed_names = {path.name for path in output.iterdir() if path.is_file()}
    if observed_names != set(_ARTIFACT_NAMES):
        raise PilotExecutionError("published P1.5 artifact set is incomplete")
    raw = (output / "data_raw.csv").read_bytes()
    performance = (output / "performance.csv").read_bytes()
    failures = (output / "failures.csv").read_bytes()
    first = derive_p1_5_artifacts(raw, performance, failures)
    second = derive_p1_5_artifacts(raw, performance, failures)
    if first != second:
        raise PilotExecutionError("P1.5 derived regeneration is not deterministic")
    for name, payload in first.items():
        if (output / name).read_bytes() != payload:
            raise PilotExecutionError(f"published {name} differs from regeneration")
    return {
        str(PILOT_OUTPUT_RELATIVE / name): sha256(
            (output / name).read_bytes()
        ).hexdigest()
        for name in _ARTIFACT_NAMES
    }
