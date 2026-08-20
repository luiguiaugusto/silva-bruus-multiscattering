"""Single-attempt, checkpointed runner for the frozen P1.6 campaign."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import resource
import tempfile
import time
from typing import Any, Callable, Mapping

import numpy as np

from .model_be import solve_model_be_nodal
from .model_e import (
    ModelENodalResult,
    evaluate_model_e_numerical_diagnostics,
    solve_model_e_nodal,
)
from .p1_campaign_artifacts import ARTIFACT_PATHS, CHANNELS
from .p1_pilot import _PilotWallTimeout, _ResourceState, _resource_limits
from .paper_pipeline import (
    P1_FROZEN_MANIFEST_SHA256,
    validate_executable_manifest_file,
    validate_manifest_file,
)
from .silva_bruus import nodal_pair_forces


CAMPAIGN_ID = "p1_dimer_confirmatory"
PILOT_ID = "p1_dimer_resource_pilot"
MANIFEST_RELATIVE = Path("campaigns/p1/campaign_manifest.yaml")
PILOT_MANIFEST_RELATIVE = Path("campaigns/p1/pilot_manifest.yaml")
DEFAULT_STATE_RELATIVE = Path("campaigns/p1/.p1_6_checkpoint")


class CampaignExecutionError(RuntimeError):
    """Raised when the P1.6 single-attempt contract would be violated."""


class CampaignCaseTimeout(TimeoutError):
    """A case-local frozen wall limit was reached."""


class CampaignGlobalTimeout(CampaignCaseTimeout):
    """The remaining frozen campaign wall budget was reached inside a case."""


@dataclass(frozen=True)
class CampaignConfiguration:
    """Validated P1.6 configuration, still containing no responses."""

    root: Path
    manifest: Mapping[str, Any]
    cases: tuple[Mapping[str, Any], ...]
    state_directory: Path


@dataclass(frozen=True)
class CampaignRunSummary:
    """Checkpoint state after one invocation of the campaign orchestrator."""

    attempted_this_run: tuple[str, ...]
    completed_count: int
    interrupted_count: int
    never_started_count: int
    accumulated_wall_seconds: float
    closed: bool
    stop_reason: str


CaseExecutor = Callable[[Mapping[str, Any], Mapping[str, Any]], Mapping[str, Any]]


def _utc_text(value: datetime) -> str:
    if value.tzinfo is None or value.utcoffset() is None:
        raise ValueError("UTC clock must return a timezone-aware datetime")
    return value.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.%fZ")


def _json_bytes(value: object) -> bytes:
    return (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(_json_bytes(value))
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _require_exact_manifest(manifest: Mapping[str, Any]) -> None:
    numerical = manifest["numerical"]
    expected_numerical = {
        "model": "E",
        "lmax_min": 2,
        "lmax_max": 21,
        "minimum_stop_lmax": 5,
        "convergence_tolerance": 1.0e-5,
        "consecutive_passes": 2,
        "required_channels": list(CHANNELS),
        "require_all_applicable_channels": True,
        "failure_policy": (
            "retain every attempt; no imputation or physics/tolerance/solver "
            "retry; record controlled stage and reason"
        ),
    }
    if numerical != expected_numerical:
        raise CampaignExecutionError("P1.6 numerical policy differs from P1.4")
    resources = manifest["resources"]
    if resources != {
        "worker_count": 1,
        "blas_threads": 1,
        "peak_rss_bytes_per_case": 4 * 1024**3,
        "wall_seconds_per_case": 1800,
        "wall_seconds_campaign": 64800,
        "limits_status": "frozen",
    }:
        raise CampaignExecutionError("P1.6 resource limits differ from P1.6A")
    expected_outputs = {
        "raw": ARTIFACT_PATHS["data_raw.csv"],
        "derived": ARTIFACT_PATHS["data_derived.csv"],
        "plot_ready": ARTIFACT_PATHS["data_plot.csv"],
        "failure_log": ARTIFACT_PATHS["failures.csv"],
    }
    if manifest["outputs"] != expected_outputs:
        raise CampaignExecutionError("P1.6 output paths differ from P1.4")
    cases = list(manifest["cases"])
    if len(cases) != 102:
        raise CampaignExecutionError("P1.6 requires exactly 102 cases")
    if [case["case_order"] for case in cases] != list(range(1, 103)):
        raise CampaignExecutionError("P1.6 case order must be exactly 1..102")
    if not all(case["enabled"] for case in cases):
        raise CampaignExecutionError("all 102 P1.6 cases must be enabled")
    primary = [
        case
        for case in cases
        if case["parameters"]["evidence_role"] == "primary"
    ]
    audits = [
        case
        for case in cases
        if case["parameters"]["evidence_role"] == "rotational_audit"
    ]
    if len(primary) != 96 or len(audits) != 6:
        raise CampaignExecutionError("P1.6 requires 96 primaries and six audits")
    if not all(case["parameters"]["include_in_scientific_tables"] for case in primary):
        raise CampaignExecutionError("all P1.6 primaries must enter scientific tables")
    if any(case["parameters"]["include_in_scientific_tables"] for case in audits):
        raise CampaignExecutionError("P1.6 audits must remain outside scientific tables")


def load_p1_6_configuration(
    root: str | Path,
    *,
    state_directory: str | Path | None = None,
) -> CampaignConfiguration:
    """Validate the enabled response-blind manifest without solving a case."""

    repository = Path(root).resolve()
    manifest = validate_executable_manifest_file(
        repository / MANIFEST_RELATIVE,
        expected_campaign_id=CAMPAIGN_ID,
    )
    _require_exact_manifest(manifest)
    if manifest["provenance"]["manifest_sha256"] != P1_FROZEN_MANIFEST_SHA256[
        CAMPAIGN_ID
    ]:
        raise CampaignExecutionError("confirmatory hash differs from public lock")
    pilot = validate_manifest_file(
        repository / PILOT_MANIFEST_RELATIVE,
        kind="campaign",
    )
    if pilot["campaign_id"] != PILOT_ID:
        raise CampaignExecutionError("unexpected P1.5 pilot identity")
    if pilot["provenance"]["manifest_sha256"] != P1_FROZEN_MANIFEST_SHA256[
        PILOT_ID
    ]:
        raise CampaignExecutionError("P1.5 pilot manifest changed after execution")
    if [case["case_id"] for case in pilot["cases"] if case["enabled"]] != [
        "p1_pilot_rigid_ka010_d0210_t000"
    ]:
        raise CampaignExecutionError("P1.5 pilot enablement changed")
    directory = (
        repository / DEFAULT_STATE_RELATIVE
        if state_directory is None
        else Path(state_directory).resolve()
    )
    return CampaignConfiguration(
        root=repository,
        manifest=manifest,
        cases=tuple(manifest["cases"]),
        state_directory=directory,
    )


def _initial_ledger(configuration: CampaignConfiguration, created_utc: str) -> dict[str, Any]:
    manifest = configuration.manifest
    return {
        "schema_version": "1.0.0",
        "campaign_id": manifest["campaign_id"],
        "manifest_sha256": manifest["provenance"]["manifest_sha256"],
        "created_utc": created_utc,
        "updated_utc": created_utc,
        "accumulated_wall_seconds": 0.0,
        "closed": False,
        "stop_reason": "not_started",
        "cases": [
            {
                "case_order": int(case["case_order"]),
                "case_id": case["case_id"],
                "parameters": case["parameters"],
                "state": "never_started",
                "attempt_count": 0,
                "started_utc": None,
                "completed_utc": None,
                "failure_stage": None,
                "failure_reason": None,
            }
            for case in configuration.cases
        ],
    }


def _validate_ledger(
    configuration: CampaignConfiguration, ledger: Mapping[str, Any]
) -> None:
    if ledger.get("campaign_id") != CAMPAIGN_ID:
        raise CampaignExecutionError("checkpoint campaign identity mismatch")
    if ledger.get("manifest_sha256") != configuration.manifest["provenance"][
        "manifest_sha256"
    ]:
        raise CampaignExecutionError("checkpoint manifest hash mismatch")
    entries = ledger.get("cases")
    if not isinstance(entries, list) or len(entries) != 102:
        raise CampaignExecutionError("checkpoint must retain all 102 cases")
    expected = [
        (case["case_order"], case["case_id"]) for case in configuration.cases
    ]
    observed = [(entry.get("case_order"), entry.get("case_id")) for entry in entries]
    if observed != expected:
        raise CampaignExecutionError("checkpoint IDs or order differ from manifest")
    for entry in entries:
        count = entry.get("attempt_count")
        state = entry.get("state")
        if count not in (0, 1):
            raise CampaignExecutionError("a P1.6 case may have at most one attempt")
        if state not in {"never_started", "started", "completed", "interrupted"}:
            raise CampaignExecutionError("checkpoint contains an unknown state")
        if state == "never_started" and count != 0:
            raise CampaignExecutionError("never-started case has an attempt")
        if state != "never_started" and count != 1:
            raise CampaignExecutionError("started case must retain one attempt")


def _read_ledger(configuration: CampaignConfiguration, now_utc: str) -> dict[str, Any]:
    ledger_path = configuration.state_directory / "campaign_ledger.json"
    if not ledger_path.exists():
        for relative in ARTIFACT_PATHS.values():
            if (configuration.root / relative).exists():
                raise FileExistsError(
                    "campaign output exists; overwrite or a second campaign is forbidden"
                )
        ledger = _initial_ledger(configuration, now_utc)
        _atomic_json(ledger_path, ledger)
        return ledger
    try:
        ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise CampaignExecutionError(f"cannot load campaign checkpoint: {exc}") from exc
    _validate_ledger(configuration, ledger)
    recovered = False
    for entry in ledger["cases"]:
        if entry["state"] == "started":
            checkpoint_path = _checkpoint_path(
                configuration, int(entry["case_order"])
            )
            checkpoint: Mapping[str, Any] | None = None
            if checkpoint_path.exists():
                try:
                    loaded = json.loads(checkpoint_path.read_text(encoding="utf-8"))
                    if isinstance(loaded, Mapping):
                        checkpoint = loaded
                except (OSError, json.JSONDecodeError):
                    checkpoint = None
            if checkpoint is not None and checkpoint.get("state") in {
                "completed",
                "interrupted",
            }:
                for field in (
                    "state",
                    "completed_utc",
                    "failure_stage",
                    "failure_reason",
                ):
                    entry[field] = checkpoint.get(field)
                if checkpoint.get("state") == "completed" and isinstance(
                    checkpoint.get("outcome"), Mapping
                ):
                    ledger["accumulated_wall_seconds"] = float(
                        ledger["accumulated_wall_seconds"]
                    ) + float(checkpoint["outcome"].get("wall_seconds", 0.0))
                elif checkpoint.get("failure_stage") == "timeout":
                    ledger["accumulated_wall_seconds"] = float(
                        ledger["accumulated_wall_seconds"]
                    ) + float(
                        configuration.manifest["resources"][
                            "wall_seconds_per_case"
                        ]
                    )
            else:
                entry["state"] = "interrupted"
                entry["completed_utc"] = now_utc
                entry["failure_stage"] = "interrupted"
                entry["failure_reason"] = "previous_process_interrupted_after_start"
                _atomic_json(
                    checkpoint_path,
                    {**entry, "outcome": None},
                )
            recovered = True
    if recovered:
        ledger["updated_utc"] = now_utc
        ledger["stop_reason"] = "resumed_after_interruption"
        global_limit = float(
            configuration.manifest["resources"]["wall_seconds_campaign"]
        )
        if float(ledger["accumulated_wall_seconds"]) >= global_limit:
            ledger["accumulated_wall_seconds"] = global_limit
            ledger["closed"] = True
            ledger["stop_reason"] = "global_wall_limit_exhausted"
            for remaining in ledger["cases"]:
                if remaining["state"] == "never_started":
                    remaining["failure_stage"] = "global_limit"
                    remaining["failure_reason"] = (
                        "wall_seconds_campaign_exhausted_before_attempt"
                    )
        _atomic_json(ledger_path, ledger)
    return ledger


def _finite_number(value: object, name: str, *, minimum: float | None = None) -> float:
    if isinstance(value, bool):
        raise CampaignExecutionError(f"{name} must be a finite number")
    result = float(value)
    if not math.isfinite(result) or (minimum is not None and result < minimum):
        raise CampaignExecutionError(f"{name} must be finite and at least {minimum}")
    return result


def _forces(value: object, name: str, *, required: bool) -> list[list[float]] | None:
    if value is None and not required:
        return None
    array = np.asarray(value, dtype=float)
    if array.shape != (2, 3) or not np.all(np.isfinite(array)):
        raise CampaignExecutionError(f"{name} must be a finite (2, 3) array")
    return array.tolist()


def _normalize_outcome(
    value: Mapping[str, Any], manifest: Mapping[str, Any]
) -> dict[str, Any]:
    attempted = tuple(int(order) for order in value.get("attempted_lmax", ()))
    evaluated = tuple(int(order) for order in value.get("evaluated_lmax", ()))
    allowed = tuple(range(2, 22))
    if not attempted or attempted != tuple(sorted(set(attempted))):
        raise CampaignExecutionError("attempted_lmax must be ordered and unique")
    if any(order not in allowed for order in attempted):
        raise CampaignExecutionError("attempted_lmax is outside L=2..21")
    if evaluated != tuple(sorted(set(evaluated))) or any(
        order not in attempted for order in evaluated
    ):
        raise CampaignExecutionError("evaluated_lmax must be an ordered subset")
    final_lmax = value.get("final_lmax")
    if final_lmax is not None:
        final_lmax = int(final_lmax)
    if final_lmax != (evaluated[-1] if evaluated else None):
        raise CampaignExecutionError("final_lmax must equal the last evaluated order")
    solve_count = int(value.get("model_e_solve_count", -1))
    if solve_count != len(evaluated):
        raise CampaignExecutionError(
            "each evaluated dimer order must be solved exactly once"
        )
    orders = list(value.get("orders", ()))
    if [int(order.get("lmax", -1)) for order in orders] != list(evaluated):
        raise CampaignExecutionError("order records must match evaluated_lmax")
    normalized_orders: list[dict[str, Any]] = []
    for order in orders:
        channels = order.get("channels")
        if not isinstance(channels, Mapping) or set(channels) != set(CHANNELS):
            raise CampaignExecutionError("each order must store the four frozen channels")
        normalized_channels: dict[str, Any] = {}
        for name in CHANNELS:
            channel = channels[name]
            normalized_channels[name] = {
                "forces_xyz": _forces(
                    channel.get("forces_xyz"), f"{name}.forces_xyz", required=True
                ),
                "successive_change": _finite_number(
                    channel.get("successive_change"),
                    f"{name}.successive_change",
                    minimum=0.0,
                ),
                "absolute_change": _finite_number(
                    channel.get("absolute_change"),
                    f"{name}.absolute_change",
                    minimum=0.0,
                ),
                "applicable": bool(channel.get("applicable")),
                "confirmed": bool(channel.get("confirmed")),
                "confirmation_lmax": channel.get("confirmation_lmax"),
            }
        diagnostics = dict(order.get("diagnostics", {}))
        normalized_orders.append(
            {
                "lmax": int(order["lmax"]),
                "wall_seconds": _finite_number(
                    order.get("wall_seconds", 0.0),
                    "order.wall_seconds",
                    minimum=0.0,
                ),
                "peak_rss_bytes": int(order.get("peak_rss_bytes", 0)),
                "diagnostics": diagnostics,
                "channels": normalized_channels,
            }
        )
    converged = bool(value.get("converged"))
    eligible = bool(value.get("eligible"))
    model_a = _forces(value.get("model_a_forces_xyz"), "Model A forces", required=True)
    model_be = _forces(
        value.get("model_be_forces_xyz"), "Model B_E forces", required=eligible
    )
    model_e = _forces(
        value.get("model_e_forces_xyz"), "Model E forces", required=eligible
    )
    failure_stage = value.get("failure_stage")
    failure_reason = value.get("failure_reason")
    if eligible and (
        not converged
        or failure_stage is not None
        or failure_reason is not None
        or model_be is None
        or model_e is None
    ):
        raise CampaignExecutionError("eligible outcome violates convergence contract")
    if not eligible and not failure_reason:
        raise CampaignExecutionError("ineligible outcome requires a reason")
    normalized = {
        "attempted_lmax": list(attempted),
        "evaluated_lmax": list(evaluated),
        "final_lmax": final_lmax,
        "model_e_solve_count": solve_count,
        "orders": normalized_orders,
        "model_a_forces_xyz": model_a,
        "model_be_forces_xyz": model_be,
        "model_e_forces_xyz": model_e,
        "converged": converged,
        "eligible": eligible,
        "failure_stage": failure_stage,
        "failure_reason": failure_reason,
        "wall_seconds": _finite_number(
            value.get("wall_seconds"), "wall_seconds", minimum=0.0
        ),
        "peak_rss_bytes": int(value.get("peak_rss_bytes")),
    }
    if normalized["peak_rss_bytes"] < 0:
        raise CampaignExecutionError("peak_rss_bytes must be non-negative")
    resources = manifest["resources"]
    if normalized["wall_seconds"] > float(resources["wall_seconds_per_case"]):
        normalized["eligible"] = False
        normalized["failure_stage"] = "timeout"
        normalized["failure_reason"] = "wall_seconds_per_case_exceeded"
    if normalized["peak_rss_bytes"] > int(resources["peak_rss_bytes_per_case"]):
        normalized["eligible"] = False
        normalized["failure_stage"] = "memory"
        normalized["failure_reason"] = "peak_rss_bytes_per_case_exceeded"
    _json_bytes(normalized)
    return normalized


def _checkpoint_path(configuration: CampaignConfiguration, order: int) -> Path:
    return configuration.state_directory / "cases" / f"{order:03d}.json"


def run_p1_6_campaign(
    root: str | Path,
    *,
    executor: CaseExecutor,
    state_directory: str | Path | None = None,
    utc_now: Callable[[], datetime] = lambda: datetime.now(timezone.utc),
    monotonic: Callable[[], float] = time.monotonic,
    max_new_cases: int | None = None,
) -> CampaignRunSummary:
    """Attempt each enabled case once, checkpointing before and after it.

    ``max_new_cases`` exists only to make interruption/resumption tests cheap;
    the production CLI never exposes it.
    """

    if not callable(executor):
        raise TypeError("executor must be callable")
    if max_new_cases is not None and max_new_cases < 1:
        raise ValueError("max_new_cases must be positive")
    configuration = load_p1_6_configuration(root, state_directory=state_directory)
    now = _utc_text(utc_now())
    ledger = _read_ledger(configuration, now)
    _validate_ledger(configuration, ledger)
    never_started = [entry for entry in ledger["cases"] if entry["state"] == "never_started"]
    if ledger["closed"] or not never_started:
        raise CampaignExecutionError(
            "campaign is closed; overwrite, retry and second execution are forbidden"
        )
    attempted_this_run: list[str] = []
    global_limit = float(configuration.manifest["resources"]["wall_seconds_campaign"])
    ledger_path = configuration.state_directory / "campaign_ledger.json"

    for case, entry in zip(configuration.cases, ledger["cases"]):
        if entry["state"] != "never_started":
            continue
        if max_new_cases is not None and len(attempted_this_run) >= max_new_cases:
            ledger["stop_reason"] = "invocation_case_limit"
            break
        if float(ledger["accumulated_wall_seconds"]) >= global_limit:
            ledger["closed"] = True
            ledger["stop_reason"] = "global_wall_limit_exhausted"
            for remaining in ledger["cases"]:
                if remaining["state"] == "never_started":
                    remaining["failure_stage"] = "global_limit"
                    remaining["failure_reason"] = (
                        "wall_seconds_campaign_exhausted_before_attempt"
                    )
            break
        global_remaining = global_limit - float(
            ledger["accumulated_wall_seconds"]
        )
        effective_wall = min(
            float(configuration.manifest["resources"]["wall_seconds_per_case"]),
            global_remaining,
        )
        runtime_manifest = {
            **configuration.manifest,
            "resources": {
                **configuration.manifest["resources"],
                "wall_seconds_per_case": effective_wall,
            },
        }
        started_utc = _utc_text(utc_now())
        case_wall_started = monotonic()
        entry.update(
            {
                "state": "started",
                "attempt_count": 1,
                "started_utc": started_utc,
                "completed_utc": None,
                "failure_stage": None,
                "failure_reason": None,
            }
        )
        ledger["updated_utc"] = started_utc
        ledger["stop_reason"] = "case_started"
        _atomic_json(ledger_path, ledger)
        started_record = {**entry, "outcome": None}
        _atomic_json(
            _checkpoint_path(configuration, int(entry["case_order"])),
            started_record,
        )
        attempted_this_run.append(entry["case_id"])

        try:
            raw_outcome = executor(case, runtime_manifest)
            outcome = _normalize_outcome(raw_outcome, configuration.manifest)
            outcome["wall_seconds"] = max(
                float(outcome["wall_seconds"]),
                max(0.0, monotonic() - case_wall_started),
            )
            global_exhausted = outcome["wall_seconds"] > global_remaining
            if global_exhausted:
                outcome["eligible"] = False
                outcome["failure_stage"] = "global_limit"
                outcome["failure_reason"] = (
                    "wall_seconds_campaign_exhausted_during_attempt"
                )
        except Exception as exc:
            completed_utc = _utc_text(utc_now())
            measured_wall = max(0.0, monotonic() - case_wall_started)
            if isinstance(exc, CampaignGlobalTimeout):
                failure_stage = "global_limit"
            elif isinstance(exc, CampaignCaseTimeout):
                failure_stage = "timeout"
            elif isinstance(exc, MemoryError):
                failure_stage = "memory"
            elif isinstance(exc, CampaignExecutionError):
                failure_stage = "contract"
            else:
                failure_stage = "interrupted"
            entry.update(
                {
                    "state": "interrupted",
                    "completed_utc": completed_utc,
                    "failure_stage": failure_stage,
                    "failure_reason": f"{type(exc).__name__}: {exc}",
                }
            )
            if failure_stage == "global_limit":
                ledger["accumulated_wall_seconds"] = global_limit
            elif failure_stage == "timeout":
                ledger["accumulated_wall_seconds"] = float(
                    ledger["accumulated_wall_seconds"]
                ) + float(
                    configuration.manifest["resources"]["wall_seconds_per_case"]
                )
            else:
                ledger["accumulated_wall_seconds"] = float(
                    ledger["accumulated_wall_seconds"]
                ) + measured_wall
            _atomic_json(
                _checkpoint_path(configuration, int(entry["case_order"])),
                {**entry, "outcome": None},
            )
            ledger["updated_utc"] = completed_utc
            ledger["stop_reason"] = (
                "global_wall_limit_exhausted"
                if failure_stage == "global_limit"
                else "case_failure_continued"
            )
            if failure_stage == "global_limit":
                ledger["closed"] = True
                for remaining in ledger["cases"]:
                    if remaining["state"] == "never_started":
                        remaining["failure_stage"] = "global_limit"
                        remaining["failure_reason"] = (
                            "wall_seconds_campaign_exhausted_before_attempt"
                        )
            _atomic_json(ledger_path, ledger)
            if failure_stage == "global_limit":
                break
            continue

        completed_utc = _utc_text(utc_now())
        entry.update(
            {
                "state": "completed",
                "completed_utc": completed_utc,
                "failure_stage": outcome["failure_stage"],
                "failure_reason": outcome["failure_reason"],
            }
        )
        ledger["accumulated_wall_seconds"] = min(
            global_limit,
            float(ledger["accumulated_wall_seconds"])
            + float(outcome["wall_seconds"]),
        )
        _atomic_json(
            _checkpoint_path(configuration, int(entry["case_order"])),
            {**entry, "outcome": outcome},
        )
        ledger["updated_utc"] = completed_utc
        ledger["stop_reason"] = (
            "global_wall_limit_exhausted"
            if global_exhausted
            else "case_completed"
        )
        if global_exhausted:
            ledger["closed"] = True
            for remaining in ledger["cases"]:
                if remaining["state"] == "never_started":
                    remaining["failure_stage"] = "global_limit"
                    remaining["failure_reason"] = (
                        "wall_seconds_campaign_exhausted_before_attempt"
                    )
        _atomic_json(ledger_path, ledger)
        if global_exhausted:
            break

    terminal = all(entry["state"] != "never_started" for entry in ledger["cases"])
    if terminal:
        ledger["closed"] = True
        ledger["stop_reason"] = "all_cases_attempted"
    ledger["updated_utc"] = _utc_text(utc_now())
    _atomic_json(ledger_path, ledger)
    return CampaignRunSummary(
        attempted_this_run=tuple(attempted_this_run),
        completed_count=sum(entry["state"] == "completed" for entry in ledger["cases"]),
        interrupted_count=sum(entry["state"] == "interrupted" for entry in ledger["cases"]),
        never_started_count=sum(entry["state"] == "never_started" for entry in ledger["cases"]),
        accumulated_wall_seconds=float(ledger["accumulated_wall_seconds"]),
        closed=bool(ledger["closed"]),
        stop_reason=str(ledger["stop_reason"]),
    )


def _peak_rss_bytes() -> int:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value if os.uname().sysname == "Darwin" else value * 1024


def execute_model_e_case(
    case: Mapping[str, Any],
    manifest: Mapping[str, Any],
    *,
    solver: Callable[..., ModelENodalResult] = solve_model_e_nodal,
    clock: Callable[[], float] = time.perf_counter,
    rss_reader: Callable[[], int] = _peak_rss_bytes,
) -> Mapping[str, Any]:
    """Execute one dimer with one Model-E solve per attempted order.

    P1.6A does not call this adapter.  It is the production adapter frozen for
    P1.6B and is testable by injecting a fake Model-E solver.
    """

    parameters = case["parameters"]
    radius = float(manifest["physical"]["radius_m"])
    energy = float(manifest["physical"]["energy_density_j_m3"])
    distance = radius * float(parameters["distance_ratio"])
    theta = float(parameters["theta_rad"])
    displacement = np.array(
        [distance * np.cos(theta), distance * np.sin(theta), 0.0], dtype=float
    )
    positions = np.vstack((-0.5 * displacement, 0.5 * displacement))
    results: list[ModelENodalResult] = []
    order_wall: list[float] = []
    order_rss: list[int] = []

    def collecting_solver(*args: object, **kwargs: object) -> ModelENodalResult:
        started = clock()
        result = solver(*args, **kwargs)
        order_wall.append(max(0.0, clock() - started))
        order_rss.append(int(rss_reader()))
        results.append(result)
        return result

    started = clock()
    be_result = solve_model_be_nodal(
        positions,
        parameters["k_rad_m"],
        radius,
        energy,
        parameters["f0"],
        parameters["f1"],
        lmax_min=manifest["numerical"]["lmax_min"],
        lmax_max=manifest["numerical"]["lmax_max"],
        minimum_stop_lmax=manifest["numerical"]["minimum_stop_lmax"],
        convergence_tolerance=manifest["numerical"]["convergence_tolerance"],
        solver=collecting_solver,
    )
    wall_seconds = max(0.0, clock() - started)
    pair = be_result.pair_ledger[0]
    convergence = {summary.channel: summary for summary in pair.convergence}
    orders: list[dict[str, Any]] = []
    scale = radius**2 * energy
    if scale <= 0.0:
        raise CampaignExecutionError("P1.6 force normalization scale must be positive")
    channel_attributes = {
        "total": "total_forces_xyz",
        "interaction": "interaction_forces_xyz",
        "external_scattered": "external_scattered_forces_xyz",
        "scattered_scattered": "scattered_scattered_forces_xyz",
    }
    for index, result in enumerate(results):
        diagnostics = evaluate_model_e_numerical_diagnostics(result)
        diagnostic_payload = asdict(diagnostics)
        diagnostic_payload.update(
            {
                "full_modes_per_particle": len(result.solution.modes),
                "active_modes_per_particle": len(result.solution.active_modes),
                "system_dimension": int(
                    result.solution.balanced_system_matrix.shape[0]
                ),
            }
        )
        channel_payload: dict[str, Any] = {}
        for channel in CHANNELS:
            summary = convergence[channel]
            step = summary.history[index]
            final_window_confirmed = bool(
                index >= 2
                and summary.history[index - 1].applicable
                and step.applicable
                and summary.history[index - 1].successive_change
                <= manifest["numerical"]["convergence_tolerance"]
                and step.successive_change
                <= manifest["numerical"]["convergence_tolerance"]
            )
            channel_payload[channel] = {
                "forces_xyz": (
                    np.asarray(getattr(result, channel_attributes[channel]), dtype=float)
                    / scale
                ).tolist(),
                "successive_change": step.successive_change,
                "absolute_change": step.absolute_change,
                "applicable": step.applicable,
                "confirmed": final_window_confirmed,
                "confirmation_lmax": summary.confirmation_lmax,
            }
        orders.append(
            {
                "lmax": int(result.lmax),
                "wall_seconds": order_wall[index],
                "peak_rss_bytes": order_rss[index],
                "diagnostics": diagnostic_payload,
                "channels": channel_payload,
            }
        )
    force_a_xy = nodal_pair_forces(
        positions[0, :2],
        positions[1, :2],
        parameters["k_rad_m"],
        radius,
        energy,
        parameters["f1"],
    )
    force_a = np.column_stack((np.asarray(force_a_xy) / scale, np.zeros(2)))
    model_be = None if be_result.forces_xyz is None else (be_result.forces_xyz / scale).tolist()
    model_e = (
        None
        if pair.interaction_forces_xyz is None
        else (pair.interaction_forces_xyz / scale).tolist()
    )
    return {
        "attempted_lmax": list(pair.attempted_lmax),
        "evaluated_lmax": list(pair.evaluated_lmax),
        "final_lmax": pair.final_lmax,
        "model_e_solve_count": len(results),
        "orders": orders,
        "model_a_forces_xyz": force_a.tolist(),
        "model_be_forces_xyz": model_be,
        "model_e_forces_xyz": model_e,
        "converged": pair.converged,
        "eligible": be_result.eligible,
        "failure_stage": be_result.failure_stage,
        "failure_reason": be_result.failure_reason,
        "wall_seconds": wall_seconds,
        "peak_rss_bytes": max([int(rss_reader()), *order_rss]),
    }


def execute_model_e_case_with_limits(
    case: Mapping[str, Any],
    manifest: Mapping[str, Any],
    *,
    executor: CaseExecutor = execute_model_e_case,
) -> Mapping[str, Any]:
    """Apply the frozen P1.5-tested wall/RLIMIT_AS guard to one P1.6 case."""

    resources = manifest["resources"]
    state = _ResourceState()
    try:
        with _resource_limits(
            int(resources["wall_seconds_per_case"]),
            int(resources["peak_rss_bytes_per_case"]),
            state,
        ):
            return executor(case, manifest)
    except _PilotWallTimeout as exc:
        if float(resources["wall_seconds_per_case"]) < 1800.0:
            raise CampaignGlobalTimeout(str(exc)) from exc
        raise CampaignCaseTimeout(str(exc)) from exc
