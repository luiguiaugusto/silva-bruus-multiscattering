"""Fake-only tests for the P1.6 runner, checkpoints, artifacts and G1."""

from __future__ import annotations

import ast
import csv
from datetime import datetime, timezone
import io
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import numpy as np
import pytest

from acoustic_ms.p1_campaign import (
    CampaignCaseTimeout,
    CampaignExecutionError,
    execute_model_e_case,
    run_p1_6_campaign,
)
from acoustic_ms.p1_campaign_artifacts import (
    G1_BUDGET,
    build_campaign_artifacts,
    evaluate_g1,
    load_checkpoint_records,
    publish_campaign_artifacts,
)


ROOT = Path(__file__).resolve().parents[1]
UTC = datetime(2026, 8, 20, 12, 0, tzinfo=timezone.utc)


def _root(tmp_path: Path, name: str) -> Path:
    root = tmp_path / name
    destination = root / "campaigns" / "p1"
    destination.mkdir(parents=True)
    for manifest in ("campaign_manifest.yaml", "pilot_manifest.yaml"):
        shutil.copy2(ROOT / "campaigns" / "p1" / manifest, destination / manifest)
    return root


def _radial(case: dict[str, object], scale: float = 1.0) -> list[list[float]]:
    theta = float(case["parameters"]["theta_rad"])
    vector = [scale * np.cos(theta), scale * np.sin(theta), 0.0]
    return [vector, [-component for component in vector]]


def _outcome(
    case: dict[str, object],
    *,
    wall_seconds: float = 1.0,
    peak_rss_bytes: int = 128 * 1024**2,
) -> dict[str, object]:
    orders = []
    changes = {2: (0.0, False), 3: (1.0e-4, True), 4: (1.0e-6, True), 5: (1.0e-7, True)}
    for order in range(2, 6):
        change, applicable = changes[order]
        channels = {}
        for index, channel in enumerate((
            "total", "interaction", "external_scattered", "scattered_scattered"
        ), start=1):
            channels[channel] = {
                "forces_xyz": _radial(case, float(index)),
                "successive_change": change,
                "absolute_change": change,
                "applicable": applicable,
                "confirmed": order == 5,
                "confirmation_lmax": 5,
            }
        orders.append(
            {
                "lmax": order,
                "wall_seconds": 0.25,
                "peak_rss_bytes": peak_rss_bytes,
                "diagnostics": {"passed": True, "finite": True},
                "channels": channels,
            }
        )
    interaction = orders[-1]["channels"]["interaction"]["forces_xyz"]
    return {
        "attempted_lmax": [2, 3, 4, 5],
        "evaluated_lmax": [2, 3, 4, 5],
        "final_lmax": 5,
        "model_e_solve_count": 4,
        "orders": orders,
        "model_a_forces_xyz": _radial(case, 1.0),
        "model_be_forces_xyz": interaction,
        "model_e_forces_xyz": interaction,
        "converged": True,
        "eligible": True,
        "failure_stage": None,
        "failure_reason": None,
        "wall_seconds": wall_seconds,
        "peak_rss_bytes": peak_rss_bytes,
    }


class FakeExecutor:
    def __init__(self, *, wall_seconds: float = 1.0) -> None:
        self.case_ids: list[str] = []
        self.wall_limits: list[float] = []
        self.wall_seconds = wall_seconds

    def __call__(self, case, manifest):
        self.case_ids.append(case["case_id"])
        self.wall_limits.append(
            float(manifest["resources"]["wall_seconds_per_case"])
        )
        return _outcome(case, wall_seconds=self.wall_seconds)


def test_102_case_order_single_attempt_g1_and_no_solver_regeneration(
    tmp_path: Path,
) -> None:
    root = _root(tmp_path, "complete")
    state = tmp_path / "state-complete"
    executor = FakeExecutor()
    summary = run_p1_6_campaign(
        root,
        executor=executor,
        state_directory=state,
        utc_now=lambda: UTC,
    )
    records = load_checkpoint_records(state)
    manifest = json.loads(
        (root / "campaigns" / "p1" / "campaign_manifest.yaml").read_text()
    )
    calls_before = list(executor.case_ids)
    first, first_gate = build_campaign_artifacts(manifest, records)
    second, second_gate = build_campaign_artifacts(manifest, records)

    assert summary.closed and summary.stop_reason == "all_cases_attempted"
    assert summary.completed_count == 102
    assert summary.interrupted_count == summary.never_started_count == 0
    assert executor.case_ids == [case["case_id"] for case in manifest["cases"]]
    assert all(record["attempt_count"] == 1 for record in records)
    assert all(record["outcome"]["model_e_solve_count"] == 4 for record in records)
    assert first == second
    assert first_gate == second_gate
    assert executor.case_ids == calls_before
    assert first_gate.gate_status == "PASS_G1"
    assert first_gate.decision == "GO_P2"
    assert first_gate.attempted_count == 102
    assert first_gate.eligible_primary_count == 96
    assert first_gate.eligible_audit_count == 6
    assert first_gate.eligible_audit_twin_pairs == 6
    assert len(first_gate.covered_strata) == 12
    assert first_gate.identity_error_max == 0.0
    assert first_gate.rotation_error_max is not None
    assert first_gate.rotation_error_max <= G1_BUDGET

    raw = list(csv.DictReader(io.StringIO(first["data_raw.csv"].decode("utf-8"))))
    derived = list(
        csv.DictReader(io.StringIO(first["data_derived.csv"].decode("utf-8")))
    )
    plot = list(csv.DictReader(io.StringIO(first["data_plot.csv"].decode("utf-8"))))
    failures = list(
        csv.DictReader(io.StringIO(first["failures.csv"].decode("utf-8")))
    )
    performance = list(
        csv.DictReader(io.StringIO(first["performance.csv"].decode("utf-8")))
    )
    assert len(raw) == 102 * (4 * 4 * 2 + 2 * 2)
    assert {row["force_channel"] for row in raw} == {
        "total", "interaction", "external_scattered", "scattered_scattered"
    }
    assert all(row["diagnostics_json"] for row in raw)
    assert sum(row["metric"] == "rotational_covariance_error" for row in derived) == 6
    assert len(plot) == 96
    assert failures == []
    assert len(performance) == 102
    assert all(row["model_e_solve_count"] == "4" for row in performance)

    published = publish_campaign_artifacts(root, first)
    assert set(published.values())
    with pytest.raises(FileExistsError, match="second publication"):
        publish_campaign_artifacts(root, first)
    with pytest.raises(CampaignExecutionError, match="second execution"):
        run_p1_6_campaign(
            root,
            executor=executor,
            state_directory=state,
            utc_now=lambda: UTC,
        )


def test_interrupted_case_is_never_retried_and_resume_uses_next_case(
    tmp_path: Path,
) -> None:
    root = _root(tmp_path, "resume")
    state = tmp_path / "state-resume"
    calls: list[int] = []

    def interrupted(case, manifest):
        del manifest
        order = int(case["case_order"])
        calls.append(order)
        if order == 2:
            raise KeyboardInterrupt("simulated process loss")
        return _outcome(case)

    with pytest.raises(KeyboardInterrupt, match="simulated process loss"):
        run_p1_6_campaign(
            root,
            executor=interrupted,
            state_directory=state,
            utc_now=lambda: UTC,
        )
    resumed_calls: list[int] = []

    def resumed(case, manifest):
        del manifest
        resumed_calls.append(int(case["case_order"]))
        return _outcome(case)

    summary = run_p1_6_campaign(
        root,
        executor=resumed,
        state_directory=state,
        utc_now=lambda: UTC,
        max_new_cases=1,
    )
    records = load_checkpoint_records(state)

    assert calls == [1, 2]
    assert resumed_calls == [3]
    assert records[1]["state"] == "interrupted"
    assert records[1]["attempt_count"] == 1
    assert records[1]["failure_reason"] == "previous_process_interrupted_after_start"
    assert summary.completed_count == 2
    assert summary.interrupted_count == 1
    assert summary.never_started_count == 99


def test_timeout_memory_and_local_failure_continue_without_retry(
    tmp_path: Path,
) -> None:
    root = _root(tmp_path, "local-failures")
    state = tmp_path / "state-local"
    calls: list[int] = []

    def executor(case, manifest):
        del manifest
        order = int(case["case_order"])
        calls.append(order)
        if order == 1:
            raise CampaignCaseTimeout("1800 second limit")
        if order == 2:
            raise MemoryError("4 GiB limit")
        return _outcome(case)

    summary = run_p1_6_campaign(
        root,
        executor=executor,
        state_directory=state,
        utc_now=lambda: UTC,
        max_new_cases=3,
    )
    records = load_checkpoint_records(state)

    assert calls == [1, 2, 3]
    assert [records[index]["state"] for index in range(3)] == [
        "interrupted", "interrupted", "completed"
    ]
    assert records[0]["failure_stage"] == "timeout"
    assert records[1]["failure_stage"] == "memory"
    assert summary.interrupted_count == 2
    assert summary.completed_count == 1


def test_returned_local_limits_and_global_limit_are_controlled(
    tmp_path: Path,
) -> None:
    root = _root(tmp_path, "limits")
    local_state = tmp_path / "state-returned-local"

    def local_executor(case, manifest):
        del manifest
        if case["case_order"] == 1:
            return _outcome(case, wall_seconds=1800.5)
        return _outcome(case, peak_rss_bytes=4 * 1024**3 + 1)

    run_p1_6_campaign(
        root,
        executor=local_executor,
        state_directory=local_state,
        utc_now=lambda: UTC,
        max_new_cases=2,
    )
    local_records = load_checkpoint_records(local_state)
    assert local_records[0]["state"] == "completed"
    assert local_records[0]["outcome"]["failure_stage"] == "timeout"
    assert local_records[1]["outcome"]["failure_stage"] == "memory"
    assert not local_records[0]["outcome"]["eligible"]
    assert not local_records[1]["outcome"]["eligible"]

    global_state = tmp_path / "state-global"
    global_executor = FakeExecutor(wall_seconds=1000.0)
    global_summary = run_p1_6_campaign(
        root,
        executor=global_executor,
        state_directory=global_state,
        utc_now=lambda: UTC,
    )
    global_records = load_checkpoint_records(global_state)
    global_gate = evaluate_g1(
        json.loads(
            (root / "campaigns" / "p1" / "campaign_manifest.yaml").read_text()
        ),
        global_records,
    )

    assert global_summary.closed
    assert global_summary.stop_reason == "global_wall_limit_exhausted"
    assert len(global_executor.case_ids) == 65
    assert global_executor.wall_limits[:64] == [1800.0] * 64
    assert global_executor.wall_limits[64] == 800.0
    assert global_summary.accumulated_wall_seconds == 64800.0
    assert global_summary.never_started_count == 37
    assert global_gate.decision == "INCONCLUSIVE_P1"
    assert "not_all_102_cases_attempted" in global_gate.reasons


def test_g1_identity_symmetry_and_coverage_classifications(tmp_path: Path) -> None:
    root = _root(tmp_path, "g1")
    state = tmp_path / "state-g1"
    run_p1_6_campaign(
        root,
        executor=FakeExecutor(),
        state_directory=state,
        utc_now=lambda: UTC,
    )
    records = list(load_checkpoint_records(state))
    manifest = json.loads(
        (root / "campaigns" / "p1" / "campaign_manifest.yaml").read_text()
    )
    broken_identity = json.loads(json.dumps(records))
    broken_identity[0]["outcome"]["model_be_forces_xyz"][0][0] += 4.0e-12
    assert evaluate_g1(manifest, broken_identity).decision == "NO_GO_P2"

    ineligible_audit = json.loads(json.dumps(records))
    ineligible_audit[96]["outcome"]["eligible"] = False
    ineligible_audit[96]["failure_stage"] = "convergence"
    ineligible_audit[96]["failure_reason"] = "unconfirmed_at_21"
    assert evaluate_g1(manifest, ineligible_audit).decision == "INCONCLUSIVE_P1"

    broken_contract = evaluate_g1(manifest, records, contract_valid=False)
    assert broken_contract.decision == "NO_GO_P2"

    reordered = list(records)
    reordered[0], reordered[1] = reordered[1], reordered[0]
    assert evaluate_g1(manifest, reordered).decision == "NO_GO_P2"


def test_artifact_module_has_no_scientific_or_solver_imports() -> None:
    path = ROOT / "src" / "acoustic_ms" / "p1_campaign_artifacts.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    imports = {
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    }
    imports.update(
        node.module or ""
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom)
    )
    assert not any(
        name.startswith(("acoustic_ms", "numpy", "scipy")) for name in imports
    )


class TickClock:
    def __init__(self) -> None:
        self.value = 0.0

    def __call__(self) -> float:
        self.value += 0.1
        return self.value


class FakeModelESolver:
    def __init__(self) -> None:
        self.orders: list[int] = []

    def __call__(self, positions, k, radius, energy, f0, f1, order):
        del k, radius, energy, f0, f1
        self.orders.append(order)
        particle_count = len(np.asarray(positions))
        modes = tuple(range((order + 1) ** 2))
        active_modes = ((1, -1), (1, 0))
        base = np.array([[1.0, 0.25, 0.0], [-1.0, -0.25, 0.0]])
        solution = SimpleNamespace(
            balanced_condition_number=1.5,
            balanced_backward_error=1.0e-15,
            effective_incident_closure_error=2.0e-15,
            scattering_closure_error=3.0e-15,
            scattered_coefficients=np.zeros((particle_count, len(modes)), dtype=complex),
            active_modes=active_modes,
            modes=modes,
            balanced_system_matrix=np.eye(particle_count * len(active_modes), dtype=complex),
            production_solver="balanced_sqrt",
        )
        return SimpleNamespace(
            solution=solution,
            total_forces_xyz=4.0 * base,
            external_forces_xyz=base,
            interaction_forces_xyz=3.0 * base,
            external_scattered_forces_xyz=2.0 * base,
            scattered_scattered_forces_xyz=base,
            scattered_incident_coefficients=np.zeros_like(solution.scattered_coefficients),
            decomposition_residual=4.0e-15,
            lmax=order,
        )


def test_real_adapter_with_fake_solver_reuses_each_dimer_order_once() -> None:
    manifest = json.loads(
        (ROOT / "campaigns" / "p1" / "campaign_manifest.yaml").read_text()
    )
    case = manifest["cases"][0]
    solver = FakeModelESolver()
    outcome = execute_model_e_case(
        case,
        manifest,
        solver=solver,
        clock=TickClock(),
        rss_reader=lambda: 128 * 1024**2,
    )

    assert solver.orders == [2, 3, 4, 5]
    assert outcome["model_e_solve_count"] == 4
    assert outcome["evaluated_lmax"] == [2, 3, 4, 5]
    np.testing.assert_allclose(
        outcome["model_be_forces_xyz"], outcome["model_e_forces_xyz"]
    )
    assert tuple(outcome["orders"][0]["channels"]) == (
        "total", "interaction", "external_scattered", "scattered_scattered"
    )
