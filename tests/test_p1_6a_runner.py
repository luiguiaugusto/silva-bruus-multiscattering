"""Fake-only tests for the P1.6 runner, checkpoints, artifacts and G1."""

from __future__ import annotations

import ast
import csv
from contextlib import contextmanager
from datetime import datetime, timezone
import io
import json
from pathlib import Path
import shutil
import subprocess
from types import SimpleNamespace

import numpy as np
import pytest

import acoustic_ms.p1_campaign as p1_campaign_module
import acoustic_ms.p1_campaign_artifacts as p1_campaign_artifacts_module

from acoustic_ms.p1_campaign import (
    CampaignCaseTimeout,
    CampaignExecutionError,
    CampaignInfrastructureError,
    CampaignSerializationError,
    NUMERIC_ENVIRONMENT_KEYS,
    capture_p1_6_execution_provenance,
    execute_model_e_case,
    execute_model_e_case_with_limits,
    run_p1_6_campaign,
)
from acoustic_ms.p1_campaign_artifacts import (
    ARTIFACT_PATHS,
    G1_BUDGET,
    build_campaign_artifacts,
    evaluate_g1,
    load_campaign_checkpoint,
    load_checkpoint_records,
    normalized_rms_error_xyz_pure,
    publish_campaign_artifacts,
)
from acoustic_ms.model_e_comparison import normalized_rms_error_xyz


ROOT = Path(__file__).resolve().parents[1]
UTC = datetime(2026, 8, 20, 12, 0, tzinfo=timezone.utc)
EXECUTION_COMMIT = "b" * 40
MANIFEST_SHA256 = "a041e07ae93e9a858bad809427039bf593641ad1f9e341ed89b9d91f648f297d"


def _provenance(root: Path, **changes: object) -> dict[str, object]:
    value: dict[str, object] = {
        "schema_version": "1.0.0",
        "git_commit": EXECUTION_COMMIT,
        "manifest_sha256": MANIFEST_SHA256,
        "branch": "agent/p1-6b-execute",
        "directory": str(root.resolve()),
        "sys_executable": "/frozen/python",
        "sys_argv": ["scripts/run_p1_6_campaign.py", "--execute"],
        "python_version": "3.test",
        "platform": "test-platform",
        "numpy_version": "test-numpy",
        "scipy_version": "test-scipy",
        "numeric_environment": {
            key: "0" if key == "PYTHONHASHSEED" else "1"
            for key in NUMERIC_ENVIRONMENT_KEYS
        },
    }
    value.update(changes)
    return value


def _root(tmp_path: Path, name: str) -> Path:
    root = tmp_path / name
    destination = root / "campaigns" / "p1"
    destination.mkdir(parents=True)
    for manifest in ("campaign_manifest_r2.yaml", "pilot_manifest.yaml"):
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


def test_real_git_head_and_allowlisted_environment_are_captured() -> None:
    environment = {
        key: "0" if key == "PYTHONHASHSEED" else "1"
        for key in NUMERIC_ENVIRONMENT_KEYS
    }
    environment["SECRET_MUST_NOT_BE_SERIALIZED"] = "secret"
    provenance = capture_p1_6_execution_provenance(
        ROOT,
        manifest_sha256=MANIFEST_SHA256,
        environ=environment,
        argv=["runner", "--execute"],
    )
    expected_head = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()

    assert provenance["git_commit"] == expected_head
    assert provenance["manifest_sha256"] == MANIFEST_SHA256
    assert provenance["sys_argv"] == ["runner", "--execute"]
    assert provenance["numeric_environment"] == {
        key: environment[key] for key in NUMERIC_ENVIRONMENT_KEYS
    }
    assert "SECRET_MUST_NOT_BE_SERIALIZED" not in json.dumps(provenance)


@pytest.mark.parametrize("pythonhashseed", [None, "1"])
def test_pythonhashseed_must_be_present_and_zero(pythonhashseed: str | None) -> None:
    environment = {
        key: "1" for key in NUMERIC_ENVIRONMENT_KEYS if key != "PYTHONHASHSEED"
    }
    if pythonhashseed is not None:
        environment["PYTHONHASHSEED"] = pythonhashseed
    with pytest.raises(CampaignExecutionError, match="frozen numeric environment"):
        capture_p1_6_execution_provenance(
            ROOT,
            manifest_sha256=MANIFEST_SHA256,
            environ=environment,
        )


@pytest.mark.parametrize(
    ("mutation", "message"),
    (
        ({"git_commit": "c" * 40}, "execution HEAD differs"),
        ({"manifest_sha256": "0" * 64}, "manifest hash mismatch"),
        (
            {
                "numeric_environment": {
                    **{
                        key: "0" if key == "PYTHONHASHSEED" else "1"
                        for key in NUMERIC_ENVIRONMENT_KEYS
                    },
                    "OMP_NUM_THREADS": "2",
                }
            },
            "frozen numeric environment",
        ),
    ),
)
def test_resume_rejects_head_hash_or_numeric_environment_change(
    tmp_path: Path,
    mutation: dict[str, object],
    message: str,
) -> None:
    root = _root(tmp_path, f"provenance-{message.split()[0]}")
    state = tmp_path / f"state-{message.split()[0]}"
    run_p1_6_campaign(
        root,
        executor=FakeExecutor(),
        state_directory=state,
        utc_now=lambda: UTC,
        max_new_cases=1,
        execution_provenance=_provenance(root),
    )
    ledger, _ = load_campaign_checkpoint(state)
    assert ledger["execution_provenance"] == _provenance(root)
    resumed_executor = FakeExecutor()
    with pytest.raises(CampaignExecutionError, match=message):
        run_p1_6_campaign(
            root,
            executor=resumed_executor,
            state_directory=state,
            utc_now=lambda: UTC,
            max_new_cases=1,
            execution_provenance=_provenance(root, **mutation),
        )
    assert resumed_executor.case_ids == []


def test_resume_rejects_real_head_change_before_executor(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    root = _root(tmp_path, "real-head-change")
    state = tmp_path / "state-real-head-change"
    subprocess.run(["git", "init", "-b", "agent/p1-6b-execute"], cwd=root, check=True)
    subprocess.run(["git", "config", "user.name", "P1.6A Test"], cwd=root, check=True)
    subprocess.run(
        ["git", "config", "user.email", "p1-6a@example.invalid"],
        cwd=root,
        check=True,
    )
    subprocess.run(["git", "add", "campaigns"], cwd=root, check=True)
    subprocess.run(["git", "commit", "-m", "frozen manifests"], cwd=root, check=True)
    for key in NUMERIC_ENVIRONMENT_KEYS:
        monkeypatch.setenv(key, "0" if key == "PYTHONHASHSEED" else "1")

    run_p1_6_campaign(
        root,
        executor=FakeExecutor(),
        state_directory=state,
        utc_now=lambda: UTC,
        max_new_cases=1,
    )
    (root / "execution-marker.txt").write_text("new head\n", encoding="utf-8")
    subprocess.run(["git", "add", "execution-marker.txt"], cwd=root, check=True)
    subprocess.run(["git", "commit", "-m", "change head"], cwd=root, check=True)
    resumed_executor = FakeExecutor()

    with pytest.raises(CampaignExecutionError, match="execution HEAD differs"):
        run_p1_6_campaign(
            root,
            executor=resumed_executor,
            state_directory=state,
            utc_now=lambda: UTC,
            max_new_cases=1,
        )
    assert resumed_executor.case_ids == []


def test_resume_rejects_preexisting_output_before_executor(tmp_path: Path) -> None:
    root = _root(tmp_path, "resume-output")
    state = tmp_path / "state-resume-output"
    run_p1_6_campaign(
        root,
        executor=FakeExecutor(),
        state_directory=state,
        utc_now=lambda: UTC,
        max_new_cases=1,
        execution_provenance=_provenance(root),
    )
    output = root / "campaigns" / "p1" / "p1_6b_r2" / "data_raw.csv"
    output.parent.mkdir(parents=True)
    output.write_text("preexisting response\n", encoding="utf-8")
    resumed_executor = FakeExecutor()

    with pytest.raises(FileExistsError, match="campaign output exists"):
        run_p1_6_campaign(
            root,
            executor=resumed_executor,
            state_directory=state,
            utc_now=lambda: UTC,
            max_new_cases=1,
            execution_provenance=_provenance(root),
        )
    assert resumed_executor.case_ids == []


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
        execution_provenance=_provenance(root),
    )
    ledger, records = load_campaign_checkpoint(state)
    manifest = json.loads(
        (root / "campaigns" / "p1" / "campaign_manifest_r2.yaml").read_text()
    )
    calls_before = list(executor.case_ids)
    first, first_gate = build_campaign_artifacts(
        manifest, records, _provenance(root)
    )
    second, second_gate = build_campaign_artifacts(
        manifest, records, _provenance(root)
    )

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
    assert ledger["execution_provenance"] == _provenance(root)
    assert ledger["manifest_provenance"] == manifest["provenance"]

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
    assert {row["git_commit"] for row in raw} == {EXECUTION_COMMIT}
    assert {row["manifest_git_commit"] for row in raw} == {
        manifest["provenance"]["git_commit"]
    }
    assert sum(row["metric"] == "rotational_covariance_error" for row in derived) == 6
    epsilon_a = [row for row in derived if row["metric"] == "epsilon_a_e"]
    epsilon_be = [row for row in derived if row["metric"] == "epsilon_be_e"]
    absolute = [row for row in derived if row["metric"] == "be_minus_a_rms"]
    assert len(epsilon_a) == len(epsilon_be) == len(absolute) == 102
    assert {row["value"] for row in epsilon_a} == {"0.5"}
    assert {row["value"] for row in epsilon_be} == {"0"}
    assert {row["value"] for row in absolute} == {"1"}
    assert len(plot) == 96
    assert {row["y_name"] for row in plot} == {"epsilon_a_e"}
    assert {row["y_value"] for row in plot} == {"0.5"}
    assert all(row["applicable"] == "true" for row in plot)
    assert failures == []
    assert len(performance) == 102
    assert all(row["model_e_solve_count"] == "4" for row in performance)
    assert {row["git_commit"] for row in performance} == {EXECUTION_COMMIT}
    assert {row["manifest_git_commit"] for row in performance} == {
        manifest["provenance"]["git_commit"]
    }
    assert all(row["effective_wall_seconds"] == "1800" for row in performance)

    published = publish_campaign_artifacts(root, first)
    assert set(published.values())
    with pytest.raises(FileExistsError, match="second publication"):
        publish_campaign_artifacts(root, first)
    with pytest.raises(FileExistsError, match="campaign output exists"):
        run_p1_6_campaign(
            root,
            executor=executor,
            state_directory=state,
            utc_now=lambda: UTC,
            execution_provenance=_provenance(root),
        )


def test_artifact_set_is_not_visible_when_staging_write_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    root = tmp_path / "publication-failure"
    artifacts = {
        name: f"{name}\n".encode("utf-8")
        for name in ARTIFACT_PATHS
    }
    original = p1_campaign_artifacts_module._atomic_write
    calls = 0

    def fail_on_second_write(path: Path, payload: bytes) -> None:
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("injected staging failure")
        original(path, payload)

    monkeypatch.setattr(
        p1_campaign_artifacts_module,
        "_atomic_write",
        fail_on_second_write,
    )

    with pytest.raises(OSError, match="injected staging failure"):
        publish_campaign_artifacts(root, artifacts)

    output_directory = root / "campaigns" / "p1" / "p1_6b_r2"
    assert not output_directory.exists()
    assert not list(output_directory.parent.glob(".p1_6b_r2.*"))


def test_interrupted_case_is_never_retried_and_resume_uses_next_case(
    tmp_path: Path,
) -> None:
    root = _root(tmp_path, "resume")
    state = tmp_path / "state-resume"
    calls: list[int] = []
    reservations_seen_inside_executor: list[float] = []

    def interrupted(case, manifest):
        del manifest
        order = int(case["case_order"])
        calls.append(order)
        live_ledger = json.loads(
            (state / "campaign_ledger.json").read_text(encoding="utf-8")
        )
        reservations_seen_inside_executor.append(
            live_ledger["cases"][order - 1]["effective_wall_seconds"]
        )
        if order == 2:
            raise KeyboardInterrupt("simulated process loss")
        return _outcome(case)

    with pytest.raises(KeyboardInterrupt, match="simulated process loss"):
        run_p1_6_campaign(
            root,
            executor=interrupted,
            state_directory=state,
            utc_now=lambda: UTC,
            execution_provenance=_provenance(root),
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
        execution_provenance=_provenance(root),
    )
    records = load_checkpoint_records(state)

    assert calls == [1, 2]
    assert reservations_seen_inside_executor == [1800.0, 1800.0]
    assert resumed_calls == [3]
    assert records[1]["state"] == "interrupted"
    assert records[1]["attempt_count"] == 1
    assert records[1]["failure_reason"] == "previous_process_interrupted_after_start"
    assert records[1]["effective_wall_seconds"] == 1800.0
    assert records[1]["wall_seconds_debited"] == 1800.0
    assert summary.completed_count == 2
    assert summary.interrupted_count == 1
    assert summary.never_started_count == 99
    assert summary.accumulated_wall_seconds == 1802.0


def test_accumulated_abandoned_reservations_exhaust_global_budget(
    tmp_path: Path,
) -> None:
    root = _root(tmp_path, "abandoned-reservations")
    state = tmp_path / "state-abandoned-reservations"
    attempted: list[int] = []

    def interrupted(case, manifest):
        del manifest
        attempted.append(int(case["case_order"]))
        raise KeyboardInterrupt("late interruption")

    for _ in range(36):
        with pytest.raises(KeyboardInterrupt, match="late interruption"):
            run_p1_6_campaign(
                root,
                executor=interrupted,
                state_directory=state,
                utc_now=lambda: UTC,
                execution_provenance=_provenance(root),
            )

    def forbidden_executor(case, manifest):
        del case, manifest
        raise AssertionError("global exhaustion must not start another case")

    summary = run_p1_6_campaign(
        root,
        executor=forbidden_executor,
        state_directory=state,
        utc_now=lambda: UTC,
        execution_provenance=_provenance(root),
    )
    records = load_checkpoint_records(state)

    assert attempted == list(range(1, 37))
    assert summary.closed
    assert summary.campaign_decision == "INCONCLUSIVE_P1"
    assert summary.accumulated_wall_seconds == 64800.0
    assert summary.interrupted_count == 36
    assert summary.never_started_count == 66
    assert all(record["wall_seconds_debited"] == 1800.0 for record in records[:36])
    assert all(record["failure_stage"] == "global_limit" for record in records[36:])


def test_fractional_global_balance_is_reserved_and_closes_without_overdebit(
    tmp_path: Path,
) -> None:
    root = _root(tmp_path, "fractional-balance")
    state = tmp_path / "state-fractional-balance"
    run_p1_6_campaign(
        root,
        executor=FakeExecutor(),
        state_directory=state,
        utc_now=lambda: UTC,
        max_new_cases=1,
        execution_provenance=_provenance(root),
    )
    ledger_path = state / "campaign_ledger.json"
    ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    ledger["accumulated_wall_seconds"] = 64799.75
    ledger_path.write_text(
        json.dumps(ledger, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    executor = FakeExecutor(wall_seconds=1.0)
    summary = run_p1_6_campaign(
        root,
        executor=executor,
        state_directory=state,
        utc_now=lambda: UTC,
        execution_provenance=_provenance(root),
    )
    records = load_checkpoint_records(state)

    assert executor.wall_limits == [0.25]
    assert records[1]["effective_wall_seconds"] == 0.25
    assert records[1]["wall_seconds_debited"] == 0.25
    assert summary.accumulated_wall_seconds == 64800.0
    assert summary.campaign_decision == "INCONCLUSIVE_P1"


def test_fractional_effective_limit_reaches_resource_timer(monkeypatch) -> None:
    captured: list[float] = []

    @contextmanager
    def fake_limits(wall_seconds, memory_bytes, state):
        del memory_bytes, state
        captured.append(wall_seconds)
        yield

    monkeypatch.setattr(p1_campaign_module, "_resource_limits", fake_limits)
    manifest = {"resources": {
        "wall_seconds_per_case": 0.25,
        "peak_rss_bytes_per_case": 4 * 1024**3,
    }}
    result = execute_model_e_case_with_limits(
        {},
        manifest,
        executor=lambda case, runtime_manifest: {
            "case": case,
            "wall": runtime_manifest["resources"]["wall_seconds_per_case"],
        },
    )

    assert captured == [0.25]
    assert result["wall"] == 0.25


def test_measured_wall_after_normalization_reapplies_local_limit(
    tmp_path: Path,
) -> None:
    root = _root(tmp_path, "measured-wall")
    state = tmp_path / "state-measured-wall"
    ticks = iter((0.0, 1800.25))
    summary = run_p1_6_campaign(
        root,
        executor=FakeExecutor(wall_seconds=1.0),
        state_directory=state,
        utc_now=lambda: UTC,
        monotonic=lambda: next(ticks),
        max_new_cases=1,
        execution_provenance=_provenance(root),
    )
    record = load_checkpoint_records(state)[0]

    assert record["outcome"]["wall_seconds"] == 1800.25
    assert record["outcome"]["failure_stage"] == "timeout"
    assert not record["outcome"]["eligible"]
    assert record["wall_seconds_debited"] == 1800.25
    assert summary.accumulated_wall_seconds == 1800.25


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
        execution_provenance=_provenance(root),
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
        execution_provenance=_provenance(root),
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
        execution_provenance=_provenance(root),
    )
    global_records = load_checkpoint_records(global_state)
    global_gate = evaluate_g1(
        json.loads(
            (root / "campaigns" / "p1" / "campaign_manifest_r2.yaml").read_text()
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
    assert global_summary.campaign_decision == "INCONCLUSIVE_P1"
    assert all(
        record["failure_stage"] == "global_limit"
        for record in global_records[65:]
    )
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
        execution_provenance=_provenance(root),
    )
    records = list(load_checkpoint_records(state))
    manifest = json.loads(
        (root / "campaigns" / "p1" / "campaign_manifest_r2.yaml").read_text()
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

    scientific_magnitude_only = json.loads(json.dumps(records))
    scientific_magnitude_only[0]["outcome"]["model_a_forces_xyz"] = [
        [1.0e12, -2.0e12, 3.0e12],
        [-1.0e12, 2.0e12, -3.0e12],
    ]
    assert evaluate_g1(manifest, scientific_magnitude_only).decision == "GO_P2"


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


@pytest.mark.parametrize(
    ("reference", "model"),
    (
        (
            [[1.0, -2.0, 0.5], [-1.0, 2.0, -0.5]],
            [[0.9, -1.8, 0.4], [-0.9, 1.8, -0.4]],
        ),
        (
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            [[1.0e-20, 0.0, 0.0], [-1.0e-20, 0.0, 0.0]],
        ),
        (
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        ),
    ),
)
def test_pure_epsilon_matches_established_helper_exactly(reference, model) -> None:
    expected_value, expected_applicable = normalized_rms_error_xyz(reference, model)
    observed_value, observed_applicable = normalized_rms_error_xyz_pure(
        reference,
        model,
    )

    assert observed_value == pytest.approx(expected_value, rel=0.0, abs=1.0e-30)
    assert observed_applicable is expected_applicable


def test_epsilon_is_scale_invariant_without_floor_or_clipping() -> None:
    reference = [[1.0, 2.0, 0.0], [-1.0, -2.0, 0.0]]
    model = [[0.75, 1.5, 0.0], [-0.75, -1.5, 0.0]]
    baseline, applicable = normalized_rms_error_xyz_pure(reference, model)

    for scale in (1.0e-18, 1.0e18):
        scaled_reference = [
            [scale * component for component in vector] for vector in reference
        ]
        scaled_model = [
            [scale * component for component in vector] for vector in model
        ]
        observed, scaled_applicable = normalized_rms_error_xyz_pure(
            scaled_reference,
            scaled_model,
        )
        assert observed == pytest.approx(baseline, rel=2.0e-15)
        assert scaled_applicable is applicable


def test_zero_reference_is_explicitly_inapplicable_in_derived_and_plot(
    tmp_path: Path,
) -> None:
    root = _root(tmp_path, "zero-reference")
    state = tmp_path / "state-zero-reference"

    def executor(case, manifest):
        del manifest
        outcome = _outcome(case)
        if case["case_order"] == 1:
            zeros = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
            outcome["model_be_forces_xyz"] = zeros
            outcome["model_e_forces_xyz"] = zeros
        return outcome

    run_p1_6_campaign(
        root,
        executor=executor,
        state_directory=state,
        utc_now=lambda: UTC,
        execution_provenance=_provenance(root),
    )
    manifest = json.loads(
        (root / "campaigns" / "p1" / "campaign_manifest_r2.yaml").read_text()
    )
    records = load_checkpoint_records(state)
    first, _ = build_campaign_artifacts(manifest, records, _provenance(root))
    second, _ = build_campaign_artifacts(manifest, records, _provenance(root))
    derived = list(
        csv.DictReader(io.StringIO(first["data_derived.csv"].decode("utf-8")))
    )
    plot = list(csv.DictReader(io.StringIO(first["data_plot.csv"].decode("utf-8"))))
    case_id = manifest["cases"][0]["case_id"]
    epsilon_row = next(
        row
        for row in derived
        if row["case_id"] == case_id and row["metric"] == "epsilon_a_e"
    )
    absolute_row = next(
        row
        for row in derived
        if row["case_id"] == case_id and row["metric"] == "be_minus_a_rms"
    )
    plot_row = next(row for row in plot if row["case_id"] == case_id)

    assert first == second
    assert epsilon_row["applicable"] == "false"
    assert epsilon_row["reason"] == (
        "reference_rms_numerically_zero;value_is_absolute_rms"
    )
    assert epsilon_row["value"] == absolute_row["value"] == "1"
    assert plot_row["y_name"] == "epsilon_a_e"
    assert plot_row["y_value"] == ""
    assert plot_row["applicable"] == "false"
    assert plot_row["reason"] == "reference_rms_numerically_zero"


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
        (ROOT / "campaigns" / "p1" / "campaign_manifest_r2.yaml").read_text()
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
    normalized = p1_campaign_module._normalize_outcome(outcome, manifest)
    encoded = p1_campaign_module._json_bytes(normalized)
    decoded = json.loads(encoded)
    planar = decoded["orders"][0]["diagnostics"]["planar_symmetry_pass"]

    assert isinstance(planar, bool)
    assert planar is True



@pytest.mark.parametrize(
    ("value", "expected", "expected_type"),
    (
        (np.bool_(True), True, bool),
        (np.int64(7), 7, int),
        (np.float64(1.25), 1.25, float),
    ),
)
def test_json_boundary_converts_supported_numpy_scalars(
    value: np.generic,
    expected: object,
    expected_type: type[object],
) -> None:
    decoded = json.loads(p1_campaign_module._json_bytes({"value": value}))

    assert decoded["value"] == expected
    assert type(decoded["value"]) is expected_type


@pytest.mark.parametrize("value", (np.complex128(1.0 + 2.0j), object()))
def test_json_boundary_rejects_complex_and_unknown_types(value: object) -> None:
    with pytest.raises(CampaignSerializationError, match=r"\$\.value.*unsupported"):
        p1_campaign_module._json_bytes({"value": value})


def test_real_adapter_to_checkpoint_and_artifact_regeneration_without_real_solver(
    tmp_path: Path,
) -> None:
    root = _root(tmp_path, "real-adapter-end-to-end")
    state = tmp_path / "state-real-adapter-end-to-end"
    solver = FakeModelESolver()

    def executor(case, manifest):
        return execute_model_e_case(
            case,
            manifest,
            solver=solver,
            clock=TickClock(),
            rss_reader=lambda: 128 * 1024**2,
        )

    summary = run_p1_6_campaign(
        root,
        executor=executor,
        state_directory=state,
        utc_now=lambda: UTC,
        max_new_cases=1,
        execution_provenance=_provenance(root),
    )
    ledger, records = load_campaign_checkpoint(state)
    manifest = json.loads(
        (root / "campaigns" / "p1" / "campaign_manifest_r2.yaml").read_text()
    )
    checkpoint = json.loads((state / "cases" / "001.json").read_text())
    first, first_gate = build_campaign_artifacts(
        manifest,
        records,
        ledger["execution_provenance"],
    )
    second, second_gate = build_campaign_artifacts(
        manifest,
        records,
        ledger["execution_provenance"],
    )

    assert summary.completed_count == 1
    assert solver.orders == [2, 3, 4, 5]
    assert checkpoint["state"] == "completed"
    assert checkpoint["outcome"]["model_e_solve_count"] == 4
    assert isinstance(
        checkpoint["outcome"]["orders"][0]["diagnostics"][
            "planar_symmetry_pass"
        ],
        bool,
    )
    assert first == second
    assert first_gate == second_gate
    assert first["data_raw.csv"]
    assert first["performance.csv"]


def test_unexpected_serialization_failure_stops_immediately_and_publishes_nothing(
    tmp_path: Path,
) -> None:
    root = _root(tmp_path, "fatal-serialization")
    state = tmp_path / "state-fatal-serialization"
    calls: list[int] = []

    def executor(case, manifest):
        del manifest
        calls.append(int(case["case_order"]))
        outcome = _outcome(case)
        outcome["orders"][0]["diagnostics"]["bad"] = object()
        return outcome

    with pytest.raises(CampaignInfrastructureError, match="serialization"):
        run_p1_6_campaign(
            root,
            executor=executor,
            state_directory=state,
            utc_now=lambda: UTC,
            execution_provenance=_provenance(root),
        )

    ledger, records = load_campaign_checkpoint(state)
    assert calls == [1]
    assert ledger["closed"] is True
    assert ledger["stop_reason"] == "invalid_infrastructure"
    assert ledger["campaign_decision"] == "INVALID_P1.6B_R2_INFRASTRUCTURE"
    assert records[0]["state"] == "interrupted"
    assert records[0]["failure_stage"] == "serialization"
    assert records[0]["outcome"] is None
    assert all(record["state"] == "never_started" for record in records[1:])
    assert not any((root / relative).exists() for relative in ARTIFACT_PATHS.values())


@pytest.mark.parametrize(
    ("exception", "expected_stage"),
    (
        (CampaignExecutionError("synthetic contract failure"), "contract"),
        (RuntimeError("synthetic infrastructure failure"), "infrastructure"),
    ),
)
def test_unexpected_contract_and_infrastructure_failures_are_fatal(
    tmp_path: Path,
    exception: Exception,
    expected_stage: str,
) -> None:
    root = _root(tmp_path, f"fatal-{expected_stage}")
    state = tmp_path / f"state-fatal-{expected_stage}"
    calls: list[int] = []

    def executor(case, manifest):
        del manifest
        calls.append(int(case["case_order"]))
        raise exception

    with pytest.raises(CampaignInfrastructureError, match=expected_stage):
        run_p1_6_campaign(
            root,
            executor=executor,
            state_directory=state,
            utc_now=lambda: UTC,
            execution_provenance=_provenance(root),
        )

    ledger, records = load_campaign_checkpoint(state)
    assert calls == [1]
    assert ledger["campaign_decision"] == "INVALID_P1.6B_R2_INFRASTRUCTURE"
    assert records[0]["failure_stage"] == expected_stage
    assert all(record["state"] == "never_started" for record in records[1:])
