"""No-solver tests for the deterministic single-use P1.5 runner."""

from __future__ import annotations

import csv
from datetime import datetime, timezone
from hashlib import sha256
import io
from pathlib import Path
import shutil
from types import SimpleNamespace

import numpy as np
import pytest

import acoustic_ms.p1_pilot as pilot_module
from acoustic_ms.p1_pilot import (
    PILOT_CASE_ID,
    PILOT_OUTPUT_RELATIVE,
    PilotExecutionError,
    execute_p1_5_pilot,
    load_p1_5_configuration,
    verify_p1_5_derivations,
)
from acoustic_ms.paper_pipeline import P1_FROZEN_MANIFEST_SHA256


ROOT = Path(__file__).resolve().parents[1]
SOURCE_COMMIT = "a" * 40
COMMAND = "python scripts/run_p1_5_timed_pilot.py --execute"
ENVIRONMENT = {
    "python_version": "test",
    "numpy_version": "test",
    "scipy_version": "test",
    "openblas_num_threads": "1",
}
UTC = datetime(2026, 8, 19, 20, 0, tzinfo=timezone.utc)


class TickClock:
    def __init__(self) -> None:
        self.value = 0.0

    def __call__(self) -> float:
        self.value += 0.125
        return self.value


class FakeStableSolver:
    def __init__(self) -> None:
        self.orders: list[int] = []

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
        self.orders.append(order)
        particle_count = len(np.asarray(positions))
        modes = tuple(range((order + 1) ** 2))
        active_modes = ((1, -1), (1, 0))
        dimension = particle_count * len(active_modes)
        base = np.array([[1.0, 0.25, 0.0], [-1.0, -0.25, 0.0]])
        solution = SimpleNamespace(
            balanced_condition_number=1.5,
            balanced_backward_error=1.0e-15,
            effective_incident_closure_error=2.0e-15,
            scattering_closure_error=3.0e-15,
            scattered_coefficients=np.zeros(
                (particle_count, len(modes)), dtype=complex
            ),
            active_modes=active_modes,
            modes=modes,
            balanced_system_matrix=np.eye(dimension, dtype=complex),
            production_solver="balanced_sqrt",
        )
        return SimpleNamespace(
            solution=solution,
            total_forces_xyz=4.0 * base,
            external_forces_xyz=base,
            interaction_forces_xyz=3.0 * base,
            external_scattered_forces_xyz=2.0 * base,
            scattered_scattered_forces_xyz=base,
            scattered_incident_coefficients=np.zeros_like(
                solution.scattered_coefficients
            ),
            decomposition_residual=4.0e-15,
            lmax=order,
        )


class FailingSolver:
    def __init__(self, error: Exception) -> None:
        self.error = error
        self.orders: list[int] = []

    def __call__(self, *args: object, **kwargs: object) -> None:
        del kwargs
        self.orders.append(int(args[6]))
        raise self.error


class NeverConvergingSolver(FakeStableSolver):
    def __call__(self, *args: object, **kwargs: object) -> SimpleNamespace:
        result = super().__call__(*args, **kwargs)
        order = int(args[6])
        for attribute in (
            "total_forces_xyz",
            "interaction_forces_xyz",
            "external_scattered_forces_xyz",
            "scattered_scattered_forces_xyz",
        ):
            setattr(result, attribute, order * getattr(result, attribute))
        return result


def _pilot_root(tmp_path: Path, name: str) -> Path:
    root = tmp_path / name
    destination = root / "campaigns" / "p1"
    destination.mkdir(parents=True)
    for manifest in ("campaign_manifest.yaml", "pilot_manifest.yaml"):
        shutil.copy2(ROOT / "campaigns" / "p1" / manifest, destination / manifest)
    return root


def _read(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(io.StringIO(path.read_text(encoding="utf-8"))))


def _execute(root: Path, solver: object):
    return execute_p1_5_pilot(
        root,
        source_commit=SOURCE_COMMIT,
        command=COMMAND,
        solver=solver,
        clock=TickClock(),
        utc_now=lambda: UTC,
        rss_reader=lambda: 128 * 1024**2,
        environment=ENVIRONMENT,
        enforce_resource_limits=False,
    )


def test_only_frozen_pilot_is_executable_and_manifests_remain_immutable(
    tmp_path: Path,
) -> None:
    root = _pilot_root(tmp_path, "configuration")
    before = {
        path.name: path.read_bytes()
        for path in (root / "campaigns" / "p1").glob("*_manifest.yaml")
    }
    configuration = load_p1_5_configuration(root)

    assert configuration.case["case_id"] == PILOT_CASE_ID
    assert configuration.case["enabled"] is True
    assert not any(
        case["enabled"]
        for case in configuration.confirmatory_manifest["cases"]
    )
    assert len(configuration.confirmatory_manifest["cases"]) == 102
    assert configuration.pilot_manifest["provenance"]["manifest_sha256"] == (
        P1_FROZEN_MANIFEST_SHA256["p1_dimer_resource_pilot"]
    )
    assert {
        path.name: path.read_bytes()
        for path in (root / "campaigns" / "p1").glob("*_manifest.yaml")
    } == before


def test_fake_runner_is_deterministic_atomic_complete_and_single_use(
    tmp_path: Path,
) -> None:
    roots = [_pilot_root(tmp_path, name) for name in ("first", "second")]
    summaries = []
    published = []
    solvers = []
    for root in roots:
        solver = FakeStableSolver()
        solvers.append(solver)
        summaries.append(_execute(root, solver))
        output = root / PILOT_OUTPUT_RELATIVE
        published.append(
            {path.name: path.read_bytes() for path in sorted(output.iterdir())}
        )

    assert summaries[0].decision == "GO_P1.6A_BLIND_FREEZE"
    assert summaries[0].stop_reason == "all_channels_confirmed"
    assert summaries[0].attempted_lmax == (2, 3, 4, 5)
    assert summaries[0].evaluated_lmax == (2, 3, 4, 5)
    assert summaries[0].final_lmax == 5
    assert summaries[0].eligible and summaries[0].converged
    assert solvers[0].orders == [2, 3, 4, 5]
    assert published[0] == published[1]
    assert set(published[0]) == {
        "data_raw.csv",
        "data_derived.csv",
        "data_plot.csv",
        "failures.csv",
        "performance.csv",
    }

    output = roots[0] / PILOT_OUTPUT_RELATIVE
    raw = _read(output / "data_raw.csv")
    performance = _read(output / "performance.csv")
    derived = _read(output / "data_derived.csv")
    failures = _read(output / "failures.csv")
    assert len(raw) == 4 * 2 * 4
    assert {row["force_channel"] for row in raw} == {
        "total",
        "interaction",
        "external_scattered",
        "scattered_scattered",
    }
    assert {row["lmax"] for row in raw} == {"2", "3", "4", "5"}
    assert all(row["include_in_scientific_tables"] == "false" for row in raw)
    assert all(row["classification"] == "development" for row in raw)
    assert len(performance) == 4
    assert all(row["command"] == COMMAND for row in performance)
    assert all(row["worker_count"] == row["blas_threads"] == "1" for row in performance)
    assert all(row["started_utc"].endswith("Z") for row in performance)
    assert all(row["completed_utc"].endswith("Z") for row in performance)
    assert all("force" not in row["metric"] for row in derived)
    assert failures[0]["include_in_scientific_tables"] == "false"
    assert failures[0]["attempted_lmax"] == "2;3;4;5"
    hashes = verify_p1_5_derivations(roots[0])
    assert hashes == summaries[0].artifact_sha256
    assert hashes == {
        str(PILOT_OUTPUT_RELATIVE / name): sha256(payload).hexdigest()
        for name, payload in published[0].items()
    }

    calls_before = list(solvers[0].orders)
    with pytest.raises(FileExistsError, match="second execution"):
        _execute(roots[0], solvers[0])
    assert solvers[0].orders == calls_before


@pytest.mark.parametrize(
    ("error", "decision", "stop_reason"),
    [
        (MemoryError("injected memory limit"), "NO_GO_P1.6_RESOURCE_LIMIT", "memory_limit_exceeded"),
        (RuntimeError("injected infrastructure"), "INCONCLUSIVE_P1.5", "infrastructure_failure"),
    ],
)
def test_resource_and_infrastructure_failures_are_serialized(
    tmp_path: Path,
    error: Exception,
    decision: str,
    stop_reason: str,
) -> None:
    root = _pilot_root(tmp_path, stop_reason)
    solver = FailingSolver(error)
    summary = _execute(root, solver)

    assert summary.decision == decision
    assert summary.stop_reason == stop_reason
    assert summary.attempted_lmax == (2,)
    assert not summary.eligible
    failures = _read(root / PILOT_OUTPUT_RELATIVE / "failures.csv")
    assert failures[0]["decision"] == decision
    assert failures[0]["stop_reason"] == stop_reason
    assert failures[0]["include_in_scientific_tables"] == "false"


def test_unconfirmed_at_21_is_valid_resource_evidence(tmp_path: Path) -> None:
    root = _pilot_root(tmp_path, "unconfirmed")
    summary = _execute(root, NeverConvergingSolver())

    assert summary.decision == "GO_P1.6A_BLIND_FREEZE"
    assert summary.stop_reason == "unconfirmed_at_21"
    assert summary.attempted_lmax == tuple(range(2, 22))
    assert summary.evaluated_lmax == tuple(range(2, 22))
    assert summary.final_lmax == 21
    assert not summary.converged
    assert not summary.eligible
    failure = _read(root / PILOT_OUTPUT_RELATIVE / "failures.csv")[0]
    assert failure["failure_reason"] == "unconfirmed_at_21"
    assert failure["include_in_scientific_tables"] == "false"


def test_invalid_manifest_hash_blocks_runner_before_solver(tmp_path: Path) -> None:
    root = _pilot_root(tmp_path, "invalid")
    path = root / "campaigns" / "p1" / "pilot_manifest.yaml"
    path.write_bytes(path.read_bytes().replace(b'"title": ', b'"title":  ', 1))
    solver = FakeStableSolver()

    with pytest.raises(Exception, match="stored hash"):
        _execute(root, solver)
    assert solver.orders == []
    assert not (root / PILOT_OUTPUT_RELATIVE).exists()


def test_atomic_directory_is_removed_when_publication_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output = tmp_path / "pilot"
    artifacts = {name: name.encode("ascii") for name in pilot_module._ARTIFACT_NAMES}
    original = pilot_module._write_bytes_fsynced
    count = 0

    def fail_on_second(path: Path, payload: bytes) -> None:
        nonlocal count
        count += 1
        if count == 2:
            raise OSError("injected publication failure")
        original(path, payload)

    monkeypatch.setattr(pilot_module, "_write_bytes_fsynced", fail_on_second)
    with pytest.raises(OSError, match="injected publication failure"):
        pilot_module._publish_directory_atomic(output, artifacts)

    assert not output.exists()
    assert list(tmp_path.iterdir()) == []


def test_incomplete_or_mutated_derivation_is_rejected(tmp_path: Path) -> None:
    root = _pilot_root(tmp_path, "mutated")
    _execute(root, FakeStableSolver())
    path = root / PILOT_OUTPUT_RELATIVE / "data_derived.csv"
    path.write_bytes(path.read_bytes() + b"\n")

    with pytest.raises(PilotExecutionError, match="differs from regeneration"):
        verify_p1_5_derivations(root)
