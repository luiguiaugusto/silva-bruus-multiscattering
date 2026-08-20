"""No-solver tests for the P1.6A.2 command-line worktree guard."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest

import acoustic_ms.p1_campaign as campaign_module
from acoustic_ms.p1_campaign import CampaignRunSummary


ROOT = Path(__file__).resolve().parents[1]
CLI_PATH = ROOT / "scripts" / "run_p1_6_campaign.py"
NUMERIC_ENVIRONMENT = {
    "OPENBLAS_NUM_THREADS": "1",
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "NUMEXPR_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
    "PYTHONHASHSEED": "0",
}


def _git(root: Path, *arguments: str) -> str:
    return subprocess.run(
        ["git", *arguments],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def _load_cli():
    name = "p1_6_campaign_cli_under_test"
    specification = importlib.util.spec_from_file_location(name, CLI_PATH)
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


@pytest.fixture
def cli_repository(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    root = tmp_path / "repository"
    root.mkdir()
    _git(root, "init", "-b", "agent/p1-6b-execute")
    _git(root, "config", "user.name", "P1.6A Test")
    _git(root, "config", "user.email", "p1-6a@example.invalid")
    (root / "tracked.txt").write_text("frozen\n", encoding="utf-8")
    _git(root, "add", "tracked.txt")
    _git(root, "commit", "-m", "fixture")

    cli = _load_cli()
    monkeypatch.setattr(cli, "ROOT", root)
    monkeypatch.setattr(
        cli,
        "STATE",
        root / "campaigns" / "p1" / ".p1_6_checkpoint",
    )
    for key, value in NUMERIC_ENVIRONMENT.items():
        monkeypatch.setenv(key, value)
    return cli, root


def _write_started_ledger(root: Path) -> Path:
    ledger = (
        root
        / "campaigns"
        / "p1"
        / ".p1_6_checkpoint"
        / "campaign_ledger.json"
    )
    ledger.parent.mkdir(parents=True)
    ledger.write_text(
        json.dumps({"cases": [{"state": "started"}]}) + "\n",
        encoding="utf-8",
    )
    return ledger


def test_first_execution_rejects_any_dirty_worktree(cli_repository) -> None:
    cli, root = cli_repository
    (root / "untracked.txt").write_text("dirty\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="absolutely clean"):
        cli._require_execution_environment()


def test_resume_accepts_only_the_exact_checkpoint_directory(cli_repository) -> None:
    cli, root = cli_repository
    _write_started_ledger(root)
    case_path = cli.STATE / "cases" / "001.json"
    case_path.parent.mkdir()
    case_path.write_text('{"state":"started"}\n', encoding="utf-8")

    cli._require_execution_environment()
    assert cli._git_status_entries() == (
        (
            "??",
            (
                "campaigns/p1/.p1_6_checkpoint/campaign_ledger.json",
            ),
        ),
        ("??", ("campaigns/p1/.p1_6_checkpoint/cases/001.json",)),
    )


def test_started_checkpoint_reaches_resume_through_cli_without_solver(
    cli_repository,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    cli, root = cli_repository
    ledger = _write_started_ledger(root)
    reached_runner: list[Path] = []
    solver_calls: list[str] = []

    def forbidden_solver(case, manifest):
        del case, manifest
        solver_calls.append("called")
        raise AssertionError("worktree verification must not solve")

    def fake_run(repository, *, executor):
        del executor
        assert json.loads(ledger.read_text(encoding="utf-8"))["cases"][0][
            "state"
        ] == "started"
        reached_runner.append(Path(repository))
        return CampaignRunSummary(
            attempted_this_run=(),
            completed_count=0,
            interrupted_count=1,
            never_started_count=101,
            accumulated_wall_seconds=1800.0,
            closed=False,
            stop_reason="resumed_after_interruption",
            campaign_decision=None,
        )

    monkeypatch.setattr(
        campaign_module,
        "execute_model_e_case_with_limits",
        forbidden_solver,
    )
    monkeypatch.setattr(campaign_module, "run_p1_6_campaign", fake_run)
    monkeypatch.setattr(sys, "argv", [str(CLI_PATH), "--execute"])

    cli.main()

    assert reached_runner == [root]
    assert solver_calls == []
    payload = json.loads(capsys.readouterr().out)
    assert payload["run"]["attempted_this_run"] == []
    assert payload["run"]["interrupted_count"] == 1
    assert payload["run"]["never_started_count"] == 101
    assert payload["run"]["stop_reason"] == "resumed_after_interruption"


def test_resume_rejects_untracked_file_outside_checkpoint(cli_repository) -> None:
    cli, root = cli_repository
    _write_started_ledger(root)
    (root / "outside.csv").write_text("response\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="exact .*checkpoint"):
        cli._require_execution_environment()


@pytest.mark.parametrize("staged", [False, True])
def test_resume_rejects_tracked_modified_or_staged_file(
    cli_repository, staged: bool
) -> None:
    cli, root = cli_repository
    _write_started_ledger(root)
    (root / "tracked.txt").write_text("changed\n", encoding="utf-8")
    if staged:
        _git(root, "add", "tracked.txt")

    with pytest.raises(RuntimeError, match="exact .*checkpoint"):
        cli._require_execution_environment()


def test_resume_rejects_preexisting_confirmatory_csv(cli_repository) -> None:
    cli, root = cli_repository
    _write_started_ledger(root)
    output = root / "campaigns" / "p1" / "data_raw.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("response\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="preexisting confirmatory output"):
        cli._require_execution_environment()


def test_wrong_branch_remains_rejected(cli_repository) -> None:
    cli, root = cli_repository
    _git(root, "checkout", "-b", "agent/p1-6a-audit")

    with pytest.raises(RuntimeError, match="outside a P1.6B branch"):
        cli._require_execution_environment()


def test_resume_rejects_changed_numeric_environment_before_solver(
    cli_repository, monkeypatch: pytest.MonkeyPatch
) -> None:
    cli, root = cli_repository
    _write_started_ledger(root)
    monkeypatch.setenv("OMP_NUM_THREADS", "2")

    with pytest.raises(RuntimeError, match="frozen numeric environment"):
        cli._require_execution_environment()


def test_porcelain_z_parser_preserves_whitespace_in_path(cli_repository) -> None:
    cli, root = cli_repository
    unusual = "space and\nnewline.txt"
    (root / unusual).write_text("dirty\n", encoding="utf-8")

    assert cli._git_status_entries() == (("??", (unusual,)),)


def test_checkpoint_directory_is_not_ignored() -> None:
    ignore = (ROOT / ".gitignore").read_text(encoding="utf-8")
    assert ".p1_6_checkpoint" not in ignore
