#!/usr/bin/env python3
"""Read-only status gate for the detached P1.6B-R2 campaign."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
from typing import Mapping


STATE_RELATIVE = Path("campaigns/p1/.p1_6b_r2_checkpoint")
ARTIFACT_RELATIVE = (
    Path("campaigns/p1/p1_6b_r2/data_raw.csv"),
    Path("campaigns/p1/p1_6b_r2/data_derived.csv"),
    Path("campaigns/p1/p1_6b_r2/data_plot.csv"),
    Path("campaigns/p1/p1_6b_r2/failures.csv"),
    Path("campaigns/p1/p1_6b_r2/performance.csv"),
)
VALID_DECISIONS = {"GO_P2", "INCONCLUSIVE_P1", "NO_GO_P2"}
INVALID_DECISION = "INVALID_P1.6B_R2_INFRASTRUCTURE"
EXIT_CODES = {
    "READY_FOR_POSTPROCESSING": 0,
    "RUNNING": 10,
    "FAILED_INVALID": 20,
    "FAILED_OR_INCOMPLETE": 21,
}


def _json(path: Path) -> Mapping[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"{path} must contain a JSON object")
    return value


def assess_terminal_evidence(
    worktree: str | Path,
    stdout_path: str | Path,
) -> str:
    """Classify immutable terminal evidence without invoking a solver."""

    root = Path(worktree)
    ledger_path = root / STATE_RELATIVE / "campaign_ledger.json"
    try:
        ledger = _json(ledger_path)
    except (OSError, ValueError, json.JSONDecodeError):
        return "FAILED_OR_INCOMPLETE"

    entries = ledger.get("cases")
    if not isinstance(entries, list):
        return "FAILED_INVALID"
    completed = sum(
        isinstance(entry, dict) and entry.get("state") == "completed"
        for entry in entries
    )
    interrupted = sum(
        isinstance(entry, dict) and entry.get("state") == "interrupted"
        for entry in entries
    )
    never_started = sum(
        isinstance(entry, dict) and entry.get("state") == "never_started"
        for entry in entries
    )
    decision = ledger.get("campaign_decision")
    if decision == INVALID_DECISION:
        return "FAILED_INVALID"
    if (
        ledger.get("closed") is True
        and completed == 0
        and interrupted == len(entries)
        and len(entries) > 0
    ):
        return "FAILED_INVALID"
    if ledger.get("closed") is not True or decision not in VALID_DECISIONS:
        return "FAILED_OR_INCOMPLETE"
    if completed == 0:
        return "FAILED_OR_INCOMPLETE"

    try:
        payload = _json(Path(stdout_path))
    except (OSError, ValueError, json.JSONDecodeError):
        return "FAILED_OR_INCOMPLETE"
    run = payload.get("run")
    artifacts = payload.get("artifacts")
    if not isinstance(run, dict) or not isinstance(artifacts, dict):
        return "FAILED_OR_INCOMPLETE"
    expected_run = {
        "completed_count": completed,
        "interrupted_count": interrupted,
        "never_started_count": never_started,
        "closed": True,
        "campaign_decision": decision,
    }
    if any(run.get(key) != value for key, value in expected_run.items()):
        return "FAILED_INVALID"

    g1 = artifacts.get("g1")
    published_hashes = artifacts.get("artifact_sha256")
    if (
        not isinstance(g1, dict)
        or g1.get("decision") != decision
        or not isinstance(published_hashes, dict)
    ):
        return "FAILED_INVALID"
    expected_paths = {path.as_posix() for path in ARTIFACT_RELATIVE}
    if set(published_hashes) != expected_paths:
        return "FAILED_INVALID"
    for relative in ARTIFACT_RELATIVE:
        path = root / relative
        try:
            digest = sha256(path.read_bytes()).hexdigest()
        except OSError:
            return "FAILED_OR_INCOMPLETE"
        if published_hashes.get(relative.as_posix()) != digest:
            return "FAILED_INVALID"
    return "READY_FOR_POSTPROCESSING"


def _service_properties(unit: str) -> dict[str, str]:
    completed = subprocess.run(
        [
            "systemctl",
            "--user",
            "show",
            unit,
            "--property=LoadState",
            "--property=SubState",
            "--property=MainPID",
            "--property=ExecMainStatus",
            "--property=InvocationID",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        return {}
    properties: dict[str, str] = {}
    for line in completed.stdout.splitlines():
        key, separator, value = line.partition("=")
        if separator:
            properties[key] = value
    return properties


def _process_active() -> bool:
    patterns = (
        "[r]un_p1_6_campaign.py --execute",
        "[p]1_6b_r2_launcher.sh",
        "[t]imeout --signal=TERM --kill-after=30s 68400s",
    )
    return any(
        subprocess.run(
            ["pgrep", "-f", "--", pattern],
            check=False,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        ).returncode
        == 0
        for pattern in patterns
    )


def determine_status(
    *,
    unit: str,
    worktree: str | Path,
    stdout_path: str | Path,
) -> str:
    properties = _service_properties(unit)
    if not properties or properties.get("LoadState") != "loaded":
        return "FAILED_OR_INCOMPLETE"
    active = _process_active()
    if (
        properties.get("SubState") == "running"
        or properties.get("MainPID") not in (None, "0")
        or active
    ):
        return "RUNNING"

    evidence = assess_terminal_evidence(worktree, stdout_path)
    if evidence == "FAILED_INVALID":
        return evidence
    if (
        properties.get("SubState") != "exited"
        or properties.get("MainPID") != "0"
        or properties.get("ExecMainStatus") != "0"
    ):
        return "FAILED_OR_INCOMPLETE"
    return evidence


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--unit", default="p1-6b-r2.service")
    parser.add_argument("--worktree", required=True)
    parser.add_argument("--stdout", required=True)
    arguments = parser.parse_args()
    state = determine_status(
        unit=arguments.unit,
        worktree=arguments.worktree,
        stdout_path=arguments.stdout,
    )
    print(state)
    return EXIT_CODES[state]


if __name__ == "__main__":
    raise SystemExit(main())
