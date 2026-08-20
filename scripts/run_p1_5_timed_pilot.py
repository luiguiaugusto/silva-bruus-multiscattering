#!/usr/bin/env python3
"""Execute once or deterministically audit the frozen P1.5 resource pilot."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
THREAD_ENVIRONMENT = (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)


def _git(*arguments: str) -> str:
    return subprocess.run(
        ["git", *arguments],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def _require_execution_environment() -> str:
    invalid = {
        key: os.environ.get(key)
        for key in THREAD_ENVIRONMENT
        if os.environ.get(key) != "1"
    }
    if invalid:
        raise RuntimeError(
            f"P1.5 requires every BLAS/thread variable to equal 1: {invalid}"
        )
    branch = _git("branch", "--show-current")
    if branch != "agent/p1-5-timed-pilot":
        raise RuntimeError(f"P1.5 execution refused on branch {branch!r}")
    status = _git("status", "--porcelain", "--untracked-files=all")
    if status:
        raise RuntimeError(
            "P1.5 execution requires a clean pre-solve commit; dirty paths:\n"
            f"{status}"
        )
    return _git("rev-parse", "HEAD")


def main() -> None:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--execute", action="store_true")
    mode.add_argument("--verify-derived", action="store_true")
    arguments = parser.parse_args()

    sys.path.insert(0, str(ROOT / "src"))
    from acoustic_ms.p1_pilot import (  # noqa: PLC0415
        execute_p1_5_pilot,
        verify_p1_5_derivations,
    )

    if arguments.verify_derived:
        print(
            json.dumps(
                verify_p1_5_derivations(ROOT),
                indent=2,
                sort_keys=True,
            )
        )
        return

    source_commit = _require_execution_environment()
    command = shlex.join([sys.executable, *sys.argv])
    summary = execute_p1_5_pilot(
        ROOT,
        source_commit=source_commit,
        command=command,
    )
    print(json.dumps(asdict(summary), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
