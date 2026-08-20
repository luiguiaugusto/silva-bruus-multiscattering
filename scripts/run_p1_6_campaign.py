#!/usr/bin/env python3
"""Execute P1.6B once or regenerate P1.6 artifacts without solver imports."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
STATE = ROOT / "campaigns" / "p1" / ".p1_6_checkpoint"
THREAD_ENVIRONMENT = (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "PYTHONHASHSEED",
)


def _git(*arguments: str) -> str:
    return subprocess.run(
        ["git", *arguments],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def _artifact_module():
    """Load the pure derivative module without importing ``acoustic_ms``."""

    path = ROOT / "src" / "acoustic_ms" / "p1_campaign_artifacts.py"
    specification = importlib.util.spec_from_file_location(
        "p1_campaign_artifacts_standalone", path
    )
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot load standalone P1.6 artifact module")
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    specification.loader.exec_module(module)
    return module


def _regenerate(*, publish: bool) -> dict[str, object]:
    artifacts_module = _artifact_module()
    manifest = json.loads(
        (ROOT / "campaigns" / "p1" / "campaign_manifest.yaml").read_text(
            encoding="utf-8"
        )
    )
    ledger, records = artifacts_module.load_campaign_checkpoint(STATE)
    provenance = ledger["execution_provenance"]
    first, first_gate = artifacts_module.build_campaign_artifacts(
        manifest, records, provenance
    )
    second, second_gate = artifacts_module.build_campaign_artifacts(
        manifest, records, provenance
    )
    if first != second or first_gate != second_gate:
        raise RuntimeError("two no-solver P1.6 regenerations differ byte-for-byte")
    hashes = artifacts_module.artifact_sha256(first)
    if publish:
        published = artifacts_module.publish_campaign_artifacts(ROOT, first)
        if published != hashes:
            raise RuntimeError("published P1.6 hashes differ from regeneration")
    else:
        for name, relative in artifacts_module.ARTIFACT_PATHS.items():
            path = ROOT / relative
            if path.exists() and path.read_bytes() != first[name]:
                raise RuntimeError(f"published {relative} differs from regeneration")
    return {"g1": asdict(first_gate), "artifact_sha256": hashes}


def _require_execution_environment() -> None:
    expected = {
        key: "0" if key == "PYTHONHASHSEED" else "1"
        for key in THREAD_ENVIRONMENT
    }
    invalid = {
        key: os.environ.get(key)
        for key, expected_value in expected.items()
        if os.environ.get(key) != expected_value
    }
    if invalid:
        raise RuntimeError(
            f"P1.6 requires the frozen numeric environment {expected}; got {invalid}"
        )
    branch = _git("branch", "--show-current")
    if not branch.startswith("agent/p1-6b"):
        raise RuntimeError(
            f"confirmatory execution is refused outside a P1.6B branch: {branch!r}"
        )
    if _git("status", "--porcelain", "--untracked-files=all"):
        raise RuntimeError("P1.6 execution requires a clean committed worktree")


def main() -> None:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--execute", action="store_true")
    mode.add_argument("--verify-derived", action="store_true")
    arguments = parser.parse_args()

    if arguments.verify_derived:
        print(json.dumps(_regenerate(publish=False), indent=2, sort_keys=True))
        return

    _require_execution_environment()
    sys.path.insert(0, str(ROOT / "src"))
    from acoustic_ms.p1_campaign import (  # noqa: PLC0415
        execute_model_e_case_with_limits,
        run_p1_6_campaign,
    )

    summary = run_p1_6_campaign(
        ROOT,
        executor=execute_model_e_case_with_limits,
    )
    payload: dict[str, object] = {"run": asdict(summary)}
    if summary.closed:
        payload["artifacts"] = _regenerate(publish=True)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
