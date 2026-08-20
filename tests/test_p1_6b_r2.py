"""Infrastructure-only gates for the P1.6B-R2 replacement."""

from __future__ import annotations

from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import sys

from acoustic_ms.paper_pipeline import (
    P1_FROZEN_MANIFEST_SHA256,
    P1_HISTORICAL_MANIFEST_SHA256,
    manifest_file_sha256,
    verify_manifest_sha256,
)


ROOT = Path(__file__).resolve().parents[1]
ORIGINAL = ROOT / "campaigns" / "p1" / "campaign_manifest.yaml"
R2 = ROOT / "campaigns" / "p1" / "campaign_manifest_r2.yaml"
R2_LOCK = "a041e07ae93e9a858bad809427039bf593641ad1f9e341ed89b9d91f648f297d"
P1_6A_LOCK = "3a63fd66501f8a7ec967ba26fbb8a46f8219fcd65ef1aca4c3ae999803ace6fe"
STATUS_PATH = ROOT / "scripts" / "p1_6b_r2_status.py"


def _load_status_module():
    name = "p1_6b_r2_status_under_test"
    specification = importlib.util.spec_from_file_location(name, STATUS_PATH)
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    sys.modules[name] = module
    specification.loader.exec_module(module)
    return module


def test_r2_lock_and_scientific_contract_are_exact() -> None:
    original = json.loads(ORIGINAL.read_text(encoding="utf-8"))
    r2 = json.loads(R2.read_text(encoding="utf-8"))

    assert manifest_file_sha256(R2) == R2_LOCK
    assert verify_manifest_sha256(R2.read_bytes(), expected_sha256=R2_LOCK) == R2_LOCK
    assert P1_FROZEN_MANIFEST_SHA256["p1_dimer_confirmatory_r2"] == R2_LOCK
    assert P1_FROZEN_MANIFEST_SHA256["p1_dimer_confirmatory"] == P1_6A_LOCK
    assert (
        P1_HISTORICAL_MANIFEST_SHA256["p1_dimer_confirmatory_p1_6a"]
        == P1_6A_LOCK
    )
    for field in (
        "schema_version",
        "classification",
        "status",
        "physical",
        "geometry",
        "numerical",
        "resources",
        "cases",
    ):
        assert r2[field] == original[field]
    assert r2["campaign_id"] == "p1_dimer_confirmatory_r2"
    assert r2["outputs"] == {
        "raw": "campaigns/p1/p1_6b_r2/data_raw.csv",
        "derived": "campaigns/p1/p1_6b_r2/data_derived.csv",
        "plot_ready": "campaigns/p1/p1_6b_r2/data_plot.csv",
        "failure_log": "campaigns/p1/p1_6b_r2/failures.csv",
    }
    assert [case["case_order"] for case in r2["cases"]] == list(range(1, 103))
    assert len({case["case_id"] for case in r2["cases"]}) == 102
    assert sum(
        case["parameters"]["evidence_role"] == "primary" for case in r2["cases"]
    ) == 96
    assert sum(
        case["parameters"]["evidence_role"] == "rotational_audit"
        for case in r2["cases"]
    ) == 6


def test_status_rejects_zero_completed_and_102_interrupted(tmp_path: Path) -> None:
    status = _load_status_module()
    state = tmp_path / status.STATE_RELATIVE
    state.mkdir(parents=True)
    ledger = {
        "closed": True,
        "campaign_decision": None,
        "cases": [{"state": "interrupted"} for _ in range(102)],
    }
    (state / "campaign_ledger.json").write_text(
        json.dumps(ledger),
        encoding="utf-8",
    )

    assert (
        status.assess_terminal_evidence(tmp_path, tmp_path / "stdout.json")
        == "FAILED_INVALID"
    )


def test_status_validates_decision_counts_and_output_hashes(tmp_path: Path) -> None:
    status = _load_status_module()
    state = tmp_path / status.STATE_RELATIVE
    state.mkdir(parents=True)
    entries = [{"state": "completed"}] + [
        {"state": "never_started"} for _ in range(101)
    ]
    ledger = {
        "closed": True,
        "campaign_decision": "INCONCLUSIVE_P1",
        "cases": entries,
    }
    (state / "campaign_ledger.json").write_text(
        json.dumps(ledger),
        encoding="utf-8",
    )
    hashes: dict[str, str] = {}
    for index, relative in enumerate(status.ARTIFACT_RELATIVE):
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = f"artifact-{index}\n".encode()
        path.write_bytes(payload)
        hashes[relative.as_posix()] = sha256(payload).hexdigest()
    stdout_path = tmp_path / "stdout.json"
    stdout_path.write_text(
        json.dumps(
            {
                "run": {
                    "completed_count": 1,
                    "interrupted_count": 0,
                    "never_started_count": 101,
                    "closed": True,
                    "campaign_decision": "INCONCLUSIVE_P1",
                },
                "artifacts": {
                    "g1": {"decision": "INCONCLUSIVE_P1"},
                    "artifact_sha256": hashes,
                },
            }
        ),
        encoding="utf-8",
    )

    assert (
        status.assess_terminal_evidence(tmp_path, stdout_path)
        == "READY_FOR_POSTPROCESSING"
    )
    (tmp_path / status.ARTIFACT_RELATIVE[0]).write_bytes(b"mutated\n")
    assert status.assess_terminal_evidence(tmp_path, stdout_path) == "FAILED_INVALID"
