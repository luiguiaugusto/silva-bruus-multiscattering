"""Immutable, no-solver checks for the single P1.5 pilot execution."""

from __future__ import annotations

import csv
from datetime import datetime
from hashlib import sha256
import json
from pathlib import Path

from acoustic_ms.p1_pilot import verify_p1_5_derivations
from acoustic_ms.paper_pipeline import validate_manifest_file


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "campaigns" / "p1" / "pilot"
PILOT_MANIFEST = ROOT / "campaigns" / "p1" / "pilot_manifest.yaml"
CONFIRMATORY_MANIFEST = ROOT / "campaigns" / "p1" / "campaign_manifest.yaml"
SOURCE_COMMIT = "a5a2a9c58f5e65b7986e24c7c64879246d946131"
ARTIFACT_SHA256 = {
    "campaigns/p1/pilot/data_derived.csv": "ccd1f7a1aac92a25c51dfb822530fc55f2c27d77c0d85bbcc4215397a3bf2026",
    "campaigns/p1/pilot/data_plot.csv": "9a94fb1203ae122f89a3eb3f49074bea89f1fd7879224b58c6af7e9cafc38424",
    "campaigns/p1/pilot/data_raw.csv": "a4416cae58654371ddcf680ce1a8470ab227c58760b8e1d507893e91883574da",
    "campaigns/p1/pilot/failures.csv": "cd60af766e7340aa04e1b3a1fb2f4b7948f7901163ea75fb5ac42ef4e93e3e8f",
    "campaigns/p1/pilot/performance.csv": "9bb573c524a31a183856289610dc91478d14c8cacb254c9f3a85a2ad00048222",
}


def _read(name: str) -> list[dict[str, str]]:
    with (OUTPUT / name).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def test_artifact_hashes_and_two_regenerations_are_exact() -> None:
    assert {path.name for path in OUTPUT.iterdir()} == {
        "data_raw.csv",
        "data_derived.csv",
        "data_plot.csv",
        "failures.csv",
        "performance.csv",
    }
    observed = {
        str(path.relative_to(ROOT)): sha256(path.read_bytes()).hexdigest()
        for path in sorted(OUTPUT.iterdir())
    }
    assert observed == ARTIFACT_SHA256
    assert verify_p1_5_derivations(ROOT) == ARTIFACT_SHA256


def test_real_attempt_retains_all_orders_channels_and_diagnostics() -> None:
    raw = _read("data_raw.csv")
    performance = _read("performance.csv")

    assert len(raw) == 20 * 4 * 2
    assert len(performance) == 20
    assert [int(row["lmax"]) for row in performance] == list(range(2, 22))
    assert {row["force_channel"] for row in raw} == {
        "total",
        "interaction",
        "external_scattered",
        "scattered_scattered",
    }
    assert all(row["classification"] == "development" for row in raw)
    assert all(row["include_in_scientific_tables"] == "false" for row in raw)
    assert all(row["git_commit"] == SOURCE_COMMIT for row in raw)
    assert all(row["attempt_success"] == "true" for row in performance)
    assert all(row["finite"] == "true" for row in performance)
    assert all(row["diagnostics_pass"] == "true" for row in performance)
    assert all(
        float(row["assembly_diagnostics_seconds"])
        + float(row["linear_solve_seconds"])
        + float(row["force_postprocess_seconds"])
        <= float(row["order_wall_seconds"]) + 1.0e-12
        for row in performance
    )


def test_resource_budget_environment_and_final_dimension_pass() -> None:
    performance = _read("performance.csv")
    final = performance[-1]
    environment = json.loads(final["environment_json"])

    assert float(final["case_wall_seconds"]) == 494.13323493499774
    assert max(int(row["peak_rss_bytes"]) for row in performance) == 311857152
    assert float(final["case_wall_seconds"]) < 1800.0
    assert max(int(row["peak_rss_bytes"]) for row in performance) < 4 * 1024**3
    assert int(final["full_modes_per_particle"]) == 484
    assert int(final["active_modes_per_particle"]) == 231
    assert int(final["system_dimension"]) == 462
    assert final["production_solver"] == "balanced_sqrt"
    assert final["command"] == (
        "/home/luigui/Documents/silva-bruus-multiscattering/.venv/bin/python "
        "scripts/run_p1_5_timed_pilot.py --execute"
    )
    assert environment["pythonhashseed"] == "0"
    for key in (
        "openblas_num_threads",
        "omp_num_threads",
        "mkl_num_threads",
        "numexpr_num_threads",
        "veclib_maximum_threads",
    ):
        assert environment[key] == "1"
    started = datetime.fromisoformat(final["started_utc"].replace("Z", "+00:00"))
    completed = datetime.fromisoformat(final["completed_utc"].replace("Z", "+00:00"))
    assert completed > started


def test_unconfirmed_outcome_and_final_channel_window_are_preserved() -> None:
    failure = _read("failures.csv")[0]
    raw = _read("data_raw.csv")
    final = {
        row["force_channel"]: row
        for row in raw
        if row["lmax"] == "21" and row["particle_index"] == "0"
    }

    assert failure["attempted_lmax"] == ";".join(str(value) for value in range(2, 22))
    assert failure["evaluated_lmax"] == failure["attempted_lmax"]
    assert failure["final_lmax"] == "21"
    assert failure["converged"] == "false"
    assert failure["eligible"] == "false"
    assert failure["failure_stage"] == "convergence"
    assert failure["failure_reason"] == failure["stop_reason"] == "unconfirmed_at_21"
    assert failure["decision"] == "GO_P1.6A_BLIND_FREEZE"
    assert failure["include_in_scientific_tables"] == "false"
    assert final["total"]["confirmed"] == "true"
    assert final["interaction"]["confirmed"] == "true"
    assert final["external_scattered"]["confirmed"] == "true"
    assert final["scattered_scattered"]["change_applicable"] == "true"
    assert final["scattered_scattered"]["confirmed"] == "false"


def test_resource_derivatives_exclude_scientific_force_metrics() -> None:
    derived = _read("data_derived.csv")
    plot = _read("data_plot.csv")

    assert all("force" not in row["metric"] for row in derived)
    assert all(
        row["inclusion_rule"]
        == "resource_evidence_only;exclude_from_scientific_tables"
        for row in derived
    )
    assert {row["series_id"] for row in plot} == {
        "order_wall_seconds",
        "peak_rss_bytes",
    }


def test_pilot_is_immutable_after_confirmatory_manifest_activation() -> None:
    pilot = validate_manifest_file(PILOT_MANIFEST, kind="campaign")
    confirmatory = validate_manifest_file(CONFIRMATORY_MANIFEST, kind="campaign")

    assert pilot["provenance"]["manifest_sha256"] == (
        "d8f56ce20f6f0821d84fd6f36e1f76c855f63f55d809ba9a7201ba52097a43bf"
    )
    assert confirmatory["provenance"]["manifest_sha256"] == (
        "3a63fd66501f8a7ec967ba26fbb8a46f8219fcd65ef1aca4c3ae999803ace6fe"
    )
    assert len(confirmatory["cases"]) == 102
    assert all(case["enabled"] for case in confirmatory["cases"])
    assert confirmatory["resources"]["limits_status"] == "frozen"
