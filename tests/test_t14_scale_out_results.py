"""Post-revelation integrity and scientific regressions for T14."""

from __future__ import annotations

import csv
from hashlib import sha256
from pathlib import Path
import subprocess

import numpy as np

from acoustic_ms import EXPECTED_SCALE_OUT_CASE_IDS


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURE = ROOT / "results" / "figures" / "t14_scale_out_validation.png"
PHASE_A = (
    DATA / "t14_scale_manifest.csv",
    DATA / "t14_frozen_predictions.csv",
    DATA / "t14_frozen_protocol.csv",
    DATA / "t14_prior_artifact_hashes.csv",
)
DERIVED = (
    DATA / "t14_forces.csv", DATA / "t14_case_summary.csv",
    DATA / "t14_scale_predictions.csv", DATA / "t14_metrics.csv",
    DATA / "t14_threshold_audit.csv", DATA / "t14_matched_scale_pairs.csv",
    DATA / "t14_performance.csv", DATA / "t14_gate.csv", FIGURE,
)


def _read(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _finite_csv(path: Path) -> None:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert rows
    for row in rows:
        for value in row.values():
            assert value is not None
            assert value.lower() not in {"nan", "inf", "+inf", "-inf"}


def test_phase_a_artifacts_are_identical_to_first_pushed_commit() -> None:
    for path in PHASE_A:
        expected = subprocess.run(
            ["git", "show", f"HEAD:{path.relative_to(ROOT)}"], cwd=ROOT,
            check=True, capture_output=True,
        ).stdout
        assert sha256(path.read_bytes()).digest() == sha256(expected).digest()


def test_raw_campaign_has_exact_cases_orders_and_diagnostics() -> None:
    rows = _read("t14_model_e_convergence.csv")
    assert len(rows) == 162
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(row["case_id"], []).append(row)
        assert int(row["lmax"]) <= 13
        assert row["production_solver"] == "balanced_sqrt"
        assert row["mode_dimension_consistent"] == "true"
        assert row["finite"] == "true"
        assert row["diagnostics_pass"] == "true"
        assert row["planar_symmetry_pass"] == "true"
        assert row["campaign_complete"] == "true"
        assert float(row["balanced_condition_number"]) < 10.0
    assert tuple(grouped) == EXPECTED_SCALE_OUT_CASE_IDS
    for case_rows in grouped.values():
        orders = [int(row["lmax"]) for row in case_rows]
        assert orders == list(range(2, orders[-1] + 1))
        final = case_rows[-1]
        assert final["stop_reason"] == "all_channels_confirmed"
        for channel in ("total", "interaction", "external_scattered", "scattered_scattered"):
            assert final[f"{channel}_confirmed"] == "true"
            assert int(final[f"{channel}_minimum_confirmed_lmax"]) > 0


def test_summary_is_complete_eligible_and_uses_interaction_reference() -> None:
    rows = _read("t14_case_summary.csv")
    assert len(rows) == 24
    assert tuple(row["case_id"] for row in rows) == EXPECTED_SCALE_OUT_CASE_IDS
    assert all(row["eligible"] == "true" for row in rows)
    assert all(row["interaction_confirmed"] == "true" for row in rows)
    assert max(int(row["final_lmax"]) for row in rows) == 11
    assert np.allclose(
        [float(row["epsilon_a_e"]) for row in rows[:4]],
        [0.0023359609785337352, 0.009571497495310029,
         0.022747734540650624, 0.06461686393697703],
        rtol=3e-13, atol=3e-15,
    )
    revealed = _read("t14_scale_predictions.csv")
    assert len(revealed) == 48
    assert {row["response_source"] for row in revealed} == {
        "Model_E_interaction_force_only"
    }


def test_m1_and_p3_metrics_match_revealed_regressions() -> None:
    rows = _read("t14_metrics.csv")
    assert len(rows) == 20
    metrics = {(row["model"], row["scope_type"], row["scope"]): row for row in rows}
    expected = {
        ("M1", "global", "all"): (0.34301137024205131, 1.0, 0.92368507165895797),
        ("M1", "particle_count", "N=15"): (0.37098685064152775, 1.0, 0.94553159559598077),
        ("M1", "particle_count", "N=28"): (0.31254176821944923, 1.0, 0.93664294422131833),
        ("P3", "global", "all"): (0.31440958235577582, 1.0, 0.90434782608695652),
        ("P3", "particle_count", "N=15"): (0.31540615238634706, 1.0, 0.8951048951048951),
        ("P3", "particle_count", "N=28"): (0.31340984348151413, 1.0, 0.91608391608391626),
    }
    for key, values in expected.items():
        row = metrics[key]
        assert np.allclose(
            [float(row["rmse_log"]), float(row["fraction_within_factor_2"]), float(row["spearman"])],
            values, rtol=3e-13, atol=3e-15,
        )


def test_conservative_classifications_and_literal_gate() -> None:
    audits = _read("t14_threshold_audit.csv")
    assert len(audits) == 18
    m1_global = [row for row in audits if row["model"] == "M1" and row["scope"] == "all"]
    assert [int(row["predicted_safe_count"]) for row in m1_global] == [6, 12, 18]
    assert [int(row["false_safe_count"]) for row in m1_global] == [0, 0, 0]
    gate = _read("t14_gate.csv")
    assert len(gate) == 29
    assert all(row["passed"] == "true" for row in gate)
    assert gate[-1]["decision"] == "PASS_T14_SCALE_OUT_FROZEN_LAMBDA_MAX"
    assert gate[-1]["next_gate"] == "GO_T15_SYNTHESIS_AND_MANUSCRIPT"


def test_matched_scale_pairs_are_complete_and_size_effect_is_small() -> None:
    rows = _read("t14_matched_scale_pairs.csv")
    assert len(rows) == 12
    assert len({(row["family"], row["target_level"]) for row in rows}) == 12
    ratios = np.asarray([float(row["ratio_28_over_15"]) for row in rows])
    assert np.all(np.isfinite(ratios))
    assert np.isclose(np.min(ratios), 0.72932901634129876, rtol=0.0, atol=1e-16)
    assert np.isclose(np.max(ratios), 1.0028447919465278, rtol=0.0, atol=1e-16)


def test_all_t14_tables_are_finite_and_have_stable_counts() -> None:
    expected = {
        "t14_scale_manifest.csv": 24, "t14_frozen_predictions.csv": 48,
        "t14_frozen_protocol.csv": 59, "t14_prior_artifact_hashes.csv": 81,
        "t14_model_e_convergence.csv": 162, "t14_forces.csv": 516,
        "t14_case_summary.csv": 24, "t14_scale_predictions.csv": 48,
        "t14_metrics.csv": 20, "t14_threshold_audit.csv": 18,
        "t14_matched_scale_pairs.csv": 12, "t14_performance.csv": 24,
        "t14_gate.csv": 29,
    }
    for name, count in expected.items():
        rows = _read(name)
        assert len(rows) == count
        _finite_csv(DATA / name)
        header = (DATA / name).read_text(encoding="utf-8").splitlines()[0].split(",")
        assert len(header) == len(set(header))


def test_previous_artifact_hash_manifest_is_preserved() -> None:
    for row in _read("t14_prior_artifact_hashes.csv"):
        path = ROOT / row["path"]
        assert path.stat().st_size == int(row["size_bytes"])
        assert sha256(path.read_bytes()).hexdigest() == row["sha256"]


def test_analyze_only_is_model_e_free_and_byte_deterministic() -> None:
    source = (ROOT / "scripts" / "analyze_t14_scale_out.py").read_text(encoding="utf-8")
    assert "solve_model_e_nodal" not in source
    command = [str(ROOT / ".venv/bin/python"), "scripts/run_t14_scale_out.py", "--analyze-only"]
    subprocess.run(command, cwd=ROOT, check=True, capture_output=True)
    first = tuple(sha256(path.read_bytes()).hexdigest() for path in DERIVED)
    subprocess.run(command, cwd=ROOT, check=True, capture_output=True)
    second = tuple(sha256(path.read_bytes()).hexdigest() for path in DERIVED)
    assert first == second
