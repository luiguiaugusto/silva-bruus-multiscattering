"""Post-revelation integrity and regression tests for T14.1."""

from __future__ import annotations

import csv
from hashlib import sha256
from pathlib import Path

import numpy as np

from acoustic_ms import EXPECTED_LARGE_N_CASE_IDS


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURE = ROOT / "results" / "figures" / "t14_1_large_n_validation.png"
PHASE_A = (
    DATA / "t14_1_large_n_manifest.csv",
    DATA / "t14_1_local_coupling.csv",
    DATA / "t14_1_frozen_predictions.csv",
    DATA / "t14_1_frozen_protocol.csv",
    DATA / "t14_1_prior_artifact_hashes.csv",
)


def _read(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def _truth(value: str) -> bool:
    assert value in {"true", "false"}
    return value == "true"


def _assert_numeric_fields_finite(
    rows: list[dict[str, str]], excluded: set[str]
) -> None:
    for row in rows:
        for field, value in row.items():
            if field in excluded or value in {"", "true", "false"}:
                continue
            assert np.isfinite(float(value)), (field, value)


def test_raw_campaign_identity_order_and_convergence() -> None:
    rows = _read("t14_1_model_e_convergence.csv")
    assert len(rows) == 160
    case_ids = tuple(dict.fromkeys(row["case_id"] for row in rows))
    assert case_ids == EXPECTED_LARGE_N_CASE_IDS
    grouped = {
        case_id: [row for row in rows if row["case_id"] == case_id]
        for case_id in case_ids
    }
    for scale_order, (case_id, case_rows) in enumerate(grouped.items(), 1):
        orders = [int(row["lmax"]) for row in case_rows]
        assert orders == list(range(2, orders[-1] + 1))
        assert all(int(row["scale_order"]) == scale_order for row in case_rows)
        assert all(_truth(row["campaign_complete"]) for row in case_rows)
        final = case_rows[-1]
        assert int(final["final_lmax"]) == orders[-1]
        assert final["stop_reason"] == "all_channels_confirmed"
        for channel in (
            "total", "interaction", "external_scattered", "scattered_scattered"
        ):
            assert _truth(final[f"{channel}_confirmed"])
        assert _truth(final["diagnostics_pass"])


def test_final_orders_and_eligibility_regression() -> None:
    rows = _read("t14_1_case_summary.csv")
    assert len(rows) == 24
    assert tuple(row["case_id"] for row in rows) == EXPECTED_LARGE_N_CASE_IDS
    assert all(_truth(row["eligible"]) for row in rows)
    expected = {
        45: {
            "linear": (7, 8, 9, 11),
            "compact": (6, 6, 7, 8),
            "irregular": (6, 7, 8, 9),
        },
        105: {
            "linear": (7, 8, 9, 11),
            "compact": (6, 6, 7, 8),
            "irregular": (6, 7, 8, 9),
        },
    }
    for particle_count, families in expected.items():
        for family, orders in families.items():
            selected = [
                int(row["final_lmax"])
                for row in rows
                if int(row["particle_count"]) == particle_count
                and row["family"] == family
            ]
            assert tuple(selected) == orders


def test_campaign_numerical_diagnostics() -> None:
    rows = _read("t14_1_model_e_convergence.csv")
    _assert_numeric_fields_finite(
        rows,
        {
            "case_id", "family", "stop_reason", "production_solver",
            "coordinates_xyz", "coordinate_sha256", "local_coupling_sha256",
            "model_a_forces_xyz", "total_forces_xyz", "external_forces_xyz",
            "interaction_forces_xyz", "external_scattered_forces_xyz",
            "scattered_scattered_forces_xyz",
        },
    )
    assert all(row["production_solver"] == "balanced_sqrt" for row in rows)
    assert max(float(row["balanced_condition_number"]) for row in rows) < 10.0
    for field in (
        "balanced_backward_error", "effective_incident_closure_error",
        "scattering_closure_error", "force_decomposition_residual",
        "max_abs_fz",
    ):
        assert max(float(row[field]) for row in rows) < 1.0e-12
    assert np.isclose(
        max(float(row["physical_residual_relative"]) for row in rows),
        0.13760124781707805,
        rtol=2e-14,
        atol=2e-15,
    )


def test_force_table_is_complete_finite_and_ordered() -> None:
    rows = _read("t14_1_forces.csv")
    assert len(rows) == 1800
    expected_keys = [
        (scale_order, case_id, particle_index)
        for scale_order, case_id in enumerate(EXPECTED_LARGE_N_CASE_IDS, 1)
        for particle_index in range(45 if scale_order <= 12 else 105)
    ]
    observed_keys = [
        (int(row["scale_order"]), row["case_id"], int(row["particle_index"]))
        for row in rows
    ]
    assert observed_keys == expected_keys
    _assert_numeric_fields_finite(rows, {"case_id"})


def test_derived_tables_have_stable_counts_and_unique_keys() -> None:
    expected = {
        "t14_1_case_summary.csv": 24,
        "t14_1_large_n_predictions.csv": 48,
        "t14_1_metrics.csv": 20,
        "t14_1_threshold_audit.csv": 18,
        "t14_1_matched_large_n_pairs.csv": 12,
        "t14_1_combined_scale_sequence.csv": 48,
        "t14_1_performance.csv": 24,
        "t14_1_gate.csv": 30,
    }
    for name, count in expected.items():
        rows = _read(name)
        assert len(rows) == count
        assert len({tuple(row.items()) for row in rows}) == count
    assert FIGURE.is_file() and FIGURE.stat().st_size > 0


def test_frozen_m1_metrics_and_threshold_gate() -> None:
    metrics = {
        (row["model"], row["scope_type"], row["scope"]): row
        for row in _read("t14_1_metrics.csv")
    }
    expected = {
        ("M1", "global", "all"): (24, 0.30428741758101263, 1.0, 0.9493149726525656),
        ("M1", "particle_count", "N=45"): (12, 0.2998410355482452, 1.0, 0.9230769230769231),
        ("M1", "particle_count", "N=105"): (12, 0.30866975620815057, 1.0, 0.9370629370629372),
        ("P3", "global", "all"): (24, 0.3601285659000717, 0.9583333333333334, 0.957391304347826),
        ("P3", "particle_count", "N=45"): (12, 0.33540151142667074, 1.0, 0.965034965034965),
        ("P3", "particle_count", "N=105"): (12, 0.3832636091350044, 0.9166666666666666, 0.9720279720279721),
    }
    for key, (count, rmse, factor_two, spearman) in expected.items():
        row = metrics[key]
        assert int(row["point_count"]) == count
        assert np.isclose(float(row["rmse_log"]), rmse, rtol=2e-14, atol=2e-15)
        assert np.isclose(float(row["fraction_within_factor_2"]), factor_two)
        assert np.isclose(float(row["spearman"]), spearman)
    audits = [
        row for row in _read("t14_1_threshold_audit.csv")
        if row["model"] == "M1" and row["scope"] == "all"
    ]
    assert [int(row["predicted_safe_count"]) for row in audits] == [6, 12, 18]
    assert all(int(row["false_safe_count"]) == 0 for row in audits)


def test_gate_and_matched_large_n_trend_are_literal() -> None:
    final = _read("t14_1_gate.csv")[-1]
    assert final["decision"] == "PASS_T14_1_LARGE_N_FROZEN_LAMBDA_MAX"
    assert final["next_gate"] == "GO_T15_SYNTHESIS_AND_MANUSCRIPT"
    assert final["large_n_trend"] == "NO_SYSTEMATIC_DETERIORATION"
    pairs = _read("t14_1_matched_large_n_pairs.csv")
    assert len(pairs) == 12
    assert all(_truth(row["applicable"]) for row in pairs)
    ratios = np.asarray([float(row["ratio_105_over_45"]) for row in pairs])
    assert np.median(ratios) <= 1.10
    assert np.percentile(ratios, 90, method="linear") <= 1.25


def test_phase_a_and_all_prior_artifacts_remain_byte_identical() -> None:
    protocol = {row["key"]: row["value"] for row in _read("t14_1_frozen_protocol.csv")}
    for path in PHASE_A:
        key = f"artifact_sha256:{path.name}"
        if key in protocol:
            assert sha256(path.read_bytes()).hexdigest() == protocol[key]
    prior = _read("t14_1_prior_artifact_hashes.csv")
    assert len(prior) == 95
    for row in prior:
        path = ROOT / row["path"]
        assert path.stat().st_size == int(row["size_bytes"])
        assert sha256(path.read_bytes()).hexdigest() == row["sha256"]


def test_all_t14_1_csv_cells_exclude_nan_and_infinity() -> None:
    forbidden = {"nan", "inf", "+inf", "-inf", "infinity", "+infinity", "-infinity"}
    for path in sorted(DATA.glob("t14_1_*.csv")):
        with path.open(encoding="utf-8", newline="") as stream:
            rows = list(csv.reader(stream))
        assert rows and rows[0]
        assert all(value.strip().lower() not in forbidden for row in rows for value in row)


def test_analysis_is_postprocessing_only() -> None:
    source = (ROOT / "scripts" / "analyze_t14_1_large_n.py").read_text(
        encoding="utf-8"
    )
    assert "solve_model_e_nodal" not in source
    assert "run_model_e_convergence" not in source
