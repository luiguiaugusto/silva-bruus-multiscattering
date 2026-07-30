import csv
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"


def _read(name):
    with (DATA / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def test_t10_validation_table_grid_and_order_are_complete():
    rows = _read("t10_mie_rayleigh_validation.csv")
    assert len(rows) == 4 * 101 * 6
    keys = [
        (float(row["f1"]), float(row["ka"]), int(row["ell"]))
        for row in rows
    ]
    assert len(keys) == len(set(keys))
    assert keys == sorted(keys)
    assert {key[0] for key in keys} == {0.1, 0.4, 0.8, 1.0}
    assert {key[2] for key in keys} == set(range(6))
    for f1 in (0.1, 0.4, 0.8, 1.0):
        assert len({ka for contrast, ka, _ in keys if contrast == f1}) == 101


def test_t10_validation_table_finiteness_and_explicit_applicability():
    rows = _read("t10_mie_rayleigh_validation.csv")
    always_finite = (
        "mie_real",
        "mie_imag",
        "rayleigh_real",
        "rayleigh_imag",
        "mie_magnitude",
        "rayleigh_magnitude",
        "complex_relative_error",
        "magnitude_relative_error",
        "absolute_error",
        "velocity_boundary_residual",
        "unitarity_defect",
    )
    for row in rows:
        assert all(np.isfinite(float(row[field])) for field in always_finite)
        rigid = row["is_rigid_limit"] == "True"
        assert (row["material_ratios_applicable"] == "True") is not rigid
        for field in ("density_ratio", "sound_speed_ratio"):
            assert bool(np.isnan(float(row[field]))) == rigid
        phase_applicable = row["phase_difference_applicable"] == "True"
        assert bool(np.isfinite(float(row["phase_difference_rad"]))) == phase_applicable
        pressure_applicable = (
            row["pressure_boundary_residual_applicable"] == "True"
        )
        assert bool(np.isfinite(float(row["pressure_boundary_residual"]))) == pressure_applicable


def test_t10_summary_contains_all_channels_and_rayleigh_slopes():
    rows = _read("t10_mie_rayleigh_summary.csv")
    assert len(rows) == 24
    assert [
        (float(row["f1"]), int(row["ell"])) for row in rows
    ] == [
        (f1, ell) for f1 in (0.1, 0.4, 0.8, 1.0) for ell in range(6)
    ]
    for row in rows:
        assert int(row["point_count"]) == 101
        if int(row["ell"]) > 0:
            assert abs(float(row["asymptotic_relative_error_slope"]) - 2.0) < 2e-4
        else:
            assert float(row["asymptotic_relative_error_slope"]) == 0.0
    assert max(float(row["maximum_velocity_boundary_residual"]) for row in rows) < 2e-16
    assert max(float(row["maximum_unitarity_defect"]) for row in rows) < 2e-23


def test_t10_figure_exists_and_is_nonempty():
    figure = FIGURES / "t10_mie_rayleigh_error.png"
    assert figure.is_file()
    assert figure.stat().st_size > 100_000
