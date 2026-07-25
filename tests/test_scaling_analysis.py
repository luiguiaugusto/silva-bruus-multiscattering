"""T06.1 coupling predictors, scaling fits, and data regressions."""

import csv
from pathlib import Path

import numpy as np
import pytest

from acoustic_ms import (
    coupling_eta,
    fit_power_law,
    irregular_quartet,
    linear_quartet,
    maximum_geometric_coupling,
    square_quartet,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
GEOMETRIES = {
    "linear_chain": linear_quartet,
    "square": square_quartet,
    "irregular": irregular_quartet,
}
EXPECTED_CG = {
    "linear_chain": 2.125,
    "square": 2.353553390593274,
    "irregular": 1.996580257145743,
}
EXPECTED_GEOMETRY_FITS = {
    "linear_chain": ((0.9477457629207996, 0.9984681402009536), (1.9209166221002896, 0.9991340421510149)),
    "square": ((0.9487770531356282, 0.9983499381833476), (1.9198426051755009, 0.9988356442244823)),
    "irregular": ((0.9427068521180622, 0.9979957125494455), (1.9128217876764344, 0.9987739440349565)),
}
EXPECTED_GROUPED = {
    ("eta", "Y3"): (0.9464098893914966, 0.9881180078496571, 0.16270474366899418),
    ("eta", "Y4"): (1.9178603383174067, 0.9920672952487812, 0.2688671716039741),
    ("lambda_max", "Y3"): (0.9477200372149447, 0.9927091201707364, 0.12745167707014535),
    ("lambda_max", "Y4"): (1.9207265477838105, 0.9968960305154440, 0.16818445600557255),
}
EXPECTED_BODY_ROWS = (
    (2, "pair", 0.04413376694829994, 0.0, 0.0, 0.0, 0.0),
    (3, "linear_chain", 0.08319737863715175, 0.04308674386016288, 0.0, 0.04308674386016288, 0.0),
    (3, "equilateral", 0.08826749345670055, 0.04617144636179811, 0.0, 0.04617144636179811, 0.0),
    (3, "scalene", 0.06634541632888573, 0.031143703608647644, 0.0, 0.031143703608647644, 0.0),
    (4, "linear_chain", 0.09977255014752409, 0.06573916615217347, 0.0020850021070293536, 0.06427311991532472, 0.0020850021070293536),
    (4, "square", 0.10419668455946075, 0.06689686555910586, 0.0027105144543686835, 0.06418635110473715, 0.0027105144543686835),
    (4, "irregular", 0.07877328473634540, 0.047516541713439955, 0.0014053466386455625, 0.046142664298696304, 0.0014053466386455622),
)


def _rows(name):
    with (DATA / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_coupling_eta_validation_and_zero_contrast():
    assert coupling_eta(1.0, 2.1, 0.0) == 0.0
    assert np.isclose(coupling_eta(1.0, 2.1, -0.8), 0.8 / 2.1**3)
    for radius in (0.0, -1.0, np.nan, np.inf):
        with pytest.raises(ValueError):
            coupling_eta(radius, 2.1, 0.8)
    for distance in (0.0, -1.0, 1.9, np.nan, np.inf):
        with pytest.raises(ValueError):
            coupling_eta(1.0, distance, 0.8)
    for contrast in (-2.1, 1.1, np.nan, np.inf):
        with pytest.raises(ValueError):
            coupling_eta(1.0, 2.1, contrast)


def test_maximum_geometric_coupling_validation():
    good = irregular_quartet(2.1)
    invalid_positions = (
        good[:, :2],
        good[None, :, :],
        good[:1],
        np.vstack((good, [np.nan, 0.0, 0.0])),
    )
    for positions in invalid_positions:
        with pytest.raises(ValueError):
            maximum_geometric_coupling(positions, 1.0, 0.8)
    coincident = good.copy()
    coincident[1] = coincident[0]
    with pytest.raises(ValueError):
        maximum_geometric_coupling(coincident, 1.0, 0.8)
    with pytest.raises(ValueError):
        maximum_geometric_coupling(irregular_quartet(1.9), 1.0, 0.8)


def test_maximum_geometric_coupling_invariances():
    positions = irregular_quartet(2.3)
    expected = maximum_geometric_coupling(positions, 1.0, 0.8)
    order = np.array([3, 1, 0, 2])
    angle = 0.47
    rotation = np.array([
        [np.cos(angle), -np.sin(angle), 0.0],
        [np.sin(angle), np.cos(angle), 0.0],
        [0.0, 0.0, 1.0],
    ])
    assert np.isclose(maximum_geometric_coupling(positions[order], 1.0, 0.8), expected)
    assert np.isclose(maximum_geometric_coupling(positions + [4.0, -3.0, 1.2], 1.0, 0.8), expected)
    assert np.isclose(maximum_geometric_coupling(positions @ rotation.T, 1.0, 0.8), expected)
    assert np.isclose(maximum_geometric_coupling(2.7 * positions, 2.7, 0.8), expected)


@pytest.mark.parametrize("geometry", GEOMETRIES)
def test_canonical_geometric_factors(geometry):
    eta = coupling_eta(1.0, 2.1, 0.8)
    value = maximum_geometric_coupling(GEOMETRIES[geometry](2.1), 1.0, 0.8)
    assert np.isclose(value / eta, EXPECTED_CG[geometry], rtol=5e-12, atol=5e-13)


def test_fit_power_law_exact_synthetic_data():
    x = np.geomspace(1e-3, 2.0, 40)
    fit = fit_power_law(x, 2.7 * x**1.35)
    assert fit.point_count == 40
    assert np.isclose(fit.prefactor, 2.7, rtol=2e-14)
    assert np.isclose(fit.exponent, 1.35, rtol=2e-14)
    assert np.isclose(fit.r_squared_log, 1.0, atol=2e-15)
    assert fit.rmse_log < 5e-15
    assert fit.max_abs_log_residual < 1e-14


@pytest.mark.parametrize(
    "x,y",
    (
        ([[1.0, 2.0]], [1.0, 2.0]),
        ([1.0, 2.0], [[1.0, 2.0]]),
        ([1.0], [2.0]),
        ([1.0, 2.0], [1.0]),
        ([1.0, np.nan], [1.0, 2.0]),
        ([1.0, 2.0], [1.0, np.inf]),
        ([0.0, 2.0], [1.0, 2.0]),
        ([1.0, 2.0], [-1.0, 2.0]),
        ([1.0, 1.0], [1.0, 2.0]),
        ([1.0, 2.0], [3.0, 3.0]),
    ),
)
def test_fit_power_law_rejects_invalid_data(x, y):
    with pytest.raises(ValueError):
        fit_power_law(x, y)


def test_t06_sweep_integrity():
    rows = _rows("t06_quartet_sweep.csv")
    assert len(rows) == 1920
    assert {row["geometry"] for row in rows} == set(GEOMETRIES)
    assert {float(row["f1"]) for row in rows} == {0.1, 0.4, 0.8, 1.0}
    keys = [(row["geometry"], float(row["f1"]), float(row["d_min_over_a"])) for row in rows]
    assert len(keys) == len(set(keys))
    for geometry in GEOMETRIES:
        selected_geometry = [row for row in rows if row["geometry"] == geometry]
        assert len(selected_geometry) == 640
        for f1 in (0.1, 0.4, 0.8, 1.0):
            selected = [row for row in selected_geometry if float(row["f1"]) == f1]
            assert len(selected) == 160
    for row in rows:
        assert float(row["relative_three_body_sum_amplitude"]) > 0.0
        assert float(row["relative_four_body_amplitude"]) > 0.0


def test_raw_sweep_reproduces_all_fit_regressions():
    rows = _rows("t06_quartet_sweep.csv")
    augmented = []
    for row in rows:
        geometry = row["geometry"]
        distance = float(row["d_min_over_a"])
        f1 = float(row["f1"])
        augmented.append({
            "geometry": geometry,
            "eta": coupling_eta(1.0, distance, f1),
            "lambda_max": maximum_geometric_coupling(
                GEOMETRIES[geometry](distance), 1.0, f1
            ),
            "Y3": float(row["relative_three_body_sum_amplitude"]),
            "Y4": float(row["relative_four_body_amplitude"]),
        })

    for group in (*GEOMETRIES, "grouped"):
        selected = (
            augmented
            if group == "grouped"
            else [row for row in augmented if row["geometry"] == group]
        )
        for predictor in ("eta", "lambda_max"):
            x = np.array([row[predictor] for row in selected])
            for response in ("Y3", "Y4"):
                fit = fit_power_law(x, [row[response] for row in selected])
                if group == "grouped":
                    expected = EXPECTED_GROUPED[(predictor, response)]
                    assert np.allclose(
                        (fit.exponent, fit.r_squared_log, fit.rmse_log),
                        expected,
                        rtol=5e-12,
                        atol=5e-13,
                    )
                else:
                    index = 0 if response == "Y3" else 1
                    exponent, r2 = EXPECTED_GEOMETRY_FITS[group][index]
                    assert np.isclose(
                        fit.exponent, exponent, rtol=5e-12, atol=5e-13
                    )
                    assert np.isclose(
                        fit.r_squared_log, r2, rtol=5e-12, atol=5e-13
                    )


def test_scaling_fit_regressions_and_within_geometry_equivalence():
    rows = _rows("t06_1_scaling_fits.csv")
    assert len(rows) == 16
    indexed = {(row["group"], row["predictor"], row["response"]): row for row in rows}
    for geometry, expected in EXPECTED_GEOMETRY_FITS.items():
        for index, response in enumerate(("Y3", "Y4")):
            eta = indexed[(geometry, "eta", response)]
            lam = indexed[(geometry, "lambda_max", response)]
            exponent, r2 = expected[index]
            assert np.isclose(float(eta["exponent"]), exponent, rtol=5e-12, atol=5e-13)
            assert np.isclose(float(eta["r_squared_log"]), r2, rtol=5e-12, atol=5e-13)
            for field in ("exponent", "r_squared_log", "rmse_log", "max_abs_log_residual"):
                assert np.isclose(float(eta[field]), float(lam[field]), rtol=5e-12, atol=5e-13)
            assert not np.isclose(float(eta["prefactor"]), float(lam["prefactor"]), rtol=1e-3)
    for key, expected in EXPECTED_GROUPED.items():
        row = indexed[("grouped", *key)]
        actual = tuple(float(row[field]) for field in ("exponent", "r_squared_log", "rmse_log"))
        assert np.allclose(actual, expected, rtol=5e-12, atol=5e-13)


def test_collapse_reductions_and_body_order_regressions():
    collapse = _rows("t06_1_collapse_summary.csv")
    assert len(collapse) == 2
    expected_reductions = {"Y3": 0.2166689538601746, "Y4": 0.3744700961361747}
    for row in collapse:
        assert np.isclose(float(row["relative_rmse_reduction"]), expected_reductions[row["response"]], rtol=5e-12, atol=5e-13)

    body = _rows("t06_1_body_order_summary.csv")
    assert len(body) == 7
    fields = (
        "rms_a_vs_c",
        "rms_b_vs_c",
        "rms_up_to_three_vs_c",
        "relative_three_body_amplitude",
        "relative_four_body_amplitude",
    )
    for row, expected in zip(body, EXPECTED_BODY_ROWS, strict=True):
        count, geometry, *values = expected
        assert int(row["particle_count"]) == count
        assert row["geometry"] == geometry
        assert np.allclose([float(row[field]) for field in fields], values, rtol=5e-12, atol=5e-13)
        assert np.isfinite(float(row["rms_model_c"]))
        assert float(row["rms_model_c"]) > 0.0
