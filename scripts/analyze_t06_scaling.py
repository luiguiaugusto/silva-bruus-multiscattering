"""Post-process existing T05/T06 data for connected-body scaling diagnostics."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

from acoustic_ms import (
    coupling_eta,
    fit_power_law,
    irregular_quartet,
    linear_quartet,
    maximum_geometric_coupling,
    nodal_pair_forces,
    rms_relative_error,
    rms_vector_magnitude,
    solve_rayleigh_nodal_interaction_forces,
    square_quartet,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"
T05_REGRESSION = DATA / "t05_trimer_regression.csv"
T05_SWEEP = DATA / "t05_trimer_sweep.csv"
T06_REGRESSION = DATA / "t06_quartet_regression.csv"
T06_SWEEP = DATA / "t06_quartet_sweep.csv"
GEOMETRIES = {
    "linear_chain": linear_quartet,
    "square": square_quartet,
    "irregular": irregular_quartet,
}
F1_VALUES = (0.1, 0.4, 0.8, 1.0)
DISTANCES = np.linspace(2.1, 10.0, 160)
RESPONSES = {
    "Y3": ("relative_three_body_sum_amplitude", 1),
    "Y4": ("relative_four_body_amplitude", 2),
}
EXPECTED_GEOMETRY_FITS = {
    "linear_chain": ((0.9477457629207996, 0.9984681402009536), (1.9209166221002896, 0.9991340421510149)),
    "square": ((0.9487770531356282, 0.9983499381833476), (1.9198426051755009, 0.9988356442244823)),
    "irregular": ((0.9427068521180622, 0.9979957125494455), (1.9128217876764344, 0.9987739440349565)),
}
EXPECTED_GROUPED = {
    "eta": ((0.9464098893914966, 0.9881180078496571, 0.16270474366899418), (1.9178603383174067, 0.9920672952487812, 0.2688671716039741)),
    "lambda_max": ((0.9477200372149447, 0.9927091201707364, 0.12745167707014535), (1.9207265477838105, 0.9968960305154440, 0.16818445600557255)),
}
EXPECTED_CG = {
    "linear_chain": 2.125,
    "square": 2.353553390593274,
    "irregular": 1.996580257145743,
}


def _read_csv(path, expected_rows, required_fields):
    if not path.is_file():
        raise RuntimeError(f"missing input CSV: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fields = set(reader.fieldnames or ())
    if len(rows) != expected_rows:
        raise RuntimeError(f"{path} must contain {expected_rows} data rows")
    missing = set(required_fields) - fields
    if missing:
        raise RuntimeError(f"{path} is missing fields: {sorted(missing)}")
    return rows


def _finite_float(row, field):
    try:
        value = float(row[field])
    except (KeyError, ValueError) as exc:
        raise RuntimeError(f"invalid numeric field {field}") from exc
    if not np.isfinite(value):
        raise RuntimeError(f"non-finite numeric field {field}")
    return value


def _validate_inputs():
    t05_regression = _read_csv(
        T05_REGRESSION, 9,
        {"geometry", "particle", "d_min_over_a", "c_x", "c_y"},
    )
    t05_sweep = _read_csv(
        T05_SWEEP, 1920,
        {"geometry", "f1", "d_min_over_a", "rms_a_vs_c", "rms_b_vs_c", "rms_irreducible_multibody"},
    )
    t06_regression = _read_csv(
        T06_REGRESSION, 12,
        {"geometry", "particle", "d_min_over_a", "c_x", "c_y"},
    )
    t06_sweep = _read_csv(
        T06_SWEEP, 1920,
        {"geometry", "ka", "f0", "f1", "d_min_over_a", "rms_a_vs_c", "rms_b_vs_c", "rms_up_to_three_vs_c", "relative_three_body_sum_amplitude", "relative_four_body_amplitude"},
    )

    for rows, geometries, particles in (
        (t05_regression, ("linear_chain", "equilateral", "scalene"), 3),
        (t06_regression, tuple(GEOMETRIES), 4),
    ):
        actual = []
        for row in rows:
            actual.append((row["geometry"], int(row["particle"])))
            for field in ("ka", "f0", "f1", "d_min_over_a", "c_x", "c_y"):
                _finite_float(row, field)
        expected = [
            (geometry, particle)
            for geometry in geometries
            for particle in range(particles)
        ]
        if actual != expected:
            raise RuntimeError("regression CSV ordering or particle coverage is invalid")

    t05_expected_order = [
        (geometry, f1, float(distance))
        for geometry in ("linear_chain", "equilateral", "scalene")
        for f1 in F1_VALUES
        for distance in DISTANCES
    ]
    t05_actual_order = []
    t05_keys = set()
    for row in t05_sweep:
        key = (
            row["geometry"],
            _finite_float(row, "f1"),
            _finite_float(row, "d_min_over_a"),
        )
        if key in t05_keys:
            raise RuntimeError("duplicate T05 sweep row")
        t05_keys.add(key)
        t05_actual_order.append(key)
        for field in (
            "ka", "f0", "rms_a_vs_c", "rms_b_vs_c",
            "rms_irreducible_multibody",
        ):
            _finite_float(row, field)
    if any(
        actual[:2] != expected[:2]
        or not np.isclose(actual[2], expected[2], rtol=0, atol=2e-15)
        for actual, expected in zip(
            t05_actual_order, t05_expected_order, strict=True
        )
    ):
        raise RuntimeError("T05 sweep order is not deterministic")

    expected_order = [
        (geometry, f1, float(distance))
        for geometry in GEOMETRIES
        for f1 in F1_VALUES
        for distance in DISTANCES
    ]
    actual_order = []
    keys = set()
    counts = {}
    for row in t06_sweep:
        geometry = row["geometry"]
        f1 = _finite_float(row, "f1")
        distance = _finite_float(row, "d_min_over_a")
        if geometry not in GEOMETRIES or f1 not in F1_VALUES:
            raise RuntimeError("unexpected T06 geometry or contrast")
        key = (geometry, f1, distance)
        if key in keys:
            raise RuntimeError("duplicate T06 sweep row")
        keys.add(key)
        actual_order.append(key)
        counts[(geometry, f1)] = counts.get((geometry, f1), 0) + 1
        for field, _ in RESPONSES.values():
            if _finite_float(row, field) <= 0.0:
                raise RuntimeError("body-order amplitudes must be strictly positive")
        for field in ("ka", "f0", "rms_a_vs_c", "rms_b_vs_c", "rms_up_to_three_vs_c"):
            _finite_float(row, field)
    if set(GEOMETRIES) != {key[0] for key in keys}:
        raise RuntimeError("T06 sweep must contain exactly three geometries")
    if any(value != 160 for value in counts.values()) or len(counts) != 12:
        raise RuntimeError("T06 sweep must contain 160 rows per geometry-contrast pair")
    if any(a[:2] != b[:2] or not np.isclose(a[2], b[2], rtol=0, atol=2e-15) for a, b in zip(actual_order, expected_order, strict=True)):
        raise RuntimeError("T06 sweep order is not deterministic")
    return t05_regression, t05_sweep, t06_regression, t06_sweep


def _augment_t06(rows):
    augmented = []
    factors = {}
    for row in rows:
        geometry = row["geometry"]
        f1 = _finite_float(row, "f1")
        distance = _finite_float(row, "d_min_over_a")
        eta = coupling_eta(1.0, distance, f1)
        lambda_max = maximum_geometric_coupling(
            GEOMETRIES[geometry](distance), 1.0, f1
        )
        if eta <= 0.0 or lambda_max <= 0.0:
            raise RuntimeError("coupling predictors must be positive")
        factor = lambda_max / eta
        previous = factors.setdefault(geometry, factor)
        if not np.isclose(factor, previous, rtol=3e-14, atol=3e-15):
            raise RuntimeError("geometry factor must be dilation-invariant")
        augmented.append({
            **row,
            "eta": eta,
            "lambda_max": lambda_max,
            "Y3": _finite_float(row, RESPONSES["Y3"][0]),
            "Y4": _finite_float(row, RESPONSES["Y4"][0]),
        })
    for geometry, expected in EXPECTED_CG.items():
        if not np.isclose(factors[geometry], expected, rtol=5e-12, atol=5e-13):
            raise RuntimeError(f"unexpected geometric factor for {geometry}")
    return augmented, factors


def _fits(rows):
    output = []
    groups = list(GEOMETRIES) + ["grouped"]
    for group in groups:
        selected = rows if group == "grouped" else [row for row in rows if row["geometry"] == group]
        for predictor in ("eta", "lambda_max"):
            x = np.array([row[predictor] for row in selected])
            for response, (_, order) in RESPONSES.items():
                y = np.array([row[response] for row in selected])
                fit = fit_power_law(x, y)
                output.append({
                    "group": group,
                    "predictor": predictor,
                    "response": response,
                    "expected_order": order,
                    "point_count": fit.point_count,
                    "prefactor": fit.prefactor,
                    "exponent": fit.exponent,
                    "exponent_over_expected_order": fit.exponent / order,
                    "r_squared_log": fit.r_squared_log,
                    "rmse_log": fit.rmse_log,
                    "max_abs_log_residual": fit.max_abs_log_residual,
                })
    if len(output) != 16:
        raise RuntimeError("exactly 16 power-law fits are required")
    return output


def _fit_row(rows, group, predictor, response):
    return next(
        row for row in rows
        if row["group"] == group and row["predictor"] == predictor and row["response"] == response
    )


def _validate_fits(fits):
    for geometry, expected_responses in EXPECTED_GEOMETRY_FITS.items():
        for index, response in enumerate(RESPONSES):
            eta = _fit_row(fits, geometry, "eta", response)
            lam = _fit_row(fits, geometry, "lambda_max", response)
            expected_exponent, expected_r2 = expected_responses[index]
            for actual in (eta, lam):
                if not np.isclose(actual["exponent"], expected_exponent, rtol=5e-12, atol=5e-13):
                    raise RuntimeError("per-geometry exponent regression failed")
                if not np.isclose(actual["r_squared_log"], expected_r2, rtol=5e-12, atol=5e-13):
                    raise RuntimeError("per-geometry R2 regression failed")
            for field in ("exponent", "r_squared_log", "rmse_log", "max_abs_log_residual"):
                if not np.isclose(eta[field], lam[field], rtol=5e-12, atol=5e-13):
                    raise RuntimeError("within-geometry predictor invariance failed")
    for predictor, expected_responses in EXPECTED_GROUPED.items():
        for index, response in enumerate(RESPONSES):
            actual = _fit_row(fits, "grouped", predictor, response)
            expected_exponent, expected_r2, expected_rmse = expected_responses[index]
            for field, expected in (("exponent", expected_exponent), ("r_squared_log", expected_r2), ("rmse_log", expected_rmse)):
                if not np.isclose(actual[field], expected, rtol=5e-12, atol=5e-13):
                    raise RuntimeError(f"grouped {predictor} {response} regression failed")


def _collapse_summary(fits):
    rows = []
    for response, (_, order) in RESPONSES.items():
        eta = _fit_row(fits, "grouped", "eta", response)
        lam = _fit_row(fits, "grouped", "lambda_max", response)
        rows.append({
            "response": response,
            "expected_order": order,
            "eta_exponent": eta["exponent"],
            "eta_r_squared_log": eta["r_squared_log"],
            "eta_rmse_log": eta["rmse_log"],
            "eta_max_abs_log_residual": eta["max_abs_log_residual"],
            "lambda_exponent": lam["exponent"],
            "lambda_r_squared_log": lam["r_squared_log"],
            "lambda_rmse_log": lam["rmse_log"],
            "lambda_max_abs_log_residual": lam["max_abs_log_residual"],
            "relative_rmse_reduction": 1.0 - lam["rmse_log"] / eta["rmse_log"],
        })
    expected = (0.2166689538601746, 0.3744700961361747)
    for row, value in zip(rows, expected, strict=True):
        if not np.isclose(row["relative_rmse_reduction"], value, rtol=5e-12, atol=5e-13):
            raise RuntimeError("RMSE-reduction regression failed")
    return rows


def _rms_from_regression(rows, geometry, distance):
    selected = [
        row for row in rows
        if row["geometry"] == geometry
        and np.isclose(_finite_float(row, "d_min_over_a"), distance, rtol=0, atol=2e-15)
    ]
    if not selected:
        return None
    vectors = np.array([[_finite_float(row, "c_x"), _finite_float(row, "c_y")] for row in selected])
    return rms_vector_magnitude(vectors)


def _select_sweep(rows, geometry):
    selected = [
        row for row in rows
        if row["geometry"] == geometry
        and np.isclose(_finite_float(row, "f1"), .8)
        and np.isclose(_finite_float(row, "d_min_over_a"), 2.1, rtol=0, atol=2e-15)
    ]
    if len(selected) != 1:
        raise RuntimeError(f"missing canonical sweep row for {geometry}")
    return selected[0]


def _body_order_summary(t05_regression, t05_sweep, t06_regression, t06_sweep):
    pair_positions = np.array([[-1.05, 0.0, 0.0], [1.05, 0.0, 0.0]])
    pair_a = np.array(nodal_pair_forces(pair_positions[0, :2], pair_positions[1, :2], .1, 1., 1., .8))
    pair_c = solve_rayleigh_nodal_interaction_forces(pair_positions, .1, 1., 1., 0., .8).forces_xy
    rows = [{
        "particle_count": 2,
        "geometry": "pair",
        "rms_a_vs_c": rms_relative_error(pair_c, pair_a),
        "rms_b_vs_c": 0.0,
        "rms_up_to_three_vs_c": 0.0,
        "relative_three_body_amplitude": 0.0,
        "relative_four_body_amplitude": 0.0,
        "rms_model_c": rms_vector_magnitude(pair_c),
        "source": "single_public_pair_evaluation",
    }]
    for geometry in ("linear_chain", "equilateral", "scalene"):
        sweep = _select_sweep(t05_sweep, geometry)
        amplitude = _rms_from_regression(t05_regression, geometry, 2.1)
        source = "t05_regression_and_sweep"
        if amplitude is None:
            amplitude = _finite_float(sweep, "rms_irreducible_multibody") / _finite_float(sweep, "rms_b_vs_c")
            source = "t05_sweep_identity_regression_unavailable_at_2.1"
        b_vs_c = _finite_float(sweep, "rms_b_vs_c")
        rows.append({
            "particle_count": 3,
            "geometry": geometry,
            "rms_a_vs_c": _finite_float(sweep, "rms_a_vs_c"),
            "rms_b_vs_c": b_vs_c,
            "rms_up_to_three_vs_c": 0.0,
            "relative_three_body_amplitude": b_vs_c,
            "relative_four_body_amplitude": 0.0,
            "rms_model_c": amplitude,
            "source": source,
        })
    for geometry in GEOMETRIES:
        sweep = _select_sweep(t06_sweep, geometry)
        amplitude = _rms_from_regression(t06_regression, geometry, 2.1)
        if amplitude is None:
            raise RuntimeError("missing T06 canonical regression")
        rows.append({
            "particle_count": 4,
            "geometry": geometry,
            "rms_a_vs_c": _finite_float(sweep, "rms_a_vs_c"),
            "rms_b_vs_c": _finite_float(sweep, "rms_b_vs_c"),
            "rms_up_to_three_vs_c": _finite_float(sweep, "rms_up_to_three_vs_c"),
            "relative_three_body_amplitude": _finite_float(sweep, "relative_three_body_sum_amplitude"),
            "relative_four_body_amplitude": _finite_float(sweep, "relative_four_body_amplitude"),
            "rms_model_c": amplitude,
            "source": "t06_regression_and_sweep",
        })
    expected = [
        (2, "pair", .04413376694829994, 0., 0., 0., 0.),
        (3, "linear_chain", .08319737863715175, .04308674386016288, 0., .04308674386016288, 0.),
        (3, "equilateral", .08826749345670055, .04617144636179811, 0., .04617144636179811, 0.),
        (3, "scalene", .06634541632888573, .031143703608647644, 0., .031143703608647644, 0.),
        (4, "linear_chain", .09977255014752409, .06573916615217347, .0020850021070293536, .06427311991532472, .0020850021070293536),
        (4, "square", .10419668455946075, .06689686555910586, .0027105144543686835, .06418635110473715, .0027105144543686835),
        (4, "irregular", .07877328473634540, .047516541713439955, .0014053466386455625, .046142664298696304, .0014053466386455622),
    ]
    fields = ("rms_a_vs_c", "rms_b_vs_c", "rms_up_to_three_vs_c", "relative_three_body_amplitude", "relative_four_body_amplitude")
    for row, expected_row in zip(rows, expected, strict=True):
        if (row["particle_count"], row["geometry"]) != expected_row[:2]:
            raise RuntimeError("body-order table ordering failed")
        for field, value in zip(fields, expected_row[2:], strict=True):
            if not np.isclose(row[field], value, rtol=5e-12, atol=5e-13):
                raise RuntimeError(f"body-order regression failed: {row['geometry']} {field}")
    return rows


def _write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _plot(rows, fits, predictor, output):
    figure, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    colors = {"linear_chain": "C0", "square": "C1", "irregular": "C2"}
    styles = {0.1: "-", 0.4: "--", 0.8: ":", 1.0: "-."}
    for axis, (response, (_, order)) in zip(axes, RESPONSES.items(), strict=True):
        for geometry in GEOMETRIES:
            for f1 in F1_VALUES:
                selected = [row for row in rows if row["geometry"] == geometry and _finite_float(row, "f1") == f1]
                selected.sort(key=lambda row: row[predictor])
                x = np.array([row[predictor] for row in selected]) ** order
                y = 100 * np.array([row[response] for row in selected])
                axis.plot(x, y, color=colors[geometry], linestyle=styles[f1], linewidth=1.0)
        grouped = _fit_row(fits, "grouped", predictor, response)
        raw_x = np.array([row[predictor] for row in rows])
        x_line = np.geomspace(np.min(raw_x) ** order, np.max(raw_x) ** order, 300)
        y_line = 100 * grouped["prefactor"] * x_line ** (grouped["exponent"] / order)
        axis.plot(x_line, y_line, color="black", linewidth=2.5, label="grouped fit")
        axis.set_xscale("log")
        axis.set_yscale("log")
        symbol = r"\eta" if predictor == "eta" else r"\Lambda_{\max}"
        axis.set_xlabel(fr"${symbol}$" if order == 1 else fr"${symbol}^2$")
        axis.set_ylabel(fr"$100Y_{3 if response == 'Y3' else 4}$ (%)")
        axis.set_title(fr"RMS amplitude ratio $Y_{3 if response == 'Y3' else 4}$")
        axis.grid(True, which="both", alpha=.25)
    geometry_handles = [Line2D([0], [0], color=color, label=name) for name, color in colors.items()]
    contrast_handles = [Line2D([0], [0], color="0.35", linestyle=style, label=fr"$f_1={f1}$") for f1, style in styles.items()]
    fit_handle = Line2D([0], [0], color="black", linewidth=2.5, label="grouped fit")
    axes[1].legend(handles=geometry_handles + contrast_handles + [fit_handle], fontsize=8, ncol=2)
    figure.suptitle("Connected-body RMS amplitude ratios; descriptive fits only")
    figure.savefig(FIGURES / output, dpi=200)
    plt.close(figure)


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    t05_regression, t05_sweep, t06_regression, t06_sweep = _validate_inputs()
    rows, factors = _augment_t06(t06_sweep)
    fits = _fits(rows)
    _validate_fits(fits)
    collapse = _collapse_summary(fits)
    body_order = _body_order_summary(t05_regression, t05_sweep, t06_regression, t06_sweep)
    _write_csv(DATA / "t06_1_scaling_fits.csv", fits)
    _write_csv(DATA / "t06_1_collapse_summary.csv", collapse)
    _write_csv(DATA / "t06_1_body_order_summary.csv", body_order)
    _plot(rows, fits, "eta", "t06_1_eta_scaling.png")
    _plot(rows, fits, "lambda_max", "t06_1_lambda_scaling.png")
    print("T06.1: 16 fits, 2 collapse rows, 7 body-order rows")
    for geometry in GEOMETRIES:
        print(f"{geometry}: Cg={factors[geometry]:.15g}")


if __name__ == "__main__":
    main()
