#!/usr/bin/env python3
"""Audit the T11.1 square-root-balanced Model-E formulation."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np

from acoustic_ms import complete_radiation_force_from_bsc, solve_model_e_nodal
from analyze_t11_model_e import _cases, _format, _rms, _write


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"
CHANNELS = (
    ("total", "total_forces_xyz"),
    ("interaction", "interaction_forces_xyz"),
    ("external_scattered", "external_scattered_forces_xyz"),
    ("scattered_scattered", "scattered_scattered_forces_xyz"),
)
ORDERS = tuple(range(2, 10))
MP_DPS = 80


def _difference(first: np.ndarray, second: np.ndarray) -> tuple[float, float, bool]:
    first = np.asarray(first)
    second = np.asarray(second)
    absolute = float(np.linalg.norm(first - second))
    scale = max(float(np.linalg.norm(first)), float(np.linalg.norm(second)))
    tolerance = 128.0 * np.finfo(float).eps * scale
    if scale <= tolerance:
        return absolute, absolute, False
    return absolute / scale, absolute, True


def _full_from_active(solution, values: np.ndarray) -> np.ndarray:
    particle_count = len(solution.effective_incident_coefficients)
    full = np.zeros_like(solution.effective_incident_coefficients)
    full[:, np.asarray(solution.active_mode_indices)] = np.asarray(values).reshape(
        particle_count, len(solution.active_modes)
    )
    return full


def _force_channels(result, effective_active: np.ndarray, k: float) -> dict[str, np.ndarray]:
    solution = result.solution
    effective = _full_from_active(solution, effective_active)
    external_coefficients = solution.external_incident_coefficients
    other = effective - external_coefficients
    particle_count = len(effective)
    total = np.empty((particle_count, 3), dtype=float)
    external = np.empty_like(total)
    scattered_scattered = np.empty_like(total)
    for particle in range(particle_count):
        total[particle] = complete_radiation_force_from_bsc(
            effective[particle], solution.scattering_coefficients, k, 1.0
        )
        external[particle] = complete_radiation_force_from_bsc(
            external_coefficients[particle], solution.scattering_coefficients, k, 1.0
        )
        scattered_scattered[particle] = complete_radiation_force_from_bsc(
            other[particle], solution.scattering_coefficients, k, 1.0
        )
    interaction = total - external
    return {
        "total": total,
        "interaction": interaction,
        "external_scattered": interaction - scattered_scattered,
        "scattered_scattered": scattered_scattered,
    }


def _audit_case(case):
    case_id, geometry, positions, ka, f0, f1 = case
    rows = []
    results = {}
    for order in ORDERS:
        result = solve_model_e_nodal(positions, ka, 1.0, 1.0, f0, f1, order)
        solution = result.solution
        balanced_b = solution.effective_incident_coefficients.ravel()[
            np.concatenate([
                particle * len(solution.modes) + np.asarray(solution.active_mode_indices)
                for particle in range(len(positions))
            ])
        ]
        balanced_d = solution.scattered_coefficients.ravel()[
            np.concatenate([
                particle * len(solution.modes) + np.asarray(solution.active_mode_indices)
                for particle in range(len(positions))
            ])
        ]
        legacy_b = np.linalg.solve(
            solution.effective_incident_system_matrix,
            solution.effective_incident_right_hand_side,
        )
        legacy_d = solution.scattering_diagonal * legacy_b
        scattered_d = np.linalg.solve(
            solution.scattered_system_matrix,
            solution.scattered_right_hand_side,
        )
        scattered_d += np.linalg.solve(
            solution.scattered_system_matrix,
            solution.scattered_right_hand_side
            - solution.scattered_system_matrix @ scattered_d,
        )
        scattered_b = (
            solution.effective_incident_right_hand_side
            + solution.translation_matrix @ scattered_d
        )
        channels = {
            "balanced": {
                short: getattr(result, attribute) for short, attribute in CHANNELS
            },
            "legacy": _force_channels(result, legacy_b, ka),
            "scattered": _force_channels(result, scattered_b, ka),
        }
        row = {
            "case_id": case_id,
            "geometry": geometry,
            "particle_count": len(positions),
            "ka": _format(ka),
            "f0": _format(f0),
            "f1": _format(f1),
            "lmax": order,
            "system_dimension": len(solution.balanced_coefficients),
            "condition_number_effective_incident": _format(solution.condition_number),
            "condition_number_scattered": _format(solution.scattered_condition_number),
            "condition_number_balanced": _format(solution.balanced_condition_number),
            "legacy_residual_relative": _format(solution.residual_relative),
            "balanced_backward_error": _format(solution.balanced_backward_error),
            "effective_incident_closure_error": _format(
                solution.effective_incident_closure_error
            ),
            "scattering_closure_error": _format(solution.scattering_closure_error),
            "production_solver": solution.production_solver,
        }
        for label, first, second in (
            ("balanced_vs_scattered_b", balanced_b, scattered_b),
            ("balanced_vs_scattered_d", balanced_d, scattered_d),
            ("balanced_vs_legacy_b", balanced_b, legacy_b),
            ("balanced_vs_legacy_d", balanced_d, legacy_d),
        ):
            relative, absolute, applicable = _difference(first, second)
            row[f"{label}_relative"] = _format(relative)
            row[f"{label}_absolute"] = _format(absolute)
            row[f"{label}_relative_applicable"] = str(applicable).lower()
        for channel, _ in CHANNELS:
            for comparator in ("scattered", "legacy"):
                relative, absolute, applicable = _difference(
                    channels["balanced"][channel], channels[comparator][channel]
                )
                label = f"{channel}_balanced_vs_{comparator}"
                row[f"{label}_relative"] = _format(relative)
                row[f"{label}_absolute"] = _format(absolute)
                row[f"{label}_relative_applicable"] = str(applicable).lower()
        rows.append(row)
        results[order] = result
    return rows, results


def _mp_matrix(values: np.ndarray) -> mp.matrix:
    return mp.matrix([
        [mp.mpc(float(value.real), float(value.imag)) for value in row]
        for row in np.asarray(values)
    ])


def _mp_vector(values: np.ndarray) -> mp.matrix:
    return mp.matrix([
        mp.mpc(float(value.real), float(value.imag)) for value in np.asarray(values)
    ])


def _numpy_vector(values: mp.matrix) -> np.ndarray:
    return np.array([complex(values[index]) for index in range(len(values))])


def _high_precision_rows(campaign):
    mp.mp.dps = MP_DPS
    rows = []
    sentinels = ("dimer_axis", "trimer_scalene")
    case_by_id = {case[0]: case for case in _cases()}
    for case_id in sentinels:
        case = case_by_id[case_id]
        ka = case[3]
        result = campaign[case_id][9]
        solution = result.solution
        matrix = _mp_matrix(solution.balanced_system_matrix)
        rhs = _mp_vector(solution.balanced_right_hand_side)
        q_mp = mp.lu_solve(matrix, rhs)
        square_root = _mp_vector(solution.square_root_scattering_diagonal)
        d_mp = mp.matrix([
            square_root[index] * q_mp[index] for index in range(len(q_mp))
        ])
        translation = _mp_matrix(solution.translation_matrix)
        incident = _mp_vector(solution.effective_incident_right_hand_side)
        b_mp = incident + translation * d_mp
        q_reference = _numpy_vector(q_mp)
        d_reference = _numpy_vector(d_mp)
        b_reference = _numpy_vector(b_mp)
        balanced_b = solution.effective_incident_coefficients.ravel()[
            np.concatenate([
                particle * len(solution.modes) + np.asarray(solution.active_mode_indices)
                for particle in range(len(result.total_forces_xyz))
            ])
        ]
        balanced_d = solution.scattered_coefficients.ravel()[
            np.concatenate([
                particle * len(solution.modes) + np.asarray(solution.active_mode_indices)
                for particle in range(len(result.total_forces_xyz))
            ])
        ]
        for quantity, computed, reference in (
            ("q", solution.balanced_coefficients, q_reference),
            ("d", balanced_d, d_reference),
            ("b", balanced_b, b_reference),
        ):
            relative, absolute, applicable = _difference(computed, reference)
            rows.append({
                "case_id": case_id,
                "lmax": 9,
                "mpmath_dps": MP_DPS,
                "quantity": quantity,
                "relative_error_or_absolute_if_unresolved": _format(relative),
                "absolute_error": _format(absolute),
                "quantity_resolved": str(applicable).lower(),
            })
        hp_channels = _force_channels(result, b_reference, ka)
        for channel, attribute in CHANNELS:
            relative, absolute, applicable = _difference(
                getattr(result, attribute), hp_channels[channel]
            )
            rows.append({
                "case_id": case_id,
                "lmax": 9,
                "mpmath_dps": MP_DPS,
                "quantity": channel,
                "relative_error_or_absolute_if_unresolved": _format(relative),
                "absolute_error": _format(absolute),
                "quantity_resolved": str(applicable).lower(),
            })
    return rows


def _plot(rows, high_precision_rows):
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8), constrained_layout=True)
    colors = plt.cm.tab10.colors
    case_ids = list(dict.fromkeys(row["case_id"] for row in rows))
    for color, case_id in zip(colors, case_ids):
        case_rows = [row for row in rows if row["case_id"] == case_id]
        orders = [int(row["lmax"]) for row in case_rows]
        axes[0].semilogy(
            orders,
            [float(row["condition_number_effective_incident"]) for row in case_rows],
            color=color, linestyle=":", linewidth=1.1,
        )
        axes[0].semilogy(
            orders,
            [float(row["condition_number_scattered"]) for row in case_rows],
            color=color, linestyle="--", linewidth=1.1,
        )
        axes[0].semilogy(
            orders,
            [float(row["condition_number_balanced"]) for row in case_rows],
            color=color, linewidth=1.4, label=case_id,
        )
        axes[1].semilogy(
            orders,
            [max(float(row["balanced_backward_error"]), np.finfo(float).tiny)
             for row in case_rows],
            color=color, linewidth=1.3,
        )
        axes[1].semilogy(
            orders,
            [max(float(row["effective_incident_closure_error"]), np.finfo(float).tiny)
             for row in case_rows],
            color=color, linestyle="--", linewidth=1.0,
        )
        axes[1].semilogy(
            orders,
            [max(float(row["scattering_closure_error"]), np.finfo(float).tiny)
             for row in case_rows],
            color=color, linestyle=":", linewidth=1.0,
        )
        discrepancy = []
        for row in case_rows:
            values = [
                float(row[f"{channel}_balanced_vs_scattered_relative"])
                for channel, _ in CHANNELS
                if row[f"{channel}_balanced_vs_scattered_relative_applicable"] == "true"
            ]
            discrepancy.append(max(values) if values else np.finfo(float).tiny)
        axes[2].semilogy(
            orders, discrepancy, color=color, linewidth=1.3, label=case_id
        )
    hp_values = [
        float(row["relative_error_or_absolute_if_unresolved"])
        for row in high_precision_rows if row["quantity_resolved"] == "true"
    ]
    if hp_values:
        axes[2].axhline(
            max(hp_values), color="black", linestyle="-.", linewidth=1.2,
            label="max high-precision discrepancy",
        )
    axes[0].set(
        xlabel=r"$L_{\max}$", ylabel="matrix condition number",
        title=r"$A_b$ (:), $A_d$ (--), $A_q$ (-)",
    )
    axes[1].axhline(1e-12, color="black", linestyle="--", linewidth=1.0)
    axes[1].set(
        xlabel=r"$L_{\max}$", ylabel="relative diagnostic",
        title="Backward error and physical closures",
    )
    axes[2].axhline(1e-9, color="black", linestyle="--", linewidth=1.0)
    axes[2].set(
        xlabel=r"$L_{\max}$", ylabel="relative or scaled absolute discrepancy",
        title="Force agreement across formulations",
    )
    for axis in axes:
        axis.grid(True, alpha=0.25)
    axes[0].legend(fontsize=7)
    axes[2].legend(fontsize=7)
    fig.savefig(
        FIGURES / "t11_1_model_e_stability.png",
        dpi=220,
        metadata={"Software": "acoustic_ms T11.1"},
    )
    plt.close(fig)


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    rows = []
    campaign = {}
    for case in _cases():
        case_rows, results = _audit_case(case)
        rows.extend(case_rows)
        campaign[case[0]] = results
    high_precision_rows = _high_precision_rows(campaign)
    _write(DATA / "t11_1_solver_stability.csv", list(rows[0]), rows)
    _write(
        DATA / "t11_1_high_precision_oracle.csv",
        list(high_precision_rows[0]),
        high_precision_rows,
    )
    _plot(rows, high_precision_rows)
    print(
        f"T11.1: {len(rows)} stability rows, "
        f"{len(high_precision_rows)} high-precision rows"
    )


if __name__ == "__main__":
    main()
