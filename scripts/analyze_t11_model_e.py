#!/usr/bin/env python3
"""Generate the compact deterministic T11 Model-E validation campaign."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from acoustic_ms import (
    complete_radiation_force_from_bsc,
    equilateral_trimer,
    irregular_quartet,
    scalene_trimer,
    solve_model_e_nodal,
)
from t11_stress_oracle import stress_tensor_force


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"
TOLERANCE = 1e-5
BASE_ORDERS = tuple(range(2, 8))


def _pair(distance: float, angle: float = 0.0) -> np.ndarray:
    direction = np.array([np.cos(angle), np.sin(angle), 0.0])
    return np.array([-0.5 * distance * direction, 0.5 * distance * direction])


def _cases():
    return (
        ("dimer_axis", "dimer", _pair(2.5), 0.1, 0.0, 0.8),
        ("dimer_diagonal", "dimer", _pair(4.0, np.pi / 4.0), 0.05, 0.0, 0.4),
        ("dimer_rigid", "dimer", _pair(3.0), 0.1, 0.0, 1.0),
        ("trimer_equilateral", "equilateral_trimer", equilateral_trimer(3.0), 0.1, 0.0, 0.8),
        ("trimer_scalene", "scalene_trimer", scalene_trimer(2.7), 0.1, 0.0, 0.8),
        ("quartet_irregular", "irregular_quartet", irregular_quartet(2.8), 0.1, 0.0, 0.8),
    )


def _rms(vectors: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum(np.asarray(vectors) ** 2, axis=1))))


def _successive(current: np.ndarray, previous: np.ndarray) -> tuple[float, bool, float]:
    denominator = max(_rms(current), _rms(previous))
    absolute = _rms(current - previous)
    scale = max(_rms(current), _rms(previous))
    zero_tolerance = 128.0 * np.finfo(float).eps * scale
    if denominator <= zero_tolerance:
        return absolute, False, absolute
    return absolute / denominator, True, absolute


def _minimum_distance(positions: np.ndarray) -> float:
    return min(
        float(np.linalg.norm(positions[i] - positions[j]))
        for i in range(len(positions)) for j in range(i)
    )


def _format(value: float) -> str:
    return format(float(value), ".17g")


def _write(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _campaign():
    rows = []
    series = {}
    for case_id, geometry, positions, ka, f0, f1 in _cases():
        orders = list(BASE_ORDERS)
        results = [
            solve_model_e_nodal(positions, ka, 1.0, 1.0, f0, f1, order)
            for order in orders
        ]
        channels = (
            "total_forces_xyz", "interaction_forces_xyz",
            "external_scattered_forces_xyz", "scattered_scattered_forces_xyz",
        )

        def all_confirmed() -> bool:
            if len(results) < 3:
                return False
            for channel in channels:
                last, last_applicable, _ = _successive(
                    getattr(results[-1], channel), getattr(results[-2], channel)
                )
                prior, prior_applicable, _ = _successive(
                    getattr(results[-2], channel), getattr(results[-3], channel)
                )
                if not (last_applicable and prior_applicable and last <= TOLERANCE and prior <= TOLERANCE):
                    return False
            return True

        for order in (8, 9):
            if all_confirmed():
                break
            orders.append(order)
            results.append(
                solve_model_e_nodal(positions, ka, 1.0, 1.0, f0, f1, order)
            )

        confirmed = {}
        for channel in channels:
            confirmed[channel] = 0
            for index in range(2, len(results)):
                first, first_applicable, _ = _successive(
                    getattr(results[index - 1], channel), getattr(results[index - 2], channel)
                )
                second, second_applicable, _ = _successive(
                    getattr(results[index], channel), getattr(results[index - 1], channel)
                )
                if first_applicable and second_applicable and first <= TOLERANCE and second <= TOLERANCE:
                    confirmed[channel] = orders[index]
                    break

        coordinates = ";".join(
            ":".join(_format(value) for value in position) for position in positions
        )
        for index, (order, result) in enumerate(zip(orders, results)):
            previous = results[index - 1] if index else None
            metrics = {}
            for short, channel in (
                ("total", "total_forces_xyz"),
                ("interaction", "interaction_forces_xyz"),
                ("external_scattered", "external_scattered_forces_xyz"),
                ("scattered_scattered", "scattered_scattered_forces_xyz"),
            ):
                if previous is None:
                    change, applicable, absolute = 0.0, False, 0.0
                else:
                    change, applicable, absolute = _successive(
                        getattr(result, channel), getattr(previous, channel)
                    )
                metrics[f"{short}_rms"] = _rms(getattr(result, channel))
                metrics[f"{short}_successive_change"] = change
                metrics[f"{short}_absolute_change"] = absolute
                metrics[f"{short}_change_applicable"] = str(applicable).lower()
                metrics[f"{short}_minimum_confirmed_lmax"] = confirmed[channel]
            row = {
                "case_id": case_id,
                "geometry": geometry,
                "particle_count": len(positions),
                "coordinates_xyz": coordinates,
                "ka": _format(ka),
                "f0": _format(f0),
                "f1": _format(f1),
                "minimum_distance_over_radius": _format(_minimum_distance(positions)),
                "lmax": order,
                "full_modes_per_particle": (order + 1) ** 2,
                "active_modes_per_particle": len(result.solution.active_modes),
                "system_residual": _format(result.solution.residual_relative),
                "condition_number": _format(result.solution.condition_number),
                "decomposition_residual": _format(result.decomposition_residual),
                "used_planar_symmetry": str(result.solution.used_planar_symmetry).lower(),
            }
            row.update({key: _format(value) if isinstance(value, float) else value for key, value in metrics.items()})
            rows.append(row)
        series[case_id] = (orders, results, confirmed)
    return rows, series


def _oracle_rows(series):
    rows = []
    result = series["dimer_axis"][1][2]
    channels = (
        ("effective_total", result.solution.effective_incident_coefficients),
        ("external", result.solution.external_incident_coefficients),
        ("other_particles", result.scattered_incident_coefficients),
    )
    for channel_name, coefficients in channels:
        for particle in range(2):
            analytic = complete_radiation_force_from_bsc(
                coefficients[particle], result.solution.scattering_coefficients, 0.1, 1.0
            )
            analytic_norm = float(np.linalg.norm(analytic))
            resolved = analytic_norm > 128.0 * np.finfo(float).eps
            for control_radius in (1.01, 1.04):
                for theta_order, phi_count in ((24, 48), (32, 64)):
                    integrated = stress_tensor_force(
                        coefficients[particle], result.solution.scattering_coefficients,
                        0.1, 1.0, control_radius, theta_order, phi_count,
                    )
                    for component, label in enumerate(("x", "y", "z")):
                        absolute_error = abs(integrated[component] - analytic[component])
                        relative_error = absolute_error / analytic_norm if resolved else absolute_error
                        rows.append({
                            "case_id": "dimer_axis", "field_channel": channel_name,
                            "particle": particle, "component": label, "lmax": 4,
                            "control_radius_over_a": _format(control_radius),
                            "theta_order": theta_order, "phi_count": phi_count,
                            "analytic_force": _format(analytic[component]),
                            "integrated_force": _format(integrated[component]),
                            "absolute_error": _format(absolute_error),
                            "relative_error_or_absolute_if_unresolved": _format(relative_error),
                            "force_resolved": str(resolved).lower(),
                        })
    return rows


def _decomposition_rows(series):
    rows = []
    for case_id, (orders, results, _) in series.items():
        result = results[-1]
        interaction_rms = _rms(result.interaction_forces_xyz)
        ss_rms = _rms(result.scattered_scattered_forces_xyz)
        applicable = interaction_rms > 128.0 * np.finfo(float).eps * max(interaction_rms, ss_rms)
        ratio = ss_rms / interaction_rms if applicable else 0.0
        for particle in range(len(result.total_forces_xyz)):
            row = {
                "case_id": case_id, "particle": particle, "lmax": orders[-1],
                "decomposition_residual": _format(result.decomposition_residual),
                "scattered_scattered_over_interaction_rms": _format(ratio),
                "ratio_applicable": str(applicable).lower(),
            }
            for prefix, values in (
                ("total", result.total_forces_xyz),
                ("external", result.external_forces_xyz),
                ("interaction", result.interaction_forces_xyz),
                ("external_scattered", result.external_scattered_forces_xyz),
                ("scattered_scattered", result.scattered_scattered_forces_xyz),
            ):
                for component, label in enumerate(("x", "y", "z")):
                    row[f"{prefix}_{label}"] = _format(values[particle, component])
            rows.append(row)
    return rows


def _plot(series, oracle_rows):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)
    colors = plt.cm.tab10.colors
    for color, (case_id, (orders, results, confirmed)) in zip(colors, series.items()):
        changes = [
            _successive(results[index].interaction_forces_xyz, results[index - 1].interaction_forces_xyz)[0]
            for index in range(1, len(results))
        ]
        axes[0].semilogy(orders[1:], changes, marker="o", linewidth=1.2, color=color, label=case_id)
        if confirmed["interaction_forces_xyz"] == 0:
            axes[0].plot(orders[-1], changes[-1], marker="x", color=color, markersize=8)
    axes[0].axhline(TOLERANCE, color="black", linewidth=1.0, linestyle="--", label=r"$10^{-5}$")
    axes[0].set(xlabel=r"$L_{\max}$", ylabel="successive interaction-force change", title="Model-E convergence")

    resolved_rows = [row for row in oracle_rows if row["force_resolved"] == "true"]
    groups = {}
    for row in resolved_rows:
        key = (row["field_channel"], row["theta_order"], row["phi_count"])
        groups.setdefault(key, []).append(float(row["relative_error_or_absolute_if_unresolved"]))
    labels, values = [], []
    for (field, theta_order, phi_count), errors in groups.items():
        labels.append(f"{field}\n{theta_order}x{phi_count}")
        values.append(max(errors))
    axes[1].bar(np.arange(len(values)), values, color=colors[:len(values)])
    axes[1].axhline(1e-5, color="black", linewidth=1.0, linestyle="--")
    axes[1].set_yscale("log")
    axes[1].set_xticks(np.arange(len(values)), labels, rotation=35, ha="right", fontsize=7)
    axes[1].set(ylabel="maximum relative component error", title="Independent stress oracle")

    names, ext_sc, ss = [], [], []
    for case_id, (_, results, _) in series.items():
        names.append(case_id)
        ext_sc.append(_rms(results[-1].external_scattered_forces_xyz))
        ss.append(_rms(results[-1].scattered_scattered_forces_xyz))
    x = np.arange(len(names))
    axes[2].bar(x - 0.18, ext_sc, 0.36, label="external-scattered")
    axes[2].bar(x + 0.18, ss, 0.36, label="scattered-scattered")
    axes[2].set_yscale("log")
    axes[2].set_xticks(x, names, rotation=35, ha="right", fontsize=7)
    axes[2].set(ylabel=r"RMS force / $(a^2 E_0)$", title="Interaction-force channels")
    axes[2].legend(fontsize=8)
    for axis in axes:
        axis.grid(True, alpha=0.25)
    axes[0].legend(fontsize=7)
    fig.savefig(
        FIGURES / "t11_model_e_validation.png", dpi=220,
        metadata={"Software": "acoustic_ms T11"},
    )
    plt.close(fig)


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    convergence_rows, series = _campaign()
    oracle_rows = _oracle_rows(series)
    decomposition_rows = _decomposition_rows(series)
    _write(DATA / "t11_model_e_convergence.csv", list(convergence_rows[0]), convergence_rows)
    _write(DATA / "t11_force_oracle.csv", list(oracle_rows[0]), oracle_rows)
    _write(DATA / "t11_force_decomposition.csv", list(decomposition_rows[0]), decomposition_rows)
    _plot(series, oracle_rows)
    print(
        f"T11: {len(convergence_rows)} convergence rows, "
        f"{len(oracle_rows)} oracle rows, {len(decomposition_rows)} decomposition rows"
    )


if __name__ == "__main__":
    main()
