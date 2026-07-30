#!/usr/bin/env python3
"""Audit the analytical rho_1 operator and illustrate its Neumann expansion."""

import csv
import os
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parents[1] / "tmp" / "matplotlib"),
)

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

from acoustic_ms.cluster_families import enumerate_transferability_configurations
from acoustic_ms.rho_foundation import (
    dipolar_balanced_coupling_matrix,
    dipolar_coupling_diagnostics,
    neumann_partial_solutions,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"
K = 0.1
RADIUS = 1.0
NEUMANN_CASES = (
    "n4_irregular_f0.4_d4.0",
    "n6_irregular_f0.8_d2.5",
    "n10_compact_f1.0_d2.1",
)


def _symbolic_checks() -> tuple[bool, bool]:
    x = sp.symbols("x", positive=True, real=True)
    h0 = -sp.I * sp.exp(sp.I * x) / x
    h1 = -sp.exp(sp.I * x) * (x + sp.I) / x**2
    h2 = 3 * h1 / x - h0
    closed = -3 * sp.exp(sp.I * x) * (x + sp.I) / x**3
    hankel_identity = sp.simplify(h0 + h2 - closed) == 0
    series = sp.series(sp.exp(sp.I * x) * (1 - sp.I * x), x, 0, 5)
    expected = 1 + x**2 / 2 + sp.I * x**3 / 3 - x**4 / 8
    near_field_series = sp.simplify(series.removeO() - expected) == 0
    return bool(hankel_identity), bool(near_field_series)


def _read_frozen_rho() -> dict[str, float]:
    with (DATA / "t08_cases.csv").open(encoding="utf-8", newline="") as stream:
        return {
            row["case_id"]: float(row["rho_l1"])
            for row in csv.DictReader(stream)
        }


def _format_rows(rows: list[dict]) -> list[dict]:
    return [
        {
            key: format(value, ".17g") if isinstance(value, float) else value
            for key, value in row.items()
        }
        for row in rows
    ]


def _write(path: Path, rows: list[dict]) -> None:
    rows = _format_rows(rows)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def _operator_rows() -> list[dict]:
    frozen = _read_frozen_rho()
    rows = []
    for configuration in enumerate_transferability_configurations():
        diagnostics = dipolar_coupling_diagnostics(
            configuration.positions_xyz, K, RADIUS, configuration.f1
        )
        rho = diagnostics.spectral_radius
        near_field_relative_difference = (
            abs(diagnostics.near_field_spectral_radius - rho) / rho
        )
        rows.append(
            {
                "case_id": configuration.case_id,
                "split": configuration.split,
                "particle_count": configuration.particle_count,
                "family": configuration.family,
                "f1": configuration.f1,
                "distance_ratio": configuration.distance_ratio,
                "rho_frozen": frozen[configuration.case_id],
                "rho_analytic": rho,
                "rho_absolute_difference": abs(
                    frozen[configuration.case_id] - rho
                ),
                "rho_near_field": diagnostics.near_field_spectral_radius,
                "near_field_relative_difference": near_field_relative_difference,
                "spectral_norm": diagnostics.spectral_norm,
                "infinity_norm": diagnostics.infinity_norm,
                "spectral_norm_to_rho": diagnostics.spectral_norm / rho,
                "infinity_norm_to_rho": diagnostics.infinity_norm / rho,
                "normalized_commutator": diagnostics.normalized_commutator,
            }
        )
    return rows


def _neumann_rows() -> list[dict]:
    configurations = {
        item.case_id: item
        for item in enumerate_transferability_configurations()
    }
    rows = []
    for case_id in NEUMANN_CASES:
        configuration = configurations[case_id]
        coupling = dipolar_balanced_coupling_matrix(
            configuration.positions_xyz, K, RADIUS, configuration.f1
        )
        rho = float(np.max(np.abs(np.linalg.eigvals(coupling))))
        source = np.ones(configuration.particle_count, dtype=complex)
        exact = np.linalg.solve(
            np.eye(configuration.particle_count) - coupling, source
        )
        partials = neumann_partial_solutions(coupling, source, 12)
        errors = np.linalg.norm(partials - exact, axis=1) / np.linalg.norm(exact)
        for order, error in enumerate(errors):
            rows.append(
                {
                    "case_id": case_id,
                    "particle_count": configuration.particle_count,
                    "family": configuration.family,
                    "f1": configuration.f1,
                    "distance_ratio": configuration.distance_ratio,
                    "rho_analytic": rho,
                    "partial_order": order,
                    "relative_solution_error": float(error),
                    "rho_reference": float(errors[0] * rho**order),
                }
            )
    return rows


def _summary_rows(
    operator_rows: list[dict], symbolic_checks: tuple[bool, bool]
) -> list[dict]:
    contact = [
        row for row in operator_rows if row["distance_ratio"] == 2.1
    ]
    metrics = (
        ("symbolic_hankel_identity", float(symbolic_checks[0])),
        ("symbolic_near_field_series", float(symbolic_checks[1])),
        (
            "maximum_frozen_rho_absolute_difference",
            max(row["rho_absolute_difference"] for row in operator_rows),
        ),
        ("maximum_rho", max(row["rho_analytic"] for row in operator_rows)),
        (
            "maximum_spectral_norm_to_rho",
            max(row["spectral_norm_to_rho"] for row in operator_rows),
        ),
        (
            "maximum_infinity_norm_to_rho",
            max(row["infinity_norm_to_rho"] for row in operator_rows),
        ),
        (
            "maximum_normalized_commutator",
            max(row["normalized_commutator"] for row in operator_rows),
        ),
        (
            "maximum_near_field_relative_difference",
            max(
                row["near_field_relative_difference"]
                for row in operator_rows
            ),
        ),
        (
            "maximum_contact_near_field_relative_difference",
            max(
                row["near_field_relative_difference"]
                for row in contact
            ),
        ),
    )
    return [{"metric": name, "value": value} for name, value in metrics]


def _plot(operator_rows: list[dict], neumann_rows: list[dict]) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(10.8, 4.4), constrained_layout=True)

    exact = np.asarray([row["rho_analytic"] for row in operator_rows])
    near = np.asarray([row["rho_near_field"] for row in operator_rows])
    distance = np.asarray([row["distance_ratio"] for row in operator_rows])
    scatter = axes[0].scatter(
        exact, near, c=distance, cmap="viridis", s=20, alpha=0.78,
        linewidths=0.0,
    )
    limits = (0.8 * min(exact.min(), near.min()), 1.2 * max(exact.max(), near.max()))
    axes[0].plot(limits, limits, color="black", linewidth=1.0, linestyle="--")
    axes[0].set(
        xscale="log", yscale="log", xlim=limits, ylim=limits,
        xlabel=r"exact dipolar $\rho_1$",
        ylabel=r"near-field $\rho_1^{\mathrm{nf}}$",
        title="Inverse-cube limit",
    )
    colorbar = figure.colorbar(scatter, ax=axes[0], pad=0.02)
    colorbar.set_label(r"$d_{\min}/a$")

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for color, case_id in zip(colors, NEUMANN_CASES):
        selected = [row for row in neumann_rows if row["case_id"] == case_id]
        order = np.asarray([row["partial_order"] for row in selected])
        error = np.asarray([row["relative_solution_error"] for row in selected])
        reference = np.asarray([row["rho_reference"] for row in selected])
        rho = selected[0]["rho_analytic"]
        label = (
            rf"$N={selected[0]['particle_count']}$, "
            rf"$\rho_1={rho:.3g}$"
        )
        axes[1].plot(
            order, np.maximum(error, np.finfo(float).eps),
            marker="o", markersize=3.4, color=color, label=label,
        )
        axes[1].plot(
            order, np.maximum(reference, np.finfo(float).eps),
            linestyle="--", linewidth=1.0, color=color, alpha=0.7,
        )
    axes[1].set(
        yscale="log", xlabel="highest rescattering order $P$",
        ylabel=r"$\|q-q^{(P)}\|_2/\|q\|_2$",
        title="Neumann-series convergence",
    )
    axes[1].legend(fontsize=8)
    figure.savefig(
        FIGURES / "t09_rho_foundation.png",
        dpi=220,
        metadata={"Software": "acoustic_ms T09"},
    )
    plt.close(figure)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    symbolic_checks = _symbolic_checks()
    if not all(symbolic_checks):
        raise RuntimeError("T09 symbolic verification failed")
    operator_rows = _operator_rows()
    neumann_rows = _neumann_rows()
    summary_rows = _summary_rows(operator_rows, symbolic_checks)
    _write(DATA / "t09_operator_audit.csv", operator_rows)
    _write(DATA / "t09_neumann_convergence.csv", neumann_rows)
    _write(DATA / "t09_analytic_summary.csv", summary_rows)
    _plot(operator_rows, neumann_rows)
    print(
        "T09 analytical audit passed for "
        f"{len(operator_rows)} frozen configurations"
    )


if __name__ == "__main__":
    main()
