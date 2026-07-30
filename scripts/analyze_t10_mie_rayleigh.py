#!/usr/bin/env python3
"""Compare exact isolated-sphere Mie coefficients with Rayleigh limits."""

import csv
import os
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    "/tmp/acoustic_ms_t10_matplotlib",
)

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import spherical_jn, spherical_yn

from acoustic_ms.mie_scattering import (
    material_ratios_from_contrasts,
    mie_scattering_coefficients_from_contrasts,
)
from acoustic_ms.multipolar_scattering import (
    rayleigh_multipolar_scattering_coefficients,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "results" / "data"
FIGURES = ROOT / "results" / "figures"
F0 = 0.0
F1_VALUES = (0.1, 0.4, 0.8, 1.0)
KA_VALUES = np.geomspace(1e-3, 1e-1, 101)
LMAX = 5


def _hankel(ell: int, value: float, derivative: bool = False) -> complex:
    return spherical_jn(ell, value, derivative=derivative) + 1j * spherical_yn(
        ell, value, derivative=derivative
    )


def _boundary_residuals(
    ka: float,
    density_ratio: float,
    compressibility_ratio: float,
    ell: int,
    coefficient: complex,
) -> tuple[float, float]:
    internal_ka = ka * np.sqrt(density_ratio * compressibility_ratio)
    beta = np.sqrt(compressibility_ratio / density_ratio)
    exterior_pressure = spherical_jn(ell, ka) + coefficient * _hankel(
        ell, ka
    )
    interior_value = spherical_jn(ell, internal_ka)
    internal_amplitude = exterior_pressure / interior_value
    pressure_residual = abs(
        exterior_pressure - internal_amplitude * interior_value
    )
    velocity_residual = abs(
        spherical_jn(ell, ka, derivative=True)
        + coefficient * _hankel(ell, ka, derivative=True)
        - beta
        * internal_amplitude
        * spherical_jn(ell, internal_ka, derivative=True)
    )
    return float(pressure_residual), float(velocity_residual)


def _rigid_boundary_residual(
    ka: float, ell: int, coefficient: complex
) -> float:
    return float(
        abs(
            spherical_jn(ell, ka, derivative=True)
            + coefficient * _hankel(ell, ka, derivative=True)
        )
    )


def _phase_difference(left: complex, right: complex) -> float:
    return float(np.angle(np.exp(1j * (np.angle(left) - np.angle(right)))))


def _format_value(value: object) -> object:
    if isinstance(value, (float, np.floating)):
        return format(float(value), ".17g")
    return value


def _write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError("cannot write an empty T10 table")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(rows[0]), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(
            {key: _format_value(value) for key, value in row.items()}
            for row in rows
        )


def _validation_rows() -> list[dict]:
    rows = []
    for f1 in F1_VALUES:
        rigid = f1 == 1.0
        if rigid:
            density_ratio = np.nan
            compressibility_ratio = 1.0 - F0
            sound_speed_ratio = np.nan
        else:
            (
                density_ratio,
                compressibility_ratio,
                sound_speed_ratio,
            ) = material_ratios_from_contrasts(F0, f1)
        for ka in KA_VALUES:
            mie = mie_scattering_coefficients_from_contrasts(
                ka, F0, f1, LMAX
            )
            rayleigh = rayleigh_multipolar_scattering_coefficients(
                ka, F0, f1, LMAX
            )
            for ell in range(LMAX + 1):
                mie_magnitude = abs(mie[ell])
                rayleigh_magnitude = abs(rayleigh[ell])
                relative_applicable = mie_magnitude > 0.0
                phase_applicable = (
                    mie_magnitude > 0.0 and rayleigh_magnitude > 0.0
                )
                if rigid:
                    pressure_residual = np.nan
                    velocity_residual = _rigid_boundary_residual(
                        ka, ell, mie[ell]
                    )
                else:
                    pressure_residual, velocity_residual = _boundary_residuals(
                        ka,
                        density_ratio,
                        compressibility_ratio,
                        ell,
                        mie[ell],
                    )
                absolute_error = abs(mie[ell] - rayleigh[ell])
                rows.append(
                    {
                        "ka": ka,
                        "f0": F0,
                        "f1": f1,
                        "density_ratio": density_ratio,
                        "compressibility_ratio": compressibility_ratio,
                        "sound_speed_ratio": sound_speed_ratio,
                        "is_rigid_limit": rigid,
                        "material_ratios_applicable": not rigid,
                        "ell": ell,
                        "mie_real": mie[ell].real,
                        "mie_imag": mie[ell].imag,
                        "rayleigh_real": rayleigh[ell].real,
                        "rayleigh_imag": rayleigh[ell].imag,
                        "mie_magnitude": mie_magnitude,
                        "rayleigh_magnitude": rayleigh_magnitude,
                        "complex_relative_error": (
                            absolute_error / mie_magnitude
                            if relative_applicable
                            else np.nan
                        ),
                        "complex_relative_error_applicable": relative_applicable,
                        "magnitude_relative_error": (
                            abs(mie_magnitude - rayleigh_magnitude)
                            / mie_magnitude
                            if relative_applicable
                            else np.nan
                        ),
                        "magnitude_relative_error_applicable": relative_applicable,
                        "phase_difference_rad": (
                            _phase_difference(mie[ell], rayleigh[ell])
                            if phase_applicable
                            else np.nan
                        ),
                        "phase_difference_applicable": phase_applicable,
                        "absolute_error": absolute_error,
                        "pressure_boundary_residual": pressure_residual,
                        "pressure_boundary_residual_applicable": not rigid,
                        "velocity_boundary_residual": velocity_residual,
                        "unitarity_defect": abs(
                            mie[ell].real + mie_magnitude**2
                        ),
                        "relative_channel_physically_active_in_nodal_problem": ell > 0,
                    }
                )
    return rows


def _summary_rows(rows: list[dict]) -> list[dict]:
    summaries = []
    for f1 in F1_VALUES:
        for ell in range(LMAX + 1):
            selected = [
                row
                for row in rows
                if row["f1"] == f1 and row["ell"] == ell
            ]
            asymptotic = selected[:31]
            relative = np.asarray(
                [row["complex_relative_error"] for row in asymptotic]
            )
            slope = np.polyfit(
                np.log([row["ka"] for row in asymptotic]),
                np.log(relative),
                1,
            )[0]
            at_upper = selected[-1]
            pressure = [
                row["pressure_boundary_residual"]
                for row in selected
                if row["pressure_boundary_residual_applicable"]
            ]
            summaries.append(
                {
                    "f1": f1,
                    "ell": ell,
                    "is_rigid_limit": f1 == 1.0,
                    "point_count": len(selected),
                    "maximum_complex_relative_error": max(
                        row["complex_relative_error"] for row in selected
                    ),
                    "complex_relative_error_at_ka_0p1": at_upper[
                        "complex_relative_error"
                    ],
                    "magnitude_relative_error_at_ka_0p1": at_upper[
                        "magnitude_relative_error"
                    ],
                    "absolute_error_at_ka_0p1": at_upper["absolute_error"],
                    "asymptotic_relative_error_slope": float(slope),
                    "relative_channel_physically_active_in_nodal_problem": ell > 0,
                    "maximum_pressure_boundary_residual": (
                        max(pressure) if pressure else np.nan
                    ),
                    "maximum_velocity_boundary_residual": max(
                        row["velocity_boundary_residual"] for row in selected
                    ),
                    "maximum_unitarity_defect": max(
                        row["unitarity_defect"] for row in selected
                    ),
                }
            )
    return summaries


def _plot(rows: list[dict]) -> None:
    figure, axes = plt.subplots(
        1, 2, figsize=(10.8, 4.4), constrained_layout=True
    )
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for color, f1 in zip(colors, F1_VALUES):
        dipole = [
            row for row in rows if row["f1"] == f1 and row["ell"] == 1
        ]
        axes[0].plot(
            [row["ka"] for row in dipole],
            [100 * row["complex_relative_error"] for row in dipole],
            color=color,
            linewidth=1.6,
            label=rf"$f_1={f1:g}$",
        )
        upper = [
            row
            for row in rows
            if row["f1"] == f1
            and row["ell"] > 0
            and np.isclose(row["ka"], 0.1)
        ]
        axes[1].plot(
            [row["ell"] for row in upper],
            [100 * row["complex_relative_error"] for row in upper],
            marker="o",
            markersize=4.0,
            color=color,
            linewidth=1.4,
            label=rf"$f_1={f1:g}$",
        )
    axes[0].set(
        xscale="log",
        yscale="log",
        xlabel=r"size parameter $ka$",
        ylabel=r"coefficient error $100\,\varepsilon_{s_1}$ (%)",
        title=r"Dominant dipole ($\ell=1$)",
    )
    axes[1].set(
        yscale="log",
        xlabel=r"multipole order $\ell$",
        ylabel=r"coefficient error at $ka=0.1$ (%)",
        title="Rayleigh truncation by order",
        xticks=range(1, LMAX + 1),
    )
    for axis in axes:
        axis.grid(True, which="both", linewidth=0.45, alpha=0.35)
        axis.legend(fontsize=8)
    figure.savefig(
        FIGURES / "t10_mie_rayleigh_error.png",
        dpi=220,
        metadata={"Software": "acoustic_ms T10"},
    )
    plt.close(figure)


def _audit(rows: list[dict], summaries: list[dict]) -> None:
    if len(rows) != len(F1_VALUES) * len(KA_VALUES) * (LMAX + 1):
        raise RuntimeError("unexpected T10 validation row count")
    if len(summaries) != len(F1_VALUES) * (LMAX + 1):
        raise RuntimeError("unexpected T10 summary row count")
    if any(not np.isfinite(row["velocity_boundary_residual"]) for row in rows):
        raise RuntimeError("non-finite boundary residual")
    if any(not np.isfinite(row["unitarity_defect"]) for row in rows):
        raise RuntimeError("non-finite unitarity defect")
    for summary in summaries:
        if summary["ell"] > 0 and not 1.85 < summary[
            "asymptotic_relative_error_slope"
        ] < 2.15:
            raise RuntimeError("unexpected Rayleigh correction slope")


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    rows = _validation_rows()
    summaries = _summary_rows(rows)
    _audit(rows, summaries)
    _write_csv(DATA / "t10_mie_rayleigh_validation.csv", rows)
    _write_csv(DATA / "t10_mie_rayleigh_summary.csv", summaries)
    _plot(rows)
    maximum_boundary = max(
        row["velocity_boundary_residual"] for row in rows
    )
    maximum_unitarity = max(row["unitarity_defect"] for row in rows)
    print(
        f"T10: {len(rows)} validation rows, {len(summaries)} summary rows; "
        f"max boundary residual={maximum_boundary:.3e}; "
        f"max unitarity defect={maximum_unitarity:.3e}"
    )


if __name__ == "__main__":
    main()
