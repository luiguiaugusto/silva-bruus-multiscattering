"""Reproduce Figure 2 relative-error data for the corrected pair benchmark."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from acoustic_ms import (
    corrected_nodal_pair_force_magnitude,
    nodal_pair_force_magnitude,
)


KA = 0.1
RADIUS_M = 1.0
ENERGY_DENSITY_J_M3 = 1.0
KD_VALUES = np.linspace(0.2, 0.3, 501)
F1_VALUES = (0.1, 0.4, 0.8, 1.0)
ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "results" / "data" / "figure_2_relative_error.csv"
FIGURE_PATH = ROOT / "results" / "figures" / "figure_2_relative_error.png"


def compute_rows() -> list[dict[str, float]]:
    """Calculate Figure 2 data through the package public force APIs."""
    k = KA / RADIUS_M
    rows: list[dict[str, float]] = []
    for f1 in F1_VALUES:
        for kd in KD_VALUES:
            distance = kd / k
            corrected = corrected_nodal_pair_force_magnitude(k, RADIUS_M, distance, ENERGY_DENSITY_J_M3, f1)
            silva_bruus = nodal_pair_force_magnitude(k, RADIUS_M, distance, ENERGY_DENSITY_J_M3, f1)
            if corrected == 0.0:
                raise ZeroDivisionError("relative error is undefined for zero corrected force")
            rows.append({
                "ka": KA,
                "kd": float(kd),
                "separation_ratio": RADIUS_M / distance,
                "f1": f1,
                "corrected_force": corrected,
                "silva_bruus_force": silva_bruus,
                "relative_error_percent": 100.0 * abs(corrected - silva_bruus) / abs(corrected),
            })
    return rows


def write_csv(rows: list[dict[str, float]]) -> None:
    """Write reproducible numerical data, including the contact-limit point."""
    CSV_PATH.parent.mkdir(parents=True, exist_ok=True)
    with CSV_PATH.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot(rows: list[dict[str, float]]) -> None:
    """Render the four Figure 2 relative-error curves from calculated rows."""
    FIGURE_PATH.parent.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(7.0, 4.5), constrained_layout=True)
    for f1 in F1_VALUES:
        subset = [row for row in rows if row["f1"] == f1]
        axis.plot([row["kd"] for row in subset], [row["relative_error_percent"] for row in subset], label=fr"$f_1={f1}$")
    axis.set(xlabel=r"$kd$", ylabel="Relative error (%)", xlim=(0.20, 0.30))
    axis.grid(True, alpha=0.3)
    axis.legend(title="Contact-limit point $kd=2ka=0.2$ included")
    figure.savefig(FIGURE_PATH, dpi=200)
    plt.close(figure)


def main() -> None:
    rows = compute_rows()
    write_csv(rows)
    plot(rows)
    print(f"Wrote {len(rows)} rows to {CSV_PATH.relative_to(ROOT)}")
    print(f"Wrote figure to {FIGURE_PATH.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
