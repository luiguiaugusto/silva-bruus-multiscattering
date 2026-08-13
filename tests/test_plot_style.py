from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

from acoustic_ms.plot_style import (
    ACCESSIBLE_COLORS,
    ONE_COLUMN_MM,
    REDUNDANT_MARKERS,
    TWO_COLUMN_MM,
    diagnostic_rc_context,
    diagnostic_style,
    figure_size,
    millimetres_to_inches,
    save_diagnostic_figure,
)


def test_physical_figure_dimensions_are_converted_exactly() -> None:
    one = figure_size("one", height_mm=50.8)
    two = figure_size("two", height_mm=25.4)
    assert one == pytest.approx((ONE_COLUMN_MM / 25.4, 2.0))
    assert two == pytest.approx((TWO_COLUMN_MM / 25.4, 1.0))
    with pytest.raises(ValueError):
        millimetres_to_inches(0.0)


def test_diagnostic_style_is_accessible_and_reversible() -> None:
    before = plt.rcParams["axes.labelsize"]
    style = diagnostic_style()
    assert style["axes.labelsize"] == 9.0
    assert style["xtick.labelsize"] == 8.0
    assert len(ACCESSIBLE_COLORS) == len(REDUNDANT_MARKERS)
    with diagnostic_rc_context():
        assert plt.rcParams["axes.labelsize"] == 9.0
        assert plt.rcParams["legend.fontsize"] == 8.0
    assert plt.rcParams["axes.labelsize"] == before


def test_noninteractive_smoke_exports_pdf_svg_png_deterministically(tmp_path: Path) -> None:
    with diagnostic_rc_context():
        figure, axis = plt.subplots(figsize=figure_size("one"))
        axis.plot([0.0, 1.0], [0.0, 1.0], label="diagnostic")
        axis.set(xlabel="x", ylabel="y")
        axis.legend()
        paths = save_diagnostic_figure(figure, tmp_path / "smoke")
        first = {path.suffix: path.read_bytes() for path in paths}
        repeated = save_diagnostic_figure(figure, tmp_path / "smoke")
        second = {path.suffix: path.read_bytes() for path in repeated}
        plt.close(figure)
    assert {path.suffix for path in paths} == {".pdf", ".svg", ".png"}
    assert all(path.stat().st_size > 0 for path in paths)
    assert first == second
