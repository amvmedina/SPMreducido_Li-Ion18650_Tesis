from __future__ import annotations

from pathlib import Path

import pytest

from analysis.battery_audit_visualization import (
    load_panasonic_plot,
    load_trc_plot,
    plot_is_finite,
)


def test_panasonic_c20_plot_uses_finite_measured_and_model_curves() -> None:
    plot = load_panasonic_plot("C/20", "Voltaje")
    assert plot_is_finite(plot)
    assert [series.label for series in plot.panels[0].series] == ["Datos experimentales", "Modelo"]


def test_trc_plot_uses_capacity_and_accepted_condition_mapping() -> None:
    source_directory = Path(__file__).resolve().parents[1] / "data" / "trc18650"
    if not source_directory.is_dir():
        pytest.skip("Configured TRC files are not installed.")
    plot = load_trc_plot("Rionegro · Aire 22–24 °C", "Descarga", "Referencia", "Voltaje")
    assert plot_is_finite(plot)
    assert plot.panels[0].x_label == "Capacidad descargada (Ah)"
    assert "cycle_cell" not in plot.panels[0].title
