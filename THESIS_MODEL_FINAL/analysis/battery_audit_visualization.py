"""Plot-data adapters for the battery audit application.

All curves come from accepted loaders, saved measurements, or the existing
forward-model functions.  This module performs no fitting and creates no
synthetic curves.
"""

from __future__ import annotations

import csv
import tempfile
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Sequence

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "thesis_data"


@dataclass(frozen=True)
class PlotSeries:
    label: str
    values: np.ndarray
    style: str = "line"


@dataclass(frozen=True)
class PlotPanel:
    title: str
    x: np.ndarray
    x_label: str
    y_label: str
    series: tuple[PlotSeries, ...]


@dataclass(frozen=True)
class ExperimentPlot:
    experiment: str
    variable: str
    role: str
    panels: tuple[PlotPanel, ...]
    note: str = ""


PANASONIC_CATALOG = {
    "C/20": {
        "role": "Identificación de la representación de voltaje de equilibrio",
        "variables": ("Voltaje",),
    },
    "HPPC": {
        "role": "Calibración dinámica y validación reservada a bajo SOC",
        "variables": ("Voltaje",),
    },
    "Descarga a 1C": {
        "role": "Validación eléctrica y calibración térmica",
        "variables": ("Voltaje", "Temperatura"),
    },
    "Carga CC–CV": {
        "role": "Validación eléctrica y térmica en carga",
        "variables": ("Voltaje", "Temperatura"),
    },
    "UDDS": {
        "role": "Validación independiente con perfil dinámico",
        "variables": ("Voltaje", "Temperatura"),
    },
}


def panasonic_catalog_rows() -> tuple[dict[str, object], ...]:
    return tuple(
        {
            "experimento": name,
            "función": item["role"],
            "variables": ", ".join(item["variables"]),
        }
        for name, item in PANASONIC_CATALOG.items()
    )


def _panasonic_c20() -> ExperimentPlot:
    from models.ocp_panasonic_qi_frozen import V_eq
    from scripts.verify_panasonic_c20_qi_frozen import _load_c20

    z, measured = _load_c20()
    predicted = np.asarray(V_eq(z), dtype=float)
    return ExperimentPlot(
        experiment="C/20",
        variable="Voltaje",
        role=str(PANASONIC_CATALOG["C/20"]["role"]),
        panels=(
            PlotPanel(
                title="Voltaje durante la descarga de baja tasa",
                x=np.asarray(z),
                x_label="Capacidad descargada normalizada, z",
                y_label="Voltaje (V)",
                series=(
                    PlotSeries("Datos experimentales", np.asarray(measured), "points"),
                    PlotSeries("Modelo", predicted, "dashed"),
                ),
            ),
        ),
    )


def _panasonic_hppc() -> ExperimentPlot:
    from scripts.audit_panasonic_hppc_identifiability import load_hppc_data, load_pulse_records
    from scripts.calibrate_panasonic_hppc_primary import evaluate, prepare_groups

    parameters = {
        "R0": 0.0259594900586,
        "alpha_p": 3204.22765511,
        "d_n": 0.000183029692927,
    }
    data = load_hppc_data()
    records = load_pulse_records()
    high_groups = prepare_groups(data, records, minimum_pulse_start_soc=1.0 / 3.0)
    low_groups = prepare_groups(data, records, maximum_pulse_start_soc=1.0 / 3.0)
    evaluations = (evaluate(high_groups, parameters), evaluate(low_groups, parameters))

    def select(groups: Sequence[object], evaluation: object, target_soc: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
        candidates: list[tuple[float, int, object]] = []
        for group_index, group in enumerate(groups):
            for pulse in group.pulses:
                if pulse.record.nominal_rate == "2C":
                    candidates.append((abs(pulse.record.soc - target_soc), group_index, pulse))
        _, group_index, pulse = min(candidates, key=lambda item: item[0])
        group = groups[group_index]
        model = evaluation.group_voltage_v[group_index]
        start = int(pulse.phase_local_indices["immediate"][0])
        stop = int(pulse.phase_local_indices["post_relaxation"][-1])
        indices = np.arange(pulse.pre_local_index, stop + 1)
        return (
            group.time_s[indices] - group.time_s[start],
            group.measured_voltage_v[indices],
            model[indices],
            float(pulse.record.soc),
        )

    panels: list[PlotPanel] = []
    for groups, evaluation, target, label in (
        (high_groups, evaluations[0], 0.55, "Pulso utilizado en calibración"),
        (low_groups, evaluations[1], 0.20, "Pulso reservado para validación"),
    ):
        time, measured, predicted, soc = select(groups, evaluation, target)
        panels.append(
            PlotPanel(
                title=f"{label} · SOC = {soc:.2f}",
                x=time,
                x_label="Tiempo desde el inicio del pulso (s)",
                y_label="Voltaje (V)",
                series=(
                    PlotSeries("Datos experimentales", measured, "line"),
                    PlotSeries("Modelo", predicted, "dashed"),
                ),
            )
        )
    return ExperimentPlot(
        experiment="HPPC",
        variable="Voltaje",
        role=str(PANASONIC_CATALOG["HPPC"]["role"]),
        panels=tuple(panels),
    )


@lru_cache(maxsize=1)
def _electrical_reference_arrays() -> dict[str, dict[str, np.ndarray]]:
    from scripts.build_final_thesis_package import electrical_forward_arrays

    return electrical_forward_arrays()


@lru_cache(maxsize=1)
def _thermal_reference_arrays() -> dict[str, dict[str, np.ndarray]]:
    from scripts.build_final_thesis_package import thermal_forward_arrays

    return thermal_forward_arrays()


def load_panasonic_plot(experiment: str, variable: str) -> ExperimentPlot:
    """Load a real accepted Panasonic experiment/model comparison."""

    if experiment not in PANASONIC_CATALOG:
        raise ValueError(f"Experimento Panasonic desconocido: {experiment}")
    if variable not in PANASONIC_CATALOG[experiment]["variables"]:
        raise ValueError(f"{variable} no está disponible para {experiment}.")
    if experiment == "C/20":
        return _panasonic_c20()
    if experiment == "HPPC":
        return _panasonic_hppc()

    key = {"Descarga a 1C": "one_c", "Carga CC–CV": "cccv", "UDDS": "udds"}[experiment]
    role = str(PANASONIC_CATALOG[experiment]["role"])
    if variable == "Voltaje":
        values = _electrical_reference_arrays()[key]
        scale = 3600.0 if experiment == "UDDS" else 60.0
        unit = "h" if experiment == "UDDS" else "min"
        x = np.asarray(values["time"]) / scale
        panel = PlotPanel(
            title=f"{experiment} · comparación eléctrica",
            x=x,
            x_label=f"Tiempo transcurrido ({unit})",
            y_label="Voltaje (V)",
            series=(
                PlotSeries("Datos experimentales", np.asarray(values["measured"]), "line"),
                PlotSeries("Modelo", np.asarray(values["predicted"]), "dashed"),
            ),
        )
    else:
        values = _thermal_reference_arrays()[key]
        scale = 3600.0 if experiment == "UDDS" else 60.0
        unit = "h" if experiment == "UDDS" else "min"
        model_key = "predicted" if experiment == "Descarga a 1C" else "b"
        panel = PlotPanel(
            title=f"{experiment} · comparación térmica",
            x=np.asarray(values["time"]) / scale,
            x_label=f"Tiempo transcurrido ({unit})",
            y_label="Temperatura (°C)",
            series=(
                PlotSeries("Datos experimentales", np.asarray(values["measured"]), "points"),
                PlotSeries("Modelo", np.asarray(values[model_key]), "dashed"),
            ),
        )
    return ExperimentPlot(experiment, variable, role, (panel,))


TRC_OPTIONS = (
    ("Rionegro · Aire 22–24 °C", "Descarga", "Referencia", ("Voltaje", "Corriente", "Temperatura")),
    ("Rionegro · Aire 22–24 °C", "Descarga", "Resistencia externa", ("Voltaje", "Corriente")),
    ("Rionegro · Aire 22–24 °C", "Carga CC–CV", "Referencia", ("Voltaje", "Corriente", "Temperatura")),
    ("Rionegro · Aire 22–24 °C", "Carga CC–CV", "Resistencia externa", ("Voltaje", "Corriente")),
    ("Cali · Aire 29–31 °C", "Descarga", "Resistencia externa", ("Voltaje", "Corriente", "Temperatura")),
    ("Cali · Aire 29–31 °C", "Carga CC–CV", "Resistencia externa", ("Voltaje", "Corriente", "Temperatura")),
    ("Valledupar · Aire 33–37 °C", "Descarga", "Resistencia externa", ("Voltaje", "Temperatura")),
    ("Valledupar · Aire 33–37 °C", "Carga CC–CV", "Resistencia externa", ("Voltaje", "Corriente", "Temperatura")),
)


def trc_catalog_rows() -> tuple[dict[str, object], ...]:
    return tuple(
        {
            "condición": condition,
            "ensayo": test,
            "montaje": setup,
            "variables": ", ".join(variables),
        }
        for condition, test, setup, variables in TRC_OPTIONS
    )


def _trc_source_ids(condition: str, test: str, setup: str) -> tuple[str, str | None]:
    electrical = {
        ("Rionegro · Aire 22–24 °C", "Descarga", "Referencia"): "cycle_cell1_ambient",
        ("Rionegro · Aire 22–24 °C", "Descarga", "Resistencia externa"): "cycle_cell2_undeclared",
        ("Rionegro · Aire 22–24 °C", "Carga CC–CV", "Referencia"): "cycle_cell1_ambient",
        ("Rionegro · Aire 22–24 °C", "Carga CC–CV", "Resistencia externa"): "cycle_cell2_undeclared",
        ("Cali · Aire 29–31 °C", "Descarga", "Resistencia externa"): "cycle_cell2_30c",
        ("Cali · Aire 29–31 °C", "Carga CC–CV", "Resistencia externa"): "cycle_cell2_30c",
        ("Valledupar · Aire 33–37 °C", "Descarga", "Resistencia externa"): "cycle_cell1_30c",
        ("Valledupar · Aire 33–37 °C", "Carga CC–CV", "Resistencia externa"): "cycle_cell1_35c",
    }
    thermal = {
        "Rionegro · Aire 22–24 °C": "cycle_cell1_ambient",
        "Cali · Aire 29–31 °C": "cycle_cell1_30c",
        "Valledupar · Aire 33–37 °C": "cycle_cell1_35c",
    }
    return electrical[(condition, test, setup)], thermal.get(condition)


def load_trc_plot(condition: str, test: str, setup: str, variable: str) -> ExperimentPlot:
    """Load one accepted TRC presentation curve without cross-file merging."""

    option = next(
        (item for item in TRC_OPTIONS if item[:3] == (condition, test, setup)),
        None,
    )
    if option is None or variable not in option[3]:
        raise ValueError("La combinación seleccionada no tiene correspondencia respaldada.")
    from analysis.own_case_study_reference import segment_complete_cycle
    from scripts.audit_own_trc18650_data import load_all
    from scripts.build_uniandes_charge_temperature_figures import _charge_record

    experiments = {item.metadata.dataset_id: item for item in load_all()}
    electrical_id, thermal_id = _trc_source_ids(condition, test, setup)
    note = ""
    if variable == "Temperatura":
        if thermal_id is None:
            raise ValueError("No existe correspondencia térmica documentada.")
        experiment = experiments[thermal_id]
        if thermal_id != electrical_id:
            note = "La temperatura procede de un registro documentado separado y no se combina muestra a muestra con la curva eléctrica."
        if test == "Descarga":
            segment = next(item for item in segment_complete_cycle(experiment) if item.name == "main_1c_discharge")
            rows = np.arange(segment.start, segment.end + 1)
        else:
            rows = _charge_record(experiment).rows
        raw_capacity = np.asarray(experiment.numeric["Capacity(Ah)"][rows], dtype=float)
        x = np.abs(raw_capacity - raw_capacity[0])
        series = tuple(
            PlotSeries(label, np.asarray(experiment.numeric[name][rows]), "line")
            for name, label in (
                ("Temp1", "Superficie central"),
                ("Temp2", "Superficie cerca del extremo"),
                ("Temp7", "Aire de la cámara"),
            )
        )
        panel = PlotPanel(
            title=f"{condition} · {test}",
            x=x,
            x_label="Capacidad descargada (Ah)" if test == "Descarga" else "Capacidad cargada (Ah)",
            y_label="Temperatura (°C)",
            series=series,
        )
    else:
        experiment = experiments[electrical_id]
        if test == "Descarga":
            segment = next(item for item in segment_complete_cycle(experiment) if item.name == "main_1c_discharge")
            rows = np.arange(segment.start, segment.end + 1)
            raw_capacity = np.asarray(experiment.numeric["Capacity(Ah)"][rows], dtype=float)
            x = np.abs(raw_capacity - raw_capacity[0])
        else:
            charge = _charge_record(experiment)
            rows = charge.rows
            x = charge.charged_capacity_ah
        if variable == "Voltaje":
            values = np.asarray(experiment.numeric["Voltage(V)"][rows])
            y_label = "Voltaje terminal (V)"
        else:
            values = np.asarray(experiment.numeric["Current(A)"][rows])
            y_label = "Corriente (A)"
        panel = PlotPanel(
            title=f"{condition} · {test} · {setup}",
            x=np.asarray(x),
            x_label="Capacidad descargada (Ah)" if test == "Descarga" else "Capacidad cargada (Ah)",
            y_label=y_label,
            series=(PlotSeries("Datos experimentales", values, "line"),),
        )
    return ExperimentPlot(
        experiment=f"{condition} · {test}",
        variable=variable,
        role="Caracterización experimental",
        panels=(panel,),
        note=note,
    )


def plot_is_finite(plot: ExperimentPlot) -> bool:
    return all(
        panel.x.size > 1
        and np.all(np.isfinite(panel.x))
        and all(series.values.size == panel.x.size and np.all(np.isfinite(series.values)) for series in panel.series)
        for panel in plot.panels
    )


def _load_uploaded_experiment(spec: object) -> object:
    """Materialize uploaded bytes through the accepted Chroma reader."""

    from analysis.experimental_data_audit import ExperimentMetadata, load_chroma_experiment

    with tempfile.TemporaryDirectory(prefix="battery_plot_upload_") as directory:
        path = Path(directory) / Path(str(spec.filename)).name
        path.write_bytes(spec.content)
        metadata = ExperimentMetadata(
            dataset_id=str(spec.dataset_id),
            path=path,
            cell_id=str(spec.cell_id),
            nominal_temperature=str(spec.condition),
            encoded_date="no suministrada",
            temperature_roles=dict(spec.temperature_roles or {}),
        )
        return load_chroma_experiment(metadata, spec.fallback_columns)


def uploaded_plot_options(spec: object) -> tuple[tuple[str, ...], tuple[str, ...]]:
    """Return supported tests and user-facing variables for one uploaded file."""

    experiment = _load_uploaded_experiment(spec)
    current = np.asarray(experiment.current_a, dtype=float)
    tests: list[str] = []
    if np.any(current < -1e-6):
        tests.append("Descarga")
    if np.any(current > 1e-6):
        tests.append("Carga")
    if not tests:
        tests.append("Registro completo")
    variables = ["Voltaje", "Corriente"]
    if any(name.startswith("Temp") for name in experiment.numeric):
        variables.append("Temperatura")
    return tuple(tests), tuple(variables)


def load_uploaded_plot(spec: object, test: str, variable: str) -> ExperimentPlot:
    """Load one uploaded experiment without presenting reconstructed time."""

    experiment = _load_uploaded_experiment(spec)
    current = np.asarray(experiment.current_a, dtype=float)
    if test == "Descarga":
        indices = np.flatnonzero(current < -1e-6)
    elif test == "Carga":
        indices = np.flatnonzero(current > 1e-6)
    elif test == "Registro completo":
        indices = np.arange(current.size)
    else:
        raise ValueError("Tipo de ensayo no disponible en el registro.")
    if indices.size < 2:
        raise ValueError("El tramo seleccionado no contiene suficientes muestras.")

    capacity = np.asarray(experiment.numeric["Capacity(Ah)"][indices], dtype=float)
    x = np.abs(capacity - capacity[0])
    x_label = "Capacidad recorrida desde el inicio del tramo (Ah)"
    if variable == "Voltaje":
        series = (PlotSeries("Datos experimentales", np.asarray(experiment.voltage_v[indices]), "line"),)
        y_label = "Voltaje (V)"
    elif variable == "Corriente":
        series = (PlotSeries("Datos experimentales", current[indices], "line"),)
        y_label = "Corriente (A)"
    elif variable == "Temperatura":
        labels = {
            "cell surface middle": "Superficie central",
            "cell surface end": "Superficie cerca del extremo",
            "chamber/ambient": "Aire de la cámara",
            "temperature location unknown": "Temperatura medida",
        }
        available: list[PlotSeries] = []
        for name, role in dict(spec.temperature_roles or {}).items():
            if name in experiment.numeric:
                available.append(PlotSeries(labels.get(role, "Temperatura medida"), np.asarray(experiment.numeric[name][indices]), "line"))
        if not available:
            for name in experiment.numeric:
                if name.startswith("Temp"):
                    available.append(PlotSeries("Temperatura medida", np.asarray(experiment.numeric[name][indices]), "line"))
        if not available:
            raise ValueError("El registro no contiene temperatura disponible.")
        series = tuple(available)
        y_label = "Temperatura (°C)"
    else:
        raise ValueError("Variable no disponible.")
    panel = PlotPanel(
        title=f"{spec.condition} · {test}",
        x=x,
        x_label=x_label,
        y_label=y_label,
        series=series,
    )
    return ExperimentPlot(
        experiment=f"{spec.condition} · {test}",
        variable=variable,
        role="Datos experimentales cargados",
        panels=(panel,),
        note="El eje usa el contador de capacidad; no se presenta tiempo reconstruido como tiempo medido.",
    )


__all__ = [
    "ExperimentPlot",
    "PANASONIC_CATALOG",
    "PlotPanel",
    "PlotSeries",
    "TRC_OPTIONS",
    "load_panasonic_plot",
    "load_trc_plot",
    "load_uploaded_plot",
    "panasonic_catalog_rows",
    "plot_is_finite",
    "trc_catalog_rows",
    "uploaded_plot_options",
]
