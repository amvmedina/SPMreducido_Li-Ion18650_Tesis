"""Build own-data CC-CV charge and temperature figures from original Chroma CSVs.

This module is presentation and data-audit only.  It does not run or calibrate
the battery model.  Every plotted point comes from one original source file;
no cross-file electrical/thermal alignment or extrapolation is performed.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.own_case_study_reference import (  # noqa: E402
    PhysicalSegment,
    segment_complete_cycle,
)
from scripts.audit_own_trc18650_data import load_all  # noqa: E402


DATA = ROOT / "thesis_data"
THESIS = ROOT / "thesis_figures" / "final"
PRESENTATION = ROOT / "presentation_assets" / "final"
FINAL_ES = ROOT / "thesis_figures" / "final_es"
LOCALIZED_THESIS = ROOT / "figuras_finales_es" / "tesis"
LOCALIZED_PRESENTATION = ROOT / "figuras_finales_es" / "presentacion"

INK = "#172033"
PURPLE = "#7656A5"
BLUE = "#2563A6"
ORANGE = "#E07A2D"
GREEN = "#2E8B57"
RED = "#C7473D"
TEAL = "#16858C"
GRAY = "#5F6775"

SENSOR_COLORS = {"Temp1": RED, "Temp2": TEAL, "Temp7": GRAY}
SENSOR_LABELS = {
    "Temp1": "Superficie central",
    "Temp2": "Superficie cerca del extremo",
    "Temp7": "Aire de la cámara",
}
VALLEDUPAR_AIR_MEDIAN_WINDOW_SAMPLES = 11

CHARGE_EXTERNAL = (
    {
        "key": "rionegro_external",
        "dataset_id": "cycle_cell2_undeclared",
        "label": "Rionegro · Aire 22–24 °C",
        "color": PURPLE,
        "assignment": "same-file charge voltage/current/Capacity(Ah); discharge-only reassignment not applied",
    },
    {
        "key": "cali_external",
        "dataset_id": "cycle_cell2_30c",
        "label": "Cali · Aire 29–31 °C",
        "color": BLUE,
        "assignment": "same-file charge voltage/current/Capacity(Ah); discharge-only reassignment not applied",
    },
    {
        "key": "valledupar_external",
        "dataset_id": "cycle_cell1_35c",
        "label": "Valledupar · Aire 33–37 °C",
        "color": ORANGE,
        "assignment": "same-file charge voltage/current/Capacity(Ah); discharge-only reassignment not applied",
    },
)

CHARGE_CONTACT = (
    {
        "key": "rionegro_reference",
        "dataset_id": "cycle_cell1_ambient",
        "label": "Montaje de referencia",
        "color": GREEN,
        "assignment": "same-file charge voltage/current/Capacity(Ah)",
    },
    {
        "key": "rionegro_external",
        "dataset_id": "cycle_cell2_undeclared",
        "label": "Montaje con resistencia externa",
        "color": PURPLE,
        "assignment": "same-file charge voltage/current/Capacity(Ah)",
    },
)

THERMAL_CONDITIONS = (
    {
        "key": "rionegro_reference",
        "dataset_id": "cycle_cell1_ambient",
        "label": "Rionegro · Aire 22–24 °C",
        "condition_color": PURPLE,
        "setup": "montaje de referencia",
    },
    {
        "key": "cali_cell1",
        "dataset_id": "cycle_cell1_30c",
        "label": "Cali · Aire 29–31 °C",
        "condition_color": BLUE,
        "setup": "montaje con resistencia externa",
    },
    {
        "key": "valledupar_cell1",
        "dataset_id": "cycle_cell1_35c",
        "label": "Valledupar · Aire 33–37 °C",
        "condition_color": ORANGE,
        "setup": "montaje con resistencia externa",
    },
)

OUTPUTS = (
    (
        "uniandes_cccv_external_conditions_final",
        "16_carga_cccv_resistencia_externa_uniandes",
    ),
    (
        "uniandes_cccv_contact_comparison_final",
        "17_comparacion_carga_cccv_montajes_uniandes",
    ),
    (
        "uniandes_temperature_discharge_final",
        "18_temperaturas_descarga_uniandes",
    ),
    (
        "uniandes_temperature_charge_final",
        "19_temperaturas_carga_uniandes",
    ),
)


@dataclass(frozen=True)
class ChargeRecord:
    rows: np.ndarray
    charged_capacity_ah: np.ndarray
    voltage_v: np.ndarray
    current_a: np.ndarray
    cc: PhysicalSegment
    cv: PhysicalSegment
    cv_local_index: int
    counter_gap_cc_to_cv_ah: float


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11.5,
            "axes.titlesize": 13.5,
            "axes.labelsize": 12.5,
            "axes.titleweight": "semibold",
            "legend.fontsize": 9.5,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.18,
            "grid.linewidth": 0.65,
            "savefig.facecolor": "white",
        }
    )


def _segment(experiment: object, name: str) -> PhysicalSegment:
    return next(item for item in segment_complete_cycle(experiment) if item.name == name)


def _charge_record(experiment: object) -> ChargeRecord:
    cc = _segment(experiment, "recharge_cc")
    cv = _segment(experiment, "recharge_cv")
    if cc.end + 1 != cv.start or cc.step_number != cv.step_number:
        raise ValueError(
            f"CC/CV continuity is not supported for {experiment.metadata.dataset_id}."
        )
    rows = np.arange(cc.start, cv.end + 1)
    capacity = np.asarray(experiment.numeric["Capacity(Ah)"][rows], dtype=float)
    if np.any(np.diff(capacity) < -1.0e-12):
        raise ValueError(
            f"Capacity(Ah) resets inside CC-CV charge for {experiment.metadata.dataset_id}."
        )
    charged = np.abs(capacity - capacity[0])
    cv_local = int(cv.start - cc.start)
    counter_gap = float(capacity[cv_local] - capacity[cv_local - 1])
    if counter_gap < -1.0e-12:
        raise ValueError(
            f"Capacity(Ah) resets at CC/CV boundary for {experiment.metadata.dataset_id}."
        )
    return ChargeRecord(
        rows=rows,
        charged_capacity_ah=charged,
        voltage_v=np.asarray(experiment.voltage_v[rows], dtype=float),
        current_a=np.asarray(experiment.raw_current_a[rows], dtype=float),
        cc=cc,
        cv=cv,
        cv_local_index=cv_local,
        counter_gap_cc_to_cv_ah=counter_gap,
    )


def _termination_summary(record: ChargeRecord) -> tuple[str, str]:
    end_voltage = float(record.voltage_v[-1])
    end_current = float(record.current_a[-1])
    end_capacity = float(record.charged_capacity_ah[-1])
    if abs(end_capacity - 2.9) <= 0.005:
        evidence = (
            "cambio de carga CV a reposo cerca de 2.900 Ah; compatible con un "
            "límite de capacidad"
        )
    else:
        evidence = (
            f"cambio de carga CV a reposo a {end_voltage:.4f} V y {end_current:.4f} A; "
            "la corriente venía disminuyendo, pero no alcanza un umbral documentado"
        )
    conclusion = (
        "el CSV no contiene receta, código de terminación, evento manual ni bandera "
        "de error; la causa definitiva no puede establecerse"
    )
    return evidence, conclusion


def _plot_stage_curve(axis: plt.Axes, record: ChargeRecord, values: np.ndarray, color: str) -> None:
    boundary = record.cv_local_index
    axis.plot(
        record.charged_capacity_ah[: boundary + 1],
        values[: boundary + 1],
        color=color,
        lw=2.0,
        linestyle="-",
    )
    axis.plot(
        record.charged_capacity_ah[boundary:],
        values[boundary:],
        color=color,
        lw=2.0,
        linestyle="--",
    )


def _charge_figure(
    experiments: dict[str, object],
    mappings: tuple[dict[str, object], ...],
    *,
    title: str,
    subtitle: str | None = None,
) -> tuple[plt.Figure, list[dict[str, object]]]:
    figure, axes = plt.subplots(2, 1, figsize=(12, 7.5), sharex=True)
    rows: list[dict[str, object]] = []
    condition_handles: list[Line2D] = []
    for mapping in mappings:
        experiment = experiments[str(mapping["dataset_id"])]
        record = _charge_record(experiment)
        color = str(mapping["color"])
        _plot_stage_curve(axes[0], record, record.voltage_v, color)
        _plot_stage_curve(axes[1], record, record.current_a, color)
        condition_handles.append(
            Line2D([0], [0], color=color, lw=2.4, label=str(mapping["label"]))
        )
        end_index = int(record.rows[-1])
        next_index = end_index + 1
        evidence, conclusion = _termination_summary(record)
        raw_capacity = np.asarray(
            experiment.numeric["Capacity(Ah)"][record.rows], dtype=float
        )
        rows.append(
            {
                "curve_key": mapping["key"],
                "figure_label": mapping["label"],
                "dataset_id": mapping["dataset_id"],
                "source_file": str(experiment.metadata.path),
                "columns_used": "Voltage(V);Current(A);Capacity(Ah);Step No.;Mode",
                "assignment_basis": mapping["assignment"],
                "charge_start_row_zero_based": int(record.rows[0]),
                "charge_end_row_zero_based": end_index,
                "cc_end_row_zero_based": int(record.cc.end),
                "cv_start_row_zero_based": int(record.cv.start),
                "cc_cv_same_step": record.cc.step_number == record.cv.step_number,
                "capacity_counter_reset_cc_to_cv": False,
                "capacity_counter_gap_cc_to_cv_ah": record.counter_gap_cc_to_cv_ah,
                "charged_capacity_definition": "Q_chg,i=abs(Capacity(Ah)_i-Capacity(Ah)_charge_start)",
                "charge_start_raw_capacity_ah": float(raw_capacity[0]),
                "charge_end_raw_capacity_ah": float(raw_capacity[-1]),
                "recorded_charged_capacity_ah": float(record.charged_capacity_ah[-1]),
                "final_active_voltage_v": float(record.voltage_v[-1]),
                "final_active_current_a": float(record.current_a[-1]),
                "next_record_step_number": int(experiment.numeric["Step No."][next_index]),
                "next_record_current_a": float(experiment.raw_current_a[next_index]),
                "termination_signal_evidence": evidence,
                "termination_cause_conclusion": conclusion,
                "processing": "original samples; no smoothing or extrapolation",
            }
        )

    axes[0].set_ylabel("Voltaje terminal (V)")
    axes[0].set_title("(a) Voltaje durante la carga", loc="left")
    axes[0].legend(handles=condition_handles, loc="lower right", frameon=True)
    axes[1].set_ylabel("Corriente de carga (A)")
    axes[1].set_xlabel("Capacidad cargada (Ah)")
    axes[1].set_title("(b) Corriente durante la carga", loc="left")
    axes[1].legend(
        handles=[
            Line2D([0], [0], color=INK, lw=2.0, linestyle="-", label="Etapa CC"),
            Line2D([0], [0], color=INK, lw=2.0, linestyle="--", label="Etapa CV"),
        ],
        loc="upper right",
        frameon=True,
    )
    for axis in axes:
        axis.grid(True, color=GRAY, alpha=0.17, linewidth=0.65)
    figure.suptitle(title, fontsize=18, weight="bold", color=INK, y=0.985)
    if subtitle:
        figure.text(
            0.5,
            0.945,
            subtitle,
            ha="center",
            va="top",
            fontsize=12.5,
            weight="semibold",
            color=INK,
        )
        top = 0.89
    else:
        top = 0.92
    figure.subplots_adjust(left=0.10, right=0.98, top=top, bottom=0.10, hspace=0.28)
    return figure, rows


def _thermal_segment_record(
    experiment: object,
    mode: str,
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray], float | None]:
    if mode == "discharge":
        segment = _segment(experiment, "main_1c_discharge")
        rows = np.arange(segment.start, segment.end + 1)
        cv_capacity = None
    elif mode == "charge":
        charge = _charge_record(experiment)
        rows = charge.rows
        cv_capacity = float(charge.charged_capacity_ah[charge.cv_local_index])
    else:
        raise ValueError(mode)
    raw_capacity = np.asarray(experiment.numeric["Capacity(Ah)"][rows], dtype=float)
    capacity = np.abs(raw_capacity - raw_capacity[0])
    temperatures = {
        name: np.asarray(experiment.numeric[name][rows], dtype=float)
        for name in SENSOR_COLORS
    }
    return rows, capacity, temperatures, cv_capacity


def _centered_running_median(values: np.ndarray, window: int) -> np.ndarray:
    """Visual-only median that preserves segment edges and cannot overshoot."""

    if window < 3 or window % 2 != 1:
        raise ValueError("The running-median window must be odd and at least 3.")
    source = np.asarray(values, dtype=float)
    result = source.copy()
    half = window // 2
    for index in range(half, source.size - half):
        result[index] = float(np.median(source[index - half : index + half + 1]))
    return result


def _valledupar_air_visual_curve(
    values: np.ndarray,
    capacity_ah: np.ndarray,
    cv_capacity_ah: float | None,
) -> np.ndarray:
    """Smooth Temp7 only, without allowing a charge-stage boundary crossing."""

    if cv_capacity_ah is None:
        return _centered_running_median(
            values, VALLEDUPAR_AIR_MEDIAN_WINDOW_SAMPLES
        )
    boundary = int(np.searchsorted(capacity_ah, cv_capacity_ah, side="left"))
    result = np.asarray(values, dtype=float).copy()
    result[:boundary] = _centered_running_median(
        result[:boundary], VALLEDUPAR_AIR_MEDIAN_WINDOW_SAMPLES
    )
    result[boundary:] = _centered_running_median(
        result[boundary:], VALLEDUPAR_AIR_MEDIAN_WINDOW_SAMPLES
    )
    return result


def _thermal_figure(
    experiments: dict[str, object],
    *,
    mode: str,
) -> tuple[plt.Figure, list[dict[str, object]]]:
    is_charge = mode == "charge"
    figure, axes = plt.subplots(3, 1, figsize=(12, 9), sharey=True)
    prepared: list[tuple[dict[str, object], np.ndarray, dict[str, np.ndarray], float | None]] = []
    all_temperatures: list[np.ndarray] = []
    rows_out: list[dict[str, object]] = []
    for mapping in THERMAL_CONDITIONS:
        experiment = experiments[str(mapping["dataset_id"])]
        rows, capacity, temperatures, cv_capacity = _thermal_segment_record(
            experiment, mode
        )
        prepared.append((mapping, capacity, temperatures, cv_capacity))
        all_temperatures.extend(temperatures.values())
        for sensor, values in temperatures.items():
            visual_processing = "none; original samples plotted"
            if mapping["key"] == "valledupar_cell1" and sensor == "Temp7":
                visual_processing = (
                    "centered running median, 11 samples; charge processed separately "
                    "inside CC and CV; first/last 5 samples of each processed segment unchanged"
                )
            rows_out.append(
                {
                    "phase": "charge" if is_charge else "discharge",
                    "condition_key": mapping["key"],
                    "condition_label": mapping["label"],
                    "setup": mapping["setup"],
                    "dataset_id": mapping["dataset_id"],
                    "source_file": str(experiment.metadata.path),
                    "columns_used": f"{sensor};Capacity(Ah);Step No.;Mode",
                    "sensor": sensor,
                    "sensor_role": SENSOR_LABELS[sensor],
                    "row_start_zero_based": int(rows[0]),
                    "row_end_zero_based": int(rows[-1]),
                    "capacity_extent_ah": float(capacity[-1]),
                    "temperature_min_c": float(np.nanmin(values)),
                    "temperature_max_c": float(np.nanmax(values)),
                    "metrics_input": "original unsmoothed samples",
                    "visual_processing": visual_processing,
                }
            )

    global_min = min(float(np.nanmin(values)) for values in all_temperatures)
    global_max = max(float(np.nanmax(values)) for values in all_temperatures)
    padding = max(0.4, 0.05 * (global_max - global_min))
    for axis, (mapping, capacity, temperatures, cv_capacity) in zip(
        axes, prepared, strict=True
    ):
        for sensor, values in temperatures.items():
            plotted_values = values
            if mapping["key"] == "valledupar_cell1" and sensor == "Temp7":
                plotted_values = _valledupar_air_visual_curve(
                    values, capacity, cv_capacity
                )
            axis.plot(
                capacity,
                plotted_values,
                color=SENSOR_COLORS[sensor],
                lw=1.7,
                label=SENSOR_LABELS[sensor],
            )
        axis.set_ylim(global_min - padding, global_max + padding)
        axis.set_ylabel("Temperatura (°C)")
        axis.set_xlabel(
            "Capacidad cargada (Ah)" if is_charge else "Capacidad descargada (Ah)"
        )
        axis.set_title(
            str(mapping["label"]),
            loc="left",
            color=str(mapping["condition_color"]),
            pad=8,
        )
        axis.grid(True, color=GRAY, alpha=0.17, linewidth=0.65)
        if is_charge and cv_capacity is not None:
            axis.axvline(cv_capacity, color=INK, lw=1.0, linestyle=":")
            ymin, ymax = axis.get_ylim()
            axis.text(
                0.5 * cv_capacity,
                ymax - 0.06 * (ymax - ymin),
                "CC",
                ha="center",
                va="top",
                color=INK,
                fontsize=9.5,
                weight="semibold",
            )
            axis.text(
                cv_capacity + 0.5 * (float(capacity[-1]) - cv_capacity),
                ymax - 0.06 * (ymax - ymin),
                "CV",
                ha="center",
                va="top",
                color=INK,
                fontsize=9.5,
                weight="semibold",
            )
    sensor_handles = [
        Line2D(
            [0],
            [0],
            color=SENSOR_COLORS[sensor],
            lw=2.0,
            label=SENSOR_LABELS[sensor],
        )
        for sensor in ("Temp1", "Temp2", "Temp7")
    ]
    title = (
        "Temperaturas registradas durante la carga CC–CV"
        if is_charge
        else "Temperaturas registradas durante la descarga"
    )
    figure.suptitle(title, fontsize=18, weight="bold", color=INK, y=0.995)
    figure.legend(
        handles=sensor_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.945),
        ncol=3,
        frameon=True,
        fontsize=9.0,
    )
    figure.subplots_adjust(left=0.10, right=0.98, top=0.84, bottom=0.08, hspace=0.46)
    return figure, rows_out


def _save_figure(figure: plt.Figure, source_stem: str, spanish_stem: str) -> None:
    for directory in (
        THESIS,
        PRESENTATION,
        FINAL_ES,
        LOCALIZED_THESIS,
        LOCALIZED_PRESENTATION,
    ):
        directory.mkdir(parents=True, exist_ok=True)
    destinations = (
        (THESIS, source_stem, 320),
        (PRESENTATION, source_stem, 240),
        (FINAL_ES, spanish_stem, 320),
        (LOCALIZED_THESIS, spanish_stem, 320),
        (LOCALIZED_PRESENTATION, spanish_stem, 240),
    )
    for directory, stem, dpi in destinations:
        figure.savefig(directory / f"{stem}.pdf", bbox_inches="tight")
        figure.savefig(directory / f"{stem}.png", dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _write_report(
    charge_rows: list[dict[str, object]],
    temperature_rows: list[dict[str, object]],
) -> None:
    lines = [
        "# Figuras de carga y temperatura de los ensayos propios",
        "",
        "## Correspondencia de registros",
        "",
        "La reasignación eléctrica anterior se mantiene limitada a voltaje y "
        "capacidad de **descarga**. No se trasladó a la carga, la corriente ni las "
        "temperaturas. Cada curva de carga usa voltaje, corriente y `Capacity(Ah)` "
        "del mismo CSV.",
        "",
        "Las figuras térmicas usan solamente archivos de celda 1 que contienen "
        "Temp1, Temp2 y Temp7 en el mismo registro. Temp3 y Temp4 de celda 2 no se "
        "muestran porque sus ubicaciones no están documentadas. El Temp7 compartido "
        "por los dos ensayos simultáneos de Aire 22–24 °C no se alinea muestra a "
        "muestra con el archivo de celda 2. En los paneles térmicos, `Capacity(Ah)` "
        "se usa solo como coordenada del mismo archivo para las temperaturas; esto "
        "no modifica la reasignación de la curva eléctrica de descarga.",
        "",
        "## Carga CC–CV",
        "",
        "La capacidad cargada se calcula como "
        "`Q_chg,i = abs(Capacity(Ah)_i-Capacity(Ah)_inicio_carga)`. En todos los "
        "registros seleccionados, CC y CV pertenecen al mismo paso 5, las filas son "
        "contiguas y el contador es monótono; no se detecta reinicio entre etapas.",
        "",
        "| Curva | Capacidad cargada (Ah) | Voltaje final (V) | Corriente final (A) | Terminación |",
        "|---|---:|---:|---:|---|",
    ]
    unique_charge: dict[str, dict[str, object]] = {}
    for row in charge_rows:
        unique_charge.setdefault(str(row["curve_key"]), row)
    for row in unique_charge.values():
        lines.append(
            f"| {row['figure_label']} | {float(row['recorded_charged_capacity_ah']):.4f} | "
            f"{float(row['final_active_voltage_v']):.4f} | "
            f"{float(row['final_active_current_a']):.4f} | "
            f"{row['termination_signal_evidence']}; {row['termination_cause_conclusion']} |"
        )
    lines.extend(
        [
            "",
            "## Temperaturas por fase",
            "",
            "| Fase | Condición | Sensor | Mínimo (°C) | Máximo (°C) |",
            "|---|---|---|---:|---:|",
        ]
    )
    for row in temperature_rows:
        lines.append(
            f"| {row['phase']} | {row['condition_label']} | {row['sensor']} | "
            f"{float(row['temperature_min_c']):.4f} | "
            f"{float(row['temperature_max_c']):.4f} |"
        )
    lines.extend(
        [
            "",
            "## Procesamiento gráfico",
            "",
            "Las cuatro figuras se construyen a partir de las muestras originales. No "
            "se usan tiempos reconstruidos, no se incluyen reposos sobre ejes de "
            "capacidad y no se extrapolan registros. Como única excepción visual, "
            "Temp7 de Valledupar usa una mediana móvil centrada de 11 muestras para "
            "reducir fluctuaciones rápidas sin sobreimpulso. En carga, el filtro se "
            "aplica por separado dentro de CC y CV; las primeras y últimas cinco "
            "muestras de cada tramo permanecen sin alterar. Temp1 y Temp2 no se "
            "suavizan. Todos los rangos y resultados se calculan con las muestras "
            "originales. Las ciudades se emplean solo como referencias de las "
            "condiciones térmicas.",
            "",
        ]
    )
    (DATA / "uniandes_charge_temperature_figures_report.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )


def build_central_outputs() -> None:
    configure_style()
    experiments = {item.metadata.dataset_id: item for item in load_all()}
    charge_external_figure, charge_external_rows = _charge_figure(
        experiments,
        CHARGE_EXTERNAL,
        title="Carga CC–CV con resistencia externa del montaje",
    )
    charge_contact_figure, charge_contact_rows = _charge_figure(
        experiments,
        CHARGE_CONTACT,
        title="Comparación de carga CC–CV entre montajes",
        subtitle="Rionegro · Aire 22–24 °C",
    )
    thermal_discharge_figure, discharge_temperature_rows = _thermal_figure(
        experiments, mode="discharge"
    )
    thermal_charge_figure, charge_temperature_rows = _thermal_figure(
        experiments, mode="charge"
    )
    figures = (
        charge_external_figure,
        charge_contact_figure,
        thermal_discharge_figure,
        thermal_charge_figure,
    )
    for figure, (source_stem, spanish_stem) in zip(figures, OUTPUTS, strict=True):
        _save_figure(figure, source_stem, spanish_stem)

    charge_rows_by_key: dict[str, dict[str, object]] = {}
    for row in charge_external_rows + charge_contact_rows:
        charge_rows_by_key.setdefault(str(row["curve_key"]), row)
    charge_rows = list(charge_rows_by_key.values())
    temperature_rows = discharge_temperature_rows + charge_temperature_rows
    _write_csv(DATA / "uniandes_charge_audit.csv", charge_rows)
    _write_csv(DATA / "uniandes_temperature_ranges.csv", temperature_rows)

    correspondence_rows = [
        {
            "scope": "charge electrical",
            "condition": item["label"],
            "dataset_id": item["dataset_id"],
            "source_file": str(experiments[str(item["dataset_id"])].metadata.path),
            "columns": "Voltage(V);Current(A);Capacity(Ah);Step No.;Mode",
            "status": "confirmed same-file correspondence",
        }
        for item in CHARGE_EXTERNAL
    ]
    correspondence_rows.extend(
        {
            "scope": "thermal charge and discharge",
            "condition": item["label"],
            "dataset_id": item["dataset_id"],
            "source_file": str(experiments[str(item["dataset_id"])].metadata.path),
            "columns": "Temp1;Temp2;Temp7;Capacity(Ah);Step No.;Mode",
            "status": "confirmed same-file correspondence and documented sensor roles",
        }
        for item in THERMAL_CONDITIONS
    )
    correspondence_rows.extend(
        [
            {
                "scope": "thermal external-resistance cell 2",
                "condition": "Rionegro · Aire 22–24 °C",
                "dataset_id": "cycle_cell2_undeclared",
                "source_file": str(experiments["cycle_cell2_undeclared"].metadata.path),
                "columns": "Temp3;Temp4 (Temp7 shared in separate cell-1 file)",
                "status": "not plotted: Temp3/Temp4 locations undocumented and cross-file Temp7 rows not synchronized",
            },
            {
                "scope": "thermal external-resistance cell 2",
                "condition": "Cali · Aire 29–31 °C",
                "dataset_id": "cycle_cell2_30c",
                "source_file": str(experiments["cycle_cell2_30c"].metadata.path),
                "columns": "Temp3;Temp4",
                "status": "not plotted: Temp3/Temp4 locations undocumented and no same-file Temp7",
            },
        ]
    )
    _write_csv(DATA / "uniandes_charge_temperature_correspondence.csv", correspondence_rows)
    _write_report(charge_rows, temperature_rows)


def main() -> None:
    build_central_outputs()
    print("Uniandes charge and temperature figures generated from original Chroma records.")


if __name__ == "__main__":
    main()
