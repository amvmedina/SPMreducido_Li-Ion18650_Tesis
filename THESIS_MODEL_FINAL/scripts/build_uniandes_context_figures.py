"""Build the Uniandes voltage and geographical-context figures.

This module is presentation-only.  It does not fit a model or alter the raw
measurements used by any diagnostic.  A segmented running median is used only
for the displayed voltage lines; all numerical summaries are calculated from
the original samples.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analysis.own_case_study_reference import discharge_curve, segment_complete_cycle  # noqa: E402


DATA = ROOT / "thesis_data"
BOUNDARY_PATH = ROOT / "data" / "geography" / "ne_110m_admin_0_countries.geojson"

INK = "#172033"
BLUE = "#2563A6"
ORANGE = "#E07A2D"
GREEN = "#2E8B57"
RED = "#C7473D"
PURPLE = "#7656A5"
GRAY = "#6B7280"

MEDIAN_WINDOW_SAMPLES = 7
PROTECTED_VOLTAGE_JUMP_V = 0.010
PROTECTED_CURRENT_JUMP_A = 0.05
FINAL_DISCHARGE_MAPPING = (
    {
        "curve_id": "rionegro_reference",
        "legend": "Montaje de referencia",
        "defective_legend": "Rionegro · Aire 22–24 °C",
        "city_context": "Rionegro",
        "condition": "Aire 22–24 °C",
        "electrical_dataset_id": "cycle_cell1_ambient",
        "electrical_cell_id": "cell 1",
        "temperature_dataset_id": "cycle_cell1_ambient",
        "air_source_dataset": "cycle_cell1_ambient",
        "contact": "montaje de referencia",
        "color": GREEN,
        "linewidth": 2.7,
        "alignment": "temperatura y voltaje están en el mismo archivo de celda 1; no se usa en cálculos entre archivos",
    },
    {
        "curve_id": "rionegro_defective",
        "legend": "Montaje con resistencia externa",
        "defective_legend": "Rionegro · Aire 22–24 °C",
        "city_context": "Rionegro",
        "condition": "Aire 22–24 °C",
        "electrical_dataset_id": "cycle_cell2_undeclared",
        "electrical_cell_id": "cell 2",
        "temperature_dataset_id": "cycle_cell2_undeclared",
        "air_source_dataset": "cycle_cell1_ambient",
        "contact": "contactos defectuosos",
        "color": PURPLE,
        "linewidth": 1.9,
        "alignment": "Temp7 es compartido y medido, pero la alineación muestra a muestra entre archivos no está establecida",
    },
    {
        "curve_id": "cali_defective",
        "legend": "Cali · Aire 29–31 °C",
        "defective_legend": "Cali · Aire 29–31 °C",
        "city_context": "Cali",
        "condition": "Aire 29–31 °C",
        "electrical_dataset_id": "cycle_cell2_30c",
        "electrical_cell_id": "cell 2",
        "temperature_dataset_id": "cycle_cell2_30c",
        "air_source_dataset": "cycle_cell1_30c",
        "contact": "contactos defectuosos",
        "color": BLUE,
        "linewidth": 1.9,
        "alignment": "la correspondencia temporal entre el registro eléctrico y el registro Temp7 no está establecida",
    },
    {
        "curve_id": "valledupar_defective",
        "legend": "Valledupar · Aire 33–37 °C",
        "defective_legend": "Valledupar · Aire 33–37 °C",
        "city_context": "Valledupar",
        "condition": "Aire 33–37 °C",
        "electrical_dataset_id": "cycle_cell1_30c",
        "electrical_cell_id": "cell 1",
        "temperature_dataset_id": "cycle_cell1_35c",
        "air_source_dataset": "cycle_cell1_35c",
        "contact": "contactos defectuosos",
        "color": ORANGE,
        "linewidth": 1.9,
        "alignment": "mapeo eléctrico y térmico proceden de registros distintos; no se combinan muestra a muestra",
    },
)

MAPPING_BY_ID = {item["curve_id"]: item for item in FINAL_DISCHARGE_MAPPING}
DEFECTIVE_CURVE_IDS = (
    "rionegro_defective",
    "cali_defective",
    "valledupar_defective",
)
CONTACT_COMPARISON_CURVE_IDS = ("rionegro_reference", "rionegro_defective")
EXCLUDED_ELECTRICAL_DATASET = "cycle_cell1_35c"

CITY_CONTEXT = (
    {
        "city": "Rionegro",
        "condition": "Aire 22–24 °C",
        "latitude": 6.1490188,
        "longitude": -75.3789523,
        "color": PURPLE,
        "source": "Ministerio de Vivienda, Ciudad y Territorio",
        "url": "https://www.minvivienda.gov.co/node/38738",
        "offset": (0.45, 0.25),
    },
    {
        "city": "Cali",
        "condition": "Aire 29–31 °C",
        "latitude": 3.4518914,
        "longitude": -76.5241715,
        "color": BLUE,
        "source": "Ministerio de Vivienda, Ciudad y Territorio",
        "url": "https://ventanilla.minvivienda.gov.co/node/39659",
        "offset": (0.45, -0.35),
    },
    {
        "city": "Valledupar",
        "condition": "Aire 33–37 °C",
        "latitude": 10.4742449,
        "longitude": -73.2436335,
        "color": ORANGE,
        "source": "Ministerio de Vivienda, Ciudad y Territorio",
        "url": "https://www.minvivienda.gov.co/node/39057",
        "offset": (0.42, 0.18),
    },
)


def _main_discharge(experiment: object) -> dict[str, np.ndarray | float]:
    segment = next(
        item
        for item in segment_complete_cycle(experiment)
        if item.name == "main_1c_discharge"
    )
    return discharge_curve(experiment, segment)


def segmented_running_median(
    voltage_v: np.ndarray,
    current_a: np.ndarray,
    *,
    window: int = MEDIAN_WINDOW_SAMPLES,
    voltage_jump_v: float = PROTECTED_VOLTAGE_JUMP_V,
    current_jump_a: float = PROTECTED_CURRENT_JUMP_A,
) -> tuple[np.ndarray, np.ndarray]:
    """Median-filter only windows that do not cross protected physical edges."""

    if window < 3 or window % 2 != 1:
        raise ValueError("The visualization median window must be odd and at least 3.")
    voltage = np.asarray(voltage_v, dtype=float)
    current = np.asarray(current_a, dtype=float)
    if voltage.shape != current.shape:
        raise ValueError("Voltage and current must have the same shape.")
    protected_edges = (np.abs(np.diff(voltage)) >= voltage_jump_v) | (
        np.abs(np.diff(current)) >= current_jump_a
    )
    filtered = voltage.copy()
    half = window // 2
    for index in range(half, voltage.size - half):
        start = index - half
        stop = index + half + 1
        if np.any(protected_edges[start : stop - 1]):
            continue
        filtered[index] = float(np.median(voltage[start:stop]))
    return filtered, protected_edges


def _linear_visual_curve(
    discharged_capacity_ah: np.ndarray,
    voltage_v: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Linear interpolation after duplicate collapse; it cannot overshoot."""

    order = np.argsort(discharged_capacity_ah, kind="stable")
    sorted_capacity = np.asarray(discharged_capacity_ah, dtype=float)[order]
    sorted_voltage = np.asarray(voltage_v, dtype=float)[order]
    unique_capacity, first, counts = np.unique(
        sorted_capacity, return_index=True, return_counts=True
    )
    collapsed = np.asarray(
        [
            np.mean(sorted_voltage[index : index + count])
            for index, count in zip(first, counts, strict=True)
        ]
    )
    line_capacity = np.linspace(
        float(unique_capacity[0]), float(unique_capacity[-1]), 1001
    )
    return line_capacity, np.interp(line_capacity, unique_capacity, collapsed)


def _termination_evidence(
    voltage_v: float,
    denominator_ah: float,
) -> tuple[str, str]:
    """Describe signal evidence without assigning an undocumented stop cause."""

    if abs(voltage_v - 2.6) <= 0.005:
        evidence = (
            "cambio de paso desde descarga activa a corriente cero a aproximadamente "
            "2.600 V; el patrón de señal es compatible con un umbral de voltaje"
        )
    elif abs(denominator_ah - 2.9) <= 0.005:
        evidence = (
            "cambio de paso desde descarga activa a corriente cero después de "
            "aproximadamente 2.900 Ah, con voltaje cargado de 2.8342 V; el patrón "
            "de señal es compatible con un umbral de capacidad"
        )
    else:
        evidence = (
            "cambio de paso desde descarga activa a corriente cero sin un campo "
            "explícito de terminación"
        )
    conclusion = (
        "no distinguible a partir del CSV: no se registra el límite de la receta, "
        "un código de terminación, un evento del operador ni una bandera de error"
    )
    return evidence, conclusion


def prepare_voltage_visualization(
    experiments: dict[str, object],
) -> tuple[dict[str, dict[str, np.ndarray]], list[dict[str, object]]]:
    visual: dict[str, dict[str, np.ndarray]] = {}
    summary: list[dict[str, object]] = []
    for mapping in FINAL_DISCHARGE_MAPPING:
        curve_id = str(mapping["curve_id"])
        dataset_id = str(mapping["electrical_dataset_id"])
        experiment = experiments[dataset_id]
        curve = _main_discharge(experiment)
        rows = np.asarray(curve["rows"], dtype=int)
        capacity = np.asarray(experiment.numeric["Capacity(Ah)"][rows], dtype=float)
        raw_discharged_capacity = np.abs(capacity - capacity[0])
        raw_voltage = np.asarray(curve["voltage_v"], dtype=float)
        raw_current = np.asarray(experiment.raw_current_a[rows], dtype=float)
        filtered_voltage, protected_edges = segmented_running_median(
            raw_voltage, raw_current
        )
        line_capacity, line_voltage = _linear_visual_curve(
            raw_discharged_capacity, filtered_voltage
        )
        visual[curve_id] = {
            "raw_discharged_capacity_ah": raw_discharged_capacity,
            "raw_voltage_v": raw_voltage,
            "line_discharged_capacity_ah": line_capacity,
            "line_voltage_v": line_voltage,
        }

        air_source_id = str(mapping["air_source_dataset"])
        air_source = experiments[air_source_id]
        temp7 = np.asarray(air_source.numeric["Temp7"], dtype=float)
        temperature_dataset_id = str(mapping["temperature_dataset_id"])
        temperature_experiment = experiments[temperature_dataset_id]
        temperature_channels = [
            name
            for name in temperature_experiment.numeric
            if name.startswith("Temp") and name != "Temp7"
        ]
        unique_voltage = np.unique(raw_voltage)
        positive_voltage_steps = np.diff(unique_voltage)
        positive_voltage_steps = positive_voltage_steps[positive_voltage_steps > 0.0]
        adjustment_mv = 1000.0 * np.abs(filtered_voltage - raw_voltage)
        denominator_ah = abs(float(capacity[-1] - capacity[0]))
        end_index = int(rows[-1])
        next_index = end_index + 1 if end_index + 1 < experiment.voltage_v.size else None
        termination_evidence, termination_conclusion = _termination_evidence(
            float(raw_voltage[-1]), denominator_ah
        )
        summary.append(
            {
                "curve_id": curve_id,
                "audit_curve_label": (
                    f"{mapping['city_context']} · {mapping['condition']} · "
                    + (
                        "montaje de referencia"
                        if curve_id == "rionegro_reference"
                        else "montaje con resistencia externa"
                    )
                ),
                "figure_legend": mapping["defective_legend"],
                "city_context": mapping["city_context"],
                "measured_air_condition_label": mapping["condition"],
                "mapping_scope": "discharge Voltage(V) and discharged capacity in Ah over the recorded interval only",
                "electrical_dataset_id": dataset_id,
                "electrical_cell_id_internal": mapping["electrical_cell_id"],
                "electrical_source_file": str(experiment.metadata.path),
                "voltage_source_column": "Voltage(V)",
                "discharged_capacity_source_column": "Capacity(Ah)",
                "discharge_source_row_start_zero_based": int(rows[0]),
                "discharge_source_row_end_zero_based": int(rows[-1]),
                "discharged_capacity_definition": "Q_dis,i=abs(Capacity(Ah)_i-Capacity(Ah)_start) within the selected source main_1c_discharge",
                "recorded_capacity_extent_ah": denominator_ah,
                "recorded_capacity_extent_provenance": "absolute difference between first and last Capacity(Ah) values in the selected source main_1c_discharge",
                "axis_normalization_denominator": "none",
                "prior_normalization_denominator_ah": denominator_ah,
                "axis_start_ah": float(raw_discharged_capacity[0]),
                "axis_end_ah": float(raw_discharged_capacity[-1]),
                "normalized_between_zero_and_one": False,
                "physical_soc_claim": "none; absolute discharged-capacity coordinate over recorded samples only",
                "discharge_start_capacity_ah": float(capacity[0]),
                "discharge_end_capacity_ah": float(capacity[-1]),
                "discharge_end_voltage_v": float(raw_voltage[-1]),
                "discharge_end_current_a": float(raw_current[-1]),
                "discharge_end_step_number": int(experiment.numeric["Step No."][end_index]),
                "discharge_end_mode": int(experiment.numeric["Mode"][end_index]),
                "next_record_step_number": int(experiment.numeric["Step No."][next_index]) if next_index is not None else "",
                "next_record_current_a": float(experiment.raw_current_a[next_index]) if next_index is not None else "",
                "termination_signal_evidence": termination_evidence,
                "termination_cause_conclusion": termination_conclusion,
                "temperature_dataset_id": temperature_dataset_id,
                "temperature_source_file": str(temperature_experiment.metadata.path),
                "temperature_columns_retained": ",".join(temperature_channels),
                "temp7_source_dataset": air_source_id,
                "temp7_min_c": float(np.nanmin(temp7)),
                "temp7_max_c": float(np.nanmax(temp7)),
                "cross_file_alignment_status": mapping["alignment"],
                "contact_arrangement": mapping["contact"],
                "charge_current_capacity_metadata_action": "not reassigned; remains with each source file pending provenance verification",
                "raw_discharge_samples": int(raw_voltage.size),
                "unique_capacity_coordinates": int(np.unique(raw_discharged_capacity).size),
                "duplicate_capacity_coordinates": int(raw_discharged_capacity.size - np.unique(raw_discharged_capacity).size),
                "minimum_positive_voltage_increment_mv": float(
                    1000.0 * np.min(positive_voltage_steps)
                ),
                "median_abs_adjacent_voltage_change_mv": float(
                    1000.0 * np.median(np.abs(np.diff(raw_voltage)))
                ),
                "maximum_abs_adjacent_voltage_change_mv": float(
                    1000.0 * np.max(np.abs(np.diff(raw_voltage)))
                ),
                "visualization_method": "segmented centered running median followed by linear interpolation",
                "visualization_window_samples": MEDIAN_WINDOW_SAMPLES,
                "protected_voltage_jump_threshold_mv": 1000.0 * PROTECTED_VOLTAGE_JUMP_V,
                "protected_current_jump_threshold_a": PROTECTED_CURRENT_JUMP_A,
                "protected_edge_count": int(np.count_nonzero(protected_edges)),
                "median_abs_visual_adjustment_mv": float(np.median(adjustment_mv)),
                "maximum_abs_visual_adjustment_mv": float(np.max(adjustment_mv)),
                "metrics_input": "original measurements only",
            }
        )
    return visual, summary


def _make_voltage_figure(
    visual: dict[str, dict[str, np.ndarray]],
    curve_ids: tuple[str, ...],
    *,
    title: str,
    subtitle: str | None = None,
) -> plt.Figure:
    figure, axis = plt.subplots(constrained_layout=False)
    for curve_id in curve_ids:
        values = visual[curve_id]
        mapping = MAPPING_BY_ID[curve_id]
        axis.scatter(
            values["raw_discharged_capacity_ah"],
            values["raw_voltage_v"],
            s=6,
            color=mapping["color"],
            alpha=0.14,
            edgecolors="none",
            zorder=2,
        )
        axis.plot(
            values["line_discharged_capacity_ah"],
            values["line_voltage_v"],
            color=mapping["color"],
            linewidth=mapping["linewidth"],
            label=(mapping["legend"] if subtitle else mapping["defective_legend"]),
            zorder=3,
        )
    axis.set_xlabel("Capacidad descargada (Ah)")
    axis.set_ylabel("Voltaje terminal durante la descarga (V)")
    if subtitle:
        axis.set_title(subtitle, pad=10, fontsize=13.0, color=INK)
    axis.grid(True, color=GRAY, alpha=0.17, linewidth=0.65)
    axis.legend(
        loc="upper right",
        fontsize=9.6,
        frameon=True,
        framealpha=0.92,
        handlelength=3.0,
    )
    figure.suptitle(
        title,
        fontsize=18,
        weight="bold",
        color=INK,
        y=0.975,
    )
    figure.subplots_adjust(left=0.10, right=0.98, top=0.86, bottom=0.13)
    return figure


def make_discharge_figures(
    experiments: dict[str, object],
) -> tuple[plt.Figure, plt.Figure, list[dict[str, object]]]:
    visual, summary = prepare_voltage_visualization(experiments)
    defective = _make_voltage_figure(
        visual,
        DEFECTIVE_CURVE_IDS,
        title="Voltaje en descarga con resistencia externa del montaje",
    )
    contact_comparison = _make_voltage_figure(
        visual,
        CONTACT_COMPARISON_CURVE_IDS,
        title="Comparación del voltaje en descarga entre montajes",
        subtitle="Rionegro · Aire 22–24 °C",
    )
    return defective, contact_comparison, summary


def _colombia_rings() -> list[np.ndarray]:
    payload = BOUNDARY_PATH.read_bytes()
    expected_hash = "6866c877d39cba9c357620878839b336d569f8c662d3cfab4cb1dbe2d39c977f"
    if hashlib.sha256(payload).hexdigest() != expected_hash:
        raise ValueError("The vendored Natural Earth boundary hash changed.")
    geojson = json.loads(payload)
    matches = [
        feature
        for feature in geojson["features"]
        if feature["properties"].get("ADM0_A3") == "COL"
        or feature["properties"].get("ADMIN") == "Colombia"
    ]
    if len(matches) != 1:
        raise ValueError("Expected exactly one Colombia feature in Natural Earth.")
    geometry = matches[0]["geometry"]
    polygons = (
        [geometry["coordinates"]]
        if geometry["type"] == "Polygon"
        else geometry["coordinates"]
    )
    return [np.asarray(polygon[0], dtype=float) for polygon in polygons]


def make_geographical_figure() -> plt.Figure:
    figure, axis = plt.subplots(constrained_layout=False)
    for ring in _colombia_rings():
        axis.add_patch(
            Polygon(
                ring,
                closed=True,
                facecolor="#EEF3F7",
                edgecolor="#44546A",
                linewidth=1.2,
                zorder=1,
            )
        )
    for city in CITY_CONTEXT:
        longitude = float(city["longitude"])
        latitude = float(city["latitude"])
        dx, dy = city["offset"]
        axis.scatter(
            longitude,
            latitude,
            s=78,
            color=city["color"],
            edgecolor="white",
            linewidth=1.1,
            zorder=4,
        )
        axis.annotate(
            f"{city['city']}\n{city['condition']}",
            xy=(longitude, latitude),
            xytext=(longitude + dx, latitude + dy),
            ha="left",
            va="center",
            fontsize=10.2,
            color=INK,
            weight="semibold",
            arrowprops={"arrowstyle": "-", "color": GRAY, "lw": 0.8},
            bbox={
                "boxstyle": "round,pad=0.28",
                "facecolor": "white",
                "edgecolor": "#D5DDE7",
                "alpha": 0.96,
            },
            zorder=5,
        )
    axis.set_xlim(-79.5, -66.3)
    axis.set_ylim(-4.8, 13.0)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("Longitud (°)")
    axis.set_ylabel("Latitud (°)")
    axis.grid(True, color=GRAY, alpha=0.14, linewidth=0.6)
    figure.suptitle(
        "Localidades de referencia para las condiciones térmicas",
        fontsize=18,
        weight="bold",
        color=INK,
        y=0.975,
    )
    figure.subplots_adjust(left=0.10, right=0.98, top=0.89, bottom=0.10)
    return figure


def write_context_documentation(summary: list[dict[str, object]]) -> None:
    summary_path = DATA / "uniandes_experimental_context_summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)

    mapping_fields = [
        "status",
        "curve_id",
        "figure_legend",
        "city_context",
        "measured_air_condition_label",
        "mapping_scope",
        "electrical_dataset_id",
        "electrical_cell_id_internal",
        "electrical_source_file",
        "voltage_source_column",
        "discharged_capacity_source_column",
        "discharged_capacity_definition",
        "temperature_dataset_id",
        "temperature_source_file",
        "temperature_columns_retained",
        "temp7_source_dataset",
        "cross_file_alignment_status",
        "contact_arrangement",
        "charge_current_capacity_metadata_action",
    ]
    mapping_rows = [
        {"status": "included", **{key: row.get(key, "") for key in mapping_fields if key != "status"}}
        for row in summary
    ]
    excluded_experiment = Path(
        r"C:\Users\amvel\Documents\Andes\Tesis\TRC18650experimental data\amvTRC18650LiIonCycle1C T35C cell1.csv"
    )
    mapping_rows.append(
        {
            "status": "excluded from corrected discharge comparisons",
            "curve_id": "excluded_cycle_cell1_35c_electrical",
            "figure_legend": "",
            "city_context": "",
            "measured_air_condition_label": "Aire 33–37 °C temperature record retained",
            "mapping_scope": "electrical voltage/discharged-capacity curve excluded; source preserved",
            "electrical_dataset_id": EXCLUDED_ELECTRICAL_DATASET,
            "electrical_cell_id_internal": "cell 1",
            "electrical_source_file": str(excluded_experiment),
            "voltage_source_column": "Voltage(V) not plotted",
            "discharged_capacity_source_column": "Capacity(Ah) not plotted",
            "discharged_capacity_definition": "not used in corrected figures",
            "temperature_dataset_id": "cycle_cell1_35c",
            "temperature_source_file": str(excluded_experiment),
            "temperature_columns_retained": "Temp1,Temp2,Temp7",
            "temp7_source_dataset": "cycle_cell1_35c",
            "cross_file_alignment_status": "temperature association retained; no transfer to electrical source rows",
            "contact_arrangement": "reported problematic; electrical curve excluded by user-confirmed mapping",
            "charge_current_capacity_metadata_action": "not reassigned; source file preserved",
        }
    )
    with (DATA / "uniandes_corrected_discharge_mapping_audit.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=mapping_fields)
        writer.writeheader()
        writer.writerows(mapping_rows)

    axis_fields = [
        "curve_id",
        "audit_curve_label",
        "electrical_dataset_id",
        "electrical_source_file",
        "discharge_source_row_start_zero_based",
        "discharge_source_row_end_zero_based",
        "discharged_capacity_definition",
        "recorded_capacity_extent_ah",
        "recorded_capacity_extent_provenance",
        "axis_normalization_denominator",
        "prior_normalization_denominator_ah",
        "axis_start_ah",
        "axis_end_ah",
        "normalized_between_zero_and_one",
        "physical_soc_claim",
        "discharge_start_capacity_ah",
        "discharge_end_capacity_ah",
        "discharge_end_voltage_v",
        "discharge_end_current_a",
        "discharge_end_step_number",
        "discharge_end_mode",
        "next_record_step_number",
        "next_record_current_a",
        "termination_signal_evidence",
        "termination_cause_conclusion",
    ]
    with (DATA / "uniandes_discharge_axis_audit.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=axis_fields)
        writer.writeheader()
        writer.writerows(
            {key: row.get(key, "") for key in axis_fields} for row in summary
        )

    audit_lines = [
        "# Auditoría del eje de capacidad de las descargas Uniandes",
        "",
        "El eje utilizado anteriormente no era un SOC físico. Se calculaba como "
        "`SOC_ref = 1 - abs(Q-Q_inicio)/abs(Q_fin-Q_inicio)`, por lo que el primer "
        "registro se forzaba a 1 y el último a 0. Los CSV no establecen que esos "
        "extremos representen una celda completamente cargada y completamente descargada.",
        "",
        "Las figuras corregidas usan `Q_desc,i = abs(Q_i-Q_inicio)` y el eje "
        "**Capacidad descargada (Ah)**. Cada curva comienza en 0 Ah y termina en la "
        "capacidad realmente recorrida dentro de su registro; no se normaliza entre "
        "0 y 1, no utiliza denominador, no se extrapola y no afirma una descarga completa.",
        "",
        "| Curva | Capacidad recorrida (Ah) | Último voltaje (V) | Última corriente (A) | "
        "Capacidad final (Ah) | Evidencia de terminación | Conclusión |",
        "|---|---:|---:|---:|---:|---|---|",
    ]
    for row in summary:
        audit_lines.append(
            f"| {row['audit_curve_label']} | {float(row['recorded_capacity_extent_ah']):.4f} | "
            f"{float(row['discharge_end_voltage_v']):.4f} | "
            f"{float(row['discharge_end_current_a']):.4f} | "
            f"{float(row['discharge_end_capacity_ah']):.4f} | "
            f"{row['termination_signal_evidence']} | "
            f"{row['termination_cause_conclusion']} |"
        )
    (DATA / "uniandes_discharge_axis_audit.md").write_text(
        "\n".join(audit_lines) + "\n", encoding="utf-8"
    )

    affected_artifacts = [
        {
            "artifact": "thesis_data/own_data_quality_audit.csv",
            "status": "unchanged",
            "reason": "dataset-level raw diagnostics remain tied to source files",
        },
        {
            "artifact": "thesis_data/own_data_apparent_resistance_transitions.csv",
            "status": "unchanged",
            "reason": "R_app remains calculated from original adjacent source records",
        },
        {
            "artifact": "thesis_data/own_case_study_cross_experiment_metrics.csv",
            "status": "not relabelled by city/temperature",
            "reason": "existing values remain dataset-ID diagnostics; temperature interpretation is not transferred",
        },
        {
            "artifact": "thesis_data/uniandes_experimental_context_summary.csv",
            "status": "regenerated",
            "reason": "now records the corrected discharge-only electrical mapping and separate thermal provenance",
        },
        {
            "artifact": "thesis_data/caption_figura_uniandes_voltaje.md",
            "status": "superseded",
            "reason": "combined discharge figure replaced by two figures",
        },
    ]
    with (DATA / "uniandes_mapping_affected_artifacts.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(affected_artifacts[0]))
        writer.writeheader()
        writer.writerows(affected_artifacts)

    city_path = DATA / "uniandes_geographical_context_coordinates.csv"
    city_fields = ["city", "condition", "latitude", "longitude", "source", "url"]
    with city_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=city_fields)
        writer.writeheader()
        writer.writerows(
            {key: city[key] for key in city_fields} for city in CITY_CONTEXT
        )

    (DATA / "caption_figura_uniandes_voltaje.md").write_text(
        "# Figura combinada reemplazada\n\n"
        "La figura combinada anterior fue reemplazada por las figuras separadas de "
        "contactos defectuosos y comparación de montajes. No usar este archivo como "
        "pie de figura activo.\n",
        encoding="utf-8",
    )
    (DATA / "caption_figura_uniandes_contactos_defectuosos.md").write_text(
        "# Pie de figura — resistencia externa del montaje\n\n"
        "Voltaje terminal medido durante las descargas asociadas con Rionegro "
        "(Aire 22–24 °C), Cali (Aire 29–31 °C) y Valledupar (Aire 33–37 °C), "
        "todas realizadas con el montaje reportado con resistencia externa. "
        "Las ciudades proporcionan únicamente contexto climático para los intervalos "
        "experimentales; los ensayos no se realizaron necesariamente en esas ciudades. "
        "Las diferencias entre curvas no se interpretan como un efecto aislado de la "
        "temperatura. Los puntos tenues son las mediciones originales y las líneas usan "
        "una mediana móvil segmentada de siete muestras solo para visualización. "
        "El eje horizontal muestra la capacidad descargada acumulada en Ah dentro de "
        "cada registro; no representa SOC ni demuestra una descarga completa.\n",
        encoding="utf-8",
    )
    (DATA / "caption_figura_uniandes_montajes_contacto.md").write_text(
        "# Pie de figura — comparación de montajes\n\n"
        "Comparación del voltaje de descarga para dos montajes de contacto asociados "
        "con Aire 22–24 °C. El montaje de referencia corresponde a la celda 1 y el "
        "montaje con resistencia externa a la celda 2; por tanto, no constituye una "
        "comparación antes-después sobre la misma celda ni una estimación aislada de "
        "resistencia de contacto. Los puntos tenues son los datos originales y las "
        "líneas usan una mediana móvil segmentada de siete muestras, seguida de "
        "interpolación lineal, únicamente para visualización. El eje horizontal "
        "muestra la capacidad descargada acumulada en Ah dentro del intervalo "
        "registrado; no identifica SOC ni demuestra una descarga completa.\n",
        encoding="utf-8",
    )
    (DATA / "caption_figura_contexto_geografico_colombia.md").write_text(
        "# Pie de figura — contexto geográfico\n\n"
        "Localidades colombianas empleadas como contexto geográfico para los "
        "intervalos de temperatura del aire medidos durante los ensayos. Los "
        "intervalos indicados corresponden exclusivamente a las condiciones "
        "experimentales asociadas y no representan rangos climáticos permanentes "
        "de las ciudades. Los experimentos tampoco reproducen ciclos meteorológicos "
        "diarios completos. Límite cartográfico: Natural Earth 1:110m, versión "
        "5.1.2. Coordenadas urbanas: Ministerio de Vivienda, Ciudad y Territorio.\n",
        encoding="utf-8",
    )

    report_lines = [
        "# Contexto experimental y presentación de datos Uniandes",
        "",
        "## Canal de aire compartido",
        "",
        "`cycle_cell1_ambient` y `cycle_cell2_undeclared` se realizaron simultáneamente y compartieron el canal medido de aire `Temp7` (22.4896–23.4919 °C). El canal no se copia dentro de las señales de la celda 2. `Temp3` y `Temp4` permanecen como mediciones independientes de ubicación no documentada.",
        "",
        "Los archivos tienen diferentes cantidades de registros y no incluyen tiempo medido. Por ello, la simultaneidad no demuestra alineación muestra a muestra y no se calculan diferencias de temperatura entre archivos.",
        "",
        "## Etiquetas de condición",
        "",
        "Las etiquetas `Aire 22–24 °C`, `Aire 29–31 °C` y `Aire 33–37 °C` resumen rangos medidos de `Temp7`; no son referencias programadas ni afirmaciones climáticas permanentes.",
        "",
        "## Mapeo eléctrico corregido para las figuras de descarga",
        "",
        "La corrección se aplica exclusivamente a `Voltage(V)` y a la capacidad descargada acumulada del intervalo seleccionado. Rionegro/montaje con resistencia externa usa `cycle_cell2_undeclared`; Cali/montaje con resistencia externa usa `cycle_cell2_30c`; Valledupar/montaje con resistencia externa usa el voltaje y la capacidad de `cycle_cell1_30c`; Rionegro/referencia usa `cycle_cell1_ambient`. La curva eléctrica de `cycle_cell1_35c` queda excluida.",
        "",
        "Para Valledupar, las temperaturas permanecen en `cycle_cell1_35c` (`Temp1`, `Temp2`, `Temp7`). No se transfieren columnas térmicas hacia `cycle_cell1_30c` ni se combinan registros muestra a muestra. Carga, corriente y capacidades totales permanecen asociadas con sus archivos originales hasta que su procedencia sea verificada.",
        "",
        "La auditoría confirmó que el eje anterior forzaba el primer registro a `SOC_ref=1` y el último a `SOC_ref=0` usando la propia capacidad final de cada archivo como denominador. Los CSV no demuestran que esos extremos sean 100 % y 0 % de SOC. Las figuras corregidas muestran `Q_desc = abs(Q-Q_inicio)` en Ah. Cada curva comienza en 0 Ah y termina en su capacidad descargada realmente registrada, sin normalización, sin extrapolación y sin afirmar una descarga completa.",
        "",
        "## Auditoría de apariencia y método visual",
        "",
        "Las descargas seleccionadas no contienen coordenadas de capacidad descargada duplicadas. El aspecto dentado procede principalmente de la cuantización/variación de voltaje y de saltos reales, no del trazado. Para mejorar la lectura se usa una mediana móvil centrada de 7 muestras, solamente dentro de la descarga continua. No se filtran ventanas que crucen saltos de voltaje ≥10 mV, cambios de corriente ≥0.05 A, extremos ni el punto de corte. Después se usa interpolación lineal; no se emplean splines y no puede aparecer sobreimpulso.",
        "",
        "Los puntos originales se dibujan debajo de las líneas. Todos los indicadores de resistencia aparente, capacidad y desviación continúan calculándose con las mediciones originales; el panel de desviación se retira únicamente de la figura de presentación.",
        "",
        "## Fuentes geográficas",
        "",
        "- Límite de Colombia: Natural Earth, `ne_110m_admin_0_countries.geojson`, versión 5.1.2, SHA-256 `6866c877d39cba9c357620878839b336d569f8c662d3cfab4cb1dbe2d39c977f`.",
        "- Rionegro: https://www.minvivienda.gov.co/node/38738",
        "- Cali: https://ventanilla.minvivienda.gov.co/node/39659",
        "- Valledupar: https://www.minvivienda.gov.co/node/39057",
        "",
        "Las ciudades proporcionan contexto geográfico. No se incorporan referencias meteorológicas ni se afirma que un intervalo experimental sea habitual a partir de un solo día.",
        "",
    ]
    (DATA / "uniandes_experimental_context_report.md").write_text(
        "\n".join(report_lines), encoding="utf-8"
    )


if __name__ == "__main__":
    from scripts import build_final_thesis_package as package
    from scripts.audit_own_trc18650_data import load_all

    package.configure_style()
    experiments = {item.metadata.dataset_id: item for item in load_all()}
    defective, contact_comparison, rows = make_discharge_figures(experiments)
    package.save_figure(defective, "uniandes_defective_contacts_discharge_final")
    package.save_figure(
        contact_comparison, "uniandes_contact_arrangements_comparison_final"
    )
    package.save_figure(
        make_geographical_figure(), "uniandes_geographical_context_final"
    )
    write_context_documentation(rows)
