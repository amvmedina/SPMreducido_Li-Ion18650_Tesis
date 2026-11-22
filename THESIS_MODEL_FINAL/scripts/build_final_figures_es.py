"""Generate Spanish copies of the accepted final thesis figures.

No optimization or parameter change is performed. Pointwise arrays are read or
reproduced only through the same frozen forward-only helpers used to create the
accepted English figures. English source assets are hashed before and after the
build and must remain byte-identical.
"""

from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import build_final_thesis_package as source  # noqa: E402


DATA = ROOT / "thesis_data"
OUTPUT = ROOT / "figuras_finales_es"
TESIS = OUTPUT / "tesis"
PRESENTACION = OUTPUT / "presentacion"

BLUE = source.BLUE
ORANGE = source.ORANGE
GREEN = source.GREEN
RED = source.RED
PURPLE = source.PURPLE
GRAY = source.GRAY
INK = source.INK

FIGURES = [
    (1, "01_ajuste_ocp_panasonic", "panasonic_ocp_fit_final", "Ajuste del voltaje de celda a C/20"),
    (2, "02_calibracion_hppc_transferencia_bajo_soc", "panasonic_hppc_transfer_final", "Calibración con HPPC y validación a bajo SOC"),
    (3, "03_validacion_electrica_1c_panasonic", "panasonic_1c_voltage_validation_final", "Validación eléctrica en descarga 1C"),
    (4, "04_validacion_electrica_udds_panasonic", "panasonic_udds_voltage_validation_final", "Validación eléctrica con perfil UDDS"),
    (5, "05_transferencia_electrica_carga_cccv", "panasonic_cccv_voltage_transfer_final", "Validación eléctrica durante carga CC-CV"),
    (6, "06_separacion_direccional_carga_descarga", "panasonic_directional_voltage_separation_final", "Voltaje a baja tasa durante carga y descarga"),
    (7, "07_calibracion_termica_1c", "panasonic_1c_thermal_calibration_final", "Calibración térmica de un nodo en descarga 1C"),
    (8, "08_transferencia_termica_udds", "panasonic_udds_thermal_transfer_final", "Transferencia térmica bajo perfil UDDS"),
    (9, "09_transferencia_termica_carga_cccv", "panasonic_cccv_thermal_transfer_final", "Transferencia electro-térmica durante carga CC-CV"),
    (10, "10_seleccion_modelo_termico", "thermal_model_selection_summary", "Selección de la arquitectura del modelo térmico"),
    (11, "11_respuesta_serie_aparente_uniandes", "uniandes_apparent_series_path_final", "Auditoría de respuesta aparente de trayectoria serie — Uniandes"),
    (12, "12_1_voltaje_contactos_defectuosos_uniandes", "uniandes_defective_contacts_discharge_final", "Voltaje en descarga con resistencia externa del montaje"),
    (12, "12_2_comparacion_montajes_contacto_uniandes", "uniandes_contact_arrangements_comparison_final", "Comparación del voltaje en descarga entre montajes"),
    (13, "13_transferencia_panasonic_a_trc", "panasonic_to_trc_transfer_final", "Transferencia del modelo Panasonic a la celda TRC18650"),
    (14, "14_flujo_evaluacion_dataset_bateria", "battery_dataset_assessment_flow_final", "Metodología para evaluación de nuevos conjuntos de datos"),
    (15, "15_contexto_geografico_intervalos_uniandes", "uniandes_geographical_context_final", "Localidades de referencia para las condiciones térmicas"),
]


def save_es(fig: plt.Figure, name: str) -> None:
    fig.set_size_inches(12.0, 6.75)
    fig.savefig(TESIS / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(TESIS / f"{name}.png", dpi=320, bbox_inches="tight")
    fig.savefig(PRESENTACION / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(PRESENTACION / f"{name}.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def metric_box(axis: plt.Axes, text: str, loc: tuple[float, float] = (0.02, 0.04)) -> None:
    axis.text(
        *loc,
        text,
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=10.5,
        color=INK,
        bbox={"boxstyle": "round,pad=0.42", "facecolor": "white", "edgecolor": "#C7D2E2", "alpha": 0.94},
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def original_paths() -> list[Path]:
    rows = source.csv_rows(DATA / "final_figure_manifest.csv")
    paths = []
    for row in rows:
        for key in ("thesis_pdf", "thesis_png", "presentation_pdf", "presentation_png"):
            paths.append(ROOT / row[key])
    return paths


def figure_01() -> None:
    from models.ocp_panasonic_qi_frozen import V_eq
    from scripts.verify_panasonic_c20_qi_frozen import _load_c20

    z, measured = _load_c20()
    predicted = np.asarray(V_eq(z))
    fit = (z >= 0.05) & (z <= 0.95)
    residual_mv = 1000.0 * (measured - predicted)
    fig, axes = plt.subplots(
        2,
        1,
        sharex=False,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [7.0, 3.0]},
    )

    voltage_axis, residual_axis = axes
    voltage_axis.plot(z, measured, color=INK, linewidth=2.3, label="Datos experimentales")
    voltage_axis.plot(z, predicted, color=BLUE, linewidth=2.2, linestyle="--", label="Modelo")
    for boundary in (0.05, 0.95):
        voltage_axis.axvline(boundary, color=GRAY, linewidth=1.0, linestyle=(0, (2, 3)), alpha=0.75)
    voltage_axis.set_xlim(0.0, 1.0)
    voltage_axis.set_xlabel("Capacidad descargada normalizada, z")
    voltage_axis.set_ylabel("Voltaje de celda (V)")
    voltage_axis.set_title("(a)", fontsize=12, loc="left", pad=5)
    voltage_axis.legend(loc="upper right", frameon=True)
    voltage_axis.text(
        0.025,
        0.075,
        "RMSE (5–95%) = 1.086 mV",
        transform=voltage_axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=10.2,
        color=INK,
        bbox={"boxstyle": "round,pad=0.32", "facecolor": "white", "edgecolor": "#C7D2E2", "alpha": 0.94},
    )

    residual_axis.plot(z[fit], residual_mv[fit], color=ORANGE, linewidth=1.8)
    residual_axis.axhline(0.0, color=INK, linewidth=0.9)
    residual_axis.set_xlim(0.05, 0.95)
    residual_axis.set_xlabel("Capacidad descargada normalizada, z")
    residual_axis.set_ylabel("Residuo (mV)")
    residual_axis.set_title("(b)  Residuo en la región de ajuste (5–95%)", fontsize=11.5, loc="left", pad=5)

    fig.patch.set_facecolor("white")
    fig.suptitle(FIGURES[0][3], fontsize=18, weight="bold", color=INK)
    save_es(fig, FIGURES[0][1])


def figure_02() -> None:
    from scripts.audit_panasonic_hppc_identifiability import load_hppc_data, load_pulse_records
    from scripts.calibrate_panasonic_hppc_primary import evaluate, prepare_groups

    parameters = {"R0": 0.0259594900586, "alpha_p": 3204.22765511, "d_n": 0.000183029692927}
    data = load_hppc_data()
    records = load_pulse_records()
    high_groups = prepare_groups(data, records, minimum_pulse_start_soc=1.0 / 3.0)
    low_groups = prepare_groups(data, records, maximum_pulse_start_soc=1.0 / 3.0)
    high_eval = evaluate(high_groups, parameters)
    low_eval = evaluate(low_groups, parameters)

    def select(groups, evaluation, target_soc):
        candidates = []
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
        return group.time_s[indices] - group.time_s[start], group.measured_voltage_v[indices], model[indices], pulse.record.soc

    high_t, high_m, high_p, high_soc = select(high_groups, high_eval, 0.55)
    low_t, low_m, low_p, low_soc = select(low_groups, low_eval, 0.20)
    fig, axes = plt.subplots(1, 3, constrained_layout=True, gridspec_kw={"width_ratios": [1.15, 1.15, 0.9]})
    for ax, time, measured, predicted, soc, title in (
        (axes[0], high_t, high_m, high_p, high_soc, "(a) Pulso 2C utilizado en la\ncalibración"),
        (axes[1], low_t, low_m, low_p, low_soc, "(b) Pulso 2C reservado para\nvalidación"),
    ):
        ax.plot(time, measured, color=INK, label="Datos experimentales")
        ax.plot(time, predicted, color=BLUE, linestyle="--", label="Modelo")
        ax.axvline(0.0, color=GRAY, lw=1, linestyle=":")
        ax.set(xlabel="Tiempo desde el inicio del pulso (s)", ylabel="Voltaje (V)", title=f"{title} — SOC = {soc:.2f}")
    axes[0].legend(loc="lower right", fontsize=8.8)
    labels = [
        "Respuesta al inicio\ndel pulso",
        "Durante el\npulso",
        "Relajación",
        "RMSE\nbalanceado",
    ]
    high = [43.115, 44.262, 15.890, 36.835197639]
    low = [22.531, 145.725, 17.466, 85.729432921]
    xpos = np.arange(len(labels))
    width = 0.36
    axes[2].bar(xpos - width / 2, high, width, color=BLUE, label="Calibración: SOC ≥ 1/3")
    axes[2].bar(xpos + width / 2, low, width, color=ORANGE, label="Validación: SOC < 1/3")
    axes[2].set_xticks(xpos, labels, fontsize=8.2)
    axes[2].tick_params(axis="x", pad=6)
    axes[2].set(ylabel="RMSE (mV)", title="(c) Error según la fase del pulso")
    axes[2].legend(fontsize=7.8, loc="upper right")
    axes[2].text(
        0.04,
        0.59,
        "24 pulsos reservados para\nvalidación, sin reajuste",
        transform=axes[2].transAxes,
        ha="left",
        va="top",
        fontsize=8.2,
        color=INK,
    )
    fig.suptitle(FIGURES[1][3], fontsize=18, weight="bold", color=INK)
    save_es(fig, FIGURES[1][1])


def figure_03(data: dict[str, np.ndarray]) -> None:
    """Export the shared authoritative Spanish 1C validation figure."""

    fig = source.make_1c_validation_figure(data)
    save_es(fig, FIGURES[2][1])


def figure_04(data: dict[str, np.ndarray]) -> None:
    fig, _ = source.make_udds_validation_figure(data)
    save_es(fig, FIGURES[3][1])


def figure_05(data: dict[str, np.ndarray]) -> None:
    fig, _ = source.make_cccv_validation_figure(data)
    save_es(fig, FIGURES[4][1])


def figure_06() -> None:
    fig = source.make_c20_charge_discharge_figure()
    save_es(fig, FIGURES[5][1])


def figure_07(data: dict[str, np.ndarray]) -> None:
    minutes = data["time"] / 60.0
    residual = data["predicted"] - data["measured"]
    fig, axes = plt.subplots(2, 1, sharex=True, constrained_layout=True, gridspec_kw={"height_ratios": [2.2, 1.0]})
    axes[0].plot(minutes, data["measured"], color=INK, label="Temperatura superficial medida")
    axes[0].plot(minutes, data["predicted"], color=RED, linestyle="--", label="Predicción de un nodo")
    axes[0].set(ylabel="Temperatura (°C)", title="Calibración térmica primaria 1C, SOC 0.90–0.10")
    axes[0].legend()
    metric_box(axes[0], "Cth = 117.523 J/K\nGth = 0.109112 W/K\nRMSE 0.143622 K; R² 0.964702")
    axes[1].plot(minutes, residual, color=ORANGE)
    axes[1].axhline(0.0, color=INK, lw=1)
    axes[1].set(xlabel="Tiempo transcurrido (min)", ylabel="Predicho − medido (K)", title="Residuo de temperatura")
    fig.suptitle(FIGURES[6][3], fontsize=18, weight="bold", color=INK)
    save_es(fig, FIGURES[6][1])


def figure_08(data: dict[str, np.ndarray]) -> None:
    fig = source.make_thermal_udds_figure(data)
    save_es(fig, FIGURES[7][1])


def figure_09(data: dict[str, np.ndarray]) -> None:
    fig = source.make_thermal_cccv_figure(data)
    save_es(fig, FIGURES[8][1])


def figure_10() -> None:
    categories = ["Ajuste 1C", "Transferencia UDDS"]
    one_node = [0.143622065, 0.526697652]
    two_node = [0.138221456, 0.532525609]
    positions = np.arange(2)
    width = 0.34
    fig, (axis, note) = plt.subplots(1, 2, constrained_layout=True, gridspec_kw={"width_ratios": [1.25, 1.0]})
    axis.bar(positions - width / 2, one_node, width, color=BLUE, label="Un nodo")
    axis.bar(positions + width / 2, two_node, width, color=ORANGE, label="Dos nodos limpio")
    axis.set_xticks(positions, categories)
    axis.set(ylabel="RMSE de temperatura (K)", title="Comparación de precisión")
    axis.legend()
    for x_values, values in zip((positions - width / 2, positions + width / 2), (one_node, two_node), strict=True):
        for x_value, value in zip(x_values, values, strict=True):
            axis.text(x_value, value + 0.012, f"{value:.3f}", ha="center", fontsize=11)
    note.axis("off")
    note.text(0.02, 0.92, "Auditoría de recuperabilidad de dos nodos", fontsize=15, weight="semibold", color=INK, va="top")
    note.text(
        0.02,
        0.78,
        "• la capacitancia superficial alcanzó el límite inferior\n"
        "• modo rápido ≈ 1.00×10⁻⁴ s\n"
        "• condición del Jacobiano local ≈ 2.56×10⁷\n"
        "• sin mejora en la validación UDDS independiente\n\n"
        "Decisión: conservar el modelo parsimonioso de un nodo.\n"
        "Los parámetros de dos nodos no son recuperables con el\n"
        "experimento disponible de temperatura superficial.",
        fontsize=12.5,
        color=INK,
        va="top",
        linespacing=1.5,
    )
    fig.suptitle(FIGURES[9][3], fontsize=18, weight="bold", color=INK)
    save_es(fig, FIGURES[9][1])


def figure_11() -> None:
    rows = source.csv_rows(DATA / "own_data_apparent_resistance_transitions.csv")
    order = ["cycle_cell1_ambient", "cycle_cell2_undeclared", "cycle_cell2_30c", "cycle_cell1_35c", "cycle_cell1_30c", "pulse_cell1_part1", "pulse_cell1_part2"]
    values = [[float(row["r_app_ohm"]) for row in rows if row["dataset_id"] == name] for name in order]
    labels = ["Ambiente\ncelda 1", "Celda 2\nsin declarar", "Celda 2\n30 °C", "Celda 1\n35 °C", "Celda 1\n30 °C", "Pulsos\nparte 1", "Pulsos\nparte 2"]
    fig, axis = plt.subplots(constrained_layout=True)
    boxes = axis.boxplot(values, patch_artist=True, tick_labels=labels, widths=0.62, showmeans=True)
    for index, box in enumerate(boxes["boxes"]):
        box.set_facecolor(GREEN if index == 0 else ORANGE)
        box.set_alpha(0.75)
    for index, series in enumerate(values, start=1):
        jitter = np.linspace(-0.08, 0.08, len(series)) if len(series) > 1 else np.zeros(1)
        axis.scatter(index + jitter, series, color=INK, s=25, zorder=3)
    axis.axhline(0.0358755369, color=GREEN, linestyle="--", label="Mediana de referencia limpia 0.03588 Ω")
    axis.set(ylabel=r"Resistencia aparente en transición, $R_{app}$ (Ω)", title="Diagnóstico adyacente entre reposo y carga")
    axis.legend(loc="upper left")
    metric_box(axis, "Medianas afectadas: 0.178–0.232 Ω\n≈5–6.5× la referencia limpia\nRapp no es la resistencia interna ni R0", (0.57, 0.05))
    fig.suptitle(FIGURES[10][3], fontsize=18, weight="bold", color=INK)
    save_es(fig, FIGURES[10][1])


def figure_12(experiments: dict[str, object]) -> None:
    defective_fig, contact_fig, _ = source.make_uniandes_discharge_figures(experiments)
    save_es(defective_fig, FIGURES[11][1])
    save_es(contact_fig, FIGURES[12][1])


def figure_13(soc: np.ndarray, curves: dict[str, np.ndarray], experiments: dict[str, object]) -> None:
    from analysis.own_case_study_reference import discharge_curve, segment_complete_cycle
    from models.ocp_panasonic_qi_frozen import V_eq
    from scripts.audit_panasonic_to_trc18650_transfer import forward_frozen_panasonic

    experiment = experiments["cycle_cell1_ambient"]
    segments = segment_complete_cycle(experiment)
    discharge_segment = next(item for item in segments if item.name == "main_1c_discharge")
    discharge = discharge_curve(experiment, discharge_segment)
    rows = np.arange(discharge_segment.start, discharge_segment.end + 1)
    time = (rows - rows[0]) * experiment.estimated_dt_s
    model_current = -experiment.raw_current_a[rows]
    temperature = experiment.numeric["Temp1"][rows] + 273.15
    forward = forward_frozen_panasonic(time, model_current, temperature)
    measured_raw = np.asarray(discharge["voltage_v"])
    soc_raw = np.asarray(discharge["soc"])
    backbone = np.asarray(V_eq(1.0 - soc))
    fig, axes = plt.subplots(1, 2, constrained_layout=True)
    axes[0].plot(soc, curves["cycle_cell1_ambient"], color=INK, label="TRC medido en descarga 1C")
    axes[0].plot(soc, backbone, color=ORANGE, linestyle="--", label="Base de voltaje de equilibrio Panasonic")
    axes[0].plot(soc_raw[forward.ocp_valid], forward.voltage_v[forward.ocp_valid], color=BLUE, alpha=0.9, label="Transferencia completa literal congelada")
    axes[0].invert_xaxis()
    axes[0].set(xlabel="SOC de referencia TRC", ylabel="Voltaje (V)")
    axes[0].set_title("Transferencia directa\nPanasonic a TRC", fontsize=12.5, pad=6)
    axes[0].legend(fontsize=9)
    metric_box(axes[0], "RMSE de la base 188.882 mV\nRMSE de transferencia completa 82.134 mV")
    residual_backbone = 1000.0 * (backbone - curves["cycle_cell1_ambient"])
    axes[1].plot(soc, residual_backbone, color=ORANGE, label="Base − medido")
    axes[1].plot(soc_raw[forward.ocp_valid], 1000.0 * (forward.voltage_v[forward.ocp_valid] - measured_raw[forward.ocp_valid]), color=BLUE, label="Transferencia completa − medido")
    axes[1].axhline(0.0, color=INK, lw=1)
    axes[1].invert_xaxis()
    axes[1].set(xlabel="SOC de referencia TRC", ylabel="Residuo (mV)")
    axes[1].set_title("Desajuste dependiente del SOC\ny pérdida de soporte a bajo SOC", fontsize=12.5, pad=6)
    axes[1].legend(fontsize=9)
    metric_box(axes[1], "DIAGNÓSTICO CON TIEMPO PROVISIONAL\n22 violaciones del soporte superficial positivo\nSin afirmación de precisión para TRC", (0.04, 0.62))
    fig.suptitle(FIGURES[13][3] + " — TRANSFERENCIA DIRECTA INCOMPATIBLE", fontsize=16, weight="bold", color=INK)
    save_es(fig, FIGURES[13][1])


def figure_14() -> None:
    fig, axis = plt.subplots(constrained_layout=True)
    axis.set_xlim(0, 15.8)
    axis.set_ylim(0, 6.7)
    axis.axis("off")
    nodes = {
        "new": (0.25, 5.2, "NUEVOS DATOS\nEXPERIMENTALES"),
        "audit": (3.25, 5.2, "AUDITORÍA DE\nCALIDAD DE DATOS"),
        "class": (6.25, 5.2, "CLASIFICACIÓN DE\nIDONEIDAD"),
        "diag": (6.25, 2.4, "CASO DE ESTUDIO\nDIAGNÓSTICO"),
        "transfer": (9.25, 5.2, "PRUEBA DE\nTRANSFERENCIA\nDEL MODELO"),
        "compat": (12.25, 5.2, "AUDITORÍA DE\nCOMPATIBILIDAD\nDE PARÁMETROS"),
        "limit": (9.75, 0.85, "REPORTAR\nLIMITACIÓN"),
        "valid": (12.25, 2.4, "VALIDACIÓN\nCUANTITATIVA"),
    }
    box_width, box_height = 2.55, 1.02
    for key, (x, y, label) in nodes.items():
        color = GREEN if key == "valid" else (ORANGE if key in {"diag", "limit"} else BLUE)
        axis.add_patch(plt.Rectangle((x, y), box_width, box_height, facecolor=color, edgecolor="white", linewidth=2, alpha=0.92))
        axis.text(x + box_width / 2, y + box_height / 2, label, ha="center", va="center", color="white", fontsize=8.1, weight="bold", linespacing=1.05)

    def arrow(origin, destination, label=""):
        sx, sy, _ = nodes[origin]
        tx, ty, _ = nodes[destination]
        start = (sx + box_width, sy + box_height / 2) if tx > sx else (sx + box_width / 2, sy)
        end = (tx, ty + box_height / 2) if tx > sx else (tx + box_width / 2, ty + box_height)
        axis.annotate("", xy=end, xytext=start, arrowprops={"arrowstyle": "->", "lw": 2, "color": INK})
        if label:
            vertical_offset = 0.55 if abs(start[1] - end[1]) < 0.05 else 0.17
            axis.text((start[0] + end[0]) / 2, (start[1] + end[1]) / 2 + vertical_offset, label, ha="center", fontsize=9.5, color=INK)

    arrow("new", "audit")
    arrow("audit", "class")
    arrow("class", "transfer", "evidencia idónea")
    arrow("transfer", "compat")
    arrow("class", "diag", "no ideal / insuficiente")
    arrow("diag", "limit")
    arrow("compat", "valid", "requisitos sustentados")
    arrow("compat", "limit", "adaptación no sustentada")
    axis.text(0.75, 0.38, "Sin ajuste automático • sin corrección automática de voltaje • la validación es condicional, no el destino predeterminado", color=INK, fontsize=11.5, weight="semibold")
    fig.suptitle(FIGURES[14][3], fontsize=18, weight="bold", color=INK)
    save_es(fig, FIGURES[14][1])


def figure_15() -> None:
    fig = source.make_uniandes_geographical_figure()
    save_es(fig, FIGURES[15][1])


def documentation_rows() -> list[dict[str, str]]:
    return [
        {"numero": "01", "muestra": "Curva C/20 medida, representación OCP congelada, residuo 5–95% y trayectorias nominales de electrodos.", "resultado": "RMSE 5–95% = 1.086 mV.", "interpretacion": "La representación de voltaje fue aceptada, pero las ventanas estequiométricas son débilmente identificables.", "uso": "PRESENTACION + TESIS"},
        {"numero": "02", "muestra": "Pulsos HPPC representativos y comparación de errores entre calibración y bajo SOC reservado.", "resultado": "36.835 mV en calibración; 85.729 mV a bajo SOC.", "interpretacion": "La transferencia a bajo SOC es limitada y el mayor deterioro está en el transitorio bajo carga.", "uso": "PRESENTACION + TESIS"},
        {"numero": "03", "muestra": "Validación directa de voltaje durante descarga continua 1C y residuo temporal.", "resultado": "RMSE global = 114.362 mV.", "interpretacion": "El error crece a bajo SOC y cerca del corte.", "uso": "TESIS"},
        {"numero": "04", "muestra": "Voltaje medido y predicho durante el perfil UDDS completo, con residuo.", "resultado": "RMSE global = 30.988 mV.", "interpretacion": "Buen resultado global con deterioro estructurado a bajo SOC y en cambios grandes de corriente.", "uso": "PRESENTACION + TESIS"},
        {"numero": "05", "muestra": "Transferencia eléctrica del modelo derivado de descarga a una carga CC-CV independiente.", "resultado": "RMSE = 61.034 mV; sesgo = −47.526 mV.", "interpretacion": "La transferencia entre direcciones de corriente es limitada.", "uso": "BACKUP"},
        {"numero": "06", "muestra": "Comparación direccional de carga y descarga C/20 sobre un eje común de capacidad.", "resultado": "Carga hasta 87.277% de Qcell.", "interpretacion": "Es separación direccional de voltaje con cautela de mapeo; no se interpreta como histéresis termodinámica.", "uso": "TESIS"},
        {"numero": "07", "muestra": "Temperatura superficial medida y predicha en la calibración térmica primaria 1C.", "resultado": "RMSE = 0.143622 K; R² = 0.964702.", "interpretacion": "El modelo efectivo de un nodo reproduce bien el ensayo de calibración.", "uso": "PRESENTACION + TESIS"},
        {"numero": "08", "muestra": "Transferencia térmica sin reajuste al perfil UDDS.", "resultado": "RMSE = 0.526698 K; sesgo = −0.512715 K.", "interpretacion": "Existe subpredicción térmica sistemática fuera de la calibración.", "uso": "TESIS"},
        {"numero": "09", "muestra": "Transferencia térmica en carga para la ruta de calor con voltaje medido y la predicción electro-térmica completa.", "resultado": "Validación A = 0.427781 K; B = 0.864391 K.", "interpretacion": "El déficit de calor por error eléctrico explica la mayor parte de la degradación acoplada.", "uso": "PRESENTACION + TESIS"},
        {"numero": "10", "muestra": "Comparación de un nodo y dos nodos en 1C y UDDS, con auditoría de recuperabilidad.", "resultado": "Dos nodos: 0.138221 K en 1C y 0.532526 K en UDDS; condición ≈2.56e7.", "interpretacion": "El segundo nodo no es recuperable y no mejora la transferencia independiente.", "uso": "BACKUP"},
        {"numero": "11", "muestra": "Distribución de resistencia aparente en transición para la referencia limpia y los ensayos afectados.", "resultado": "Referencia ≈0.03588 Ω; afectados ≈0.178–0.232 Ω.", "interpretacion": "Hay evidencia reproducible de respuesta aparente de trayectoria serie; Rapp no es resistencia interna ni R0.", "uso": "PRESENTACION + TESIS"},
        {"numero": "12", "muestra": "Dos comparaciones de descarga: tres condiciones con contactos defectuosos y dos montajes de contacto bajo Aire 22–24 °C.", "resultado": "Mapeo eléctrico corregido y trazable por archivo; sin panel de diferencia de voltaje.", "interpretacion": "Las ciudades son contexto climático y las diferencias no se atribuyen aisladamente a la temperatura ni a la resistencia de contacto.", "uso": "TESIS"},
        {"numero": "13", "muestra": "Transferencia literal de la base de voltaje y del modelo Panasonic completo a la celda TRC18650.", "resultado": "RMSE de base = 188.882 mV; transferencia completa = 82.134 mV.", "interpretacion": "La transferencia directa es incompatible y mezcla desajustes estáticos, de mapeo y dinámicos.", "uso": "PRESENTACION + TESIS"},
        {"numero": "14", "muestra": "Flujo reutilizable desde auditoría de datos hasta validación o reporte de limitación.", "resultado": "Sin ajuste ni corrección automática.", "interpretacion": "La validación cuantitativa es una salida condicionada por la evidencia disponible.", "uso": "PRESENTACION + TESIS"},
        {"numero": "15", "muestra": "Mapa de Colombia con Rionegro, Cali y Valledupar asociados a los intervalos experimentales medidos.", "resultado": "Contexto geográfico, sin atribución de rangos climáticos permanentes.", "interpretacion": "Las ciudades contextualizan las condiciones; los ensayos no reproducen ciclos meteorológicos completos.", "uso": "TESIS"},
    ]


def write_index_and_manifest() -> None:
    docs = {row["numero"]: row for row in documentation_rows()}
    lines = [
        "# Índice de figuras finales en español",
        "",
        "Las figuras son copias localizadas de los resultados aceptados. No se realizó ajuste, recalibración ni modificación de parámetros.",
        "",
        "| N.º | Nombre en español | Original en inglés | Título | Qué muestra | Resultado principal | Interpretación científica | Uso recomendado |",
        "|---:|---|---|---|---|---|---|---|",
    ]
    manifest_rows = []
    source_manifest = {row["stem"]: row for row in source.csv_rows(DATA / "final_figure_manifest.csv")}
    for number, name_es, original, title in FIGURES:
        item = docs[f"{number:02d}"]
        lines.append(f"| {number:02d} | `{name_es}` | `{original}` | {title} | {item['muestra']} | {item['resultado']} | {item['interpretacion']} | {item['uso']} |")
        manifest_rows.append(
            {
                "numero": f"{number:02d}",
                "nombre_es": name_es,
                "nombre_original": original,
                "titulo_es": title,
                "archivo_tesis_pdf": f"figuras_finales_es/tesis/{name_es}.pdf",
                "archivo_tesis_png": f"figuras_finales_es/tesis/{name_es}.png",
                "archivo_presentacion_png": f"figuras_finales_es/presentacion/{name_es}.png",
                "resultado_principal": item["resultado"],
                "uso_recomendado": item["uso"],
                "fuente_original": source_manifest[original]["thesis_png"],
            }
        )
    (OUTPUT / "INDICE_FIGURAS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    source.write_csv(OUTPUT / "manifest_figuras_es.csv", manifest_rows)


def write_shortlist() -> None:
    text = """# Figuras principales para la presentación de 20 minutos

| Diapositiva / tema | Figura | Mensaje en una frase | Número para decir oralmente | No explicar salvo que pregunten |
|---|---|---|---|---|
| 3. Compatibilidad OCP | `presentacion/01_ajuste_ocp_panasonic.png` | La curva C/20 queda representada con alta precisión, aunque las ventanas de electrodo no son únicas. | RMSE 5–95% = **1.086 mV**. | La historia completa de multistarts y las alternativas OCP rechazadas. |
| 4. Calibración eléctrica condicional | `presentacion/02_calibracion_hppc_transferencia_bajo_soc.png` | El ajuste HPPC funciona en su región declarada, pero no transfiere igual a bajo SOC. | **36.835 mV** en calibración frente a **85.729 mV** a bajo SOC. | El detalle de cada pulso, el Jacobiano y todas las correlaciones de parámetros. |
| 5. Validación dinámica independiente | `presentacion/04_validacion_electrica_udds_panasonic.png` | El desempeño global en UDDS es fuerte, con deterioro localizado a bajo SOC y en cambios grandes de corriente. | RMSE global = **30.988 mV**. | Todas las métricas por cuartil de corriente o los diez eventos individuales. |
| 7. Calibración térmica primaria | `presentacion/07_calibracion_termica_1c.png` | El modelo de un nodo reproduce con precisión el ensayo térmico usado para calibración. | RMSE = **0.143622 K**, R² = **0.964702**. | El detalle de multistarts y la trayectoria secundaria con calor calculado desde voltaje de modelo. |
| 8. Transferencia electro-térmica | `presentacion/09_transferencia_termica_carga_cccv.png` | El error eléctrico reduce el calor predicho y domina la falla de la predicción acoplada en carga. | A = **0.427781 K**, B = **0.864391 K**; déficit = **436.438 J**. | Todos los balances por fase CC/CV y las fracciones de calor reversible. |
| 10. Lección experimental Uniandes | `presentacion/11_respuesta_serie_aparente_uniandes.png` | La auditoría separa una referencia limpia de ensayos con respuesta serie aparente elevada antes de intentar identificar parámetros. | **0.03588 Ω** frente a **0.178–0.232 Ω**. | Una atribución física única; Rapp no es resistencia interna ni R0. |
| 11. Transferencia entre celdas | `presentacion/13_transferencia_panasonic_a_trc.png` | El modelo Panasonic no puede trasladarse literalmente a la TRC18650 sin adaptación de estado y referencia de voltaje. | **188.882 mV** de base y **82.134 mV** de transferencia completa. | La integración con tiempo provisional y la descomposición completa del error, salvo pregunta. |

La figura `14_flujo_evaluacion_dataset_bateria.png` queda como alternativa para la conclusión o como respaldo metodológico si se dispone de una diapositiva adicional.
"""
    (OUTPUT / "FIGURAS_PRESENTACION.md").write_text(text, encoding="utf-8")


def write_objectives() -> None:
    text = """# Objetivos aprobados y evidencia producida

## OBJETIVO GENERAL

> Estudiar la implementación, calibración y validación de un modelo electroquímico–térmico reducido para una batería de ion-litio tipo 18650, considerando la identificación de parámetros a partir de bases de datos reportadas y mediciones experimentales propias, con el propósito de representar la respuesta eléctrica y térmica de la celda bajo diferentes condiciones de operación.

### Evidencia producida

- Implementación documentada del modelo agrupado SPM, generación de calor y dinámica térmica de un nodo.
- Calibración eléctrica y térmica en la base Panasonic NCR18650PF.
- Validaciones independientes 1C, UDDS y CC-CV, con limitaciones declaradas.
- Auditoría experimental y transferencia al caso TRC18650 de Uniandes.

## OBJETIVOS ESPECÍFICOS

### 1.

> Implementar un modelo electroquímico–térmico reducido, basado en el Single Particle Model (SPM) o en una variante afín, para representar variables relevantes de una celda 18650, tales como voltaje terminal, estado de carga, generación de calor y temperatura.

#### Evidencia producida

- Ecuaciones implementadas y auditadas en `thesis_data/final_model_equations.md`.
- Figuras 01, 07 y 10: representación OCP, respuesta térmica y selección de arquitectura.
- Salidas explícitas de voltaje, SOC, estados de electrodo, calor y temperatura sin recorte silencioso.

### 2.

> Formular una metodología de selección y calibración de parámetros del modelo, considerando parámetros eléctricos, electroquímicos y térmicos obtenidos de la literatura, especificaciones del fabricante, bases de datos reportadas o mediciones propias.

#### Evidencia producida

- Figura 14: flujo reproducible de auditoría, clasificación, transferencia y validación.
- Figuras 01, 02 y 07: identificación OCP, calibración HPPC condicional y calibración térmica.
- Registros trazables de procedencia, dominio de soporte, identificabilidad y rol de cada conjunto de datos.

### 3.

> Contrastar el desempeño del modelo implementado frente a bases de datos experimentales reportadas y mediciones propias bajo condiciones controladas de operación, utilizando perfiles de carga, descarga, pulsos de corriente o perfiles dinámicos disponibles.

#### Evidencia producida

- Figuras 02–09: contraste con pulsos HPPC, descarga 1C, UDDS y carga CC-CV.
- Figuras 11 y 12: mediciones propias, referencia experimental y ensayos afectados.
- Métricas RMSE, MAE, sesgo y R² reportadas según el rol de evidencia de cada ensayo.

### 4.

> Analizar la transferencia de los parámetros calibrados entre diferentes conjuntos de datos o celdas 18650, evaluando cuándo los parámetros identificados son suficientes para representar la respuesta del sistema y cuándo es necesario considerar un reajuste o adaptación parcial de los mismos.

#### Evidencia producida

- Figuras 02–05, 08 y 09: transferencia entre regiones de SOC, perfiles y dirección de corriente en Panasonic.
- Figura 13: transferencia directa Panasonic a TRC18650 y diagnóstico de incompatibilidad.
- Decisión jerárquica de adaptación: mapeo de capacidad/estado, referencia de voltaje compatible y evidencia dinámica adicional cuando exista base de tiempo medida.
"""
    (OUTPUT / "objetivos_y_evidencia.md").write_text(text, encoding="utf-8")


def master_table() -> None:
    import textwrap

    rows = [
        ["C/20 OCP", "Ajuste OCP", "RMSE 5–95%", "1.086 mV", "Aceptado; ventanas débilmente identificables"],
        ["HPPC", "Calibración / validación reservada", "RMSE balanceado por fase", "36.835 / 85.729 mV", "Transferencia limitada a bajo SOC"],
        ["Descarga 1C", "Validación independiente", "RMSE de voltaje", "114.362 mV", "Limitada a bajo SOC y corte"],
        ["UDDS", "Validación independiente", "RMSE de voltaje", "30.988 mV", "Limitada a bajo SOC y cambios grandes"],
        ["Carga CC-CV", "Validación direccional", "RMSE / sesgo", "61.034 / −47.526 mV", "Transferencia direccional limitada"],
        ["Térmica 1C", "Calibración primaria", "RMSE de temperatura", "0.143622 K", "Ajuste primario de un nodo"],
        ["Térmica UDDS", "Validación A", "RMSE / sesgo", "0.526698 / −0.512715 K", "Subpredicción sistemática"],
        ["Térmica CC-CV", "Validación A / B", "RMSE de temperatura", "0.427781 / 0.864391 K", "Falla acoplada dominada por error eléctrico"],
        ["Panasonic → TRC", "Diagnóstico de transferencia", "RMSE base / completo", "188.882 / 82.134 mV", "Transferencia directa incompatible"],
    ]
    csv_rows = [
        {"Ensayo": row[0], "Rol de la evidencia": row[1], "Métrica principal": row[2], "Resultado": row[3], "Interpretación / limitación": row[4]}
        for row in rows
    ]
    source.write_csv(OUTPUT / "tabla_maestra_resultados.csv", csv_rows)
    wrapped_rows = [[textwrap.fill(value, width=width) for value, width in zip(row, [16, 25, 22, 21, 35], strict=True)] for row in rows]
    wrapped_headers = [textwrap.fill(value, width=width) for value, width in zip(
        ["Ensayo", "Rol de la evidencia", "Métrica principal", "Resultado", "Interpretación / limitación"],
        [16, 22, 20, 18, 28],
        strict=True,
    )]
    fig, axis = plt.subplots(figsize=(18, 10))
    axis.axis("off")
    table = axis.table(
        cellText=wrapped_rows,
        colLabels=wrapped_headers,
        colWidths=[0.13, 0.21, 0.18, 0.18, 0.30],
        cellLoc="left",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10.8)
    table.scale(0.98, 2.75)
    for (row, _), cell in table.get_celld().items():
        cell.set_edgecolor("white")
        if row == 0:
            cell.set_facecolor(BLUE)
            cell.get_text().set_color("white")
            cell.get_text().set_weight("bold")
        else:
            cell.set_facecolor("#F2F6FA" if row % 2 else "white")
            cell.get_text().set_color(INK)
    axis.set_title("Tabla maestra de resultados aceptados", fontsize=22, weight="bold", color=INK, pad=20)
    fig.savefig(OUTPUT / "tabla_maestra_resultados.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def contact_sheet() -> None:
    fig, axes = plt.subplots(4, 4, figsize=(20, 14), constrained_layout=True)
    short_names = [
        "Ajuste OCP", "HPPC y bajo SOC", "Validación 1C", "Validación UDDS",
        "Carga CC-CV", "Separación C/20", "Térmica 1C", "Térmica UDDS",
        "Térmica CC-CV", "Selección térmica", "Respuesta serie", "Resistencia externa",
        "Montajes de contacto", "Panasonic a TRC", "Flujo metodológico",
        "Contexto geográfico",
    ]
    for axis, (number, name, _, _), short_name in zip(axes.flat, FIGURES, short_names, strict=False):
        image = plt.imread(PRESENTACION / f"{name}.png")
        axis.imshow(image)
        axis.set_title(f"{number:02d}. {short_name}", fontsize=13, weight="semibold", color=INK)
        axis.axis("off")
    for axis in axes.flat[len(FIGURES):]:
        axis.axis("off")
    fig.suptitle("Contacto rápido — figuras finales en español", fontsize=24, weight="bold", color=INK)
    fig.savefig(OUTPUT / "contact_sheet_figuras.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def quality_check(before_hashes: dict[Path, str]) -> None:
    after_hashes = {path: sha256(path) for path in before_hashes}
    changed = [path for path in before_hashes if before_hashes[path] != after_hashes[path]]
    if changed:
        raise RuntimeError(f"Se modificaron figuras originales: {changed}")

    manifest = source.csv_rows(OUTPUT / "manifest_figuras_es.csv")
    if len(manifest) != len(FIGURES):
        raise RuntimeError(
            f"El manifiesto en español no contiene {len(FIGURES)} figuras."
        )
    missing = []
    for row in manifest:
        for key in ("archivo_tesis_pdf", "archivo_tesis_png", "archivo_presentacion_png"):
            path = ROOT / row[key]
            if not path.is_file() or path.stat().st_size == 0:
                missing.append(path)
        presentation_pdf = PRESENTACION / f"{row['nombre_es']}.pdf"
        if not presentation_pdf.is_file() or presentation_pdf.stat().st_size == 0:
            missing.append(presentation_pdf)
    if missing:
        raise RuntimeError(f"Faltan salidas en español: {missing}")

    expected_names = {item[1] for item in FIGURES}
    actual_thesis = {path.stem for path in TESIS.glob("*.png")}
    actual_presentation = {path.stem for path in PRESENTACION.glob("*.png")}
    if actual_thesis != expected_names or actual_presentation != expected_names:
        raise RuntimeError("Los nombres numerados no coinciden con la especificación.")

    report = f"""# Control de calidad de figuras localizadas

- Figuras en el manifiesto: **{len(FIGURES)}**.
- PNG de tesis: **{len(list(TESIS.glob('*.png')))}**; PDF de tesis: **{len(list(TESIS.glob('*.pdf')))}**.
- PNG de presentación: **{len(list(PRESENTACION.glob('*.png')))}**; PDF de presentación: **{len(list(PRESENTACION.glob('*.pdf')))}**.
- Nombres ASCII numerados: **verificados**.
- Archivos referenciados por el manifiesto: **verificados, no vacíos**.
- Originales ingleses: **verificados por SHA-256 antes/después; sin cambios**.
- Optimización o calibración ejecutada: **NO**.
- Parámetros, ecuaciones, datasets o métricas aceptadas modificados: **NO**.
- Objetivos oficiales: **transcritos literalmente en `objetivos_y_evidencia.md`**.
- Terminología de `R_app`: **resistencia aparente en transición / respuesta aparente de trayectoria serie; nunca resistencia interna**.
"""
    (OUTPUT / "CONTROL_CALIDAD.md").write_text(report, encoding="utf-8")


def main() -> None:
    source.configure_style()
    TESIS.mkdir(parents=True, exist_ok=True)
    PRESENTACION.mkdir(parents=True, exist_ok=True)
    originals = original_paths()
    before_hashes = {path: sha256(path) for path in originals}

    figure_01()
    figure_02()
    electrical = source.electrical_forward_arrays()
    figure_03(electrical["one_c"])
    figure_04(electrical["udds"])
    figure_05(electrical["cccv"])
    figure_06()
    thermal = source.thermal_forward_arrays()
    figure_07(thermal["one_c"])
    figure_08(thermal["udds"])
    figure_09(thermal["cccv"])
    figure_10()
    figure_11()
    soc, curves, experiments = source.own_case_arrays()
    figure_12(experiments)
    figure_13(soc, curves, experiments)
    figure_14()
    figure_15()

    write_index_and_manifest()
    write_shortlist()
    write_objectives()
    master_table()
    contact_sheet()
    quality_check(before_hashes)
    print("Paquete final en español generado sin modificar los originales.")


if __name__ == "__main__":
    main()
