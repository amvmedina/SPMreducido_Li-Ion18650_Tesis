"""Build the final HPPC selection and pulse-phase figures.

This script is visualization-only. It reuses the accepted pulse records,
forward-model parameters, representative pulses, and reported metrics. It does
not fit, optimize, resample, or modify any model or experimental data.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts import build_final_thesis_package as package  # noqa: E402
from scripts.audit_panasonic_hppc_identifiability import (  # noqa: E402
    EXCLUDED_PULSE_IDS,
    PulseRecord,
    load_hppc_data,
    load_pulse_records,
    partition_contiguous_groups,
)
from scripts.calibrate_panasonic_hppc_primary import (  # noqa: E402
    Evaluation,
    PreparedGroup,
    PreparedPulse,
    evaluate,
    prepare_groups,
)


INK = package.INK
BLUE = package.BLUE
ORANGE = package.ORANGE
GREEN = package.GREEN
GRAY = package.GRAY
LIGHT = package.LIGHT

NOMINAL_PARAMETERS = {
    "R0": 0.0259594900586,
    "alpha_p": 3204.22765511,
    "d_n": 0.000183029692927,
}

CALIBRATION_PHASE_RMSE_MV = np.asarray(
    [43.114900568558746, 44.26189689154858, 15.889782497229506, 36.835197639]
)
VALIDATION_PHASE_RMSE_MV = np.asarray(
    [22.531227929720978, 145.72537146498618, 17.466164034563598, 85.729432921]
)


def retained(records: list[PulseRecord]) -> list[PulseRecord]:
    return [
        record
        for record in records
        if record.usable and record.pulse_id not in EXCLUDED_PULSE_IDS
    ]


def select_trace(
    groups: list[PreparedGroup],
    evaluation: Evaluation,
    target_soc: float,
    nominal_rate: str = "2C",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, PreparedPulse]:
    candidates: list[tuple[float, int, PreparedPulse]] = []
    for group_index, group in enumerate(groups):
        for pulse in group.pulses:
            if pulse.record.nominal_rate == nominal_rate:
                candidates.append(
                    (abs(pulse.record.soc - target_soc), group_index, pulse)
                )
    if not candidates:
        raise ValueError(f"No retained {nominal_rate} pulse is available.")

    _, group_index, pulse = min(candidates, key=lambda item: item[0])
    group = groups[group_index]
    model_voltage = evaluation.group_voltage_v[group_index]
    immediate = int(pulse.phase_local_indices["immediate"][0])
    relaxation_end = int(pulse.phase_local_indices["post_relaxation"][-1])
    indices = np.arange(pulse.pre_local_index, relaxation_end + 1)
    return (
        group.time_s[indices] - group.time_s[immediate],
        group.measured_voltage_v[indices],
        model_voltage[indices],
        pulse,
    )


def figure_hppc_selection() -> tuple[int, int, int, int]:
    data = load_hppc_data()
    records = load_pulse_records()
    usable = retained(records)
    calibration_records = [record for record in usable if record.soc >= 1.0 / 3.0]
    validation_records = [record for record in usable if record.soc < 1.0 / 3.0]
    if (len(calibration_records), len(validation_records)) != (40, 24):
        raise ValueError("The retained HPPC calibration/validation split changed.")

    calibration_groups = prepare_groups(
        data, records, minimum_pulse_start_soc=1.0 / 3.0
    )
    validation_groups = prepare_groups(
        data, records, maximum_pulse_start_soc=1.0 / 3.0
    )
    calibration_evaluation = evaluate(calibration_groups, NOMINAL_PARAMETERS)
    validation_evaluation = evaluate(validation_groups, NOMINAL_PARAMETERS)

    calibration_time, calibration_measured, calibration_model, calibration_pulse = select_trace(
        calibration_groups, calibration_evaluation, 0.55
    )
    validation_time, validation_measured, validation_model, validation_pulse = select_trace(
        validation_groups, validation_evaluation, 0.20
    )
    if (calibration_pulse.record.pulse_id, validation_pulse.record.pulse_id) != (33, 53):
        raise ValueError("The accepted representative HPPC pulse selection changed.")

    figure, axes = plt.subplots(2, 2, constrained_layout=True)
    distribution_axis, calibration_axis, validation_axis, error_axis = axes.flat

    experiment_start = data.time_s[usable[0].pre_pulse_index]
    calibration_end = calibration_records[-1].start_index
    validation_start = validation_records[0].start_index
    split_time_s = 0.5 * (
        data.time_s[calibration_end] + data.time_s[validation_start]
    )
    elapsed_hours = (data.time_s - experiment_start) / 3600.0
    split_hours = (split_time_s - experiment_start) / 3600.0
    experiment_end_hours = (
        data.time_s[usable[-1].post_rest_index] - experiment_start
    ) / 3600.0

    distribution_axis.axvspan(
        0.0, split_hours, color=BLUE, alpha=0.045, zorder=0
    )
    distribution_axis.axvspan(
        split_hours, experiment_end_hours, color=ORANGE, alpha=0.055, zorder=0
    )
    all_groups = partition_contiguous_groups(data, records)
    first_voltage_segment = True
    for group_records in all_groups:
        usable_group = [record for record in group_records if record in usable]
        if not usable_group:
            continue
        lower = usable_group[0].pre_pulse_index
        upper = usable_group[-1].post_rest_index + 1
        distribution_axis.plot(
            elapsed_hours[lower:upper],
            data.voltage_v[lower:upper],
            color=INK,
            linewidth=1.05,
            label="Voltaje medido" if first_voltage_segment else None,
            zorder=2,
        )
        first_voltage_segment = False
    distribution_axis.scatter(
        [elapsed_hours[record.start_index] for record in calibration_records],
        [data.voltage_v[record.start_index] for record in calibration_records],
        s=20,
        color=BLUE,
        marker="o",
        label="Inicio de pulso: calibración",
        zorder=3,
    )
    distribution_axis.scatter(
        [elapsed_hours[record.start_index] for record in validation_records],
        [data.voltage_v[record.start_index] for record in validation_records],
        s=23,
        color=ORANGE,
        marker="s",
        label="Inicio de pulso: validación",
        zorder=3,
    )
    distribution_axis.axvline(
        split_hours, color=GRAY, linewidth=1.0, linestyle=(0, (4, 3))
    )
    distribution_axis.text(
        0.025,
        0.08,
        "Pulsos usados en calibración\nSOC inicial ≥ 1/3  (40 pulsos)",
        transform=distribution_axis.transAxes,
        ha="left",
        va="bottom",
        color=BLUE,
        fontsize=8.8,
        weight="semibold",
    )
    distribution_axis.text(
        0.975,
        0.92,
        "Pulsos reservados para validación\nSOC inicial < 1/3  (24 pulsos)",
        transform=distribution_axis.transAxes,
        ha="right",
        va="top",
        color=ORANGE,
        fontsize=8.8,
        weight="semibold",
    )
    distribution_axis.set(
        xlim=(0.0, experiment_end_hours),
        xlabel="Tiempo transcurrido (h)",
        ylabel="Voltaje de celda (V)",
        title="(a) Ensayo HPPC completo y selección de pulsos",
    )
    distribution_axis.legend(loc="lower right", fontsize=7.8, ncol=1)

    for axis, time, measured, model, pulse, title in (
        (
            calibration_axis,
            calibration_time,
            calibration_measured,
            calibration_model,
            calibration_pulse,
            "(b) Pulso 2C utilizado en calibración",
        ),
        (
            validation_axis,
            validation_time,
            validation_measured,
            validation_model,
            validation_pulse,
            "(c) Pulso 2C reservado para validación",
        ),
    ):
        axis.plot(time, measured, color=INK, linewidth=2.1, label="Datos experimentales")
        axis.plot(time, model, color=BLUE, linewidth=2.0, linestyle="--", label="Modelo")
        axis.axvline(0.0, color=GRAY, linewidth=0.9, linestyle=":")
        axis.set(
            xlabel="Tiempo desde el inicio del pulso (s)",
            ylabel="Voltaje (V)",
            title=f"{title} — SOC = {pulse.record.soc:.2f}",
        )
        axis.legend(loc="lower right", fontsize=8.1)

    labels = [
        "Respuesta al inicio\ndel pulso",
        "Durante el\npulso",
        "Relajación",
        "RMSE\nbalanceado",
    ]
    positions = np.arange(len(labels))
    width = 0.36
    error_axis.bar(
        positions - width / 2,
        CALIBRATION_PHASE_RMSE_MV,
        width,
        color=BLUE,
        label="Calibración: SOC ≥ 1/3",
    )
    error_axis.bar(
        positions + width / 2,
        VALIDATION_PHASE_RMSE_MV,
        width,
        color=ORANGE,
        label="Validación: SOC < 1/3",
    )
    error_axis.set_xticks(positions, labels, fontsize=8.2)
    error_axis.tick_params(axis="x", pad=5)
    error_axis.set(
        ylabel="RMSE (mV)",
        title="(d) Comparación del error por fase",
    )
    error_axis.legend(loc="upper right", fontsize=8.0)

    figure.patch.set_facecolor("white")
    figure.suptitle(
        "Calibración con HPPC y validación a bajo SOC",
        fontsize=17,
        weight="bold",
        color=INK,
    )
    package.save_figure(figure, "hppc_calibracion_validacion_final")
    return (
        len(calibration_records),
        len(validation_records),
        calibration_pulse.record.pulse_id,
        validation_pulse.record.pulse_id,
    )


def figure_hppc_phases() -> tuple[int, str, float]:
    data = load_hppc_data()
    records = retained(load_pulse_records())
    candidates = [record for record in records if record.nominal_rate == "6C"]
    pulse = min(candidates, key=lambda record: abs(record.soc - 0.50))
    if pulse.pulse_id != 35:
        raise ValueError("The explanatory HPPC pulse selection changed.")

    pulse_start_time = data.time_s[pulse.start_index]
    # Show only a short pre-pulse interval: it is a voltage-reference level,
    # not a separate transient of interest in this explanatory figure.
    window_start = int(np.searchsorted(data.time_s, pulse_start_time - 8.0, side="left"))
    window_stop = pulse.post_rest_index + 1
    indices = np.arange(window_start, window_stop)
    time = data.time_s[indices] - pulse_start_time
    current = -data.raw_current_a[indices]
    voltage = data.voltage_v[indices]

    pulse_end_time = data.time_s[pulse.end_index] - pulse_start_time
    relaxation_start_time = data.time_s[pulse.end_index + 1] - pulse_start_time
    release_boundary = 0.5 * (pulse_end_time + relaxation_start_time)
    pre_pulse_voltage = data.voltage_v[pulse.pre_pulse_index]
    immediate_voltage = data.voltage_v[pulse.start_index]

    figure, (current_axis, voltage_axis) = plt.subplots(
        2,
        1,
        sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [1.0, 1.25]},
    )

    current_axis.plot(time, current, color=INK, linewidth=2.0)
    current_axis.axvspan(time[0], 0.0, color=GRAY, alpha=0.09, label="Reposo previo")
    current_axis.axvspan(0.0, release_boundary, color=BLUE, alpha=0.09, label="Pulso de corriente")
    current_axis.axvspan(release_boundary, time[-1], color=GREEN, alpha=0.09, label="Relajación")
    current_axis.axvline(0.0, color=GRAY, linewidth=0.9, linestyle=":")
    current_axis.axvline(release_boundary, color=GRAY, linewidth=0.9, linestyle=":")
    current_axis.set(ylabel="Corriente (A)", title="(a) Corriente aplicada")
    current_axis.legend(loc="upper right", ncol=3, fontsize=8.8)

    post_pulse = time >= 0.0
    voltage_axis.axvspan(
        time[0],
        0.0,
        color=GRAY,
        alpha=0.09,
        label="Reposo previo",
    )
    voltage_axis.axvspan(
        0.0,
        release_boundary,
        color=BLUE,
        alpha=0.075,
        label="Durante el pulso",
    )
    voltage_axis.axvspan(
        release_boundary,
        time[-1],
        color=GREEN,
        alpha=0.085,
        label="Relajación",
    )
    voltage_axis.hlines(
        pre_pulse_voltage,
        time[0],
        0.0,
        color=GRAY,
        linewidth=1.4,
        linestyle=(0, (3, 3)),
    )
    voltage_axis.plot(
        time[post_pulse],
        voltage[post_pulse],
        color=INK,
        linewidth=2.0,
        label="Voltaje experimental",
    )
    voltage_axis.scatter(
        [0.0],
        [immediate_voltage],
        s=48,
        color=ORANGE,
        edgecolor="white",
        linewidth=0.8,
        zorder=4,
        label="Respuesta inmediata",
    )
    voltage_axis.annotate(
        "",
        xy=(0.0, immediate_voltage),
        xytext=(0.0, pre_pulse_voltage),
        arrowprops={"arrowstyle": "->", "color": ORANGE, "linewidth": 1.5},
        zorder=3,
    )
    voltage_axis.text(
        2.0,
        0.5 * (pre_pulse_voltage + immediate_voltage),
        "Respuesta inmediata\n(al inicio del pulso)",
        color=ORANGE,
        fontsize=8.4,
        weight="semibold",
        ha="left",
        va="center",
    )
    voltage_axis.axvline(0.0, color=GRAY, linewidth=0.9, linestyle=":")
    voltage_axis.axvline(release_boundary, color=GRAY, linewidth=0.9, linestyle=":")
    voltage_axis.set(
        xlim=(time[0], time[-1]),
        xlabel="Tiempo desde el inicio del pulso (s)",
        ylabel="Voltaje (V)",
        title="(b) Respuesta de voltaje y fases de evaluación",
    )
    handles, labels = voltage_axis.get_legend_handles_labels()
    order = [0, 4, 1, 2, 3]
    voltage_axis.legend(
        [handles[index] for index in order],
        [labels[index] for index in order],
        loc="lower right",
        ncol=3,
        fontsize=8.2,
    )

    figure.patch.set_facecolor("white")
    figure.suptitle(
        "Definición de las fases de evaluación en un pulso HPPC",
        fontsize=18,
        weight="bold",
        color=INK,
    )
    package.save_figure(figure, "hppc_fases_pulso_final")
    return pulse.pulse_id, pulse.nominal_rate, pulse.soc


def update_figure_manifest() -> None:
    manifest_path = ROOT / "thesis_data" / "final_figure_manifest.csv"
    rows = package.csv_rows(manifest_path)
    new_stems = {
        "hppc_calibracion_validacion_final",
        "hppc_fases_pulso_final",
    }
    rows = [row for row in rows if row["stem"] not in new_stems]
    rows.extend(
        [
            {
                "figure_id": "F15",
                "stem": "hppc_calibracion_validacion_final",
                "role": "Conditional dynamic calibration/withheld validation",
                "purpose": "Full HPPC voltage overview, representative pulses, and accepted phase errors",
                "thesis_pdf": "thesis_figures/final/hppc_calibracion_validacion_final.pdf",
                "thesis_png": "thesis_figures/final/hppc_calibracion_validacion_final.png",
                "presentation_pdf": "presentation_assets/final/hppc_calibracion_validacion_final.pdf",
                "presentation_png": "presentation_assets/final/hppc_calibracion_validacion_final.png",
                "generated_by": "scripts/build_hppc_final_figures.py",
                "data_policy": "accepted artifacts and frozen forward-only reproduction; no fit",
            },
            {
                "figure_id": "F16",
                "stem": "hppc_fases_pulso_final",
                "role": "HPPC metric-definition aid",
                "purpose": "Explain current and voltage phases on retained real pulse 35",
                "thesis_pdf": "thesis_figures/final/hppc_fases_pulso_final.pdf",
                "thesis_png": "thesis_figures/final/hppc_fases_pulso_final.png",
                "presentation_pdf": "presentation_assets/final/hppc_fases_pulso_final.pdf",
                "presentation_png": "presentation_assets/final/hppc_fases_pulso_final.png",
                "generated_by": "scripts/build_hppc_final_figures.py",
                "data_policy": "accepted measured HPPC pulse; explanatory figure only; no fit",
            },
        ]
    )
    package.write_csv(manifest_path, rows)


def main() -> None:
    package.configure_style()
    package.THESIS_FIGURES.mkdir(parents=True, exist_ok=True)
    package.PRESENTATION.mkdir(parents=True, exist_ok=True)
    calibration_count, validation_count, calibration_id, validation_id = (
        figure_hppc_selection()
    )
    phase_id, phase_rate, phase_soc = figure_hppc_phases()
    update_figure_manifest()
    print(
        "HPPC figures generated without fitting: "
        f"calibration={calibration_count}, validation={validation_count}, "
        f"representative pulses={calibration_id}/{validation_id}, "
        f"phase pulse={phase_id} ({phase_rate}, SOC={phase_soc:.12f})."
    )


if __name__ == "__main__":
    main()
