"""Build final thesis figures and presentation assets from accepted results.

This script performs no fitting and changes no model parameter.  It reads the
accepted result tables and, where pointwise arrays were not retained, reruns
only the already-frozen forward equations used by the accepted analyses.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
DATA = ROOT / "thesis_data"
THESIS_FIGURES = ROOT / "thesis_figures" / "final"
PRESENTATION = ROOT / "presentation_assets" / "final"

BLUE = "#2563A6"
ORANGE = "#E07A2D"
GREEN = "#2E8B57"
RED = "#C7473D"
PURPLE = "#7656A5"
GRAY = "#6B7280"
LIGHT = "#E8EEF5"
INK = "#172033"


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 12,
            "axes.titlesize": 15,
            "axes.labelsize": 13,
            "axes.titleweight": "semibold",
            "legend.fontsize": 10.5,
            "xtick.labelsize": 10.5,
            "ytick.labelsize": 10.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": True,
            "grid.alpha": 0.22,
            "grid.linewidth": 0.7,
            "lines.linewidth": 2.2,
            "savefig.facecolor": "white",
        }
    )


def save_figure(fig: plt.Figure, stem: str) -> None:
    """Save identical-data publication and 16:9 presentation copies."""

    fig.set_size_inches(12.0, 6.75)
    fig.savefig(THESIS_FIGURES / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(THESIS_FIGURES / f"{stem}.png", dpi=320, bbox_inches="tight")
    fig.savefig(PRESENTATION / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(PRESENTATION / f"{stem}.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def add_metric_box(
    axis: plt.Axes,
    text: str,
    *,
    loc: tuple[float, float] = (0.02, 0.04),
    fontsize: float = 10.5,
) -> None:
    axis.text(
        *loc,
        text,
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        fontsize=fontsize,
        color=INK,
        bbox={"boxstyle": "round,pad=0.42", "facecolor": "white", "edgecolor": "#C7D2E2", "alpha": 0.94},
    )


def csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def metric_value(path: Path, **selectors: str) -> float:
    """Read exactly one value from an authoritative metrics CSV."""

    matches = [
        row
        for row in csv_rows(path)
        if all(row.get(key) == value for key, value in selectors.items())
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one metric in {path.name} for {selectors}, found {len(matches)}."
        )
    return float(matches[0]["value"])


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def electrical_forward_arrays() -> dict[str, dict[str, np.ndarray]]:
    from models.ocp_panasonic_qi_frozen import Q_CELL_AH, X_N_START, X_P_START
    from scripts.audit_panasonic_hppc_identifiability import simulate_local_group
    from scripts.validate_panasonic_1c_forward import (
        DISCHARGE_CURRENT_THRESHOLD_A,
        NOMINAL_PARAMETERS,
        load_raw_discharge,
        longest_contiguous_true,
    )
    from scripts.validate_panasonic_udds_forward import load_udds, remove_exact_duplicate_times
    from scripts.validate_panasonic_cccv_charge import (
        RAW_SOURCE_PATH as CCCV_SOURCE_PATH,
        forward_charge,
        identify_phases,
        initial_stoichiometries,
        load_mat_channels,
        preceding_discharge_capacity_ah,
        remove_duplicate_times,
    )

    _, raw_1c = load_raw_discharge()
    idx = longest_contiguous_true(raw_1c["Current"] <= DISCHARGE_CURRENT_THRESHOLD_A)
    time_1c = raw_1c["Time"][idx]
    current_1c = -raw_1c["Current"][idx]
    measured_1c = raw_1c["Voltage"][idx]
    soc_1c = 1.0 - (raw_1c["Ah"][idx][0] - raw_1c["Ah"][idx]) / Q_CELL_AH
    result_1c = simulate_local_group(
        time_s=time_1c,
        discharge_positive_current_a=current_1c,
        temperature_k=raw_1c["Battery_Temp_degC"][idx] + 273.15,
        initial_x_n=X_N_START,
        initial_x_p=X_P_START,
        dynamic_parameters=NOMINAL_PARAMETERS,
    )

    _, raw_udds = load_udds()
    udds, _ = remove_exact_duplicate_times(raw_udds)
    time_udds = udds["Time"]
    soc_udds = 1.0 + udds["Ah"] / Q_CELL_AH
    result_udds = simulate_local_group(
        time_s=time_udds,
        discharge_positive_current_a=-udds["Current"],
        temperature_k=udds["Battery_Temp_degC"] + 273.15,
        initial_x_n=X_N_START,
        initial_x_p=X_P_START,
        dynamic_parameters=NOMINAL_PARAMETERS,
    )

    _, raw_charge = load_mat_channels(CCCV_SOURCE_PATH)
    charge, _ = remove_duplicate_times(raw_charge)
    phases = identify_phases(charge["Current"], charge["Voltage"], charge["Ah"])
    preceding_capacity, _ = preceding_discharge_capacity_ah()
    initial_soc = 1.0 - preceding_capacity / Q_CELL_AH
    x_n0, x_p0 = initial_stoichiometries(initial_soc)
    result_charge = forward_charge(
        charge["Time"],
        charge["Current"],
        charge["Ah"],
        charge["Battery_Temp_degC"],
        x_n0,
        x_p0,
    )
    return {
        "one_c": {
            "time": time_1c - time_1c[0],
            "soc": soc_1c,
            "current": current_1c,
            "measured": measured_1c,
            "predicted": result_1c.voltage_v,
        },
        "udds": {
            "time": time_udds - time_udds[0],
            "soc": soc_udds,
            "current_raw": udds["Current"],
            "measured": udds["Voltage"],
            "predicted": result_udds.voltage_v,
        },
        "cccv": {
            "time": charge["Time"] - charge["Time"][0],
            "current_raw": charge["Current"],
            "measured": charge["Voltage"],
            "predicted": result_charge.voltage_v,
            "initial_rest": phases.initial_rest,
            "cc": phases.cc,
            "cv": phases.cv,
            "final_rest": phases.final_rest,
        },
    }


def figure_ocp_fit() -> None:
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
    fig.suptitle("Ajuste del voltaje de celda a C/20", fontsize=18, weight="bold", color=INK)
    save_figure(fig, "panasonic_ocp_fit_final")


def figure_hppc_transfer() -> None:
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
        for gi, group in enumerate(groups):
            for pulse in group.pulses:
                if pulse.record.nominal_rate == "2C":
                    candidates.append((abs(pulse.record.soc - target_soc), gi, pulse))
        _, gi, pulse = min(candidates, key=lambda item: item[0])
        group = groups[gi]
        model = evaluation.group_voltage_v[gi]
        start = int(pulse.phase_local_indices["immediate"][0])
        stop = int(pulse.phase_local_indices["post_relaxation"][-1])
        indices = np.arange(pulse.pre_local_index, stop + 1)
        return group.time_s[indices] - group.time_s[start], group.measured_voltage_v[indices], model[indices], pulse.record.soc

    high_t, high_m, high_p, high_soc = select(high_groups, high_eval, 0.55)
    low_t, low_m, low_p, low_soc = select(low_groups, low_eval, 0.20)
    fig, axes = plt.subplots(1, 3, constrained_layout=True, gridspec_kw={"width_ratios": [1.15, 1.15, 0.9]})
    for ax, t, measured, predicted, soc, title in (
        (axes[0], high_t, high_m, high_p, high_soc, "(a) Pulso 2C utilizado en la\ncalibración"),
        (axes[1], low_t, low_m, low_p, low_soc, "(b) Pulso 2C reservado para\nvalidación"),
    ):
        ax.plot(t, measured, color=INK, label="Datos experimentales")
        ax.plot(t, predicted, color=BLUE, linestyle="--", label="Modelo")
        ax.axvline(0.0, color=GRAY, lw=1, linestyle=":")
        ax.set(
            xlabel="Tiempo desde el inicio del pulso (s)",
            ylabel="Voltaje (V)",
            title=f"{title} — SOC = {soc:.2f}",
        )
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
    fig.suptitle("Calibración con HPPC y validación a bajo SOC", fontsize=18, weight="bold", color=INK)
    save_figure(fig, "panasonic_hppc_transfer_final")


def make_1c_validation_figure(data: dict[str, np.ndarray]) -> plt.Figure:
    """Build the Spanish 1C validation figure from authoritative saved metrics."""

    t_min = data["time"] / 60.0
    soc = data["soc"]
    residual_mv = 1000.0 * (data["predicted"] - data["measured"])
    metrics_path = DATA / "panasonic_1c_forward_validation_metrics.csv"

    # Same whole-region limits used by validate_panasonic_1c_forward.py.
    high_middle_time = float(np.interp(2.0 / 3.0, soc[::-1], t_min[::-1]))
    middle_low_time = float(np.interp(1.0 / 3.0, soc[::-1], t_min[::-1]))
    regional_specs = (
        ("SOC alto", "soc_high", t_min[0], high_middle_time, "#DCEEFF", 0.28),
        ("SOC medio", "soc_middle", high_middle_time, middle_low_time, "#EEE5FA", 0.10),
        ("SOC bajo", "soc_low", middle_low_time, t_min[-1], "#FCE8D5", 0.84),
    )
    global_rmse_mv = 1000.0 * metric_value(
        metrics_path, scope="overall", metric="rmse_v"
    )

    fig, axes = plt.subplots(
        2,
        1,
        sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [2.2, 1.0]},
    )
    for axis in axes:
        for _, _, start, stop, facecolor, _ in regional_specs:
            axis.axvspan(start, stop, color=facecolor, alpha=0.62, zorder=0)
        for boundary in (high_middle_time, middle_low_time):
            axis.axvline(
                boundary,
                color=GRAY,
                linewidth=0.9,
                linestyle=(0, (3, 3)),
                zorder=1,
            )
        axis.grid(True, color=GRAY, alpha=0.14, linewidth=0.6)

    axes[0].plot(
        t_min,
        data["measured"],
        color=INK,
        linewidth=1.7,
        label="Datos experimentales",
        zorder=3,
    )
    axes[0].plot(
        t_min,
        data["predicted"],
        color=BLUE,
        linewidth=1.7,
        linestyle="--",
        label="Modelo",
        zorder=3,
    )
    axes[0].set(
        ylabel="Voltaje (V)",
        title="(a) Voltaje durante la descarga",
    )
    axes[0].legend(loc="lower left", fontsize=9.0)

    for name, scope, start, stop, _, vertical_position in regional_specs:
        rmse_mv = 1000.0 * metric_value(
            metrics_path, scope=scope, metric="rmse_v"
        )
        axes[0].text(
            0.5 * (start + stop),
            vertical_position,
            f"{name}\nRMSE = {rmse_mv:.2f} mV",
            transform=axes[0].get_xaxis_transform(),
            ha="center",
            va="bottom" if vertical_position < 0.5 else "top",
            fontsize=9.0,
            color=INK,
            weight="semibold",
            zorder=5,
            bbox={
                "boxstyle": "round,pad=0.26",
                "facecolor": "white",
                "edgecolor": "#C7D2E2",
                "alpha": 0.90,
            },
        )

    axes[0].text(
        0.53,
        0.91,
        f"RMSE global = {global_rmse_mv:.2f} mV",
        transform=axes[0].transAxes,
        ha="center",
        va="top",
        fontsize=9.2,
        color=INK,
        zorder=5,
        bbox={
            "boxstyle": "round,pad=0.30",
            "facecolor": "white",
            "edgecolor": "#C7D2E2",
            "alpha": 0.94,
        },
    )

    axes[1].plot(t_min, residual_mv, color=ORANGE, linewidth=1.6, zorder=3)
    axes[1].axhline(0.0, color=INK, linewidth=0.9, zorder=2)
    axes[1].set(
        xlabel="Tiempo transcurrido (min)",
        ylabel="Error de voltaje (mV)",
        title="(b) Error de voltaje",
    )
    axes[1].set_xlim(t_min[0], t_min[-1])

    soc_axis = axes[0].twiny()
    soc_axis.set_xlim(axes[0].get_xlim())
    soc_axis.set_xticks(
        [t_min[0], high_middle_time, middle_low_time, t_min[-1]],
        ["1.00", "0.67", "0.33", f"{soc[-1]:.2f}"],
    )
    soc_axis.set_xlabel("SOC")

    fig.patch.set_facecolor("white")
    fig.suptitle(
        "Validación eléctrica durante descarga a 1C",
        fontsize=18,
        weight="bold",
        color=INK,
    )
    return fig


def figure_1c(data: dict[str, np.ndarray]) -> None:
    save_figure(make_1c_validation_figure(data), "panasonic_1c_voltage_validation_final")


def make_udds_validation_figure(
    data: dict[str, np.ndarray],
) -> tuple[plt.Figure, list[dict[str, float | str]]]:
    """Build the accepted UDDS comparison without fitting or error-based selection."""

    time_s = data["time"]
    hours = time_s / 3600.0
    metrics_path = DATA / "panasonic_udds_forward_validation_metrics.csv"
    region_specs = (
        ("SOC alto", "high", 2.0 / 3.0, 1.0, 5.0 / 6.0, BLUE, "(b)"),
        ("SOC medio", "middle", 1.0 / 3.0, 2.0 / 3.0, 0.5, PURPLE, "(c)"),
        ("SOC bajo", "low", 0.0, 1.0 / 3.0, 1.0 / 6.0, ORANGE, "(d)"),
    )
    window_half_width_s = 300.0
    windows: list[dict[str, float | str]] = []
    for (
        name,
        metric_subset,
        lower_soc,
        upper_soc,
        target_soc,
        color,
        panel,
    ) in region_specs:
        rmse_mv = 1000.0 * metric_value(
            metrics_path,
            scope="soc_region",
            subset=metric_subset,
            metric="rmse_v",
        )
        if upper_soc == 1.0:
            region_mask = (data["soc"] >= lower_soc) & (data["soc"] <= upper_soc)
        else:
            region_mask = (data["soc"] >= lower_soc) & (data["soc"] < upper_soc)
        region_indices = np.flatnonzero(region_mask)
        if region_indices.size == 0:
            raise ValueError(f"No UDDS samples found for {name}.")
        center_index = int(
            region_indices[np.argmin(np.abs(data["soc"][region_indices] - target_soc))]
        )
        start_s = max(time_s[center_index] - window_half_width_s, time_s[region_indices[0]])
        stop_s = min(time_s[center_index] + window_half_width_s, time_s[region_indices[-1]])
        mask = region_mask & (time_s >= start_s) & (time_s <= stop_s)
        current_span = float(np.ptp(data["current_raw"][mask]))
        if stop_s - start_s < 540.0 or current_span < 0.5:
            raise ValueError(f"The deterministic UDDS window for {name} lacks coverage.")
        windows.append(
            {
                "name": name,
                "panel": panel,
                "rmse_mv": rmse_mv,
                "color": color,
                "start_s": float(start_s),
                "stop_s": float(stop_s),
                "soc_min": float(np.min(data["soc"][mask])),
                "soc_max": float(np.max(data["soc"][mask])),
            }
        )

    figure = plt.figure(constrained_layout=True)
    grid = figure.add_gridspec(2, 3, height_ratios=[1.55, 1.0])
    full_axis = figure.add_subplot(grid[0, :])
    zoom_axes = [figure.add_subplot(grid[1, index]) for index in range(3)]

    step = max(1, time_s.size // 18000)
    selected = slice(None, None, step)
    full_axis.plot(
        hours[selected],
        data["measured"][selected],
        color=INK,
        linewidth=1.25,
        label="Datos experimentales",
        zorder=3,
    )
    full_axis.plot(
        hours[selected],
        data["predicted"][selected],
        color=BLUE,
        linewidth=1.15,
        linestyle="--",
        label="Modelo",
        zorder=3,
    )
    for window in windows:
        start_h = float(window["start_s"]) / 3600.0
        stop_h = float(window["stop_s"]) / 3600.0
        color = str(window["color"])
        full_axis.axvspan(
            start_h,
            stop_h,
            facecolor=color,
            edgecolor=color,
            alpha=0.075,
            linewidth=1.1,
            zorder=1,
        )
        full_axis.text(
            0.5 * (start_h + stop_h),
            0.93,
            str(window["panel"]),
            transform=full_axis.get_xaxis_transform(),
            ha="center",
            va="top",
            color=color,
            fontsize=10.0,
            weight="bold",
        )
    full_axis.set(
        xlim=(hours[0], hours[-1]),
        xlabel="Tiempo transcurrido (h)",
        ylabel="Voltaje (V)",
        title="(a) Perfil UDDS completo",
    )
    full_axis.legend(loc="lower left", ncol=2, fontsize=9.0)
    global_rmse_mv = 1000.0 * metric_value(
        metrics_path,
        scope="overall",
        subset="complete_usable",
        metric="rmse_v",
    )
    add_metric_box(
        full_axis,
        f"RMSE global = {global_rmse_mv:.2f} mV",
        loc=(0.80, 0.05),
    )

    for axis, window in zip(zoom_axes, windows):
        mask = (time_s >= float(window["start_s"])) & (
            time_s <= float(window["stop_s"])
        )
        local_hours = hours[mask]
        axis.plot(
            local_hours,
            data["measured"][mask],
            color=INK,
            linewidth=1.55,
            zorder=3,
        )
        axis.plot(
            local_hours,
            data["predicted"][mask],
            color=BLUE,
            linewidth=1.45,
            linestyle="--",
            zorder=3,
        )
        axis.set(
            xlim=(local_hours[0], local_hours[-1]),
            xlabel="Tiempo transcurrido (h)",
            ylabel="Voltaje (V)",
            title=f"{window['panel']} {window['name']}",
        )
        axis.text(
            0.04,
            0.94,
            f"RMSE de la regi\u00f3n de SOC = {float(window['rmse_mv']):.2f} mV",
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=8.8,
            color=INK,
            bbox={
                "boxstyle": "round,pad=0.28",
                "facecolor": "white",
                "edgecolor": "#D9DEE7",
                "alpha": 0.92,
            },
        )

    figure.patch.set_facecolor("white")
    figure.suptitle(
        "Validación eléctrica con perfil UDDS",
        fontsize=18,
        weight="bold",
        color=INK,
    )
    return figure, windows


def figure_udds(data: dict[str, np.ndarray]) -> None:
    figure, _ = make_udds_validation_figure(data)
    save_figure(figure, "panasonic_udds_voltage_validation_final")


def make_cccv_validation_figure(
    data: dict[str, np.ndarray],
) -> tuple[plt.Figure, float]:
    """Build the accepted CC-CV electrical validation without recalculation."""

    minutes = data["time"] / 60.0
    phase_indices = {
        phase: np.flatnonzero(data[phase])
        for phase in ("initial_rest", "cc", "cv", "final_rest")
    }
    if any(indices.size == 0 for indices in phase_indices.values()):
        raise ValueError("The accepted CC-CV phase masks must all be non-empty.")
    cc_start_min = float(minutes[phase_indices["cc"][0]])
    transition_min = float(minutes[phase_indices["cv"][0]])
    final_rest_start_min = float(minutes[phase_indices["final_rest"][0]])

    fig, (ax_v, ax_i) = plt.subplots(
        2,
        1,
        sharex=True,
        constrained_layout=True,
        gridspec_kw={"height_ratios": [2.2, 1.0]},
    )
    phase_spans = (
        (float(minutes[0]), cc_start_min, "Reposo inicial", GRAY, 0.055),
        (cc_start_min, transition_min, "CC\nCorriente constante", BLUE, 0.045),
        (transition_min, final_rest_start_min, "CV\nVoltaje constante", ORANGE, 0.055),
        (final_rest_start_min, float(minutes[-1]), "Reposo final", GRAY, 0.055),
    )
    for axis in (ax_v, ax_i):
        for start, stop, _, color, alpha in phase_spans:
            axis.axvspan(start, stop, color=color, alpha=alpha, zorder=0)
        for boundary in (cc_start_min, transition_min, final_rest_start_min):
            axis.axvline(
                boundary,
                color=GRAY,
                linewidth=0.9,
                linestyle=(0, (3, 3)),
                zorder=1,
            )
        axis.grid(True, color=GRAY, alpha=0.14, linewidth=0.6)

    ax_v.plot(
        minutes,
        data["measured"],
        color=INK,
        linewidth=1.8,
        label="Datos experimentales",
        zorder=3,
    )
    ax_v.plot(
        minutes,
        data["predicted"],
        color=BLUE,
        linewidth=1.8,
        linestyle="--",
        label="Modelo",
        zorder=3,
    )
    ax_v.set(ylabel="Voltaje (V)", title="(a) Voltaje durante la carga CC-CV")
    ax_v.legend(loc="lower right", ncol=2, fontsize=9.0)
    metrics_path = DATA / "panasonic_cccv_charge_validation_metrics.csv"
    global_rmse_mv = 1000.0 * metric_value(
        metrics_path,
        scope="overall",
        subset="complete_usable_charge",
        metric="rmse_v",
    )
    global_bias_mv = 1000.0 * metric_value(
        metrics_path,
        scope="overall",
        subset="complete_usable_charge",
        metric="mean_signed_error_v",
    )
    add_metric_box(
        ax_v,
        f"RMSE global (incluye reposos) = {global_rmse_mv:.2f} mV\n"
        f"Error medio (modelo − medición) = {global_bias_mv:.2f} mV",
        loc=(0.30, 0.05),
    )
    for start, stop, label, color, _ in phase_spans:
        ax_v.text(
            0.5 * (start + stop),
            0.84,
            label,
            transform=ax_v.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=9.0,
            color=color,
            weight="semibold",
        )

    ax_i.plot(
        minutes,
        data["current_raw"],
        color=PURPLE,
        linewidth=1.8,
        zorder=3,
    )
    ax_i.axhline(0.0, color=INK, linewidth=0.9, zorder=2)
    ax_i.set(
        xlim=(minutes[0], minutes[-1]),
        xlabel="Tiempo transcurrido (min)",
        ylabel="Corriente de carga (A)",
        title="(b) Corriente aplicada durante la carga",
    )

    fig.patch.set_facecolor("white")
    fig.suptitle(
        "Validación eléctrica durante carga CC-CV",
        fontsize=18,
        weight="bold",
        color=INK,
    )
    return fig, transition_min


def figure_cccv(data: dict[str, np.ndarray]) -> None:
    figure, _ = make_cccv_validation_figure(data)
    save_figure(figure, "panasonic_cccv_voltage_transfer_final")


def make_c20_charge_discharge_figure() -> plt.Figure:
    """Build the accepted low-rate comparison on the existing common-SOC map."""

    rows = csv_rows(DATA / "panasonic_c20_charge_discharge_comparison.csv")
    soc = np.asarray([float(row["soc_common_capacity"]) for row in rows])
    discharge = np.asarray([float(row["voltage_c20_discharge_v"]) for row in rows])
    charge = np.asarray([float(row["voltage_c20_charge_v"]) for row in rows])

    fig, axis = plt.subplots(constrained_layout=True)
    axis.plot(
        soc,
        discharge,
        color=BLUE,
        linewidth=2.2,
        label="Descarga C/20",
    )
    axis.plot(
        soc,
        charge,
        color=ORANGE,
        linewidth=2.2,
        label="Carga C/20",
    )
    axis.set(
        xlim=(float(np.min(soc)), float(np.max(soc))),
        xlabel="SOC normalizado para comparación",
        ylabel="Voltaje terminal (V)",
        title="Voltaje a baja tasa durante carga y descarga",
    )
    axis.legend(loc="lower right", fontsize=10.0)
    axis.grid(True, color=GRAY, alpha=0.16, linewidth=0.6)
    fig.patch.set_facecolor("white")
    return fig


def figure_directional_separation() -> None:
    figure = make_c20_charge_discharge_figure()
    save_figure(figure, "panasonic_directional_voltage_separation_final")


def thermal_forward_arrays() -> dict[str, dict[str, np.ndarray]]:
    from models.ocp_panasonic_qi_frozen import U_n, U_p
    from models.panasonic_entropic import dUeq_dT_panasonic
    from models.panasonic_thermal_parameters import C_TH_J_PER_K, G_TH_W_PER_K
    from scripts.calibrate_panasonic_1c_thermal import exact_thermal_trajectory, prepare_data as prepare_1c
    from scripts.calibrate_panasonic_joint_1c_hppc_thermal import exact_one_node
    from scripts.validate_panasonic_udds_thermal import (
        exact_thermal_trajectory as exact_udds,
        fully_predictive_trajectory,
        prepare_data as prepare_udds,
    )
    from scripts.validate_panasonic_cccv_thermal import prepare_data as prepare_cccv, validation_b

    one_c = prepare_1c()
    predicted_1c = exact_thermal_trajectory(
        one_c.time_s,
        one_c.qdot_cal_w,
        one_c.ambient_temperature_k,
        one_c.measured_temperature_k[0],
        C_TH_J_PER_K,
        G_TH_W_PER_K,
    )

    udds, _ = prepare_udds()
    qirr_u = udds.current_a * (udds.ueq_average_v - udds.measured_voltage_v)
    qrev_u = -udds.current_a * udds.measured_temperature_k * udds.dudt_v_per_k
    predicted_udds_a = exact_udds(
        udds.time_s, qirr_u + qrev_u, udds.ambient_temperature_k, udds.measured_temperature_k[0]
    )
    predicted_udds_b, _, _, _, _ = fully_predictive_trajectory(udds)

    cccv = prepare_cccv()
    ueq_c = np.asarray(U_p(cccv.x_p_avg)) - np.asarray(U_n(cccv.x_n_avg))
    dudt_c = np.asarray(dUeq_dT_panasonic(cccv.soc))
    qirr_c = cccv.model_current_a * (ueq_c - cccv.measured_voltage_v)
    qrev_c = -cccv.model_current_a * cccv.measured_temperature_k * dudt_c
    predicted_cccv_a = exact_one_node(
        cccv.time_s,
        qirr_c + qrev_c,
        cccv.ambient_temperature_k,
        cccv.measured_temperature_k[0],
        C_TH_J_PER_K,
        G_TH_W_PER_K,
    )
    predicted_cccv_b = validation_b(cccv, ueq_c, dudt_c).temperature_k
    return {
        "one_c": {"time": one_c.time_s - one_c.time_s[0], "soc": one_c.soc, "measured": one_c.measured_temperature_k - 273.15, "predicted": predicted_1c - 273.15},
        "udds": {"time": udds.time_s, "soc": udds.soc, "measured": udds.measured_temperature_k - 273.15, "a": predicted_udds_a - 273.15, "b": predicted_udds_b - 273.15},
        "cccv": {"time": cccv.time_s, "soc": cccv.soc, "measured": cccv.measured_temperature_k - 273.15, "a": predicted_cccv_a - 273.15, "b": predicted_cccv_b - 273.15, "cc": cccv.cc_mask, "cv": cccv.cv_mask},
    }


def figure_thermal_1c(data: dict[str, np.ndarray]) -> None:
    minutes = data["time"] / 60.0
    residual = data["predicted"] - data["measured"]
    fig, axes = plt.subplots(
        2,
        1,
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0]},
    )

    # Present the measurements as discrete observations and the one-node
    # prediction as a continuous model trajectory.  The underlying accepted
    # arrays and the residual definition are unchanged.
    axes[0].plot(
        minutes,
        data["measured"],
        linestyle="none",
        marker="o",
        markersize=3.0,
        markeredgewidth=0.0,
        color="#20242B",
        alpha=0.78,
        label="Datos experimentales",
        zorder=3,
    )
    axes[0].plot(
        minutes,
        data["predicted"],
        color=RED,
        linewidth=2.0,
        linestyle="--",
        label="Modelo de un nodo",
        zorder=2,
    )
    axes[0].set_ylabel("Temperatura (°C)")
    axes[0].set_title("(a) Temperatura durante la descarga", loc="left", pad=8)
    axes[0].legend(loc="upper left", frameon=False)
    add_metric_box(
        axes[0],
        "Cth = 117.52 J/K\nGth = 0.1091 W/K\nRMSE = 0.144 K\nR² = 0.965",
        loc=(0.80, 0.04),
    )

    axes[1].scatter(
        minutes,
        residual,
        s=13,
        color=ORANGE,
        edgecolors="none",
        alpha=0.82,
        zorder=3,
    )
    axes[1].axhline(0.0, color=INK, lw=1.0, alpha=0.85, zorder=2)
    axes[1].set_xlabel("Tiempo transcurrido (min)")
    axes[1].set_ylabel("Residuo (K)")
    axes[1].set_title("(b) Residuo de calibración", loc="left", pad=8)

    fig.suptitle("Calibración térmica en descarga 1C", fontsize=18, weight="bold", color=INK, y=0.985)
    fig.text(
        0.5,
        0.935,
        "Intervalo de calibración: SOC 0.90–0.10",
        ha="center",
        va="top",
        fontsize=11.5,
        color=GRAY,
    )
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.11, top=0.84, hspace=0.36)
    save_figure(fig, "panasonic_1c_thermal_calibration_final")


def make_thermal_udds_figure(data: dict[str, np.ndarray]) -> plt.Figure:
    hours = data["time"] / 3600.0
    model_step = max(1, data["time"].size // 16000)
    measured_step = max(1, data["time"].size // 5000)
    model_sl = slice(None, None, model_step)
    measured_sl = slice(None, None, measured_step)
    residual_a = data["a"] - data["measured"]
    residual_b = data["b"] - data["measured"]
    rmse_a = float(np.sqrt(np.mean(residual_a**2)))
    rmse_b = float(np.sqrt(np.mean(residual_b**2)))
    bias_a = float(np.mean(residual_a))
    bias_b = float(np.mean(residual_b))

    metrics_path = DATA / "panasonic_udds_thermal_validation_metrics.csv"
    authoritative = {
        "rmse_a": metric_value(
            metrics_path,
            validation="validation_a",
            subset="overall",
            metric="rmse_k",
        ),
        "bias_a": metric_value(
            metrics_path,
            validation="validation_a",
            subset="overall",
            metric="mean_signed_error_k",
        ),
        "rmse_b": metric_value(
            metrics_path,
            validation="validation_b",
            subset="overall",
            metric="rmse_k",
        ),
        "bias_b": metric_value(
            metrics_path,
            validation="validation_b",
            subset="overall",
            metric="mean_signed_error_k",
        ),
    }
    computed = {"rmse_a": rmse_a, "bias_a": bias_a, "rmse_b": rmse_b, "bias_b": bias_b}
    for name, value in computed.items():
        if not np.isclose(value, authoritative[name], rtol=0.0, atol=5.0e-12):
            raise ValueError(
                f"UDDS thermal curve metric {name}={value:.15g} does not match "
                f"the authoritative result {authoritative[name]:.15g}."
            )

    def format_metric(value: float) -> str:
        return f"{value:.3f}".replace("-", "−")

    # Validation A uses measured voltage for irreversible heat and measured
    # temperature for reversible heat. Validation B uses model voltage and
    # predicted temperature (the coupled electro-thermal prediction).
    fig, axes = plt.subplots(
        2,
        1,
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0]},
    )
    axes[0].scatter(
        hours[measured_sl],
        data["measured"][measured_sl],
        s=5,
        color="#20242B",
        edgecolors="none",
        alpha=0.62,
        label="Datos experimentales",
        zorder=3,
    )
    axes[0].plot(
        hours[model_sl],
        data["a"][model_sl],
        color=RED,
        lw=1.35,
        linestyle="--",
        label="Validación A",
        zorder=2,
    )
    axes[0].plot(
        hours[model_sl],
        data["b"][model_sl],
        color=PURPLE,
        lw=1.15,
        linestyle=":",
        label="Validación B",
        zorder=2,
    )
    axes[0].set_ylabel("Temperatura (°C)")
    axes[0].set_title("(a) Temperatura durante el perfil UDDS", loc="left", pad=8)
    axes[0].legend(loc="upper left", ncol=3, frameon=False)
    add_metric_box(
        axes[0],
        "Validación A                 Validación B\n"
        f"RMSE = {authoritative['rmse_a']:.3f} K         RMSE = {authoritative['rmse_b']:.3f} K\n"
        f"Sesgo = {format_metric(authoritative['bias_a'])} K        "
        f"Sesgo = {format_metric(authoritative['bias_b'])} K",
        loc=(0.40, 0.62),
        fontsize=9.4,
    )

    axes[1].plot(
        hours[model_sl],
        residual_a[model_sl],
        color=ORANGE,
        lw=0.65,
        alpha=0.85,
        label="Validación A",
    )
    axes[1].plot(
        hours[model_sl],
        residual_b[model_sl],
        color=PURPLE,
        lw=0.60,
        alpha=0.78,
        label="Validación B",
    )
    axes[1].axhline(0.0, color=INK, lw=1.0, alpha=0.85)
    axes[1].set_xlabel("Tiempo transcurrido (h)")
    axes[1].set_ylabel("Residuo\n(modelo − medición)\n[K]")
    axes[1].set_title("(b) Residuo de validación", loc="left", pad=8)
    axes[1].legend(loc="lower left", frameon=False)

    fig.suptitle("Validación térmica con perfil UDDS", fontsize=18, weight="bold", color=INK, y=0.985)
    fig.text(
        0.5,
        0.935,
        "Intervalo analizado: SOC 0.90–0.10",
        ha="center",
        va="top",
        fontsize=11.5,
        color=GRAY,
    )
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.11, top=0.84, hspace=0.36)
    return fig


def figure_thermal_udds(data: dict[str, np.ndarray]) -> None:
    fig = make_thermal_udds_figure(data)
    save_figure(fig, "panasonic_udds_thermal_transfer_final")


def make_thermal_cccv_figure(data: dict[str, np.ndarray]) -> plt.Figure:
    minutes = data["time"] / 60.0
    residual_a = data["a"] - data["measured"]
    residual_b = data["b"] - data["measured"]
    metrics_path = DATA / "panasonic_cccv_thermal_validation_metrics.csv"
    authoritative = {
        "rmse_a": metric_value(
            metrics_path,
            scope="temperature_metrics",
            subset="validation_a",
            metric="rmse_k",
        ),
        "bias_a": metric_value(
            metrics_path,
            scope="temperature_metrics",
            subset="validation_a",
            metric="mean_signed_error_k",
        ),
        "rmse_b": metric_value(
            metrics_path,
            scope="temperature_metrics",
            subset="validation_b",
            metric="rmse_k",
        ),
        "bias_b": metric_value(
            metrics_path,
            scope="temperature_metrics",
            subset="validation_b",
            metric="mean_signed_error_k",
        ),
    }
    computed = {
        "rmse_a": float(np.sqrt(np.mean(residual_a**2))),
        "bias_a": float(np.mean(residual_a)),
        "rmse_b": float(np.sqrt(np.mean(residual_b**2))),
        "bias_b": float(np.mean(residual_b)),
    }
    for name, value in computed.items():
        if not np.isclose(value, authoritative[name], rtol=0.0, atol=5.0e-12):
            raise ValueError(
                f"CC-CV thermal curve metric {name}={value:.15g} does not match "
                f"the authoritative result {authoritative[name]:.15g}."
            )
    cc = np.flatnonzero(data["cc"])
    cv = np.flatnonzero(data["cv"])
    if not cc.size or not cv.size:
        raise ValueError("The accepted CC and CV masks are required for the final thermal figure.")

    cc_start = minutes[cc[0]]
    transition = minutes[cv[0]]
    cv_end = minutes[cv[-1]]
    cc_label_time = cc_start + 0.72 * (transition - cc_start)
    cv_label_time = 0.5 * (transition + cv_end)

    # Validation A uses measured voltage for irreversible heat and measured
    # temperature for reversible heat. Validation B uses model voltage and
    # predicted temperature (the coupled electro-thermal prediction).
    fig, axes = plt.subplots(
        2,
        1,
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1.0]},
    )

    for axis in axes:
        axis.axvspan(cc_start, transition, color="#DCEAF5", alpha=0.34, lw=0, zorder=0)
        axis.axvspan(transition, cv_end, color="#F8E6D1", alpha=0.42, lw=0, zorder=0)
        axis.axvline(transition, color=GRAY, lw=1.0, linestyle="--", alpha=0.85, zorder=1)
        axis.text(
            cc_label_time,
            1.01,
            "CC\ncorriente constante",
            transform=axis.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=9.0,
            color=BLUE,
            linespacing=1.05,
            zorder=5,
        )
        axis.text(
            cv_label_time,
            1.01,
            "CV\nvoltaje constante",
            transform=axis.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=9.0,
            color=ORANGE,
            linespacing=1.05,
            zorder=5,
        )

    axes[0].scatter(
        minutes,
        data["measured"],
        s=16,
        color="#20242B",
        edgecolors="none",
        alpha=0.82,
        label="Datos experimentales",
        zorder=4,
    )
    axes[0].plot(
        minutes,
        data["a"],
        color=RED,
        lw=1.7,
        linestyle="--",
        label="Validación A",
        zorder=3,
    )
    axes[0].plot(
        minutes,
        data["b"],
        color=PURPLE,
        lw=1.45,
        linestyle=":",
        label="Validación B",
        zorder=3,
    )
    axes[0].set_ylabel("Temperatura (°C)")
    axes[0].set_title("(a) Temperatura durante la carga", loc="left", pad=8)
    axes[0].legend(loc="upper left", ncol=3, frameon=False)
    add_metric_box(
        axes[0],
        "Validación A                 Validación B\n"
        f"RMSE = {authoritative['rmse_a']:.3f} K         "
        f"RMSE = {authoritative['rmse_b']:.3f} K\n"
        f"Sesgo = {authoritative['bias_a']:+.3f} K       "
        f"Sesgo = {authoritative['bias_b']:+.3f} K",
        loc=(0.03, 0.54),
        fontsize=9.4,
    )

    axes[1].plot(minutes, residual_a, color=ORANGE, lw=1.15, label="Validación A", zorder=3)
    axes[1].plot(minutes, residual_b, color=PURPLE, lw=1.05, label="Validación B", zorder=3)
    axes[1].axhline(0.0, color=INK, lw=1.0, alpha=0.85, zorder=2)
    axes[1].set_xlabel("Tiempo transcurrido (min)")
    axes[1].set_ylabel("Residuo\n(modelo − medición)\n[K]")
    axes[1].set_title("(b) Residuo de validación", loc="left", pad=8)
    axes[1].legend(loc="lower left", frameon=False)

    fig.suptitle("Validación térmica durante carga CC-CV", fontsize=18, weight="bold", color=INK, y=0.985)
    fig.text(
        0.5,
        0.935,
        "Intervalo analizado: SOC 0.10–0.90",
        ha="center",
        va="top",
        fontsize=11.5,
        color=GRAY,
    )
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.11, top=0.84, hspace=0.36)
    return fig


def figure_thermal_cccv(data: dict[str, np.ndarray]) -> None:
    fig = make_thermal_cccv_figure(data)
    save_figure(fig, "panasonic_cccv_thermal_transfer_final")


def figure_thermal_selection() -> None:
    categories = ["1C fit", "UDDS transfer"]
    one = [0.143622065, 0.526697652]
    two = [0.138221456, 0.532525609]
    x = np.arange(2)
    width = 0.34
    fig, (ax, note) = plt.subplots(1, 2, constrained_layout=True, gridspec_kw={"width_ratios": [1.25, 1.0]})
    ax.bar(x - width / 2, one, width, color=BLUE, label="One node")
    ax.bar(x + width / 2, two, width, color=ORANGE, label="Clean two node")
    ax.set_xticks(x, categories)
    ax.set(ylabel="Temperature RMSE (K)", title="Accuracy comparison")
    ax.legend()
    for xpos, values in zip((x - width / 2, x + width / 2), (one, two), strict=True):
        for xx, value in zip(xpos, values, strict=True):
            ax.text(xx, value + 0.012, f"{value:.3f}", ha="center", fontsize=11)
    note.axis("off")
    note.text(
        0.02,
        0.92,
        "Two-node recoverability audit",
        fontsize=16,
        weight="semibold",
        color=INK,
        va="top",
    )
    note.text(
        0.02,
        0.78,
        "• surface capacitance reached lower guard\n"
        "• fast mode ≈ 1.00×10⁻⁴ s\n"
        "• local Jacobian condition ≈ 2.56×10⁷\n"
        "• no independent UDDS improvement\n\n"
        "Decision: retain the parsimonious one-node model.\n"
        "Two-node parameters are not recoverable from the\n"
        "available surface-temperature experiment.",
        fontsize=13,
        color=INK,
        va="top",
        linespacing=1.5,
    )
    fig.suptitle("Thermal architecture selection", fontsize=18, weight="bold", color=INK)
    save_figure(fig, "thermal_model_selection_summary")


def figure_uniandes_rapp() -> None:
    rows = csv_rows(DATA / "own_data_apparent_resistance_transitions.csv")
    order = [
        "cycle_cell1_ambient",
        "cycle_cell2_undeclared",
        "cycle_cell2_30c",
        "cycle_cell1_35c",
        "cycle_cell1_30c",
        "pulse_cell1_part1",
        "pulse_cell1_part2",
    ]
    values = [[float(row["r_app_ohm"]) for row in rows if row["dataset_id"] == name] for name in order]
    labels = ["Ambient\ncell 1", "Cell 2\nundeclared", "Cell 2\n30 °C", "Cell 1\n35 °C", "Cell 1\n30 °C", "Pulse\npart 1", "Pulse\npart 2"]
    fig, ax = plt.subplots(constrained_layout=True)
    box = ax.boxplot(values, patch_artist=True, tick_labels=labels, widths=0.62, showmeans=True)
    for index, patch in enumerate(box["boxes"]):
        patch.set_facecolor(GREEN if index == 0 else ORANGE)
        patch.set_alpha(0.75)
    for index, series in enumerate(values, start=1):
        jitter = np.linspace(-0.08, 0.08, len(series)) if len(series) > 1 else np.zeros(1)
        ax.scatter(index + jitter, series, color=INK, s=25, zorder=3)
    ax.axhline(0.0358755369, color=GREEN, linestyle="--", label="Clean-reference median 0.03588 Ω")
    ax.set(ylabel=r"Apparent transition resistance, $R_{app}$ (Ω)", title="Adjacent rest-to-load diagnostic")
    ax.legend(loc="upper left")
    add_metric_box(ax, "Affected medians: 0.178–0.232 Ω\n≈5–6.5× clean reference\nRapp is not intrinsic cell R0", loc=(0.58, 0.05))
    fig.suptitle("Uniandes data: evidence of elevated series-path response", fontsize=18, weight="bold", color=INK)
    save_figure(fig, "uniandes_apparent_series_path_final")


def own_case_arrays() -> tuple[np.ndarray, dict[str, np.ndarray], dict[str, object]]:
    from analysis.own_case_study_reference import SOC_GRID, discharge_curve, interpolate_on_soc, segment_complete_cycle
    from scripts.audit_own_trc18650_data import load_all

    experiments = {item.metadata.dataset_id: item for item in load_all()}
    ids = ["cycle_cell1_ambient", "cycle_cell1_30c", "cycle_cell1_35c", "cycle_cell2_undeclared", "cycle_cell2_30c"]
    curves: dict[str, np.ndarray] = {}
    for dataset_id in ids:
        segments = segment_complete_cycle(experiments[dataset_id])
        discharge = next(segment for segment in segments if segment.name == "main_1c_discharge")
        curves[dataset_id] = np.asarray(interpolate_on_soc(discharge_curve(experiments[dataset_id], discharge))["voltage_v"])
    return SOC_GRID, curves, experiments


def make_uniandes_discharge_figures(
    experiments: dict[str, object],
) -> tuple[plt.Figure, plt.Figure]:
    from scripts.build_uniandes_context_figures import (
        make_discharge_figures,
        write_context_documentation,
    )

    defective, contact_comparison, summary = make_discharge_figures(experiments)
    write_context_documentation(summary)
    return defective, contact_comparison


def figure_uniandes_discharge_comparisons(experiments: dict[str, object]) -> None:
    defective, contact_comparison = make_uniandes_discharge_figures(experiments)
    save_figure(
        defective,
        "uniandes_defective_contacts_discharge_final",
    )
    save_figure(
        contact_comparison,
        "uniandes_contact_arrangements_comparison_final",
    )


def make_uniandes_geographical_figure() -> plt.Figure:
    from scripts.build_uniandes_context_figures import make_geographical_figure

    return make_geographical_figure()


def figure_uniandes_geographical_context() -> None:
    save_figure(
        make_uniandes_geographical_figure(),
        "uniandes_geographical_context_final",
    )


def figure_uniandes_charge_temperature() -> None:
    from scripts.build_uniandes_charge_temperature_figures import build_central_outputs

    build_central_outputs()


def figure_transfer(soc: np.ndarray, curves: dict[str, np.ndarray], experiments: dict[str, object]) -> None:
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
    axes[0].plot(soc, curves["cycle_cell1_ambient"], color=INK, label="TRC measured loaded 1C")
    axes[0].plot(soc, backbone, color=ORANGE, linestyle="--", label="Panasonic equilibrium backbone")
    axes[0].plot(soc_raw[forward.ocp_valid], forward.voltage_v[forward.ocp_valid], color=BLUE, alpha=0.9, label="Literal frozen full transfer")
    axes[0].invert_xaxis()
    axes[0].set(xlabel="TRC reference SOC", ylabel="Voltage (V)", title="Direct Panasonic-to-TRC transfer")
    axes[0].legend(fontsize=9)
    add_metric_box(axes[0], "Backbone RMSE 188.882 mV\nFull-transfer RMSE 82.134 mV")
    residual_backbone = 1000.0 * (backbone - curves["cycle_cell1_ambient"])
    axes[1].plot(soc, residual_backbone, color=ORANGE, label="Backbone − measured")
    axes[1].plot(soc_raw[forward.ocp_valid], 1000.0 * (forward.voltage_v[forward.ocp_valid] - measured_raw[forward.ocp_valid]), color=BLUE, label="Full transfer − measured")
    axes[1].axhline(0.0, color=INK, lw=1)
    axes[1].invert_xaxis()
    axes[1].set(xlabel="TRC reference SOC", ylabel="Residual (mV)", title="SOC-shaped mismatch and low-SOC support loss")
    axes[1].legend(fontsize=9)
    add_metric_box(axes[1], "PROVISIONAL-TIME DIAGNOSTIC\n22 positive-surface support violations\nNo TRC accuracy claim", loc=(0.04, 0.62))
    fig.suptitle("Panasonic benchmark transfer to Uniandes TRC18650 — INCOMPATIBLE DIRECT TRANSFER", fontsize=17, weight="bold", color=INK)
    save_figure(fig, "panasonic_to_trc_transfer_final")


def figure_assessment_flow() -> None:
    fig, ax = plt.subplots(constrained_layout=True)
    ax.set_xlim(0, 12.6)
    ax.set_ylim(0, 6.7)
    ax.axis("off")
    nodes = {
        "new": (0.25, 5.2, "NEW\nEXPERIMENTAL\nDATA"),
        "audit": (2.75, 5.2, "DATA QUALITY\nAUDIT"),
        "class": (5.25, 5.2, "SUITABILITY\nCLASSIFICATION"),
        "diag": (5.25, 2.4, "DIAGNOSTIC\nCASE STUDY"),
        "transfer": (7.75, 5.2, "MODEL TRANSFER\nTEST"),
        "compat": (10.25, 5.2, "PARAMETER-\nCOMPATIBILITY\nAUDIT"),
        "limit": (8.25, 0.85, "REPORT\nLIMITATION"),
        "valid": (10.25, 2.4, "QUANTITATIVE\nVALIDATION"),
    }
    box_w, box_h = 2.05, 1.02
    for key, (x, y, label) in nodes.items():
        color = GREEN if key == "valid" else (ORANGE if key in {"diag", "limit"} else BLUE)
        ax.add_patch(plt.Rectangle((x, y), box_w, box_h, facecolor=color, edgecolor="white", linewidth=2, alpha=0.92))
        ax.text(x + box_w / 2, y + box_h / 2, label, ha="center", va="center", color="white", fontsize=9.3, weight="bold", linespacing=1.05)

    def arrow(source, target, label=""):
        sx, sy, _ = nodes[source]
        tx, ty, _ = nodes[target]
        start = (sx + box_w, sy + box_h / 2) if tx > sx else (sx + box_w / 2, sy)
        end = (tx, ty + box_h / 2) if tx > sx else (tx + box_w / 2, ty + box_h)
        ax.annotate("", xy=end, xytext=start, arrowprops={"arrowstyle": "->", "lw": 2, "color": INK})
        if label:
            ax.text((start[0] + end[0]) / 2, (start[1] + end[1]) / 2 + 0.17, label, ha="center", fontsize=9.5, color=INK)

    arrow("new", "audit")
    arrow("audit", "class")
    arrow("class", "transfer", "suitable evidence")
    arrow("transfer", "compat")
    arrow("class", "diag", "nonideal / insufficient")
    arrow("diag", "limit")
    arrow("compat", "valid", "requirements supported")
    arrow("compat", "limit", "adaptation unsupported")
    ax.text(0.8, 0.4, "No automatic fitting • no automatic voltage correction • validation is conditional, not the default endpoint", color=INK, fontsize=13, weight="semibold")
    fig.suptitle("Reusable battery-dataset assessment and transfer methodology", fontsize=18, weight="bold", color=INK)
    save_figure(fig, "battery_dataset_assessment_flow_final")


def result_register() -> list[dict[str, object]]:
    return [
        {"result_id": "E1", "domain": "electrical", "experiment": "Panasonic C/20 discharge", "role": "OCP identification", "metric": "5–95% voltage RMSE", "value": 1.085749161, "unit": "mV", "status": "ACCEPTED REPRESENTATION", "limitation": "Stoichiometric windows weakly identifiable; endpoints are nominal effective windows", "source_file": "thesis_data/panasonic_ocp_identification_closure.md"},
        {"result_id": "E2", "domain": "electrical", "experiment": "Panasonic HPPC SOC >= 1/3", "role": "conditional dynamic calibration", "metric": "phase-balanced RMSE", "value": 36.835197639, "unit": "mV", "status": "CALIBRATION RESULT", "limitation": "Conditional three-parameter fit; parameters not globally identified", "source_file": "thesis_data/panasonic_hppc_primary_fit_no_low_soc_report.md"},
        {"result_id": "E3", "domain": "electrical", "experiment": "Panasonic HPPC SOC < 1/3", "role": "withheld forward validation", "metric": "phase-balanced RMSE", "value": 85.729432921, "unit": "mV", "status": "LOW-SOC TRANSFER LIMITED", "limitation": "Dominant loaded-transient mismatch; no refit", "source_file": "thesis_data/panasonic_hppc_low_soc_validation_report.md"},
        {"result_id": "E4", "domain": "electrical", "experiment": "Panasonic 1C discharge", "role": "independent forward validation", "metric": "overall voltage RMSE", "value": 114.362188001, "unit": "mV", "status": "1C TRANSFER LIMITED", "limitation": "Low-SOC RMSE 195.954 mV and cutoff deterioration", "source_file": "thesis_data/panasonic_1c_forward_validation_report.md"},
        {"result_id": "E5", "domain": "electrical", "experiment": "Panasonic UDDS", "role": "independent forward validation", "metric": "overall voltage RMSE", "value": 30.988191694, "unit": "mV", "status": "UDDS TRANSFER LIMITED", "limitation": "High/middle/low SOC RMSE 14.361/21.523/52.785 mV; large-step mismatch", "source_file": "thesis_data/panasonic_udds_forward_validation_report.md"},
        {"result_id": "E6", "domain": "electrical", "experiment": "Panasonic 25 C CC–CV charge", "role": "independent directional validation", "metric": "overall voltage RMSE", "value": 61.033526348, "unit": "mV", "status": "CC-CV CHARGE TRANSFER LIMITED", "limitation": "Bias −47.526 mV; discharge-derived representation; regulated CV is measured reference", "source_file": "thesis_data/panasonic_cccv_charge_validation_report.md"},
        {"result_id": "T1", "domain": "thermal", "experiment": "Panasonic 1C discharge", "role": "primary one-node calibration", "metric": "temperature RMSE", "value": 0.143622065, "unit": "K", "status": "PRIMARY THERMAL CALIBRATION", "limitation": "Measured-voltage heat; SOC 0.10–0.90; effective parameters", "source_file": "thesis_data/panasonic_1c_thermal_calibration.md"},
        {"result_id": "T2", "domain": "thermal", "experiment": "Panasonic UDDS", "role": "independent one-node transfer A", "metric": "temperature RMSE", "value": 0.526697652, "unit": "K", "status": "THERMAL TRANSFER LIMITED", "limitation": "Bias −0.512715 K; measured-voltage heat path", "source_file": "thesis_data/panasonic_udds_thermal_validation.md"},
        {"result_id": "T3", "domain": "thermal", "experiment": "Panasonic CC–CV charge", "role": "independent one-node transfer A", "metric": "temperature RMSE", "value": 0.427780949, "unit": "K", "status": "CHARGE-SIDE THERMAL TRANSFER LIMITED", "limitation": "Near-zero overall bias masks CC/CV bias reversal", "source_file": "thesis_data/panasonic_cccv_thermal_validation.md"},
        {"result_id": "T4", "domain": "electrothermal", "experiment": "Panasonic CC–CV charge", "role": "fully predictive transfer B", "metric": "temperature RMSE", "value": 0.864391266, "unit": "K", "status": "COUPLED TRANSFER FAILED", "limitation": "Electrical heat deficit 436.438 J dominates added degradation", "source_file": "thesis_data/panasonic_cccv_thermal_validation.md"},
        {"result_id": "T5", "domain": "thermal", "experiment": "Panasonic 1C / UDDS", "role": "architecture selection", "metric": "two-node 1C / UDDS RMSE", "value": "0.138221 / 0.532526", "unit": "K", "status": "TWO-NODE NOT RECOVERABLE", "limitation": "Surface capacitance at lower guard; fast mode ~1e-4 s; condition ~2.56e7", "source_file": "thesis_data/panasonic_two_node_physical_calibration.md"},
        {"result_id": "U1", "domain": "data quality", "experiment": "Uniandes TRC18650", "role": "series-path audit", "metric": "median R_app clean / affected", "value": "0.03588 / 0.178–0.232", "unit": "ohm", "status": "ELEVATED SERIES-PATH RESPONSE", "limitation": "R_app is a setup diagnostic, not intrinsic R0; affected data rejected for electrical identification", "source_file": "thesis_data/uniandes_trc18650_case_study_final.md"},
        {"result_id": "U2", "domain": "transfer", "experiment": "Panasonic model to clean Uniandes TRC cycle", "role": "hierarchical transfer diagnostic", "metric": "backbone / full-transfer RMSE", "value": "188.882 / 82.134", "unit": "mV", "status": "INCOMPATIBLE DIRECT TRANSFER", "limitation": "Provisional time; static, mapping, and dynamic mismatch; positive-surface support loss", "source_file": "thesis_data/own_model_transfer_audit.md"},
    ]


def write_results_register(rows: list[dict[str, object]]) -> None:
    write_csv(DATA / "final_results_register.csv", rows)
    lines = [
        "# Final results register",
        "",
        "Accepted values only. Every row states whether it is a calibration result, independent validation, diagnostic, or limitation.",
        "",
        "| ID | Domain | Experiment | Role | Metric | Value | Status | Limitation | Source |",
        "|---|---|---|---|---|---:|---|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['result_id']} | {row['domain']} | {row['experiment']} | {row['role']} | "
            f"{row['metric']} | {row['value']} {row['unit']} | {row['status']} | {row['limitation']} | `{row['source_file']}` |"
        )
    lines.extend(
        [
            "",
            "## Frozen parameter context",
            "",
            "- Panasonic electrical: `R0=0.0259594900586 ohm`, `alpha_p=3204.22765511 s`, `d_n=0.000183029692927 1/s`; fixed `alpha_n=3000 s`, `d_p=0.0007 1/s`.",
            "- Panasonic capacities: `Q_cell=2.994979138519 Ah`, `Q_n=23008.732269258333 C`, `Q_p=14780.898652565113 C`.",
            "- Panasonic thermal: `Cth=117.52317800099644 J/K`, `Gth=0.10911208974639978 W/K`.",
            "- Uniandes TRC own-cell capacity reference: `2.8999 Ah`; it is never mixed with the Panasonic capacity.",
        ]
    )
    (DATA / "final_results_register.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_equations() -> None:
    content = r"""# Final implemented model equations

This document reports the equations as implemented in the accepted code. It distinguishes direct code transcription, conceptual interpretation, literature constitutive input, and project-specific calibration. The repository-wide internal sign convention is **positive current for discharge** (`I > 0` discharge, `I < 0` charge). Raw experimental signs are converted only at loader boundaries.

## Symbols and units

| Symbol | Meaning | Unit |
|---|---|---:|
| `I` | internal discharge-positive cell current | A |
| `t`, `Delta t` | time, interval length | s |
| `x_k,avg`, `x_k,aux`, `x_k,surf` | normalized solid concentration states for electrode `k in {n,p}` | 1 |
| `alpha_k` | grouped solid-diffusion time scale (`R_k^2/D_k` interpretation) | s |
| `Q_k` | grouped electrode capacity scale | C |
| `d_k` | grouped kinetic-rate parameter | 1/s |
| `U_n`, `U_p`, `U_eq`, `V` | electrode OCPs, cell equilibrium voltage, terminal voltage | V |
| `eta_n`, `eta_p` | grouped electrode overpotentials | V |
| `R0` | lumped series resistance | ohm |
| `R`, `F` | gas constant, Faraday constant | J/(mol K), C/mol |
| `SOC_cell` | independent cell-level coulomb-counted SOC | 1 |
| `Q_cell,Ah` | cell capacity used for SOC | Ah |
| `Qdot` | heat rate, positive into the cell | W |
| `T`, `T_amb` | predicted cell and prescribed ambient temperature | K |
| `Cth`, `Gth` | effective thermal capacitance and conductance | J/K, W/K |

## Grouped solid-state dynamics

State order in `models/grouped_spm.py::continuous_state_matrices` is `[x_n,avg, x_n,aux, x_p,avg, x_p,aux]`. With discharge-positive current:

```math
\dot{x}_{n,avg}=-I/Q_n
```

```math
\dot{x}_{n,aux}=\frac{30}{\alpha_n}(x_{n,avg}-x_{n,aux})-\frac{19I}{7Q_n}
```

```math
\dot{x}_{p,avg}=+I/Q_p
```

```math
\dot{x}_{p,aux}=\frac{30}{\alpha_p}(x_{p,avg}-x_{p,aux})+\frac{19I}{7Q_p}
```

Surface outputs are algebraic:

```math
x_{n,surf}=x_{n,aux}-\frac{\alpha_n I}{105Q_n},\qquad
x_{p,surf}=x_{p,aux}+\frac{\alpha_p I}{105Q_p}.
```

These formulas are direct code transcription. Conceptually, the code describes them as a two-state-per-electrode parabolic diffusion reduction, not a radial finite-difference discretization. **REDUCTION METHOD PROVENANCE REQUIRES MANUAL LITERATURE VERIFICATION**. A final thesis citation is required specifically for the reduction that produces the numerical factors `30`, `19/7`, and `1/105`, and for the interpretation `alpha_k=R_k^2/D_k`; the local repository establishes software provenance but not a complete literature derivation.

The linear subsystem is propagated exactly for left-sample zero-order-held current on every actual interval. `models/grouped_spm.py::exact_zoh_matrices` evaluates

```math
\begin{bmatrix}x_{k+1}\\1\end{bmatrix}
=\exp\!\left(\begin{bmatrix}A&B\\0&0\end{bmatrix}\Delta t_k\right)
\begin{bmatrix}x_k\\I_k\end{bmatrix}.
```

## Grouped kinetics and terminal voltage

`models/grouped_spm.py::grouped_overpotentials` implements

```math
\eta_k=\frac{2RT}{F}\operatorname{asinh}\!\left[
\frac{I}{6Q_kd_k\sqrt{x_{k,surf}(1-x_{k,surf})}}
\right],\quad k\in\{n,p\}.
```

`models/grouped_spm.py::simulate_grouped_spm` then evaluates

```math
U_{eq}=U_p(x_{p,surf})-U_n(x_{n,surf}),
```

```math
V=U_{eq}+\eta_p-\eta_n-R_0I.
```

These are direct implemented equations. `R0`, `alpha_p`, and `d_n` are the accepted conditional HPPC calibration quantities; `alpha_n` and `d_p` are fixed. Their accepted values are project-specific effective parameters, not universally identified material properties.

## Frozen Panasonic OCP representation

The Panasonic-facing functions are in `models/ocp_panasonic_qi_frozen.py` and enforce their supported domains without clipping or extrapolation.

```math
U_{eq}(z)=U_p[x_p(z)]-U_n[x_n(z)]
```

with nominal effective trajectories

```math
x_n(z)=x_{n,start}+z(x_{n,end}-x_{n,start}),\qquad
x_p(z)=x_{p,start}+z(x_{p,end}-x_{p,start}).
```

The negative OCP is the fixed Guo–White (2005) MCMB discharge spline implemented in `models/ocp_panasonic_qi.py::U_n`, supported on `0.01 <= x_n <= 0.80`. The positive OCP is a project-specific monotonic piecewise-linear 13-knot representation implemented by `models/ocp_panasonic_qi.py::U_p`, supported on `0.18 <= x_p <= 0.96`. Its knot voltages and nominal windows were identified from the Panasonic C/20 curve and frozen after selecting multistart 3. The endpoints are **nominal effective stoichiometric windows**, not uniquely identified physical stoichiometries.

## Independent cell SOC

`models/electrothermal_spm.py::_cell_soc_trajectory` implements

```math
SOC_{cell,k+1}=SOC_{cell,k}-\frac{I_k\Delta t_k}{3600Q_{cell,Ah}}.
```

`Q_cell,Ah` is independent of `Q_n` and `Q_p`. No SOC clipping is used. Panasonic uses `2.994979138519 Ah`; the own TRC cell uses `2.8999 Ah` and the values are never mixed.

## Heat generation

For discharge-positive current, `models/thermal.py::calculate_heat_generation` implements

```math
\dot Q_{irr}=I(U_{eq}-V),
```

```math
\dot Q_{ohmic}=I(R_0I)=R_0I^2,\qquad
\dot Q_{kinetic}=I(\eta_n-\eta_p),
```

```math
\dot Q_{irr}=\dot Q_{ohmic}+\dot Q_{kinetic},
```

```math
\dot Q_{rev}=-IT\frac{dU_{eq}}{dT},\qquad
\dot Q_{gen}=\dot Q_{irr}+\dot Q_{rev}.
```

The generic coupled wrapper uses surface-state `U_eq` and predicted terminal voltage. The accepted Panasonic thermal calibration/validation scripts deliberately use average-state equilibrium voltage,

```math
U_{eq,avg}=U_p(x_{p,avg})-U_n(x_{n,avg}),
```

and retain two explicitly labeled heat paths:

```math
\dot Q_{irr,measured}=I(U_{eq,avg}-V_{measured}),\qquad
\dot Q_{irr,model}=I(U_{eq,avg}-V_{model}).
```

The full-cell entropy coefficient is a literature constitutive input implemented in `models/panasonic_entropic.py::dUeq_dT_panasonic`, linearly interpolated only on `0.10 <= SOC <= 0.90`. Its code provenance is He et al. (2024) Supplementary Information Table S2, citing Lu et al. (2019). No entropy extrapolation or clipping is used.

## One-node thermal model

`models/thermal.py` implements

```math
C_{th}\frac{dT}{dt}=\dot Q_{gen}-G_{th}(T-T_{amb}).
```

For constant left-sample heat and ambient over one interval and `Gth > 0`, `models/thermal.py::lumped_thermal_step` uses the exact analytical update

```math
T_{k+1}=T_{amb,k}+\frac{\dot Q_{gen,k}}{G_{th}}+
\left(T_k-T_{amb,k}-\frac{\dot Q_{gen,k}}{G_{th}}\right)
\exp\!\left(-\frac{G_{th}\Delta t_k}{C_{th}}\right).
```

For `Gth=0`,

```math
T_{k+1}=T_k+\dot Q_{gen,k}\Delta t_k/C_{th}.
```

The accepted Panasonic `Cth` and `Gth` are project-specific effective one-node parameters calibrated on the 25 °C 1C discharge over SOC 0.10–0.90 using measured-voltage irreversible heat plus literature reversible heat. They are not independently measured material properties.

## Evaluation metrics

For residual `e_i = y_pred,i - y_meas,i`, the reporting scripts use

```math
RMSE=\sqrt{\frac{1}{N}\sum_i e_i^2},\qquad
MAE=\frac{1}{N}\sum_i |e_i|,
```

```math
bias=\frac{1}{N}\sum_i e_i,
```

```math
R^2=1-\frac{\sum_i e_i^2}{\sum_i(y_{meas,i}-\bar y_{meas})^2}.
```

RMSE is the primary error-scale summary, MAE describes typical absolute deviation, signed bias shows systematic over/underprediction, and `R²` is reported only with the variance context of the measured subset. A negative `R²` is possible for low-variance subsets such as CV plateaus and does not replace RMSE/bias interpretation.

The data-quality diagnostic is

```math
R_{app}=\left|\frac{\Delta V}{\Delta I}\right|.
```

`R_app` is an adjacent rest-to-load experimental-chain diagnostic. It must not be interpreted as true cell internal resistance or as the model parameter `R0`.

## Source-status summary

| Component | Status | Code source | Thesis citation action |
|---|---|---|---|
| Four grouped solid states and surface outputs | implemented directly | `models/grouped_spm.py` | manually verify/cite reduction derivation |
| Butler–Volmer-style grouped kinetics and voltage | implemented directly; inherited software architecture | `models/grouped_spm.py` | cite theoretical kinetics and software provenance |
| Panasonic MCMB negative OCP | literature constitutive input | `models/ocp_panasonic_qi.py` | Guo–White (2005) |
| Panasonic positive OCP | project-specific identified constitutive curve | `models/ocp_panasonic_qi_frozen.py` | cite this thesis identification workflow |
| Entropy coefficient | literature constitutive input | `models/panasonic_entropic.py` | He et al. (2024) Table S2; Lu et al. (2019) as cited there |
| One-node energy balance and exact update | implemented model equation | `models/thermal.py` | cite lumped thermal modeling convention if required by thesis style |
| Electrical/thermal parameter values | project-specific calibration | frozen Panasonic modules and reports | cite this thesis calibration chapters |
"""
    (DATA / "final_model_equations.md").write_text(content, encoding="utf-8")


def write_narrative_and_objectives() -> None:
    narrative = """# Final results narrative

## Model contribution

The implemented framework couples a four-state grouped single-particle electrical model, independent cell-level coulomb counting, electrochemical heat generation, and an analytically propagated one-node thermal state. Current is internally discharge-positive, actual nonuniform sample intervals are retained, and invalid SOC/OCP states are reported rather than clipped. The final Panasonic parameterization combines a fixed MCMB negative-electrode OCP, a frozen 13-knot positive-electrode OCP, nominal effective electrode windows, conditional HPPC dynamic parameters, a literature full-cell entropy coefficient, and effective one-node thermal parameters.

## Calibration findings

The Panasonic C/20 representation is accepted on its declared 5–95% fit interval with 1.086 mV RMSE, while the electrode windows remain weakly identifiable. Dynamic calibration was deliberately conditional: only HPPC pulses beginning at SOC at or above one third entered the final three-parameter fit, giving 36.835 mV phase-balanced RMSE. The fitted `R0`, `alpha_p`, and `d_n` are therefore not claimed as globally identified Panasonic constants.

The primary thermal pair was calibrated on the 25 °C 1C discharge over SOC 0.10–0.90 using measured-voltage irreversible heat and literature reversible heat. It achieved 0.143622 K RMSE and `R²=0.964702`. These are effective lumped parameters under the tested boundary condition, not independently measured material properties.

## Validation findings

Withheld low-SOC HPPC transfer is limited: phase-balanced RMSE rises to 85.729 mV and the dominant error is the loaded transient. The independent 1C discharge shows the same direction of deterioration, with 114.362 mV overall RMSE and 195.954 mV at low SOC. UDDS voltage transfer is stronger overall at 30.988 mV RMSE, but its residual still grows from 14.361 mV at high SOC to 52.785 mV at low SOC and large regenerative current steps are overpredicted. CC–CV charge transfer remains limited at 61.034 mV RMSE with −47.526 mV bias, showing incomplete transfer across current direction.

Thermally, the one-node model underpredicts the independent UDDS temperature trajectory (0.526698 K RMSE, −0.512715 K bias). On CC–CV charge, the measured-voltage heat path gives 0.427781 K RMSE, whereas the fully predictive path degrades to 0.864391 K. The 436.438 J electrical heat deficit explains most of this added coupled error; the retained reversible contribution is 213.914 J. These results separate limitations of electrical heat prediction from limitations of the thermal state model.

## Why the one-node thermal model was retained

The clean two-node model marginally improves the fitted 1C RMSE (0.138221 K versus 0.143622 K) but does not improve independent UDDS transfer (0.532526 K versus 0.526698 K). Its surface capacitance reaches the numerical lower guard, its fast mode is about `1e-4 s`—far below the experiment’s temporal resolution—and its local Jacobian condition number is about `2.56e7`. The extra state is therefore not practically recoverable from the available surface-temperature record. The one-node architecture is retained because it is the more parsimonious identifiable representation, not because all thermal transfers are accurate.

## Uniandes data-quality lesson

The cleanest ambient cell-1 cycle has median `R_app≈0.03588 ohm`; affected cycles and pulse records lie around 0.178–0.232 ohm. This reproducible separation is evidence of elevated series-path response in the complete experimental chain, not a measurement of intrinsic cell resistance. The affected records are useful nonideality case studies but are rejected as electrical-identification truth. The cleanest cycle remains limited by an absent measured timebase and by its loaded 1C/CC–CV nature.

## Transfer findings

Literal Panasonic-to-TRC transfer is incompatible. The equilibrium-voltage backbone differs from the loaded TRC reference by 188.882 mV RMSE, and the provisional-time full transfer gives 82.134 mV RMSE only through compensating dynamic terms while losing positive-electrode OCP support near low SOC. The minimum defensible adaptation is cell-specific capacity/state mapping plus a TRC-compatible voltage-reference/OCP basis. Dynamic and thermal re-identification would also require measured-time evidence that is not available.

## Thesis contribution

The central contribution is not a claim of universal parameter transfer. It is a reproducible evidence hierarchy: audit the experiment, separate calibration from withheld validation, preserve model-domain checks, diagnose transfer before adaptation, and report unsupported identification as a limitation. This converts imperfect experimental data into defensible conclusions without silently correcting voltage, clipping states, or promoting diagnostics to physical parameters.

## Remaining limitations

- The grouped diffusion-reduction literature derivation still requires manual citation verification.
- Panasonic electrode windows and positive OCP knots are effective, weakly identifiable representations.
- Dynamic parameters are conditional on HPPC pulses with start SOC at or above one third.
- Low-SOC and current-direction transfer remain limited.
- Thermal calibration uses measured-voltage heat; fully predictive heat inherits electrical-model error.
- The Uniandes files have no measured timebase, preventing defensible dynamic or thermal identification.
"""
    (DATA / "final_results_narrative.md").write_text(narrative, encoding="utf-8")

    objectives = """# Objective-to-evidence matrix

| Thesis objective | Evidence produced | Supporting figures | Supporting metrics | Remaining caveat |
|---|---|---|---|---|
| 1. Formulate the electrochemical–thermal model. | Frozen grouped SPM, independent cell SOC, heat generation, reversible heat option, and exact one-node thermal propagation; equation/code provenance register. | `panasonic_ocp_fit_final`, `panasonic_1c_thermal_calibration_final`, `thermal_model_selection_summary` | OCP fit 1.086 mV; primary thermal fit 0.143622 K; accepted frozen parameter set | Diffusion-reduction derivation needs manual literature verification; parameters are effective and chemistry-specific. |
| 2. Build a methodology for assessment from experimental data. | Reproducible quality audit, suitability gate, conditional calibration/validation separation, domain checks, transfer hierarchy, and no-automatic-fit policy. | `battery_dataset_assessment_flow_final`, `uniandes_apparent_series_path_final`, `panasonic_hppc_transfer_final` | Clean `R_app≈0.03588 ohm`; affected `0.178–0.232 ohm`; HPPC calibration/withheld 36.835/85.729 mV | `R_app` is a setup diagnostic, not intrinsic `R0`; acceptance thresholds remain experiment-specific. |
| 3. Determine application-scenario conditions. | Forward validation across steady 1C discharge, bidirectional UDDS, CC–CV charge, low-SOC HPPC, and thermal transfer; limits tied to SOC, direction, and heat path. | `panasonic_1c_voltage_validation_final`, `panasonic_udds_voltage_validation_final`, `panasonic_cccv_voltage_transfer_final`, three thermal-transfer figures | 1C 114.362 mV; UDDS 30.988 mV; charge 61.034 mV; UDDS thermal 0.526698 K; charge A/B 0.427781/0.864391 K | No universal accuracy claim; low SOC, large regenerative steps, and charge direction are limited. |
| 4. Validate with Universidad de los Andes data. | Own-data quality classification, corrected discharge-only mapping, affected-data nonideality case study, and hierarchical Panasonic-to-TRC transfer audit. | `uniandes_defective_contacts_discharge_final`, `uniandes_contact_arrangements_comparison_final`, `panasonic_to_trc_transfer_final`, `uniandes_apparent_series_path_final` | Backbone/full transfer 188.882/82.134 mV; affected series-path elevation ≈5–6.5× | Own files lack measured time; the city/air labels are presentation context and do not reassign charge, capacity, current, or thermal samples. |
"""
    (DATA / "objective_evidence_matrix.md").write_text(objectives, encoding="utf-8")


def write_presentation_assets(rows: list[dict[str, object]]) -> None:
    equations = r"""# Main equations for the 20-minute presentation

Internal convention: **`I > 0` discharge; `I < 0` charge**.

```math
\dot x_{n,avg}=-I/Q_n,\qquad \dot x_{p,avg}=+I/Q_p
```

```math
\dot x_{k,aux}=\frac{30}{\alpha_k}(x_{k,avg}-x_{k,aux})+\text{signed current forcing}
```

```math
x_{n,surf}=x_{n,aux}-\frac{\alpha_n I}{105Q_n},\qquad
x_{p,surf}=x_{p,aux}+\frac{\alpha_p I}{105Q_p}
```

```math
\eta_k=\frac{2RT}{F}\operatorname{asinh}\left(\frac{I}{6Q_kd_k\sqrt{x_{k,surf}(1-x_{k,surf})}}\right)
```

```math
U_{eq}=U_p-U_n,\qquad V=U_{eq}+\eta_p-\eta_n-R_0I
```

```math
\dot Q_{irr}=I(U_{eq}-V),\qquad \dot Q_{rev}=-IT\frac{dU_{eq}}{dT}
```

```math
C_{th}\dot T=\dot Q_{irr}+\dot Q_{rev}-G_{th}(T-T_{amb})
```

The electrical state and thermal state are propagated with exact left-sample zero-order-hold updates on the actual nonuniform timestamps. No state/OCP clipping is used.
"""
    (PRESENTATION / "main_equations.md").write_text(equations, encoding="utf-8")

    shortlist = """# Presentation figure shortlist

| Slide purpose | Recommended asset | Key spoken message |
|---|---|---|
| Model and OCP basis | `panasonic_ocp_fit_final.png` | Excellent full-cell representation does not uniquely identify electrode windows. |
| Conditional electrical calibration | `panasonic_hppc_transfer_final.png` | Calibration is credible only in its declared SOC region; withheld low-SOC loaded transients degrade. |
| Independent electrical validation | `panasonic_udds_voltage_validation_final.png` | Strong overall drive-cycle accuracy coexists with structured low-SOC/current-step limitations. |
| Current-direction test | `panasonic_cccv_voltage_transfer_final.png` | Discharge-derived electrical parameters do not transfer completely to CC–CV charge. |
| Primary thermal result | `panasonic_1c_thermal_calibration_final.png` | One-node parameters reproduce the calibration discharge closely. |
| Thermal transfer and coupling | `panasonic_cccv_thermal_transfer_final.png` | Electrical voltage error becomes a heat-input error in the coupled prediction. |
| Architecture choice | `thermal_model_selection_summary.png` | The two-node model is not recoverable and does not improve independent UDDS transfer. |
| Own-data lesson | `uniandes_apparent_series_path_final.png` | Experimental-chain quality must be audited before electrical identification. |
| Cross-cell transfer | `panasonic_to_trc_transfer_final.png` | Direct transfer fails; minimum adaptation is capacity/state mapping plus cell-compatible voltage reference. |
| General methodology / conclusion | `battery_dataset_assessment_flow_final.png` | Validation is an evidence-gated endpoint, not an automatic consequence of having data. |

Suggested 20-minute allocation: 2 min problem/objectives; 3 min model/equations; 4 min Panasonic electrical calibration/validation; 3 min thermal calibration/transfer; 3 min architecture selection; 3 min Uniandes quality/transfer; 2 min contributions and limitations.
"""
    (DATA / "presentation_figure_shortlist.md").write_text(shortlist, encoding="utf-8")

    # Compact presentation table.
    table_rows = [
        ["C/20 OCP", "fit", "1.086 mV", "accepted; windows weakly identifiable"],
        ["HPPC", "cal / withheld", "36.835 / 85.729 mV", "low-SOC transfer limited"],
        ["1C discharge", "validation", "114.362 mV", "low SOC / cutoff limited"],
        ["UDDS", "validation", "30.988 mV", "low SOC & large steps limited"],
        ["CC–CV charge", "validation", "61.034 mV", "directional transfer limited"],
        ["1C thermal", "calibration", "0.143622 K", "primary one-node fit"],
        ["UDDS thermal", "validation A", "0.526698 K", "systematic underprediction"],
        ["CC–CV thermal", "validation A / B", "0.427781 / 0.864391 K", "coupled electrical-heat failure"],
        ["Uniandes transfer", "diagnostic", "188.882 / 82.134 mV", "direct transfer incompatible"],
    ]
    fig, ax = plt.subplots(figsize=(14, 7.875))
    ax.axis("off")
    table = ax.table(
        cellText=table_rows,
        colLabels=["Dataset", "Evidence role", "Headline metric", "Decision / limitation"],
        colWidths=[0.18, 0.17, 0.22, 0.43],
        cellLoc="left",
        loc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1, 2.05)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("white")
        if row == 0:
            cell.set_facecolor(BLUE)
            cell.get_text().set_color("white")
            cell.get_text().set_weight("bold")
        else:
            cell.set_facecolor("#F2F6FA" if row % 2 else "white")
            cell.get_text().set_color(INK)
    ax.set_title("Master results: accepted evidence and limitations", fontsize=22, weight="bold", color=INK, pad=20)
    fig.savefig(PRESENTATION / "master_results_table.png", dpi=240, bbox_inches="tight")
    plt.close(fig)

    master_rows = [
        {"dataset": row[0], "evidence_role": row[1], "headline_metric": row[2], "decision_or_limitation": row[3]}
        for row in table_rows
    ]
    write_csv(DATA / "master_results_table.csv", master_rows)


def write_manifest() -> None:
    stems = [
        ("panasonic_ocp_fit_final", "Panasonic OCP identification", "C/20 fit and nominal effective windows"),
        ("panasonic_hppc_transfer_final", "Conditional dynamic calibration/validation", "HPPC fitted versus withheld low-SOC pulses"),
        ("panasonic_1c_voltage_validation_final", "Electrical validation", "Independent continuous 1C voltage transfer"),
        ("panasonic_udds_voltage_validation_final", "Electrical validation", "Independent bidirectional UDDS voltage transfer"),
        ("panasonic_cccv_voltage_transfer_final", "Electrical validation", "Independent CC–CV charge transfer"),
        ("panasonic_directional_voltage_separation_final", "Directional diagnostic", "C/20 charge/discharge separation; not thermodynamic hysteresis"),
        ("panasonic_1c_thermal_calibration_final", "Thermal calibration", "Primary one-node 1C temperature fit"),
        ("panasonic_udds_thermal_transfer_final", "Thermal validation", "UDDS one-node forward transfer"),
        ("panasonic_cccv_thermal_transfer_final", "Electrothermal validation", "Charge measured-heat and fully predictive paths"),
        ("thermal_model_selection_summary", "Model selection", "One-node versus clean two-node accuracy and recoverability"),
        ("uniandes_apparent_series_path_final", "Own-data quality", "Relative apparent series-path response"),
        ("uniandes_defective_contacts_discharge_final", "Own-data case study", "Defective-contact discharge curves under three measured-air intervals"),
        ("uniandes_contact_arrangements_comparison_final", "Own-data case study", "Reference versus defective contact arrangements at Aire 22–24 °C"),
        ("uniandes_cccv_external_conditions_final", "Own-data charge comparison", "CC-CV charge under three measured-air conditions with external setup resistance"),
        ("uniandes_cccv_contact_comparison_final", "Own-data charge comparison", "Reference versus external-resistance contact arrangements during CC-CV charge"),
        ("uniandes_temperature_discharge_final", "Own-data thermal results", "Measured Temp1, Temp2, and Temp7 during discharge"),
        ("uniandes_temperature_charge_final", "Own-data thermal results", "Measured Temp1, Temp2, and Temp7 during CC-CV charge"),
        ("panasonic_to_trc_transfer_final", "Cross-cell transfer diagnostic", "Literal benchmark-to-own-cell transfer"),
        ("battery_dataset_assessment_flow_final", "Reusable methodology", "Evidence-gated assessment and transfer flow"),
        ("uniandes_geographical_context_final", "Experimental context", "Colombian localities associated with measured-air intervals"),
    ]
    rows = []
    for index, (stem, role, purpose) in enumerate(stems, start=1):
        rows.append(
            {
                "figure_id": f"F{index:02d}",
                "stem": stem,
                "role": role,
                "purpose": purpose,
                "thesis_pdf": f"thesis_figures/final/{stem}.pdf",
                "thesis_png": f"thesis_figures/final/{stem}.png",
                "presentation_pdf": f"presentation_assets/final/{stem}.pdf",
                "presentation_png": f"presentation_assets/final/{stem}.png",
                "generated_by": "scripts/build_final_thesis_package.py",
                "data_policy": "accepted artifacts and frozen forward-only reproduction; no fit",
            }
        )
    write_csv(DATA / "final_figure_manifest.csv", rows)


def verify_outputs() -> None:
    required = [
        DATA / "final_model_equations.md",
        DATA / "final_results_register.csv",
        DATA / "final_results_register.md",
        DATA / "presentation_figure_shortlist.md",
        DATA / "master_results_table.csv",
        DATA / "final_results_narrative.md",
        DATA / "objective_evidence_matrix.md",
        DATA / "final_figure_manifest.csv",
        PRESENTATION / "master_results_table.png",
        PRESENTATION / "main_equations.md",
    ]
    manifest = csv_rows(DATA / "final_figure_manifest.csv")
    for row in manifest:
        required.extend(ROOT / row[key] for key in ("thesis_pdf", "thesis_png", "presentation_pdf", "presentation_png"))
    missing = [path for path in required if not path.is_file() or path.stat().st_size == 0]
    if missing:
        raise RuntimeError(f"Missing or empty final outputs: {missing}")


def main() -> None:
    configure_style()
    THESIS_FIGURES.mkdir(parents=True, exist_ok=True)
    PRESENTATION.mkdir(parents=True, exist_ok=True)

    figure_ocp_fit()
    figure_hppc_transfer()
    electrical = electrical_forward_arrays()
    figure_1c(electrical["one_c"])
    figure_udds(electrical["udds"])
    figure_cccv(electrical["cccv"])
    figure_directional_separation()

    thermal = thermal_forward_arrays()
    figure_thermal_1c(thermal["one_c"])
    figure_thermal_udds(thermal["udds"])
    figure_thermal_cccv(thermal["cccv"])
    figure_thermal_selection()

    figure_uniandes_rapp()
    soc, curves, experiments = own_case_arrays()
    figure_uniandes_discharge_comparisons(experiments)
    figure_uniandes_charge_temperature()
    figure_transfer(soc, curves, experiments)
    figure_assessment_flow()
    figure_uniandes_geographical_context()

    rows = result_register()
    write_results_register(rows)
    write_equations()
    write_narrative_and_objectives()
    write_presentation_assets(rows)
    write_manifest()
    verify_outputs()
    print("Final thesis and presentation package built from accepted results only.")


if __name__ == "__main__":
    main()
