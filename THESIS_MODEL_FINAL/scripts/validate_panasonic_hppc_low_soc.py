"""Forward-only validation of the nominal conditional HPPC fit below 1/3 SOC."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from scripts.audit_panasonic_hppc_identifiability import (  # noqa: E402
    load_hppc_data,
    load_pulse_records,
)
from scripts.calibrate_panasonic_hppc_primary import (  # noqa: E402
    ALPHA_N_FIXED_S,
    D_P_FIXED_PER_S,
    PHASES,
    PreparedGroup,
    PreparedPulse,
    calibration_diagnostics,
    evaluate,
    positive_surface_state,
    prepare_groups,
)


OUTPUT_DIRECTORY = REPOSITORY_ROOT / "thesis_data"
SOC_CUTOFF = 1.0 / 3.0
CALIBRATION_REGION_RMSE_V = 0.036835197638870915
NOMINAL_CONDITIONAL_PARAMETERS = {
    "R0": 0.0259594900586,
    "alpha_p": 3204.22765511,
    "d_n": 0.000183029692927,
}


def surface_ranges(groups: list[PreparedGroup]) -> dict[str, float | bool]:
    x_n = np.concatenate([group.x_n_surface for group in groups])
    x_p = np.concatenate(
        [
            positive_surface_state(
                group, NOMINAL_CONDITIONAL_PARAMETERS["alpha_p"]
            )
            for group in groups
        ]
    )
    violation = bool(
        np.any((x_n < 0.01) | (x_n > 0.80))
        or np.any((x_p < 0.18) | (x_p > 0.96))
    )
    return {
        "x_n_surface_min": float(np.min(x_n)),
        "x_n_surface_max": float(np.max(x_n)),
        "x_p_surface_min": float(np.min(x_p)),
        "x_p_surface_max": float(np.max(x_p)),
        "ocp_domain_violation": violation,
    }


def write_metrics(
    path: Path,
    pulse_count: int,
    history_count: int,
    diagnostics: dict[str, object],
    ranges: dict[str, float | bool],
) -> None:
    metrics = dict(diagnostics["metrics"])
    rows: list[dict[str, object]] = [
        {"group": "data", "subset": "low_soc", "metric": "pulse_count", "value": pulse_count},
        {"group": "data", "subset": "low_soc", "metric": "history_count", "value": history_count},
    ]
    for name, value in metrics.items():
        rows.append({"group": "overall", "subset": "low_soc", "metric": name, "value": value})
    for rate, value in dict(diagnostics["by_rate"]).items():
        rows.append(
            {
                "group": "c_rate",
                "subset": rate,
                "metric": "phase_balanced_dynamic_rmse_v",
                "value": value,
            }
        )
    rows.append(
        {
            "group": "diagnostic",
            "subset": "low_soc",
            "metric": "voltage_step_bias_v",
            "value": diagnostics["step_bias_v"],
        }
    )
    for name, value in ranges.items():
        rows.append({"group": "state", "subset": "low_soc", "metric": name, "value": value})
    rows.append(
        {
            "group": "comparison",
            "subset": "calibration_soc_ge_one_third",
            "metric": "phase_balanced_dynamic_rmse_v",
            "value": CALIBRATION_REGION_RMSE_V,
        }
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_validation_plot(
    groups: list[PreparedGroup], group_voltages: tuple[np.ndarray, ...], path: Path
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    candidates: list[tuple[PreparedGroup, PreparedPulse, np.ndarray]] = []
    for group, voltage in zip(groups, group_voltages, strict=True):
        for pulse in group.pulses:
            candidates.append((group, pulse, voltage))
    targets = [
        ("0.5C", 0.30),
        ("1C", 0.25),
        ("2C", 0.20),
        ("4C", 0.17),
        ("6C", 0.30),
        ("1C", 0.08),
    ]
    selected: list[tuple[PreparedGroup, PreparedPulse, np.ndarray]] = []
    for rate, target_soc in targets:
        eligible = [item for item in candidates if item[1].record.nominal_rate == rate]
        if not eligible:
            continue
        choice = min(eligible, key=lambda item: abs(item[1].record.soc - target_soc))
        if all(choice[1].record.pulse_id != item[1].record.pulse_id for item in selected):
            selected.append(choice)

    figure, axes = plt.subplots(2, 3, figsize=(13.2, 7.6), constrained_layout=True)
    for axis, (group, pulse, voltage) in zip(axes.flat, selected, strict=False):
        pre = pulse.pre_local_index
        stop = int(
            np.searchsorted(group.time_s, group.time_s[pre] + 80.0, side="right")
        )
        indices = np.arange(pre, min(stop, group.time_s.size))
        time = group.time_s[indices] - group.time_s[pre]
        measured_delta = (
            group.measured_voltage_v[indices] - group.measured_voltage_v[pre]
        )
        model_delta = voltage[indices] - voltage[pre]
        axis.plot(time, measured_delta, label="measured")
        axis.plot(time, model_delta, label="forward model")
        pulse_start = pulse.phase_local_indices["immediate"][0]
        axis.axvspan(
            group.time_s[pulse_start] - group.time_s[pre],
            group.time_s[pulse.end_local_index] - group.time_s[pre],
            color="tab:red",
            alpha=0.08,
        )
        axis.set(
            title=(
                f"Pulse {pulse.record.pulse_id}: {pulse.record.nominal_rate}, "
                f"{100*pulse.record.soc:.1f}% SOC"
            ),
            xlabel="Time from pre-pulse [s]",
            ylabel="Delta voltage [V]",
        )
        axis.grid(alpha=0.25)
    for axis in axes.flat[len(selected) :]:
        axis.set_visible(False)
    if selected:
        axes.flat[0].legend()
    figure.suptitle("Forward-only low-SOC HPPC validation")
    figure.savefig(path, dpi=180)
    plt.close(figure)


def save_residual_plot(records: tuple[dict[str, object], ...], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 3, figsize=(14.0, 4.8), constrained_layout=True)
    for axis, phase in zip(axes, PHASES, strict=True):
        selected = [row for row in records if row["phase"] == phase]
        stride = max(1, len(selected) // 5000)
        shown = selected[::stride]
        axis.scatter(
            [100.0 * float(row["soc"]) for row in shown],
            [1000.0 * float(row["dynamic_residual_v"]) for row in shown],
            s=6,
            alpha=0.35,
        )
        axis.axhline(0.0, color="black", linewidth=0.8)
        axis.set(
            title=phase.replace("_", " "),
            xlabel="Pulse-start SOC [%]",
            ylabel="Dynamic residual [mV]",
        )
        axis.grid(alpha=0.25)
    figure.suptitle("Withheld low-SOC forward-validation residuals")
    figure.savefig(path, dpi=180)
    plt.close(figure)


def write_report(
    path: Path,
    pulse_count: int,
    history_count: int,
    diagnostics: dict[str, object],
    ranges: dict[str, float | bool],
) -> None:
    metrics = dict(diagnostics["metrics"])
    transfer_ratio = metrics["phase_balanced_rmse_v"] / CALIBRATION_REGION_RMSE_V
    phase_values = {
        "instantaneous step": metrics["rmse_immediate_v"],
        "loaded transient": metrics["rmse_within_pulse_v"],
        "relaxation": metrics["rmse_post_relaxation_v"],
        "absolute baseline voltage": metrics["absolute_rmse_v"],
    }
    largest_name = max(phase_values, key=phase_values.get)
    classification = (
        "LOW-SOC TRANSFER ACCEPTABLE"
        if transfer_ratio <= 1.25
        and not bool(ranges["ocp_domain_violation"])
        and max(phase_values.values()) <= 0.075
        else "LOW-SOC TRANSFER LIMITED"
    )
    lines = [
        "# Panasonic HPPC low-SOC forward validation",
        "",
        f"## {classification}",
        "",
        "Forward evaluation only. The nominal conditional parameters from the `SOC >= 1/3` diagnostic are used unchanged: `R0=0.0259594900586 ohm`, `alpha_p=3204.22765511 s`, and `d_n=0.000183029692927 1/s`. `alpha_n=3000 s` and `d_p=0.0007 1/s` remain fixed. No parameter was optimized.",
        "",
        f"- Withheld complete low-SOC pulses: **{pulse_count}** in {history_count} independently initialized histories.",
        "- Eligibility: pulse-start cell `SOC < 1/3`; cutoff-truncated IDs 60, 64, and 67 remain excluded.",
        "- Cumulative Ah initializes each history; omitted gaps are not integrated, raw negative current is converted to discharge-positive, and actual nonuniform time steps are retained.",
        "",
        "## Dynamic metrics",
        "",
        "The same pulse-relative, pulse-balanced, equal-phase metric used for calibration is applied without refitting.",
        "",
        f"- Phase-balanced dynamic RMSE: **{1000*metrics['phase_balanced_rmse_v']:.6f} mV**.",
        f"- Calibration-region (`SOC >= 1/3`) RMSE: {1000*CALIBRATION_REGION_RMSE_V:.6f} mV.",
        f"- Transfer/calibration RMSE ratio: {transfer_ratio:.6f}.",
        "",
        "| Response phase | RMSE [mV] |",
        "|---|---:|",
        f"| Immediate | {1000*metrics['rmse_immediate_v']:.6f} |",
        f"| Within pulse | {1000*metrics['rmse_within_pulse_v']:.6f} |",
        f"| Post relaxation | {1000*metrics['rmse_post_relaxation_v']:.6f} |",
        "",
        "### By C-rate",
        "",
        "| C-rate | Phase-balanced RMSE [mV] |",
        "|---|---:|",
        *[
            f"| {rate} | {1000*value:.6f} |"
            for rate, value in dict(diagnostics["by_rate"]).items()
        ],
        "",
        f"Voltage-step bias (`measured - predicted`): {1000*float(diagnostics['step_bias_v']):+.6f} mV.",
        "",
        "## Absolute-voltage diagnostic",
        "",
        f"- RMSE: {1000*metrics['absolute_rmse_v']:.6f} mV.",
        f"- MAE: {1000*metrics['absolute_mae_v']:.6f} mV.",
        f"- R2: {metrics['absolute_r2']:.9f}.",
        "",
        "## Surface-state and OCP-domain audit",
        "",
        f"- `x_n,surf` range: {float(ranges['x_n_surface_min']):.12f} to {float(ranges['x_n_surface_max']):.12f}.",
        f"- `x_p,surf` range: {float(ranges['x_p_surface_min']):.12f} to {float(ranges['x_p_surface_max']):.12f}.",
        f"- Frozen OCP-domain violation: **{'YES' if ranges['ocp_domain_violation'] else 'NO'}**.",
        "- No state or stoichiometry was clipped.",
        "",
        "## Interpretation",
        "",
        f"The withheld-region error is {transfer_ratio:.3f} times the calibration-region phase-balanced RMSE. The largest diagnostic discrepancy is in the **{largest_name}** ({1000*phase_values[largest_name]:.3f} mV RMSE).",
        "",
        f"**Final classification: {classification}.**",
        "",
        "No repair, refit, additional calibration, parameter release, OCP/support change, model change, UDDS run, or 1C validation was performed.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    data = load_hppc_data()
    records = load_pulse_records()
    groups = prepare_groups(
        data,
        records,
        maximum_pulse_start_soc=SOC_CUTOFF,
    )
    pulse_count = sum(len(group.pulses) for group in groups)
    evaluation = evaluate(groups, NOMINAL_CONDITIONAL_PARAMETERS)
    diagnostics = calibration_diagnostics(groups, evaluation)
    ranges = surface_ranges(groups)

    write_metrics(
        OUTPUT_DIRECTORY / "panasonic_hppc_low_soc_validation_metrics.csv",
        pulse_count,
        len(groups),
        diagnostics,
        ranges,
    )
    write_report(
        OUTPUT_DIRECTORY / "panasonic_hppc_low_soc_validation_report.md",
        pulse_count,
        len(groups),
        diagnostics,
        ranges,
    )
    save_validation_plot(
        groups,
        evaluation.group_voltage_v,
        OUTPUT_DIRECTORY / "panasonic_hppc_low_soc_validation.png",
    )
    save_residual_plot(
        evaluation.dynamic_records,
        OUTPUT_DIRECTORY / "panasonic_hppc_low_soc_validation_residuals.png",
    )
    print(f"Low-SOC validation pulses: {pulse_count}")
    print(
        "Phase-balanced RMSE: "
        f"{1000*float(diagnostics['metrics']['phase_balanced_rmse_v']):.6f} mV"
    )
    print((OUTPUT_DIRECTORY / "panasonic_hppc_low_soc_validation_report.md").resolve())


if __name__ == "__main__":
    main()
