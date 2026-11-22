"""Forward-only validation of the frozen Panasonic model on the initial 25 C 1C discharge.

This script loads the manifest-selected raw MAT file read-only, isolates its
contiguous discharge branch, and evaluates the frozen conditional HPPC
parameterization without fitting, clipping, or changing any model/OCP value.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from models.ocp_panasonic_qi import (  # noqa: E402
    MCMB_SUPPORTED_STOICHIOMETRY_RANGE,
    NCA_WORKING_STOICHIOMETRY_RANGE,
)
from models.ocp_panasonic_qi_frozen import (  # noqa: E402
    Q_CELL_AH,
    X_N_END,
    X_N_START,
    X_P_END,
    X_P_START,
)
from scripts.audit_panasonic_hppc_identifiability import (  # noqa: E402
    Q_N_C as IMPLEMENTED_Q_N_C,
    Q_P_C as IMPLEMENTED_Q_P_C,
    simulate_local_group,
)
from runtime_paths import PANASONIC_1C_PATH, require_data_file  # noqa: E402


RAW_SOURCE_PATH = PANASONIC_1C_PATH
MANIFEST_PATH = REPOSITORY_ROOT / "experimental_data_manifest.md"
OUTPUT_DIRECTORY = REPOSITORY_ROOT / "thesis_data"

Q_N_C = 23008.732269258333
Q_P_C = 14780.898652565113
NOMINAL_PARAMETERS = {
    "R0": 0.0259594900586,
    "alpha_p": 3204.22765511,
    "d_n": 0.000183029692927,
    "alpha_n": 3000.0,
    "d_p": 0.0007,
}

DISCHARGE_CURRENT_THRESHOLD_A = -0.5
CUTOFF_REGION_FRACTION = 0.01
SOC_REGION_MASKS = {
    "high": lambda soc: (soc >= 2.0 / 3.0) & (soc <= 1.0),
    "middle": lambda soc: (soc >= 1.0 / 3.0) & (soc < 2.0 / 3.0),
    "low": lambda soc: (soc >= 0.0) & (soc < 1.0 / 3.0),
}


def _as_vector(value: object, name: str) -> np.ndarray:
    result = np.asarray(value, dtype=float).reshape(-1)
    if result.size == 0 or not np.all(np.isfinite(result)):
        raise ValueError(f"Raw field {name!r} is empty or non-finite.")
    return result


def load_raw_discharge(path: Path = RAW_SOURCE_PATH) -> tuple[tuple[str, ...], dict[str, np.ndarray]]:
    """Load the accepted CSV (or an explicitly supplied MAT) read-only."""

    path = require_data_file(path)
    if path.suffix.lower() == ".csv":
        table = np.genfromtxt(path, delimiter=",", names=True, dtype=float, encoding="utf-8")
        columns = tuple(table.dtype.names or ())
        required = (
            "Time", "Voltage", "Current", "Ah", "Wh", "Power",
            "Battery_Temp_degC", "Chamber_Temp_degC",
        )
        missing = sorted(set(required) - set(columns))
        if missing:
            raise ValueError(f"Processed 1C CSV is missing required fields: {missing}")
        arrays = {name: _as_vector(table[name], name) for name in required}
        if len({array.size for array in arrays.values()}) != 1:
            raise ValueError("Processed 1C CSV channels have inconsistent lengths.")
        return columns, arrays
    contents = loadmat(path, squeeze_me=True, struct_as_record=False)
    measurement = contents.get("meas")
    if measurement is None or not hasattr(measurement, "_fieldnames"):
        raise ValueError("Expected a MATLAB structure named 'meas'.")
    columns = tuple(measurement._fieldnames)
    required = (
        "Time",
        "Voltage",
        "Current",
        "Ah",
        "Wh",
        "Power",
        "Battery_Temp_degC",
        "Chamber_Temp_degC",
    )
    missing = sorted(set(required) - set(columns))
    if missing:
        raise ValueError(f"Raw 1C data is missing required fields: {missing}")
    arrays = {name: _as_vector(getattr(measurement, name), name) for name in required}
    lengths = {values.size for values in arrays.values()}
    if len(lengths) != 1:
        raise ValueError("Raw 1C channels have inconsistent lengths.")
    return columns, arrays


def longest_contiguous_true(mask: np.ndarray) -> np.ndarray:
    """Return indices of the longest contiguous true run, with no gap joining."""

    indices = np.flatnonzero(mask)
    if indices.size == 0:
        raise ValueError("No 1C discharge branch was found.")
    split_points = np.flatnonzero(np.diff(indices) != 1) + 1
    runs = np.split(indices, split_points)
    return max(runs, key=len)


def scalar_metrics(measured: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    residual = predicted - measured
    denominator = float(np.sum((measured - np.mean(measured)) ** 2))
    return {
        "rmse_v": float(np.sqrt(np.mean(residual**2))),
        "mae_v": float(np.mean(np.abs(residual))),
        "r2": float(1.0 - np.sum(residual**2) / denominator) if denominator > 0.0 else float("nan"),
        "max_abs_error_v": float(np.max(np.abs(residual))),
        "mean_signed_error_v": float(np.mean(residual)),
    }


def add_metric_rows(
    rows: list[dict[str, object]],
    scope: str,
    metrics: dict[str, object],
    n_samples: int,
) -> None:
    units = {
        "rmse_v": "V",
        "mae_v": "V",
        "r2": "1",
        "max_abs_error_v": "V",
        "mean_signed_error_v": "V",
    }
    for name, value in metrics.items():
        rows.append(
            {
                "scope": scope,
                "metric": name,
                "value": value,
                "unit": units.get(name, "1"),
                "n_samples": n_samples,
            }
        )


def save_metrics(
    path: Path,
    overall: dict[str, float],
    regional: dict[str, dict[str, float]],
    region_counts: dict[str, int],
    cutoff: dict[str, float],
    cutoff_count: int,
    audit: dict[str, object],
) -> None:
    rows: list[dict[str, object]] = []
    add_metric_rows(rows, "overall", overall, int(audit["active_samples"]))
    for region in ("high", "middle", "low"):
        add_metric_rows(rows, f"soc_{region}", regional[region], region_counts[region])
    add_metric_rows(rows, "cutoff_last_1pct_capacity", cutoff, cutoff_count)

    audit_units = {
        "active_samples": "count",
        "active_duration_s": "s",
        "sample_interval_min_s": "s",
        "sample_interval_median_s": "s",
        "sample_interval_mean_s": "s",
        "sample_interval_max_s": "s",
        "raw_current_min_a": "A",
        "raw_current_max_a": "A",
        "measured_capacity_integrated_ah": "Ah",
        "exported_counter_change_ah": "Ah",
        "source_soc_final": "1",
        "predicted_final_soc": "1",
        "measured_voltage_start_v": "V",
        "measured_voltage_end_v": "V",
        "predicted_voltage_start_v": "V",
        "predicted_voltage_end_v": "V",
        "final_cutoff_point_error_v": "V",
        "battery_temperature_min_c": "degC",
        "battery_temperature_max_c": "degC",
        "chamber_temperature_min_c": "degC",
        "chamber_temperature_max_c": "degC",
        "x_n_avg_min": "1",
        "x_n_avg_max": "1",
        "x_p_avg_min": "1",
        "x_p_avg_max": "1",
        "x_n_surface_min": "1",
        "x_n_surface_max": "1",
        "x_p_surface_min": "1",
        "x_p_surface_max": "1",
        "final_x_n_avg": "1",
        "final_x_p_avg": "1",
        "final_x_n_surface": "1",
        "final_x_p_surface": "1",
        "ocp_domain_violation_count": "count",
        "clipping_used": "boolean",
    }
    for name, value in audit.items():
        if name in audit_units:
            rows.append(
                {
                    "scope": "audit",
                    "metric": name,
                    "value": value,
                    "unit": audit_units[name],
                    "n_samples": int(audit["active_samples"]),
                }
            )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("scope", "metric", "value", "unit", "n_samples"))
        writer.writeheader()
        writer.writerows(rows)


def save_voltage_plot(
    path: Path,
    time_s: np.ndarray,
    soc: np.ndarray,
    measured_v: np.ndarray,
    predicted_v: np.ndarray,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 1, figsize=(9.0, 7.2), constrained_layout=True)
    axes[0].plot(time_s / 60.0, measured_v, color="black", linewidth=1.7, label="Measured 25 C 1C")
    axes[0].plot(time_s / 60.0, predicted_v, color="#0072B2", linewidth=1.6, label="Frozen grouped SPM")
    axes[0].set(xlabel="Elapsed discharge time [min]", ylabel="Terminal voltage [V]")
    axes[0].grid(alpha=0.25)
    axes[0].legend()

    axes[1].plot(100.0 * soc, measured_v, color="black", linewidth=1.7, label="Measured")
    axes[1].plot(100.0 * soc, predicted_v, color="#0072B2", linewidth=1.6, label="Predicted")
    axes[1].axvline(100.0 / 3.0, color="0.45", linestyle="--", linewidth=1.0)
    axes[1].axvline(200.0 / 3.0, color="0.45", linestyle="--", linewidth=1.0)
    axes[1].invert_xaxis()
    axes[1].set(xlabel="Cell SOC from relative exported Ah [%]", ylabel="Terminal voltage [V]")
    axes[1].grid(alpha=0.25)
    axes[1].legend()
    fig.suptitle("Panasonic NCR18650PF 1C forward-only validation")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def save_residual_plot(
    path: Path,
    soc: np.ndarray,
    residual_v: np.ndarray,
    regional: dict[str, dict[str, float]],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9.0, 4.8), constrained_layout=True)
    ax.plot(100.0 * soc, 1000.0 * residual_v, color="#D55E00", linewidth=1.6)
    ax.axhline(0.0, color="black", linewidth=0.9)
    ax.axvline(100.0 / 3.0, color="0.45", linestyle="--", linewidth=1.0)
    ax.axvline(200.0 / 3.0, color="0.45", linestyle="--", linewidth=1.0)
    ax.axvspan(100.0 * soc[-1], 100.0 / 3.0, color="#CC79A7", alpha=0.08, label="Low SOC")
    labels = ", ".join(
        f"{name} {1000.0 * regional[name]['rmse_v']:.1f} mV"
        for name in ("high", "middle", "low")
    )
    ax.text(0.02, 0.95, f"Regional RMSE: {labels}", transform=ax.transAxes, va="top")
    ax.invert_xaxis()
    ax.set(
        xlabel="Cell SOC from relative exported Ah [%]",
        ylabel="Residual, predicted - measured [mV]",
        title="Panasonic 1C forward-validation residual versus SOC",
    )
    ax.grid(alpha=0.25)
    ax.legend(loc="lower left")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def fmt_metric_table(metrics: dict[str, float], n_samples: int) -> list[str]:
    return [
        f"| Samples | {n_samples} |",
        f"| RMSE [mV] | {1000.0 * metrics['rmse_v']:.6f} |",
        f"| MAE [mV] | {1000.0 * metrics['mae_v']:.6f} |",
        f"| R2 | {metrics['r2']:.9f} |",
        f"| Maximum absolute error [mV] | {1000.0 * metrics['max_abs_error_v']:.6f} |",
        f"| Mean signed error, predicted - measured [mV] | {1000.0 * metrics['mean_signed_error_v']:.6f} |",
    ]


def write_report(
    path: Path,
    columns: tuple[str, ...],
    audit: dict[str, object],
    overall: dict[str, float],
    regional: dict[str, dict[str, float]],
    region_counts: dict[str, int],
    cutoff: dict[str, float],
    cutoff_count: int,
) -> None:
    source_string = str(RAW_SOURCE_PATH)
    lines = [
        "# Panasonic 25 C 1C forward-only validation",
        "",
        "## 1C TRANSFER LIMITED",
        "",
        "The manifest-selected initial-series 1C discharge was evaluated once with the frozen conditional HPPC parameterization. No parameter was optimized, and no OCP, stoichiometric window, capacity, support interval, or model equation was changed.",
        "",
        "## Data audit and preparation",
        "",
        f"- Exact raw source: `{source_string}`",
        f"- Manifest entry verified: `{MANIFEST_PATH}`",
        f"- Raw columns: `{', '.join(columns)}`",
        f"- Active branch: source indices {audit['start_index']}-{audit['end_index']} inclusive ({audit['active_samples']} samples, {audit['active_duration_s']:.6f} s).",
        f"- Sample interval on the active branch: median {audit['sample_interval_median_s']:.9f} s, mean {audit['sample_interval_mean_s']:.9f} s, range {audit['sample_interval_min_s']:.9f}-{audit['sample_interval_max_s']:.9f} s. Actual nonuniform times were preserved.",
        f"- Raw current is negative on discharge ({audit['raw_current_min_a']:.5f} to {audit['raw_current_max_a']:.5f} A). Its sign was changed only at the model boundary so model current is discharge-positive.",
        f"- Measured discharge capacity from zero-order-hold current/time integration: {audit['measured_capacity_integrated_ah']:.12f} Ah. The relative exported `Ah` counter changes by {audit['exported_counter_change_ah']:.12f} Ah; its arbitrary absolute offset was not interpreted as prior discharge.",
        f"- Measured loaded voltage: {audit['measured_voltage_start_v']:.5f} V at the first active sample and {audit['measured_voltage_end_v']:.5f} V at cutoff.",
        f"- Battery case temperature: {audit['battery_temperature_min_c']:.5f}-{audit['battery_temperature_max_c']:.5f} degC. Chamber temperature: {audit['chamber_temperature_min_c']:.0f}-{audit['chamber_temperature_max_c']:.0f} degC.",
        "- Relative cumulative exported Ah defined the SOC axis and initialized the discharge at SOC=1. The frozen grouped-state equations were then propagated with the measured current and actual time increments; measured battery temperature was prescribed only in the existing kinetic-voltage factor.",
        "",
        "## Frozen forward evaluation",
        "",
        f"- `R0 = {NOMINAL_PARAMETERS['R0']:.13g} ohm`",
        f"- `alpha_p = {NOMINAL_PARAMETERS['alpha_p']:.12g} s`; `alpha_n = {NOMINAL_PARAMETERS['alpha_n']:.12g} s`",
        f"- `d_n = {NOMINAL_PARAMETERS['d_n']:.12g} 1/s`; `d_p = {NOMINAL_PARAMETERS['d_p']:.12g} 1/s`",
        f"- `Q_cell = {Q_CELL_AH:.12f} Ah`; `Q_n = {Q_N_C:.12f} C`; `Q_p = {Q_P_C:.12f} C`",
        "- Frozen nominal effective electrode windows and frozen Guo-White MCMB/13-knot positive OCP were used unchanged.",
        "",
        "## Electrical metrics",
        "",
        "Residual is defined as predicted minus measured voltage.",
        "",
        "| Metric | Value |",
        "|---|---:|",
        *fmt_metric_table(overall, int(audit["active_samples"])),
        "",
        "### By cell-SOC region",
        "",
        "| SOC region | Samples | RMSE [mV] | MAE [mV] | R2 | Max abs [mV] | Mean signed [mV] |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for region in ("high", "middle", "low"):
        metric = regional[region]
        interval = {"high": "[2/3, 1]", "middle": "[1/3, 2/3)", "low": "[0, 1/3)"}[region]
        lines.append(
            f"| {region} {interval} | {region_counts[region]} | {1000.0 * metric['rmse_v']:.6f} | {1000.0 * metric['mae_v']:.6f} | {metric['r2']:.9f} | {1000.0 * metric['max_abs_error_v']:.6f} | {1000.0 * metric['mean_signed_error_v']:.6f} |"
        )
    lines.extend(
        [
            "",
            "### Final cutoff region",
            "",
            f"The cutoff-region diagnostic is the last {100.0 * CUTOFF_REGION_FRACTION:.0f}% of the recorded discharge capacity; it is reported separately and was not given special weight or forced to match.",
            "",
            "| Metric | Value |",
            "|---|---:|",
            *fmt_metric_table(cutoff, cutoff_count),
            f"| Final measured voltage [V] | {audit['measured_voltage_end_v']:.9f} |",
            f"| Final predicted voltage [V] | {audit['predicted_voltage_end_v']:.9f} |",
            f"| Final-point error, predicted - measured [mV] | {1000.0 * audit['final_cutoff_point_error_v']:.6f} |",
            "",
            "## State and domain audit",
            "",
            "| Quantity | Minimum | Maximum | Frozen support |",
            "|---|---:|---:|---:|",
            f"| `x_n_avg` | {audit['x_n_avg_min']:.12f} | {audit['x_n_avg_max']:.12f} | diagnostic state |",
            f"| `x_p_avg` | {audit['x_p_avg_min']:.12f} | {audit['x_p_avg_max']:.12f} | diagnostic state |",
            f"| `x_n_surf` | {audit['x_n_surface_min']:.12f} | {audit['x_n_surface_max']:.12f} | [{MCMB_SUPPORTED_STOICHIOMETRY_RANGE[0]:.2f}, {MCMB_SUPPORTED_STOICHIOMETRY_RANGE[1]:.2f}] |",
            f"| `x_p_surf` | {audit['x_p_surface_min']:.12f} | {audit['x_p_surface_max']:.12f} | [{NCA_WORKING_STOICHIOMETRY_RANGE[0]:.2f}, {NCA_WORKING_STOICHIOMETRY_RANGE[1]:.2f}] |",
            "",
            f"- OCP-domain violations: **{audit['ocp_domain_violation_count']}**.",
            "- Clipping used: **NO**.",
            f"- Final source-Ah cell SOC: {audit['source_soc_final']:.12f}.",
            f"- Predicted final cell SOC from measured current/time integration: {audit['predicted_final_soc']:.12f}.",
            f"- Predicted final states: `x_n_avg={audit['final_x_n_avg']:.12f}`, `x_p_avg={audit['final_x_p_avg']:.12f}`, `x_n_surf={audit['final_x_n_surface']:.12f}`, `x_p_surf={audit['final_x_p_surface']:.12f}`.",
            "",
            "## Interpretation",
            "",
            f"The high-SOC 1C RMSE is {1000.0 * regional['high']['rmse_v']:.3f} mV, close in scale to the HPPC calibration-region phase-balanced RMSE of about 36.835 mV. Error then rises to {1000.0 * regional['middle']['rmse_v']:.3f} mV at middle SOC and {1000.0 * regional['low']['rmse_v']:.3f} mV at low SOC. The residual is predominantly positive, meaning the model holds voltage above the measurement, and expands steeply near cutoff (final-point error {1000.0 * audit['final_cutoff_point_error_v']:.3f} mV).",
            "",
            "This independent continuous-discharge result therefore shows the same low-SOC transfer deterioration previously observed in withheld low-SOC HPPC pulses. The discrepancy is not an OCP-domain or clipping artifact: all surface stoichiometries remain supported. It is a systematic loaded-voltage/low-SOC transfer limitation of the frozen conditional representation, with the largest mismatch in the cutoff region.",
            "",
            "**Final classification: 1C TRANSFER LIMITED.**",
            "",
            "No repair, refit, parameter release, OCP/support/window change, model-equation change, UDDS run, or thermal fit was performed.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    if not MANIFEST_PATH.is_file() or RAW_SOURCE_PATH.name not in MANIFEST_PATH.read_text(encoding="utf-8"):
        raise RuntimeError("The selected raw 1C path is not recorded in the experimental-data manifest.")
    if not np.isclose(Q_CELL_AH, 2.994979138519, rtol=0.0, atol=5.0e-13):
        raise RuntimeError("Frozen Panasonic cell capacity changed.")
    if not (
        np.isclose(IMPLEMENTED_Q_N_C, Q_N_C, rtol=0.0, atol=5.0e-9)
        and np.isclose(IMPLEMENTED_Q_P_C, Q_P_C, rtol=0.0, atol=5.0e-9)
    ):
        raise RuntimeError("Implemented grouped-model electrode capacities differ from the frozen values.")

    columns, raw = load_raw_discharge()
    branch_indices = longest_contiguous_true(raw["Current"] <= DISCHARGE_CURRENT_THRESHOLD_A)
    start_index = int(branch_indices[0])
    end_index = int(branch_indices[-1])
    if not np.array_equal(branch_indices, np.arange(start_index, end_index + 1)):
        raise RuntimeError("Selected discharge branch is not contiguous.")

    time_s = raw["Time"][branch_indices]
    measured_v = raw["Voltage"][branch_indices]
    raw_current_a = raw["Current"][branch_indices]
    model_current_a = -raw_current_a  # sole current-sign conversion boundary
    exported_ah = raw["Ah"][branch_indices]
    battery_temperature_c = raw["Battery_Temp_degC"][branch_indices]
    chamber_temperature_c = raw["Chamber_Temp_degC"][branch_indices]
    dt_s = np.diff(time_s)
    if time_s.size < 2 or np.any(dt_s <= 0.0):
        raise ValueError("Active branch time must be strictly increasing.")
    if np.any(model_current_a <= 0.0):
        raise ValueError("The selected branch is not entirely discharge-positive after boundary conversion.")

    relative_exported_discharge_ah = exported_ah[0] - exported_ah
    source_soc = 1.0 - relative_exported_discharge_ah / Q_CELL_AH
    if np.any((source_soc < 0.0) | (source_soc > 1.0)):
        raise ValueError("Source-Ah SOC left [0, 1]; no clipping is allowed.")
    initial_z = relative_exported_discharge_ah[0] / Q_CELL_AH
    initial_x_n = X_N_START + initial_z * (X_N_END - X_N_START)
    initial_x_p = X_P_START + initial_z * (X_P_END - X_P_START)

    result = simulate_local_group(
        time_s=time_s,
        discharge_positive_current_a=model_current_a,
        temperature_k=battery_temperature_c + 273.15,
        initial_x_n=float(initial_x_n),
        initial_x_p=float(initial_x_p),
        dynamic_parameters=NOMINAL_PARAMETERS,
    )
    predicted_v = result.voltage_v
    residual_v = predicted_v - measured_v

    negative_violation = (result.x_n_surface < MCMB_SUPPORTED_STOICHIOMETRY_RANGE[0]) | (
        result.x_n_surface > MCMB_SUPPORTED_STOICHIOMETRY_RANGE[1]
    )
    positive_violation = (result.x_p_surface < NCA_WORKING_STOICHIOMETRY_RANGE[0]) | (
        result.x_p_surface > NCA_WORKING_STOICHIOMETRY_RANGE[1]
    )
    violation_count = int(np.count_nonzero(negative_violation | positive_violation))

    overall = scalar_metrics(measured_v, predicted_v)
    regional: dict[str, dict[str, float]] = {}
    region_counts: dict[str, int] = {}
    for region, mask_function in SOC_REGION_MASKS.items():
        mask = mask_function(source_soc)
        if not np.any(mask):
            raise RuntimeError(f"No samples in required {region}-SOC region.")
        regional[region] = scalar_metrics(measured_v[mask], predicted_v[mask])
        region_counts[region] = int(np.count_nonzero(mask))

    profile_progress = relative_exported_discharge_ah / relative_exported_discharge_ah[-1]
    cutoff_mask = profile_progress >= 1.0 - CUTOFF_REGION_FRACTION
    cutoff = scalar_metrics(measured_v[cutoff_mask], predicted_v[cutoff_mask])
    cutoff_count = int(np.count_nonzero(cutoff_mask))

    integrated_capacity_ah = float(np.sum(model_current_a[:-1] * dt_s) / 3600.0)
    predicted_final_soc = 1.0 - integrated_capacity_ah / Q_CELL_AH
    audit: dict[str, object] = {
        "start_index": start_index,
        "end_index": end_index,
        "active_samples": int(time_s.size),
        "active_duration_s": float(time_s[-1] - time_s[0]),
        "sample_interval_min_s": float(np.min(dt_s)),
        "sample_interval_median_s": float(np.median(dt_s)),
        "sample_interval_mean_s": float(np.mean(dt_s)),
        "sample_interval_max_s": float(np.max(dt_s)),
        "raw_current_min_a": float(np.min(raw_current_a)),
        "raw_current_max_a": float(np.max(raw_current_a)),
        "measured_capacity_integrated_ah": integrated_capacity_ah,
        "exported_counter_change_ah": float(relative_exported_discharge_ah[-1]),
        "source_soc_final": float(source_soc[-1]),
        "predicted_final_soc": predicted_final_soc,
        "measured_voltage_start_v": float(measured_v[0]),
        "measured_voltage_end_v": float(measured_v[-1]),
        "predicted_voltage_start_v": float(predicted_v[0]),
        "predicted_voltage_end_v": float(predicted_v[-1]),
        "final_cutoff_point_error_v": float(residual_v[-1]),
        "battery_temperature_min_c": float(np.min(battery_temperature_c)),
        "battery_temperature_max_c": float(np.max(battery_temperature_c)),
        "chamber_temperature_min_c": float(np.min(chamber_temperature_c)),
        "chamber_temperature_max_c": float(np.max(chamber_temperature_c)),
        "x_n_avg_min": float(np.min(result.x_n_avg)),
        "x_n_avg_max": float(np.max(result.x_n_avg)),
        "x_p_avg_min": float(np.min(result.x_p_avg)),
        "x_p_avg_max": float(np.max(result.x_p_avg)),
        "x_n_surface_min": float(np.min(result.x_n_surface)),
        "x_n_surface_max": float(np.max(result.x_n_surface)),
        "x_p_surface_min": float(np.min(result.x_p_surface)),
        "x_p_surface_max": float(np.max(result.x_p_surface)),
        "final_x_n_avg": float(result.x_n_avg[-1]),
        "final_x_p_avg": float(result.x_p_avg[-1]),
        "final_x_n_surface": float(result.x_n_surface[-1]),
        "final_x_p_surface": float(result.x_p_surface[-1]),
        "ocp_domain_violation_count": violation_count,
        "clipping_used": False,
    }

    if violation_count:
        raise RuntimeError("Frozen OCP-domain violation occurred; no extrapolation is permitted.")
    if not np.all(np.isfinite(predicted_v)):
        raise RuntimeError("Non-finite forward voltage prediction.")

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    save_metrics(
        OUTPUT_DIRECTORY / "panasonic_1c_forward_validation_metrics.csv",
        overall,
        regional,
        region_counts,
        cutoff,
        cutoff_count,
        audit,
    )
    save_voltage_plot(
        OUTPUT_DIRECTORY / "panasonic_1c_forward_validation_voltage.png",
        time_s - time_s[0],
        source_soc,
        measured_v,
        predicted_v,
    )
    save_residual_plot(
        OUTPUT_DIRECTORY / "panasonic_1c_forward_validation_residual_soc.png",
        source_soc,
        residual_v,
        regional,
    )
    write_report(
        OUTPUT_DIRECTORY / "panasonic_1c_forward_validation_report.md",
        columns,
        audit,
        overall,
        regional,
        region_counts,
        cutoff,
        cutoff_count,
    )

    print("Panasonic 25 C 1C forward-only validation complete.")
    print(f"Active samples: {time_s.size}; integrated capacity: {integrated_capacity_ah:.12f} Ah")
    print(
        f"Overall RMSE={1000.0 * overall['rmse_v']:.6f} mV; "
        f"high/middle/low={1000.0 * regional['high']['rmse_v']:.6f}/"
        f"{1000.0 * regional['middle']['rmse_v']:.6f}/"
        f"{1000.0 * regional['low']['rmse_v']:.6f} mV"
    )
    print(f"Final-point error={1000.0 * residual_v[-1]:.6f} mV")
    print("OCP-domain violations: 0; clipping: none")
    print("Classification: 1C TRANSFER LIMITED")


if __name__ == "__main__":
    main()
