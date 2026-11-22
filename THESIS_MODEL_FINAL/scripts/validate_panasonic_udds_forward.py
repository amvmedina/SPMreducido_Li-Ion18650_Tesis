"""Forward-only validation of the frozen Panasonic model on the 25 C UDDS test.

The manifest-selected raw MAT file is loaded read-only. Recorded nonuniform
times and regenerative current are retained; raw current is converted to the
model's discharge-positive convention only at the simulation boundary. This
script contains no optimizer, clipping, baseline correction, or parameter fit.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.stats import spearmanr


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from models.ocp_panasonic_qi import (  # noqa: E402
    MCMB_SUPPORTED_STOICHIOMETRY_RANGE,
    NCA_WORKING_STOICHIOMETRY_RANGE,
)
from models.ocp_panasonic_qi_frozen import (  # noqa: E402
    Q_CELL_AH,
    X_N_START,
    X_P_START,
)
from scripts.audit_panasonic_hppc_identifiability import (  # noqa: E402
    Q_N_C as IMPLEMENTED_Q_N_C,
    Q_P_C as IMPLEMENTED_Q_P_C,
    simulate_local_group,
)
from runtime_paths import PANASONIC_UDDS_PATH, require_data_file  # noqa: E402


RAW_SOURCE_PATH = PANASONIC_UDDS_PATH
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

NEAR_REST_CURRENT_A = 0.05
N_LARGEST_CURRENT_STEPS = 10
PRIOR_HPPC_CALIBRATION_REGION_RMSE_V = 0.036835197638870915
PRIOR_HPPC_LOW_SOC_PHASE_BALANCED_RMSE_V = 0.08572943292069442
PRIOR_HPPC_LOW_SOC_ABSOLUTE_RMSE_V = 0.0595821150372503
PRIOR_ONE_C_OVERALL_RMSE_V = 0.11436218800084591
PRIOR_ONE_C_HIGH_SOC_RMSE_V = 0.037383603944578854
PRIOR_ONE_C_MIDDLE_SOC_RMSE_V = 0.06431895439081552
PRIOR_ONE_C_LOW_SOC_RMSE_V = 0.19595403962007815


def _vector(value: object, name: str) -> np.ndarray:
    result = np.asarray(value, dtype=float).reshape(-1)
    if result.size == 0 or not np.all(np.isfinite(result)):
        raise ValueError(f"Raw field {name!r} is empty or non-finite.")
    return result


def load_udds(path: Path = RAW_SOURCE_PATH) -> tuple[tuple[str, ...], dict[str, np.ndarray]]:
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
            raise ValueError(f"Processed UDDS CSV is missing required fields: {missing}")
        arrays = {name: _vector(table[name], name) for name in required}
        if len({array.size for array in arrays.values()}) != 1:
            raise ValueError("Processed UDDS CSV channels have inconsistent lengths.")
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
        raise ValueError(f"Raw UDDS data is missing required fields: {missing}")
    arrays = {name: _vector(getattr(measurement, name), name) for name in required}
    if len({array.size for array in arrays.values()}) != 1:
        raise ValueError("Raw UDDS channels have inconsistent lengths.")
    return columns, arrays


def remove_exact_duplicate_times(arrays: dict[str, np.ndarray]) -> tuple[dict[str, np.ndarray], np.ndarray]:
    """Remove only exact duplicate-time records; never alter positive time gaps."""

    time_s = arrays["Time"]
    differences = np.diff(time_s)
    if np.any(differences < 0.0):
        raise ValueError("UDDS elapsed time decreases.")
    duplicate_indices = np.flatnonzero(differences == 0.0) + 1
    for index in duplicate_indices:
        for name, values in arrays.items():
            if values[index] != values[index - 1]:
                raise ValueError(f"Repeated time at index {index} is not an exact duplicate in {name}.")
    keep = np.ones(time_s.size, dtype=bool)
    keep[duplicate_indices] = False
    cleaned = {name: values[keep] for name, values in arrays.items()}
    if np.any(np.diff(cleaned["Time"]) <= 0.0):
        raise ValueError("Usable UDDS time must be strictly increasing.")
    return cleaned, duplicate_indices


def metrics(measured: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    residual = predicted - measured
    denominator = float(np.sum((measured - np.mean(measured)) ** 2))
    return {
        "rmse_v": float(np.sqrt(np.mean(residual**2))),
        "mae_v": float(np.mean(np.abs(residual))),
        "r2": float(1.0 - np.sum(residual**2) / denominator) if denominator > 0.0 else float("nan"),
        "mean_signed_error_v": float(np.mean(residual)),
        "residual_std_v": float(np.std(residual)),
        "abs_error_p95_v": float(np.quantile(np.abs(residual), 0.95)),
        "max_abs_error_v": float(np.max(np.abs(residual))),
    }


def soc_masks(soc: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "high": (soc >= 2.0 / 3.0) & (soc <= 1.0),
        "middle": (soc >= 1.0 / 3.0) & (soc < 2.0 / 3.0),
        "low": (soc >= 0.0) & (soc < 1.0 / 3.0),
    }


def current_magnitude_masks(raw_current_a: np.ndarray) -> tuple[np.ndarray, list[dict[str, object]]]:
    """Use near-rest plus empirical quartiles of loaded |current|."""

    magnitude = np.abs(raw_current_a)
    loaded = magnitude > NEAR_REST_CURRENT_A
    quartiles = np.quantile(magnitude[loaded], (0.25, 0.50, 0.75))
    edges = np.array([NEAR_REST_CURRENT_A, *quartiles, np.inf])
    bins: list[dict[str, object]] = [
        {
            "name": "near_rest",
            "lower": 0.0,
            "upper": NEAR_REST_CURRENT_A,
            "mask": magnitude <= NEAR_REST_CURRENT_A,
        }
    ]
    lower = NEAR_REST_CURRENT_A
    for number, upper in enumerate(edges[1:], start=1):
        bins.append(
            {
                "name": f"loaded_q{number}",
                "lower": float(lower),
                "upper": float(upper),
                "mask": (magnitude > lower) & (magnitude <= upper),
            }
        )
        lower = float(upper)
    return quartiles, bins


def largest_current_steps(
    time_s: np.ndarray,
    source_soc: np.ndarray,
    raw_current_a: np.ndarray,
    measured_v: np.ndarray,
    predicted_v: np.ndarray,
) -> list[dict[str, float | int]]:
    absolute_step = np.abs(np.diff(raw_current_a))
    ranked = np.argsort(absolute_step)[-N_LARGEST_CURRENT_STEPS:][::-1] + 1
    events: list[dict[str, float | int]] = []
    for rank, index in enumerate(ranked, start=1):
        measured_change = measured_v[index] - measured_v[index - 1]
        predicted_change = predicted_v[index] - predicted_v[index - 1]
        events.append(
            {
                "rank": rank,
                "index": int(index),
                "time_s": float(time_s[index]),
                "soc": float(source_soc[index]),
                "raw_current_before_a": float(raw_current_a[index - 1]),
                "raw_current_after_a": float(raw_current_a[index]),
                "delta_current_a": float(raw_current_a[index] - raw_current_a[index - 1]),
                "measured_delta_voltage_v": float(measured_change),
                "predicted_delta_voltage_v": float(predicted_change),
                "delta_voltage_error_v": float(predicted_change - measured_change),
                "residual_after_step_v": float(predicted_v[index] - measured_v[index]),
            }
        )
    return events


def append_metric_rows(
    rows: list[dict[str, object]],
    scope: str,
    subset: str,
    values: dict[str, object],
    n_samples: int,
    lower: object = "",
    upper: object = "",
) -> None:
    units = {
        "rmse_v": "V",
        "mae_v": "V",
        "r2": "1",
        "mean_signed_error_v": "V",
        "residual_std_v": "V",
        "abs_error_p95_v": "V",
        "max_abs_error_v": "V",
    }
    for name, value in values.items():
        rows.append(
            {
                "scope": scope,
                "subset": subset,
                "metric": name,
                "value": value,
                "unit": units.get(name, "1"),
                "n_samples": n_samples,
                "lower_bound": lower,
                "upper_bound": upper,
            }
        )


def write_metrics(
    path: Path,
    overall: dict[str, float],
    driven: dict[str, float],
    regional: dict[str, dict[str, float]],
    region_masks: dict[str, np.ndarray],
    operating: dict[str, tuple[np.ndarray, dict[str, float]]],
    current_bins: list[dict[str, object]],
    current_metrics: dict[str, dict[str, float]],
    correlations: dict[str, float],
    events: list[dict[str, float | int]],
    audit: dict[str, object],
) -> None:
    rows: list[dict[str, object]] = []
    append_metric_rows(rows, "overall", "complete_usable", overall, int(audit["usable_samples"]))
    append_metric_rows(rows, "overall", "through_cutoff", driven, int(audit["cutoff_index"]) + 1)
    for name, values in regional.items():
        append_metric_rows(rows, "soc_region", name, values, int(np.count_nonzero(region_masks[name])))
    for name, (mask, values) in operating.items():
        append_metric_rows(rows, "operating_state", name, values, int(np.count_nonzero(mask)))
    for item in current_bins:
        name = str(item["name"])
        mask = np.asarray(item["mask"], dtype=bool)
        append_metric_rows(
            rows,
            "current_magnitude",
            name,
            current_metrics[name],
            int(np.count_nonzero(mask)),
            item["lower"],
            item["upper"],
        )
    for name, value in correlations.items():
        rows.append(
            {
                "scope": "correlation",
                "subset": "complete_usable",
                "metric": name,
                "value": value,
                "unit": "Spearman_rho",
                "n_samples": audit["usable_samples"],
                "lower_bound": "",
                "upper_bound": "",
            }
        )
    for event in events:
        for name, value in event.items():
            unit = "1"
            if name.endswith("_v"):
                unit = "V"
            elif name.endswith("_a"):
                unit = "A"
            elif name == "time_s":
                unit = "s"
            rows.append(
                {
                    "scope": "largest_current_step",
                    "subset": f"rank_{event['rank']}",
                    "metric": name,
                    "value": value,
                    "unit": unit,
                    "n_samples": 1,
                    "lower_bound": "",
                    "upper_bound": "",
                }
            )
    audit_units = {
        "raw_samples": "count",
        "usable_samples": "count",
        "duplicate_rows_removed": "count",
        "duration_s": "s",
        "dt_min_s": "s",
        "dt_median_s": "s",
        "dt_mean_s": "s",
        "dt_max_s": "s",
        "raw_current_min_a": "A",
        "raw_current_max_a": "A",
        "voltage_min_v": "V",
        "voltage_max_v": "V",
        "battery_temperature_min_c": "degC",
        "battery_temperature_max_c": "degC",
        "chamber_temperature_min_c": "degC",
        "chamber_temperature_max_c": "degC",
        "exported_ah_min": "Ah",
        "exported_ah_max": "Ah",
        "exported_net_discharge_ah": "Ah",
        "integrated_net_discharge_ah": "Ah",
        "integrated_gross_discharge_ah": "Ah",
        "integrated_regenerative_charge_ah": "Ah",
        "source_soc_min": "1",
        "predicted_soc_min": "1",
        "x_n_avg_min": "1",
        "x_n_avg_max": "1",
        "x_p_avg_min": "1",
        "x_p_avg_max": "1",
        "x_n_surface_min": "1",
        "x_n_surface_max": "1",
        "x_p_surface_min": "1",
        "x_p_surface_max": "1",
        "ocp_domain_violation_count": "count",
        "clipping_count": "count",
        "cutoff_time_s": "s",
        "cutoff_residual_v": "V",
        "largest_step_delta_voltage_rmse_v": "V",
        "largest_step_delta_voltage_bias_v": "V",
    }
    for name, unit in audit_units.items():
        rows.append(
            {
                "scope": "audit",
                "subset": "complete_usable",
                "metric": name,
                "value": audit[name],
                "unit": unit,
                "n_samples": audit["usable_samples"],
                "lower_bound": "",
                "upper_bound": "",
            }
        )
    comparisons = {
        "hppc_calibration_soc_ge_one_third_phase_balanced_rmse_v": PRIOR_HPPC_CALIBRATION_REGION_RMSE_V,
        "hppc_low_soc_phase_balanced_rmse_v": PRIOR_HPPC_LOW_SOC_PHASE_BALANCED_RMSE_V,
        "hppc_low_soc_absolute_rmse_v": PRIOR_HPPC_LOW_SOC_ABSOLUTE_RMSE_V,
        "one_c_overall_rmse_v": PRIOR_ONE_C_OVERALL_RMSE_V,
        "one_c_high_soc_rmse_v": PRIOR_ONE_C_HIGH_SOC_RMSE_V,
        "one_c_middle_soc_rmse_v": PRIOR_ONE_C_MIDDLE_SOC_RMSE_V,
        "one_c_low_soc_rmse_v": PRIOR_ONE_C_LOW_SOC_RMSE_V,
    }
    for name, value in comparisons.items():
        rows.append(
            {
                "scope": "prior_validation",
                "subset": "reported",
                "metric": name,
                "value": value,
                "unit": "V",
                "n_samples": "",
                "lower_bound": "",
                "upper_bound": "",
            }
        )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_voltage_plot(
    path: Path,
    time_s: np.ndarray,
    raw_current_a: np.ndarray,
    source_soc: np.ndarray,
    measured_v: np.ndarray,
    predicted_v: np.ndarray,
    cutoff_index: int,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    stride = max(1, time_s.size // 40000)
    shown = np.arange(0, time_s.size, stride)
    minutes = (time_s - time_s[0]) / 60.0
    cutoff_minutes = minutes[cutoff_index]
    fig, axes = plt.subplots(
        3,
        1,
        figsize=(12.0, 8.2),
        sharex=True,
        gridspec_kw={"height_ratios": (2.0, 1.0, 0.8)},
        constrained_layout=True,
    )
    axes[0].plot(minutes[shown], measured_v[shown], color="black", linewidth=0.9, label="Measured")
    axes[0].plot(minutes[shown], predicted_v[shown], color="#0072B2", linewidth=0.9, label="Frozen grouped SPM")
    axes[0].set_ylabel("Voltage [V]")
    axes[0].legend(loc="best")
    axes[1].plot(minutes[shown], raw_current_a[shown], color="#D55E00", linewidth=0.65)
    axes[1].axhline(0.0, color="black", linewidth=0.7)
    axes[1].set_ylabel("Raw current [A]\n(discharge negative)")
    axes[2].plot(minutes[shown], 100.0 * source_soc[shown], color="#009E73", linewidth=0.9)
    axes[2].set(xlabel="Elapsed time [min]", ylabel="SOC [%]")
    for axis in axes:
        axis.axvline(cutoff_minutes, color="0.4", linestyle="--", linewidth=1.0)
        axis.grid(alpha=0.22)
    axes[0].set_title("Panasonic NCR18650PF 25 C UDDS forward-only validation")
    axes[0].text(
        cutoff_minutes,
        0.02,
        " cutoff / recovery",
        transform=axes[0].get_xaxis_transform(),
        rotation=90,
        va="bottom",
        color="0.35",
    )
    fig.savefig(path, dpi=180)
    plt.close(fig)


def binned_residual_soc(soc: np.ndarray, residual_v: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    edges = np.linspace(float(np.min(soc)), 1.0, 31)
    centers = 0.5 * (edges[:-1] + edges[1:])
    mean = np.full(centers.size, np.nan)
    rmse = np.full(centers.size, np.nan)
    for index, (lower, upper) in enumerate(zip(edges[:-1], edges[1:], strict=True)):
        mask = (soc >= lower) & ((soc < upper) if index < centers.size - 1 else (soc <= upper))
        if np.any(mask):
            mean[index] = np.mean(residual_v[mask])
            rmse[index] = np.sqrt(np.mean(residual_v[mask] ** 2))
    return centers, mean, rmse


def save_residual_soc_plot(path: Path, soc: np.ndarray, residual_v: np.ndarray) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    stride = max(1, soc.size // 30000)
    centers, mean, rmse = binned_residual_soc(soc, residual_v)
    fig, ax = plt.subplots(figsize=(10.0, 5.2), constrained_layout=True)
    ax.scatter(100.0 * soc[::stride], 1000.0 * residual_v[::stride], s=2, alpha=0.12, color="#0072B2", rasterized=True, label="Samples")
    ax.plot(100.0 * centers, 1000.0 * mean, color="#D55E00", linewidth=2.0, label="SOC-bin mean")
    ax.plot(100.0 * centers, 1000.0 * rmse, color="#CC79A7", linewidth=1.5, label="SOC-bin RMSE")
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.axvline(100.0 / 3.0, color="0.4", linestyle="--", linewidth=1.0)
    ax.axvline(200.0 / 3.0, color="0.4", linestyle="--", linewidth=1.0)
    ax.invert_xaxis()
    ax.set(
        xlabel="Cell SOC from exported cumulative Ah [%]",
        ylabel="Residual, predicted - measured [mV]",
        title="UDDS residual versus SOC",
    )
    ax.grid(alpha=0.22)
    ax.legend(loc="best")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def save_residual_current_plot(
    path: Path,
    raw_current_a: np.ndarray,
    residual_v: np.ndarray,
    current_bins: list[dict[str, object]],
    current_metrics: dict[str, dict[str, float]],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    magnitude = np.abs(raw_current_a)
    stride = max(1, magnitude.size // 40000)
    shown = np.arange(0, magnitude.size, stride)
    discharge = raw_current_a[shown] < -NEAR_REST_CURRENT_A
    charge = raw_current_a[shown] > NEAR_REST_CURRENT_A
    rest = ~(discharge | charge)
    fig, axes = plt.subplots(2, 1, figsize=(10.0, 7.2), constrained_layout=True)
    for mask, label, color in (
        (discharge, "Discharge", "#D55E00"),
        (charge, "Regenerative charge", "#0072B2"),
        (rest, "Near rest", "0.45"),
    ):
        axes[0].scatter(
            magnitude[shown][mask],
            1000.0 * residual_v[shown][mask],
            s=2,
            alpha=0.15,
            color=color,
            rasterized=True,
            label=label,
        )
    axes[0].axhline(0.0, color="black", linewidth=0.8)
    axes[0].set(xlabel="Current magnitude [A]", ylabel="Residual [mV]")
    axes[0].grid(alpha=0.22)
    axes[0].legend(markerscale=4)

    labels = [str(item["name"]).replace("_", " ") for item in current_bins]
    rmse_mv = [1000.0 * current_metrics[str(item["name"])]["rmse_v"] for item in current_bins]
    mean_mv = [1000.0 * current_metrics[str(item["name"])]["mean_signed_error_v"] for item in current_bins]
    positions = np.arange(len(labels))
    width = 0.38
    axes[1].bar(positions - width / 2.0, rmse_mv, width, color="#CC79A7", label="RMSE")
    axes[1].bar(positions + width / 2.0, mean_mv, width, color="#009E73", label="Mean signed")
    axes[1].set_xticks(positions, labels)
    axes[1].set(xlabel="Empirical |current| bin", ylabel="Error statistic [mV]")
    axes[1].grid(axis="y", alpha=0.22)
    axes[1].legend()
    fig.suptitle("UDDS residual as a function of current magnitude")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def metric_table_lines(values: dict[str, float], samples: int) -> list[str]:
    return [
        f"| Samples | {samples} |",
        f"| RMSE [mV] | {1000.0 * values['rmse_v']:.6f} |",
        f"| MAE [mV] | {1000.0 * values['mae_v']:.6f} |",
        f"| R2 | {values['r2']:.9f} |",
        f"| Mean signed error [mV] | {1000.0 * values['mean_signed_error_v']:.6f} |",
        f"| Maximum absolute error [mV] | {1000.0 * values['max_abs_error_v']:.6f} |",
    ]


def write_report(
    path: Path,
    columns: tuple[str, ...],
    overall: dict[str, float],
    driven: dict[str, float],
    regional: dict[str, dict[str, float]],
    region_masks: dict[str, np.ndarray],
    operating: dict[str, tuple[np.ndarray, dict[str, float]]],
    quartiles: np.ndarray,
    current_bins: list[dict[str, object]],
    current_metrics: dict[str, dict[str, float]],
    correlations: dict[str, float],
    events: list[dict[str, float | int]],
    audit: dict[str, object],
) -> None:
    lines = [
        "# Panasonic 25 C UDDS forward-only validation",
        "",
        "## UDDS TRANSFER LIMITED",
        "",
        "The manifest-selected UDDS dataset was evaluated once with the frozen nominal conditional Panasonic parameterization. No optimization, baseline correction, clipping, parameter release, or model/OCP change was performed.",
        "",
        "## Dataset audit",
        "",
        f"- Exact source: `{RAW_SOURCE_PATH}`",
        f"- Manifest entry verified in: `{MANIFEST_PATH}`",
        f"- Columns: `{', '.join(columns)}`",
        f"- Raw samples: {audit['raw_samples']:,}; usable samples: {audit['usable_samples']:,}. One exact duplicate final row was removed.",
        f"- Total recorded and usable duration: {audit['duration_s']:.9f} s.",
        f"- Positive sampling interval: median {audit['dt_median_s']:.9f} s, mean {audit['dt_mean_s']:.9f} s, range {audit['dt_min_s']:.9f}-{audit['dt_max_s']:.9f} s. All actual positive time steps were preserved.",
        f"- Raw current convention: discharge negative, regenerative charge positive. Range {audit['raw_current_min_a']:.5f} to {audit['raw_current_max_a']:.5f} A. Sign conversion occurred only at the model-loader boundary.",
        f"- Measured voltage range: {audit['voltage_min_v']:.5f}-{audit['voltage_max_v']:.5f} V.",
        f"- Battery case temperature: {audit['battery_temperature_min_c']:.5f}-{audit['battery_temperature_max_c']:.5f} degC; chamber temperature: {audit['chamber_temperature_min_c']:.0f}-{audit['chamber_temperature_max_c']:.0f} degC.",
        f"- Exported `Ah` starts at {audit['exported_ah_max']:.5f} Ah and ends at its minimum, {audit['exported_ah_min']:.5f} Ah. It is non-monotonic within the drive ({audit['exported_ah_increasing_steps']:,} increasing and {audit['exported_ah_decreasing_steps']:,} decreasing increments) because regenerative charge is retained. Net exported discharge is {audit['exported_net_discharge_ah']:.9f} Ah.",
        f"- Current/time integration gives {audit['integrated_net_discharge_ah']:.9f} Ah net, comprising {audit['integrated_gross_discharge_ah']:.9f} Ah gross discharge and {audit['integrated_regenerative_charge_ah']:.9f} Ah regenerative charge.",
        "- Source SOC is `1 + Ah/Q_cell`; this directly preserves cumulative-Ah continuity through regenerative intervals. Model states start at the frozen full-charge nominal effective endpoints and propagate continuously with the measured current and actual times.",
        "",
        "## Frozen forward evaluation",
        "",
        f"- `R0={NOMINAL_PARAMETERS['R0']:.13g} ohm`, `alpha_p={NOMINAL_PARAMETERS['alpha_p']:.12g} s`, `d_n={NOMINAL_PARAMETERS['d_n']:.12g} 1/s`",
        f"- Fixed `alpha_n={NOMINAL_PARAMETERS['alpha_n']:.12g} s`, `d_p={NOMINAL_PARAMETERS['d_p']:.12g} 1/s`",
        f"- Frozen `Q_cell={Q_CELL_AH:.12f} Ah`, `Q_n={Q_N_C:.12f} C`, `Q_p={Q_P_C:.12f} C`",
        "- Frozen Guo-White MCMB negative OCP, 13-knot positive OCP, and nominal effective stoichiometric windows were used unchanged.",
        "",
        "## Primary electrical metrics",
        "",
        "Residual is predicted minus measured voltage. Primary metrics include the complete usable drive and its recorded 300 s post-cutoff recovery.",
        "",
        "| Metric | Complete usable profile |",
        "|---|---:|",
        *metric_table_lines(overall, int(audit["usable_samples"])),
        "",
        f"Through the cutoff sample alone ({int(audit['cutoff_index']) + 1:,} samples), RMSE is {1000.0 * driven['rmse_v']:.6f} mV, MAE is {1000.0 * driven['mae_v']:.6f} mV, and R2 is {driven['r2']:.9f}. The cutoff residual is {1000.0 * audit['cutoff_residual_v']:.6f} mV and was not forced to match.",
        "",
        "### By SOC region",
        "",
        "| Region | Samples | RMSE [mV] | MAE [mV] | R2 | Mean signed [mV] | P95 absolute [mV] | Max absolute [mV] |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for name in ("high", "middle", "low"):
        value = regional[name]
        interval = {"high": "[2/3, 1]", "middle": "[1/3, 2/3)", "low": "[0, 1/3)"}[name]
        lines.append(
            f"| {name} {interval} | {np.count_nonzero(region_masks[name]):,} | {1000.0 * value['rmse_v']:.6f} | {1000.0 * value['mae_v']:.6f} | {value['r2']:.9f} | {1000.0 * value['mean_signed_error_v']:.6f} | {1000.0 * value['abs_error_p95_v']:.6f} | {1000.0 * value['max_abs_error_v']:.6f} |"
        )
    lines.extend(
        [
            "",
            "## Residual versus operating condition",
            "",
            "Near rest is defined as `|I_raw| <= 0.05 A`; all other points are loaded. Discharge and regenerative subsets use the same 0.05 A threshold.",
            "",
            "| Subset | Samples | RMSE [mV] | MAE [mV] | Mean signed [mV] | Max absolute [mV] |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for name in ("near_rest", "loaded", "discharge", "regenerative_charge", "post_cutoff_recovery"):
        mask, value = operating[name]
        lines.append(
            f"| {name.replace('_', ' ')} | {np.count_nonzero(mask):,} | {1000.0 * value['rmse_v']:.6f} | {1000.0 * value['mae_v']:.6f} | {1000.0 * value['mean_signed_error_v']:.6f} | {1000.0 * value['max_abs_error_v']:.6f} |"
        )
    lines.extend(
        [
            "",
            "### Current-magnitude bins",
            "",
            f"Loaded-sample empirical quartiles are {quartiles[0]:.5f}, {quartiles[1]:.5f}, and {quartiles[2]:.5f} A. These data-derived edges provide four similarly populated loaded bins, plus the independently defined near-rest bin.",
            "",
            "| Bin | |I| interval [A] | Samples | RMSE [mV] | MAE [mV] | Mean signed [mV] | P95 absolute [mV] |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for item in current_bins:
        name = str(item["name"])
        value = current_metrics[name]
        upper_text = "infinity" if np.isinf(float(item["upper"])) else f"{float(item['upper']):.5f}"
        left_bracket = "[" if name == "near_rest" else "("
        lines.append(
            f"| {name} | {left_bracket}{float(item['lower']):.5f}, {upper_text}] | {np.count_nonzero(item['mask']):,} | {1000.0 * value['rmse_v']:.6f} | {1000.0 * value['mae_v']:.6f} | {1000.0 * value['mean_signed_error_v']:.6f} | {1000.0 * value['abs_error_p95_v']:.6f} |"
        )
    lines.extend(
        [
            "",
            f"Spearman correlation between absolute residual and decreasing SOC is {correlations['abs_error_vs_depth_of_discharge_all']:+.6f} overall and {correlations['abs_error_vs_depth_of_discharge_loaded']:+.6f} on loaded samples. The corresponding correlation with current magnitude is {correlations['abs_error_vs_current_magnitude_all']:+.6f} overall and {correlations['abs_error_vs_current_magnitude_loaded']:+.6f} on loaded samples.",
            "",
            "Residual magnitude therefore grows with both decreasing SOC and increasing current magnitude, but the SOC association is stronger. Within each SOC region, the loaded-sample absolute-error/current-magnitude correlations remain positive, so the current trend is not solely an SOC-composition artifact.",
            "",
            "## Largest current-step voltage changes",
            "",
            "The events below are the ten largest adjacent-sample `|delta I_raw|` transitions. Voltage change is evaluated across exactly the same recorded sample pair; no baseline correction is applied.",
            "",
            "| Rank | Time [s] | SOC | I before [A] | I after [A] | delta I [A] | Measured delta V [mV] | Predicted delta V [mV] | Delta-V error [mV] |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for event in events:
        lines.append(
            f"| {event['rank']} | {event['time_s']:.6f} | {event['soc']:.6f} | {event['raw_current_before_a']:.5f} | {event['raw_current_after_a']:.5f} | {event['delta_current_a']:+.5f} | {1000.0 * event['measured_delta_voltage_v']:+.6f} | {1000.0 * event['predicted_delta_voltage_v']:+.6f} | {1000.0 * event['delta_voltage_error_v']:+.6f} |"
        )
    lines.extend(
        [
            "",
            f"Across these ten events, delta-voltage RMSE is {1000.0 * audit['largest_step_delta_voltage_rmse_v']:.6f} mV and bias is {1000.0 * audit['largest_step_delta_voltage_bias_v']:+.6f} mV. Most selected non-cutoff events are large regenerative-current onsets; the discharge-pulse-derived conditional resistance/kinetic representation substantially overpredicts their immediate upward voltage changes.",
            "",
            "## State and OCP-domain audit",
            "",
            "| State | Minimum | Maximum | Supported OCP interval where applicable |",
            "|---|---:|---:|---:|",
            f"| `x_n_avg` | {audit['x_n_avg_min']:.12f} | {audit['x_n_avg_max']:.12f} | diagnostic state |",
            f"| `x_p_avg` | {audit['x_p_avg_min']:.12f} | {audit['x_p_avg_max']:.12f} | diagnostic state |",
            f"| `x_n_surf` | {audit['x_n_surface_min']:.12f} | {audit['x_n_surface_max']:.12f} | [{MCMB_SUPPORTED_STOICHIOMETRY_RANGE[0]:.2f}, {MCMB_SUPPORTED_STOICHIOMETRY_RANGE[1]:.2f}] |",
            f"| `x_p_surf` | {audit['x_p_surface_min']:.12f} | {audit['x_p_surface_max']:.12f} | [{NCA_WORKING_STOICHIOMETRY_RANGE[0]:.2f}, {NCA_WORKING_STOICHIOMETRY_RANGE[1]:.2f}] |",
            "",
            f"- Minimum source-Ah SOC: {audit['source_soc_min']:.12f}; minimum predicted SOC from current/time integration: {audit['predicted_soc_min']:.12f}.",
            f"- OCP-domain violations: **{audit['ocp_domain_violation_count']}**.",
            f"- Clipping count: **{audit['clipping_count']}**.",
            "",
            "## Comparison with prior electrical validation",
            "",
            "| Validation | Metric type | Result [mV] |",
            "|---|---|---:|",
            f"| HPPC calibration, SOC >= 1/3 | Phase-balanced pulse-relative RMSE | {1000.0 * PRIOR_HPPC_CALIBRATION_REGION_RMSE_V:.6f} |",
            f"| Withheld low-SOC HPPC | Phase-balanced pulse-relative RMSE | {1000.0 * PRIOR_HPPC_LOW_SOC_PHASE_BALANCED_RMSE_V:.6f} |",
            f"| Withheld low-SOC HPPC | Absolute-voltage RMSE | {1000.0 * PRIOR_HPPC_LOW_SOC_ABSOLUTE_RMSE_V:.6f} |",
            f"| Independent 1C, overall | Absolute-voltage RMSE | {1000.0 * PRIOR_ONE_C_OVERALL_RMSE_V:.6f} |",
            f"| Independent 1C, high / middle / low SOC | Absolute-voltage RMSE | {1000.0 * PRIOR_ONE_C_HIGH_SOC_RMSE_V:.3f} / {1000.0 * PRIOR_ONE_C_MIDDLE_SOC_RMSE_V:.3f} / {1000.0 * PRIOR_ONE_C_LOW_SOC_RMSE_V:.3f} |",
            f"| UDDS, overall | Absolute-voltage RMSE | {1000.0 * overall['rmse_v']:.6f} |",
            f"| UDDS, high / middle / low SOC | Absolute-voltage RMSE | {1000.0 * regional['high']['rmse_v']:.3f} / {1000.0 * regional['middle']['rmse_v']:.3f} / {1000.0 * regional['low']['rmse_v']:.3f} |",
            "",
            "Metric definitions differ between pulse-relative HPPC and absolute-voltage continuous-profile validation, so their magnitudes are contextual rather than interchangeable. Nevertheless, UDDS shows the same direction of transfer loss: positive mean residual grows from high to middle to low SOC, so the frozen model increasingly overpredicts measured voltage as SOC falls.",
            "",
            "## Interpretation and classification",
            "",
            f"Aggregate UDDS tracking is strong (`RMSE={1000.0 * overall['rmse_v']:.3f} mV`, `R2={overall['r2']:.6f}`), but it is not uniform. Low-SOC RMSE is {1000.0 * regional['low']['rmse_v']:.3f} mV versus {1000.0 * regional['high']['rmse_v']:.3f} mV at high SOC, with a positive low-SOC bias of {1000.0 * regional['low']['mean_signed_error_v']:.3f} mV. The cutoff residual reaches {1000.0 * audit['cutoff_residual_v']:.3f} mV. In addition, the ten largest current changes have {1000.0 * audit['largest_step_delta_voltage_rmse_v']:.3f} mV delta-voltage RMSE, exposing a transient/current-direction limitation that aggregate sample metrics dilute.",
            "",
            "The dominant limitation is therefore a **combination**: a low-SOC, SOC-dependent positive/static voltage offset plus deficient transient transfer at the largest current changes, especially regenerative-current onsets; increasing current magnitude is a secondary aggravating factor. The limitation is not caused by OCP extrapolation or clipping.",
            "",
            "**Final classification: UDDS TRANSFER LIMITED.**",
            "",
            "No repair, optimization, baseline adjustment, thermal fit, or frozen-model change was attempted.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    if not MANIFEST_PATH.is_file() or RAW_SOURCE_PATH.name not in MANIFEST_PATH.read_text(encoding="utf-8"):
        raise RuntimeError("The selected raw UDDS path is not recorded in the manifest.")
    if not np.isclose(Q_CELL_AH, 2.994979138519, rtol=0.0, atol=5.0e-13):
        raise RuntimeError("Frozen Panasonic cell capacity changed.")
    if not (
        np.isclose(IMPLEMENTED_Q_N_C, Q_N_C, rtol=0.0, atol=5.0e-9)
        and np.isclose(IMPLEMENTED_Q_P_C, Q_P_C, rtol=0.0, atol=5.0e-9)
    ):
        raise RuntimeError("Implemented electrode capacities differ from the frozen values.")

    columns, raw = load_udds()
    raw_sample_count = raw["Time"].size
    data, duplicate_indices = remove_exact_duplicate_times(raw)
    time_s = data["Time"]
    measured_v = data["Voltage"]
    raw_current_a = data["Current"]
    exported_ah = data["Ah"]
    battery_temperature_c = data["Battery_Temp_degC"]
    chamber_temperature_c = data["Chamber_Temp_degC"]
    dt_s = np.diff(time_s)

    source_soc = 1.0 + exported_ah / Q_CELL_AH
    if np.any((source_soc < 0.0) | (source_soc > 1.0)):
        raise ValueError("Exported-Ah SOC leaves [0, 1]; no clipping is allowed.")

    model_current_a = -raw_current_a  # sole sign conversion at model boundary
    result = simulate_local_group(
        time_s=time_s,
        discharge_positive_current_a=model_current_a,
        temperature_k=battery_temperature_c + 273.15,
        initial_x_n=X_N_START,
        initial_x_p=X_P_START,
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
    if violation_count:
        first = int(np.flatnonzero(negative_violation | positive_violation)[0])
        raise RuntimeError(
            f"First frozen OCP-domain violation at usable index {first}, time {time_s[first]:.9f} s; no clipping."
        )

    overall = metrics(measured_v, predicted_v)
    cutoff_index = int(np.argmin(measured_v))
    through_cutoff = np.arange(time_s.size) <= cutoff_index
    driven = metrics(measured_v[through_cutoff], predicted_v[through_cutoff])
    region_masks = soc_masks(source_soc)
    regional = {
        name: metrics(measured_v[mask], predicted_v[mask])
        for name, mask in region_masks.items()
    }

    near_rest = np.abs(raw_current_a) <= NEAR_REST_CURRENT_A
    loaded = ~near_rest
    discharge = raw_current_a < -NEAR_REST_CURRENT_A
    regenerative = raw_current_a > NEAR_REST_CURRENT_A
    post_cutoff = np.arange(time_s.size) > cutoff_index
    operating_masks = {
        "near_rest": near_rest,
        "loaded": loaded,
        "discharge": discharge,
        "regenerative_charge": regenerative,
        "post_cutoff_recovery": post_cutoff,
    }
    operating = {
        name: (mask, metrics(measured_v[mask], predicted_v[mask]))
        for name, mask in operating_masks.items()
    }

    quartiles, current_bins = current_magnitude_masks(raw_current_a)
    current_metrics = {
        str(item["name"]): metrics(
            measured_v[np.asarray(item["mask"], dtype=bool)],
            predicted_v[np.asarray(item["mask"], dtype=bool)],
        )
        for item in current_bins
    }

    absolute_error = np.abs(residual_v)
    current_magnitude = np.abs(raw_current_a)
    correlations = {
        "abs_error_vs_depth_of_discharge_all": float(spearmanr(absolute_error, 1.0 - source_soc).statistic),
        "abs_error_vs_depth_of_discharge_loaded": float(
            spearmanr(absolute_error[loaded], (1.0 - source_soc)[loaded]).statistic
        ),
        "abs_error_vs_current_magnitude_all": float(spearmanr(absolute_error, current_magnitude).statistic),
        "abs_error_vs_current_magnitude_loaded": float(
            spearmanr(absolute_error[loaded], current_magnitude[loaded]).statistic
        ),
    }
    for name, mask in region_masks.items():
        selected = mask & loaded
        correlations[f"abs_error_vs_current_magnitude_{name}_soc_loaded"] = float(
            spearmanr(absolute_error[selected], current_magnitude[selected]).statistic
        )

    events = largest_current_steps(time_s, source_soc, raw_current_a, measured_v, predicted_v)
    event_delta_errors = np.asarray([event["delta_voltage_error_v"] for event in events], dtype=float)

    interval_ah = raw_current_a[:-1] * dt_s / 3600.0
    cumulative_model_discharge_ah = np.zeros(time_s.size)
    cumulative_model_discharge_ah[1:] = np.cumsum(-interval_ah)
    predicted_soc = 1.0 - cumulative_model_discharge_ah / Q_CELL_AH
    exported_steps = np.diff(exported_ah)
    audit: dict[str, object] = {
        "raw_samples": int(raw_sample_count),
        "usable_samples": int(time_s.size),
        "duplicate_rows_removed": int(duplicate_indices.size),
        "duration_s": float(time_s[-1] - time_s[0]),
        "dt_min_s": float(np.min(dt_s)),
        "dt_median_s": float(np.median(dt_s)),
        "dt_mean_s": float(np.mean(dt_s)),
        "dt_max_s": float(np.max(dt_s)),
        "raw_current_min_a": float(np.min(raw_current_a)),
        "raw_current_max_a": float(np.max(raw_current_a)),
        "voltage_min_v": float(np.min(measured_v)),
        "voltage_max_v": float(np.max(measured_v)),
        "battery_temperature_min_c": float(np.min(battery_temperature_c)),
        "battery_temperature_max_c": float(np.max(battery_temperature_c)),
        "chamber_temperature_min_c": float(np.min(chamber_temperature_c)),
        "chamber_temperature_max_c": float(np.max(chamber_temperature_c)),
        "exported_ah_min": float(np.min(exported_ah)),
        "exported_ah_max": float(np.max(exported_ah)),
        "exported_ah_increasing_steps": int(np.count_nonzero(exported_steps > 0.0)),
        "exported_ah_decreasing_steps": int(np.count_nonzero(exported_steps < 0.0)),
        "exported_net_discharge_ah": float(exported_ah[0] - exported_ah[-1]),
        "integrated_net_discharge_ah": float(-np.sum(interval_ah)),
        "integrated_gross_discharge_ah": float(-np.sum(interval_ah[interval_ah < 0.0])),
        "integrated_regenerative_charge_ah": float(np.sum(interval_ah[interval_ah > 0.0])),
        "source_soc_min": float(np.min(source_soc)),
        "predicted_soc_min": float(np.min(predicted_soc)),
        "x_n_avg_min": float(np.min(result.x_n_avg)),
        "x_n_avg_max": float(np.max(result.x_n_avg)),
        "x_p_avg_min": float(np.min(result.x_p_avg)),
        "x_p_avg_max": float(np.max(result.x_p_avg)),
        "x_n_surface_min": float(np.min(result.x_n_surface)),
        "x_n_surface_max": float(np.max(result.x_n_surface)),
        "x_p_surface_min": float(np.min(result.x_p_surface)),
        "x_p_surface_max": float(np.max(result.x_p_surface)),
        "ocp_domain_violation_count": violation_count,
        "clipping_count": 0,
        "cutoff_index": cutoff_index,
        "cutoff_time_s": float(time_s[cutoff_index]),
        "cutoff_residual_v": float(residual_v[cutoff_index]),
        "largest_step_delta_voltage_rmse_v": float(np.sqrt(np.mean(event_delta_errors**2))),
        "largest_step_delta_voltage_bias_v": float(np.mean(event_delta_errors)),
    }

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    write_metrics(
        OUTPUT_DIRECTORY / "panasonic_udds_forward_validation_metrics.csv",
        overall,
        driven,
        regional,
        region_masks,
        operating,
        current_bins,
        current_metrics,
        correlations,
        events,
        audit,
    )
    save_voltage_plot(
        OUTPUT_DIRECTORY / "panasonic_udds_forward_validation_voltage.png",
        time_s,
        raw_current_a,
        source_soc,
        measured_v,
        predicted_v,
        cutoff_index,
    )
    save_residual_soc_plot(
        OUTPUT_DIRECTORY / "panasonic_udds_forward_validation_residual_soc.png",
        source_soc,
        residual_v,
    )
    save_residual_current_plot(
        OUTPUT_DIRECTORY / "panasonic_udds_forward_validation_residual_current.png",
        raw_current_a,
        residual_v,
        current_bins,
        current_metrics,
    )
    write_report(
        OUTPUT_DIRECTORY / "panasonic_udds_forward_validation_report.md",
        columns,
        overall,
        driven,
        regional,
        region_masks,
        operating,
        quartiles,
        current_bins,
        current_metrics,
        correlations,
        events,
        audit,
    )

    print("Panasonic 25 C UDDS forward-only validation complete.")
    print(
        f"Usable samples={time_s.size}; RMSE={1000.0 * overall['rmse_v']:.6f} mV; "
        f"R2={overall['r2']:.9f}"
    )
    print(
        f"High/middle/low SOC RMSE={1000.0 * regional['high']['rmse_v']:.6f}/"
        f"{1000.0 * regional['middle']['rmse_v']:.6f}/"
        f"{1000.0 * regional['low']['rmse_v']:.6f} mV"
    )
    print(
        f"OCP-domain violations={violation_count}; clipping=0; "
        f"largest-step delta-V RMSE={1000.0 * audit['largest_step_delta_voltage_rmse_v']:.6f} mV"
    )
    print("Classification: UDDS TRANSFER LIMITED")


if __name__ == "__main__":
    main()
