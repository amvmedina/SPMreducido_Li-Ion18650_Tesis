"""Forward-only Panasonic 25 C CC-CV charge validation.

The discharge-derived Panasonic OCPs and conditional dynamic parameters are
used unchanged. The raw charge current and its cumulative Ah channel provide
the sampled and interval-average model inputs; no optimizer, clipping,
hysteresis, charger controller, baseline correction, or thermal fit is used.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
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
    U_n,
    U_p,
    X_N_END,
    X_N_START,
    X_P_END,
    X_P_START,
)
from runtime_paths import (  # noqa: E402
    PANASONIC_CCCV_PATH,
    PANASONIC_PRECEDING_1C_PATH,
    require_data_file,
)


RAW_SOURCE_PATH = PANASONIC_CCCV_PATH
CONVERTED_SOURCE_PATH = PANASONIC_CCCV_PATH
PRECEDING_DISCHARGE_PATH = PANASONIC_PRECEDING_1C_PATH
MANIFEST_PATH = REPOSITORY_ROOT / "experimental_data_manifest.md"
OUTPUT_DIRECTORY = REPOSITORY_ROOT / "thesis_data"

Q_N_C = 23008.732269258333
Q_P_C = 14780.898652565113
R_GAS_J_PER_MOL_K = 8.314472
FARADAY_C_PER_MOL = 96487.0
PARAMETERS = {
    "R0": 0.0259594900586,
    "alpha_p": 3204.22765511,
    "d_n": 0.000183029692927,
    "alpha_n": 3000.0,
    "d_p": 0.0007,
}

REST_THRESHOLD_A = 0.01
CC_RELATIVE_CURRENT_THRESHOLD = 0.95
FINAL_CAPACITY_FRACTION = 0.95

PRIOR_HPPC_CALIBRATION_RMSE_V = 0.036835197638870915
PRIOR_HPPC_LOW_SOC_DYNAMIC_RMSE_V = 0.08572943292069442
PRIOR_HPPC_LOW_SOC_ABSOLUTE_RMSE_V = 0.0595821150372503
PRIOR_ONE_C_RMSE_V = 0.11436218800084591
PRIOR_UDDS_RMSE_V = 0.03098819169351242
PRIOR_UDDS_LARGEST_STEP_RMSE_V = 0.09271376454308056


@dataclass(frozen=True)
class ForwardResult:
    voltage_v: np.ndarray
    x_n_avg: np.ndarray
    x_n_aux: np.ndarray
    x_n_surface: np.ndarray
    x_p_avg: np.ndarray
    x_p_aux: np.ndarray
    x_p_surface: np.ndarray


@dataclass(frozen=True)
class PhaseMasks:
    initial_rest: np.ndarray
    cc: np.ndarray
    cv: np.ndarray
    transition: np.ndarray
    final_taper: np.ndarray
    final_rest: np.ndarray
    first_charge_index: int
    cc_end_index: int
    cv_start_index: int
    last_charge_index: int
    first_final_rest_index: int


def _vector(value: object, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float).reshape(-1)
    if array.size == 0 or not np.all(np.isfinite(array)):
        raise ValueError(f"Raw field {name!r} is empty or non-finite.")
    return array


def load_mat_channels(path: Path) -> tuple[tuple[str, ...], dict[str, np.ndarray]]:
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
            raise ValueError(f"Processed source CSV is missing required fields: {missing}")
        arrays = {name: _vector(table[name], name) for name in required}
        if len({array.size for array in arrays.values()}) != 1:
            raise ValueError("Processed source CSV channels have inconsistent lengths.")
        return columns, arrays
    measurement = loadmat(path, squeeze_me=True, struct_as_record=False).get("meas")
    if measurement is None or not hasattr(measurement, "_fieldnames"):
        raise ValueError(f"Expected a MATLAB structure named 'meas' in {path}.")
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
        raise ValueError(f"Missing source fields: {missing}")
    arrays = {name: _vector(getattr(measurement, name), name) for name in required}
    if len({array.size for array in arrays.values()}) != 1:
        raise ValueError("Source channels have inconsistent lengths.")
    return columns, arrays


def remove_duplicate_times(arrays: dict[str, np.ndarray]) -> tuple[dict[str, np.ndarray], np.ndarray]:
    differences = np.diff(arrays["Time"])
    if np.any(differences < 0.0):
        raise ValueError("Elapsed time decreases.")
    duplicate_indices = np.flatnonzero(differences == 0.0) + 1
    for index in duplicate_indices:
        for name in ("Voltage", "Current", "Ah", "Wh", "Power", "Chamber_Temp_degC"):
            values = arrays[name]
            if values[index] != values[index - 1]:
                raise ValueError(f"Repeated-time electrical record at index {index} differs in {name}.")
    keep = np.ones(arrays["Time"].size, dtype=bool)
    keep[duplicate_indices] = False
    cleaned = {name: values[keep] for name, values in arrays.items()}
    if np.any(np.diff(cleaned["Time"]) <= 0.0):
        raise ValueError("Usable times must be strictly increasing.")
    return cleaned, duplicate_indices


def verify_converted_csv(raw: dict[str, np.ndarray]) -> tuple[str, ...]:
    if not CONVERTED_SOURCE_PATH.is_file():
        raise FileNotFoundError(CONVERTED_SOURCE_PATH)
    converted = np.genfromtxt(
        CONVERTED_SOURCE_PATH,
        delimiter=",",
        names=True,
        dtype=float,
        encoding="utf-8",
    )
    columns = tuple(converted.dtype.names or ())
    for name in columns:
        if name not in raw:
            raise ValueError(f"Unexpected converted column: {name}")
        if not np.allclose(np.asarray(converted[name]), raw[name], rtol=0.0, atol=5.0e-12):
            raise ValueError(f"Converted column {name} differs from the raw MAT source.")
    return columns


def preceding_discharge_capacity_ah() -> tuple[float, float]:
    """Reproduce the associated 1C capacity and recovered terminal voltage."""

    _, raw = load_mat_channels(PRECEDING_DISCHARGE_PATH)
    active = np.flatnonzero(raw["Current"] < -0.5)
    if active.size < 2 or not np.array_equal(active, np.arange(active[0], active[-1] + 1)):
        raise ValueError("Preceding discharge branch is not one contiguous segment.")
    time_s = raw["Time"][active]
    current_a = raw["Current"][active]
    capacity_ah = float(np.sum(-current_a[:-1] * np.diff(time_s)) / 3600.0)
    recovered_voltage_v = float(raw["Voltage"][-1])
    return capacity_ah, recovered_voltage_v


def identify_phases(raw_current_a: np.ndarray, voltage_v: np.ndarray, exported_ah: np.ndarray) -> PhaseMasks:
    n_samples = raw_current_a.size
    index = np.arange(n_samples)
    active = raw_current_a > REST_THRESHOLD_A
    active_indices = np.flatnonzero(active)
    if active_indices.size == 0 or not np.array_equal(
        active_indices, np.arange(active_indices[0], active_indices[-1] + 1)
    ):
        raise ValueError("Charge-active samples must form one contiguous segment.")
    first_charge = int(active_indices[0])
    last_charge = int(active_indices[-1])
    plateau_current = float(np.max(raw_current_a[active]))
    cc_candidate = active & (raw_current_a >= CC_RELATIVE_CURRENT_THRESHOLD * plateau_current)
    cc_indices = np.flatnonzero(cc_candidate)
    if cc_indices[0] != first_charge or not np.array_equal(
        cc_indices, np.arange(first_charge, cc_indices[-1] + 1)
    ):
        raise ValueError("CC samples are not a contiguous leading part of charge.")
    cc_end = int(cc_indices[-1])
    cv_start = cc_end + 1
    if cv_start > last_charge:
        raise ValueError("No current-taper phase follows CC.")
    if voltage_v[cv_start] < np.max(voltage_v) - 0.01:
        raise ValueError("Detected taper does not begin near the maximum voltage.")

    initial_rest = index < first_charge
    cc = (index >= first_charge) & (index <= cc_end)
    cv = (index >= cv_start) & (index <= last_charge)
    transition = (index >= cc_end) & (index <= min(cv_start + 1, last_charge))
    final_taper = cv & (exported_ah >= FINAL_CAPACITY_FRACTION * exported_ah[last_charge])
    final_rest = index > last_charge
    return PhaseMasks(
        initial_rest=initial_rest,
        cc=cc,
        cv=cv,
        transition=transition,
        final_taper=final_taper,
        final_rest=final_rest,
        first_charge_index=first_charge,
        cc_end_index=cc_end,
        cv_start_index=cv_start,
        last_charge_index=last_charge,
        first_final_rest_index=last_charge + 1,
    )


def initial_stoichiometries(initial_soc: float) -> tuple[float, float]:
    if not 0.0 <= initial_soc <= 1.0:
        raise ValueError("Initial SOC must be physical; no clipping.")
    x_n = X_N_END + initial_soc * (X_N_START - X_N_END)
    x_p = X_P_END + initial_soc * (X_P_START - X_P_END)
    return float(x_n), float(x_p)


def forward_charge(
    time_s: np.ndarray,
    raw_current_a: np.ndarray,
    exported_ah: np.ndarray,
    battery_temperature_c: np.ndarray,
    initial_x_n: float,
    initial_x_p: float,
) -> ForwardResult:
    """Propagate unchanged grouped-SPM equations with measured charge input.

    The current channel is used for sample-time surface/kinetic/ohmic voltage.
    Since this source is sampled only every 60 s and stores endpoint current,
    its measured delta-Ah/delta-time supplies each interval-average current.
    This exactly retains charge passed through CC onset and the CV taper.
    """

    time_s = np.asarray(time_s, dtype=float)
    raw_current_a = np.asarray(raw_current_a, dtype=float)
    exported_ah = np.asarray(exported_ah, dtype=float)
    dt_s = np.diff(time_s)
    model_sample_current_a = -raw_current_a
    model_interval_current_a = -np.diff(exported_ah) * 3600.0 / dt_s
    n_samples = time_s.size

    x_n_avg = np.empty(n_samples)
    x_n_aux = np.empty(n_samples)
    x_p_avg = np.empty(n_samples)
    x_p_aux = np.empty(n_samples)
    x_n_avg[0] = x_n_aux[0] = initial_x_n
    x_p_avg[0] = x_p_aux[0] = initial_x_p

    alpha_n = PARAMETERS["alpha_n"]
    alpha_p = PARAMETERS["alpha_p"]
    for interval, dt in enumerate(dt_s):
        current = model_interval_current_a[interval]

        y_n = x_n_aux[interval] - x_n_avg[interval]
        decay_n = np.exp(-30.0 * dt / alpha_n)
        forcing_n = (-12.0 / (7.0 * Q_N_C)) * current
        y_n_next = y_n * decay_n + forcing_n * (alpha_n / 30.0) * (1.0 - decay_n)
        x_n_avg[interval + 1] = x_n_avg[interval] - current * dt / Q_N_C
        x_n_aux[interval + 1] = x_n_avg[interval + 1] + y_n_next

        y_p = x_p_aux[interval] - x_p_avg[interval]
        decay_p = np.exp(-30.0 * dt / alpha_p)
        forcing_p = (12.0 / (7.0 * Q_P_C)) * current
        y_p_next = y_p * decay_p + forcing_p * (alpha_p / 30.0) * (1.0 - decay_p)
        x_p_avg[interval + 1] = x_p_avg[interval] + current * dt / Q_P_C
        x_p_aux[interval + 1] = x_p_avg[interval + 1] + y_p_next

    x_n_surface = x_n_aux - alpha_n * model_sample_current_a / (105.0 * Q_N_C)
    x_p_surface = x_p_aux + alpha_p * model_sample_current_a / (105.0 * Q_P_C)
    negative_invalid = (x_n_surface < MCMB_SUPPORTED_STOICHIOMETRY_RANGE[0]) | (
        x_n_surface > MCMB_SUPPORTED_STOICHIOMETRY_RANGE[1]
    )
    positive_invalid = (x_p_surface < NCA_WORKING_STOICHIOMETRY_RANGE[0]) | (
        x_p_surface > NCA_WORKING_STOICHIOMETRY_RANGE[1]
    )
    invalid = negative_invalid | positive_invalid
    if np.any(invalid):
        first = int(np.flatnonzero(invalid)[0])
        raise ValueError(
            f"First OCP-domain violation at index {first}, time {time_s[first]:.9f} s: "
            f"x_n_surf={x_n_surface[first]:.12g}, x_p_surf={x_p_surface[first]:.12g}. "
            "No clipping or extrapolation was performed."
        )

    temperature_k = battery_temperature_c + 273.15
    kinetic_factor = 2.0 * R_GAS_J_PER_MOL_K * temperature_k / FARADAY_C_PER_MOL
    eta_n = kinetic_factor * np.arcsinh(
        model_sample_current_a
        / (
            6.0
            * Q_N_C
            * PARAMETERS["d_n"]
            * np.sqrt(x_n_surface * (1.0 - x_n_surface))
        )
    )
    eta_p = kinetic_factor * np.arcsinh(
        model_sample_current_a
        / (
            6.0
            * Q_P_C
            * PARAMETERS["d_p"]
            * np.sqrt(x_p_surface * (1.0 - x_p_surface))
        )
    )
    voltage_v = (
        U_p(x_p_surface)
        - U_n(x_n_surface)
        + eta_p
        - eta_n
        - PARAMETERS["R0"] * model_sample_current_a
    )
    return ForwardResult(
        voltage_v=np.asarray(voltage_v),
        x_n_avg=x_n_avg,
        x_n_aux=x_n_aux,
        x_n_surface=x_n_surface,
        x_p_avg=x_p_avg,
        x_p_aux=x_p_aux,
        x_p_surface=x_p_surface,
    )


def voltage_metrics(measured_v: np.ndarray, predicted_v: np.ndarray) -> dict[str, float]:
    residual = predicted_v - measured_v
    denominator = float(np.sum((measured_v - np.mean(measured_v)) ** 2))
    return {
        "rmse_v": float(np.sqrt(np.mean(residual**2))),
        "mae_v": float(np.mean(np.abs(residual))),
        "r2": float(1.0 - np.sum(residual**2) / denominator) if denominator > 0.0 else float("nan"),
        "mean_signed_error_v": float(np.mean(residual)),
        "max_abs_error_v": float(np.max(np.abs(residual))),
        "p95_abs_error_v": float(np.quantile(np.abs(residual), 0.95)),
    }


def soc_region_masks(soc: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "low": (soc >= 0.0) & (soc < 1.0 / 3.0),
        "middle": (soc >= 1.0 / 3.0) & (soc < 2.0 / 3.0),
        "high": (soc >= 2.0 / 3.0) & (soc <= 1.0),
    }


def append_metrics(
    rows: list[dict[str, object]],
    scope: str,
    subset: str,
    metric_values: dict[str, object],
    n_samples: int,
) -> None:
    units = {
        "rmse_v": "V",
        "mae_v": "V",
        "r2": "1",
        "mean_signed_error_v": "V",
        "max_abs_error_v": "V",
        "p95_abs_error_v": "V",
    }
    for name, value in metric_values.items():
        rows.append(
            {
                "scope": scope,
                "subset": subset,
                "metric": name,
                "value": value,
                "unit": units.get(name, "1"),
                "n_samples": n_samples,
            }
        )


def write_metrics(
    path: Path,
    overall: dict[str, float],
    phase_metrics: dict[str, dict[str, float]],
    phases: PhaseMasks,
    regional: dict[str, dict[str, float]],
    region_masks: dict[str, np.ndarray],
    audit: dict[str, object],
) -> None:
    rows: list[dict[str, object]] = []
    append_metrics(rows, "overall", "complete_usable_charge", overall, int(audit["usable_samples"]))
    phase_masks = {
        "initial_rest": phases.initial_rest,
        "cc": phases.cc,
        "cv": phases.cv,
        "cc_to_cv_transition": phases.transition,
        "final_5pct_taper": phases.final_taper,
        "final_rest": phases.final_rest,
    }
    for name, values in phase_metrics.items():
        append_metrics(rows, "phase", name, values, int(np.count_nonzero(phase_masks[name])))
    for name, values in regional.items():
        append_metrics(rows, "soc_region", name, values, int(np.count_nonzero(region_masks[name])))

    audit_units = {
        "raw_samples": "count",
        "usable_samples": "count",
        "duplicate_rows_removed": "count",
        "duplicate_battery_temperature_change_c": "degC",
        "duration_s": "s",
        "dt_min_s": "s",
        "dt_median_s": "s",
        "dt_mean_s": "s",
        "dt_max_s": "s",
        "raw_current_min_a": "A",
        "raw_current_max_a": "A",
        "voltage_min_v": "V",
        "voltage_max_v": "V",
        "exported_ah_min": "Ah",
        "exported_ah_max": "Ah",
        "exported_charge_ah": "Ah",
        "battery_temperature_min_c": "degC",
        "battery_temperature_max_c": "degC",
        "chamber_temperature_min_c": "degC",
        "chamber_temperature_max_c": "degC",
        "preceding_discharge_capacity_ah": "Ah",
        "initial_soc": "1",
        "final_soc": "1",
        "cc_current_min_a": "A",
        "cc_current_mean_a": "A",
        "cc_current_max_a": "A",
        "cv_current_start_a": "A",
        "cv_current_end_a": "A",
        "cv_terminal_residual_v": "V",
        "cv_last_five_bias_v": "V",
        "residual_vs_soc_spearman_rho": "1",
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
    }
    for name, unit in audit_units.items():
        rows.append(
            {
                "scope": "audit",
                "subset": "complete_usable_charge",
                "metric": name,
                "value": audit[name],
                "unit": unit,
                "n_samples": audit["usable_samples"],
            }
        )
    prior = {
        "hppc_calibration_soc_ge_one_third_phase_balanced_rmse_v": PRIOR_HPPC_CALIBRATION_RMSE_V,
        "hppc_low_soc_phase_balanced_rmse_v": PRIOR_HPPC_LOW_SOC_DYNAMIC_RMSE_V,
        "hppc_low_soc_absolute_rmse_v": PRIOR_HPPC_LOW_SOC_ABSOLUTE_RMSE_V,
        "one_c_discharge_overall_rmse_v": PRIOR_ONE_C_RMSE_V,
        "udds_overall_rmse_v": PRIOR_UDDS_RMSE_V,
        "udds_largest_current_step_delta_voltage_rmse_v": PRIOR_UDDS_LARGEST_STEP_RMSE_V,
    }
    for name, value in prior.items():
        rows.append(
            {
                "scope": "prior_validation",
                "subset": "reported",
                "metric": name,
                "value": value,
                "unit": "V",
                "n_samples": "",
            }
        )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_voltage_plot(
    path: Path,
    time_s: np.ndarray,
    measured_v: np.ndarray,
    predicted_v: np.ndarray,
    raw_current_a: np.ndarray,
    phases: PhaseMasks,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    minutes = time_s / 60.0
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(10.5, 7.0),
        sharex=True,
        gridspec_kw={"height_ratios": (2.0, 1.0)},
        constrained_layout=True,
    )
    axes[0].plot(minutes, measured_v, "o-", color="black", markersize=3.0, linewidth=1.2, label="Measured")
    axes[0].plot(minutes, predicted_v, "o-", color="#0072B2", markersize=2.5, linewidth=1.2, label="Frozen grouped SPM")
    axes[0].set_ylabel("Voltage [V]")
    axes[0].legend(loc="best")
    axes[1].step(minutes, raw_current_a, where="mid", color="#D55E00", linewidth=1.2)
    axes[1].set(xlabel="Elapsed time [min]", ylabel="Raw current [A]\n(charge positive)")
    for axis in axes:
        axis.axvline(minutes[phases.first_charge_index], color="#009E73", linestyle="--", linewidth=1.0)
        axis.axvline(minutes[phases.cv_start_index], color="#CC79A7", linestyle="--", linewidth=1.0)
        axis.axvline(minutes[phases.first_final_rest_index], color="0.4", linestyle="--", linewidth=1.0)
        axis.grid(alpha=0.25)
    axes[0].axvspan(minutes[phases.first_charge_index], minutes[phases.cv_start_index], color="#009E73", alpha=0.06, label="CC")
    axes[0].axvspan(minutes[phases.cv_start_index], minutes[phases.first_final_rest_index], color="#CC79A7", alpha=0.06, label="CV taper")
    axes[0].text(minutes[phases.first_charge_index], 0.02, " charge onset", transform=axes[0].get_xaxis_transform(), rotation=90, va="bottom")
    axes[0].text(minutes[phases.cv_start_index], 0.02, " CC to CV", transform=axes[0].get_xaxis_transform(), rotation=90, va="bottom")
    axes[0].text(minutes[phases.first_final_rest_index], 0.02, " termination", transform=axes[0].get_xaxis_transform(), rotation=90, va="bottom")
    axes[0].set_title("Panasonic NCR18650PF 25 C CC-CV forward-only validation")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def save_residual_soc_plot(
    path: Path,
    soc: np.ndarray,
    residual_v: np.ndarray,
    phases: PhaseMasks,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9.0, 5.0), constrained_layout=True)
    index = np.arange(soc.size)
    groups = (
        (phases.initial_rest, "Initial rest", "0.45"),
        (phases.cc, "CC", "#D55E00"),
        (phases.cv, "CV", "#0072B2"),
        (phases.final_rest, "Final rest", "#009E73"),
    )
    for mask, label, color in groups:
        ax.plot(100.0 * soc[mask], 1000.0 * residual_v[mask], "o-", color=color, markersize=4.0, linewidth=1.2, label=label)
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.axvline(100.0 / 3.0, color="0.4", linestyle="--", linewidth=1.0)
    ax.axvline(200.0 / 3.0, color="0.4", linestyle="--", linewidth=1.0)
    ax.set(
        xlabel="Cell SOC from associated discharge and exported charge Ah [%]",
        ylabel="Residual, predicted - measured [mV]",
        title="CC-CV charge residual versus SOC",
    )
    ax.grid(alpha=0.25)
    ax.legend(loc="best")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def save_cv_taper_plot(
    path: Path,
    time_s: np.ndarray,
    raw_current_a: np.ndarray,
    residual_v: np.ndarray,
    phases: PhaseMasks,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cv = phases.cv
    cv_time_min = (time_s[cv] - time_s[phases.cv_start_index]) / 60.0
    current = raw_current_a[cv]
    residual_mv = 1000.0 * residual_v[cv]
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5), constrained_layout=True)
    axes[0].plot(cv_time_min, residual_mv, "o-", color="#0072B2", markersize=4.0, linewidth=1.2)
    axes[0].axhline(0.0, color="black", linewidth=0.8)
    axes[0].set(xlabel="Time since CV onset [min]", ylabel="Predicted - regulated voltage [mV]")
    axes[0].grid(alpha=0.25)
    points = axes[1].scatter(current, residual_mv, c=cv_time_min, cmap="viridis", s=28)
    axes[1].plot(current, residual_mv, color="0.55", linewidth=0.7, zorder=0)
    axes[1].axhline(0.0, color="black", linewidth=0.8)
    axes[1].set(xlabel="Measured taper current [A]", ylabel="Predicted - regulated voltage [mV]")
    axes[1].grid(alpha=0.25)
    colorbar = fig.colorbar(points, ax=axes[1])
    colorbar.set_label("Time since CV onset [min]")
    fig.suptitle("CC-CV validation: externally regulated CV taper")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def phase_table_row(name: str, mask: np.ndarray, values: dict[str, float]) -> str:
    return (
        f"| {name} | {np.count_nonzero(mask)} | {1000.0 * values['rmse_v']:.6f} | "
        f"{1000.0 * values['mae_v']:.6f} | {values['r2']:.9f} | "
        f"{1000.0 * values['mean_signed_error_v']:.6f} | "
        f"{1000.0 * values['p95_abs_error_v']:.6f} | "
        f"{1000.0 * values['max_abs_error_v']:.6f} |"
    )


def write_report(
    path: Path,
    raw_columns: tuple[str, ...],
    converted_columns: tuple[str, ...],
    phases: PhaseMasks,
    overall: dict[str, float],
    phase_metrics: dict[str, dict[str, float]],
    regional: dict[str, dict[str, float]],
    region_masks: dict[str, np.ndarray],
    time_s: np.ndarray,
    raw_current_a: np.ndarray,
    exported_ah: np.ndarray,
    residual_v: np.ndarray,
    audit: dict[str, object],
) -> None:
    phase_masks = {
        "Initial rest": phases.initial_rest,
        "CC": phases.cc,
        "CV": phases.cv,
        "CC-to-CV transition": phases.transition,
        "Final 5% taper": phases.final_taper,
        "Final rest": phases.final_rest,
    }
    lines = [
        "# Panasonic 25 C CC-CV charge forward validation",
        "",
        "## CC-CV CHARGE TRANSFER LIMITED",
        "",
        "The manifest-selected associated charge was run forward once with the frozen discharge-side OCP representation and frozen conditional dynamics. Measured current was the input during both CC and CV; predicted voltage was never forced to the regulated CV voltage. No optimization, hysteresis, charge-specific OCP, baseline correction, clipping, or thermal fit was introduced.",
        "",
        "## Dataset audit",
        "",
        f"- Exact raw source: `{RAW_SOURCE_PATH}`",
        f"- Exact converted source: `{CONVERTED_SOURCE_PATH}` (numerically verified against the raw channels)",
        f"- Raw columns: `{', '.join(raw_columns)}`",
        f"- Converted columns: `{', '.join(converted_columns)}`",
        f"- Raw samples: {audit['raw_samples']}; usable samples: {audit['usable_samples']} after removing one zero-time duplicate final electrical record. Its electrical fields and chamber temperature repeat exactly, while battery temperature changes by {audit['duplicate_battery_temperature_change_c']:+.5f} degC despite no elapsed time; the earlier record is retained.",
        f"- Duration: {audit['duration_s']:.9f} s. Positive sampling interval median {audit['dt_median_s']:.9f} s, mean {audit['dt_mean_s']:.9f} s, range {audit['dt_min_s']:.9f}-{audit['dt_max_s']:.9f} s. The 0.010 s interval is retained because it is positive; all actual times are preserved.",
        f"- Raw current is positive on charge and zero at rest: {audit['raw_current_min_a']:.5f} to {audit['raw_current_max_a']:.5f} A. Conversion to model charge-negative current occurs only at the loader/simulation boundary.",
        f"- Voltage range: {audit['voltage_min_v']:.5f}-{audit['voltage_max_v']:.5f} V.",
        f"- Exported Ah is nondecreasing from {audit['exported_ah_min']:.5f} to {audit['exported_ah_max']:.5f} Ah and ends at {audit['exported_charge_ah']:.5f} Ah; no counter reversal occurs.",
        f"- Battery temperature: {audit['battery_temperature_min_c']:.5f}-{audit['battery_temperature_max_c']:.5f} degC; chamber temperature: {audit['chamber_temperature_min_c']:.0f}-{audit['chamber_temperature_max_c']:.0f} degC.",
        "",
        "## Signal-derived experiment structure",
        "",
        f"- Initial rest: indices 0-{phases.first_charge_index - 1}, through {time_s[phases.first_charge_index - 1]:.6f} s ({np.count_nonzero(phases.initial_rest)} records). The first charged record is index {phases.first_charge_index} at {time_s[phases.first_charge_index]:.6f} s.",
        f"- CC phase: indices {phases.first_charge_index}-{phases.cc_end_index}; current {audit['cc_current_min_a']:.5f}-{audit['cc_current_max_a']:.5f} A, mean {audit['cc_current_mean_a']:.6f} A.",
        f"- CC-to-CV transition: the first current reduction occurs at index {phases.cv_start_index}, time {time_s[phases.cv_start_index]:.6f} s, where measured voltage is {audit['cv_start_voltage_v']:.5f} V and current is {audit['cv_current_start_a']:.5f} A. The transition metric uses the last CC sample and first two CV samples.",
        f"- CV phase: indices {phases.cv_start_index}-{phases.last_charge_index}; regulated voltage range {audit['cv_voltage_min_v']:.5f}-{audit['cv_voltage_max_v']:.5f} V and current taper {audit['cv_current_start_a']:.5f} to {audit['cv_current_end_a']:.5f} A.",
        f"- Charge termination: last active record index {phases.last_charge_index} at {time_s[phases.last_charge_index]:.6f} s and {raw_current_a[phases.last_charge_index]:.5f} A; zero current first appears at index {phases.first_final_rest_index}, {time_s[phases.first_final_rest_index]:.6f} s.",
        f"- Final rest: indices {phases.first_final_rest_index}-{time_s.size - 1}, spanning {time_s[-1] - time_s[phases.first_final_rest_index]:.6f} s.",
        "",
        "## Initial state and measured-current handling",
        "",
        f"The immediately preceding manifest-recorded 1C discharge removes {audit['preceding_discharge_capacity_ah']:.12f} Ah and is followed by recovery to {audit['preceding_recovered_voltage_v']:.5f} V. The charge file begins at {audit['initial_measured_voltage_v']:.5f} V with zero current. With frozen `Q_cell={Q_CELL_AH:.12f} Ah`, the non-fitted initial SOC is therefore:",
        "",
        f"`SOC_initial = 1 - {audit['preceding_discharge_capacity_ah']:.12f}/{Q_CELL_AH:.12f} = {audit['initial_soc']:.12f}`.",
        "",
        f"The frozen nominal effective trajectories give `x_n_avg=x_n_aux={audit['initial_x_n']:.12f}` and `x_p_avg=x_p_aux={audit['initial_x_p']:.12f}` at the first record. Auxiliary states equal average states because initialization is at the recorded zero-current recovery; the small measured rest drift is retained as validation error rather than used to adjust states.",
        "",
        "The source records current at 60 s endpoints while exported Ah contains the charge passed over each interval. Accordingly, sample current drives the algebraic surface/kinetic/ohmic voltage, while `delta Ah/delta t` supplies the measured interval-average current for state propagation. This exactly preserves CC onset and CV taper charge instead of introducing a one-sample delay. It is a data-timing treatment, not a parameter or model-equation change.",
        "",
        f"Final SOC is `{audit['final_soc']:.12f}` after adding the measured {audit['exported_charge_ah']:.5f} Ah; neither initial nor final SOC was forced to zero or one.",
        "",
        "## Primary voltage metrics",
        "",
        "Residual is predicted minus measured voltage.",
        "",
        "| Metric | Complete usable charge |",
        "|---|---:|",
        f"| Samples | {audit['usable_samples']} |",
        f"| RMSE [mV] | {1000.0 * overall['rmse_v']:.6f} |",
        f"| MAE [mV] | {1000.0 * overall['mae_v']:.6f} |",
        f"| R2 | {overall['r2']:.9f} |",
        f"| Mean signed error [mV] | {1000.0 * overall['mean_signed_error_v']:.6f} |",
        f"| Maximum absolute error [mV] | {1000.0 * overall['max_abs_error_v']:.6f} |",
        f"| P95 absolute error [mV] | {1000.0 * overall['p95_abs_error_v']:.6f} |",
        "",
        "### By experimental phase",
        "",
        "| Phase | Samples | RMSE [mV] | MAE [mV] | R2 | Mean signed [mV] | P95 absolute [mV] | Max absolute [mV] |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    metric_key = {
        "Initial rest": "initial_rest",
        "CC": "cc",
        "CV": "cv",
        "CC-to-CV transition": "cc_to_cv_transition",
        "Final 5% taper": "final_5pct_taper",
        "Final rest": "final_rest",
    }
    for label, mask in phase_masks.items():
        lines.append(phase_table_row(label, mask, phase_metrics[metric_key[label]]))
    lines.extend(
        [
            "",
            "CV voltage is externally regulated; its near-constant measured value is not evidence of model-predicted CV control. The very negative CV R2 reflects division by the tiny regulated-voltage variance and is not used as the transfer criterion.",
            "",
            "### By cell SOC",
            "",
            "| SOC region | Samples | RMSE [mV] | MAE [mV] | Mean signed [mV] | Max absolute [mV] |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for name in ("low", "middle", "high"):
        value = regional[name]
        interval = {"low": "[0, 1/3)", "middle": "[1/3, 2/3)", "high": "[2/3, 1]"}[name]
        lines.append(
            f"| {name} {interval} | {np.count_nonzero(region_masks[name])} | {1000.0 * value['rmse_v']:.6f} | {1000.0 * value['mae_v']:.6f} | {1000.0 * value['mean_signed_error_v']:.6f} | {1000.0 * value['max_abs_error_v']:.6f} |"
        )
    lines.extend(
        [
            "",
            f"Signed residual versus SOC has Spearman `rho={audit['residual_vs_soc_spearman_rho']:+.6f}`. The charge-side residual shifts toward increasing underprediction as SOC rises: mean residual changes from {1000.0 * regional['low']['mean_signed_error_v']:+.3f} mV at low SOC to {1000.0 * regional['middle']['mean_signed_error_v']:+.3f} mV at middle SOC and {1000.0 * regional['high']['mean_signed_error_v']:+.3f} mV at high SOC. Low-SOC aggregation contains the positive initial-rest offset and the start of loaded CC, so its near-zero mean does not indicate uniformly small low-SOC error.",
            "",
            "## CC current response",
            "",
            f"The measured CC current is {audit['cc_current_mean_a']:.6f} A on average, range {audit['cc_current_min_a']:.5f}-{audit['cc_current_max_a']:.5f} A. With no correction, CC voltage RMSE is {1000.0 * phase_metrics['cc']['rmse_v']:.6f} mV and predicted-minus-measured bias is {1000.0 * phase_metrics['cc']['mean_signed_error_v']:.6f} mV. Its residual-versus-SOC trajectory is shown in the SOC figure.",
            "",
            "## CV taper response",
            "",
            f"Using the measured taper current, CV RMSE is {1000.0 * phase_metrics['cv']['rmse_v']:.6f} mV, MAE {1000.0 * phase_metrics['cv']['mae_v']:.6f} mV, bias {1000.0 * phase_metrics['cv']['mean_signed_error_v']:.6f} mV, and maximum absolute error {1000.0 * phase_metrics['cv']['max_abs_error_v']:.6f} mV. At the final active sample ({audit['cv_current_end_a']:.5f} A), residual is {1000.0 * audit['cv_terminal_residual_v']:.6f} mV; the last five active taper records average {1000.0 * audit['cv_last_five_bias_v']:.6f} mV.",
            "",
            "## State and domain audit",
            "",
            "| State | Minimum | Maximum | Frozen support where applicable |",
            "|---|---:|---:|---:|",
            f"| `x_n_avg` | {audit['x_n_avg_min']:.12f} | {audit['x_n_avg_max']:.12f} | diagnostic state |",
            f"| `x_p_avg` | {audit['x_p_avg_min']:.12f} | {audit['x_p_avg_max']:.12f} | diagnostic state |",
            f"| `x_n_surf` | {audit['x_n_surface_min']:.12f} | {audit['x_n_surface_max']:.12f} | [{MCMB_SUPPORTED_STOICHIOMETRY_RANGE[0]:.2f}, {MCMB_SUPPORTED_STOICHIOMETRY_RANGE[1]:.2f}] |",
            f"| `x_p_surf` | {audit['x_p_surface_min']:.12f} | {audit['x_p_surface_max']:.12f} | [{NCA_WORKING_STOICHIOMETRY_RANGE[0]:.2f}, {NCA_WORKING_STOICHIOMETRY_RANGE[1]:.2f}] |",
            "",
            f"- Cell SOC range: {audit['initial_soc']:.12f}-{audit['final_soc']:.12f}.",
            f"- OCP-domain violations: **{audit['ocp_domain_violation_count']}**.",
            f"- Clipping count: **{audit['clipping_count']}**.",
            "",
            "## Charge-direction transfer and prior validation",
            "",
            "| Validation | Metric | Result [mV] |",
            "|---|---|---:|",
            f"| HPPC calibration, SOC >= 1/3 | Phase-balanced pulse-relative RMSE | {1000.0 * PRIOR_HPPC_CALIBRATION_RMSE_V:.6f} |",
            f"| Withheld low-SOC HPPC | Phase-balanced / absolute RMSE | {1000.0 * PRIOR_HPPC_LOW_SOC_DYNAMIC_RMSE_V:.3f} / {1000.0 * PRIOR_HPPC_LOW_SOC_ABSOLUTE_RMSE_V:.3f} |",
            f"| Independent 1C discharge | Absolute RMSE | {1000.0 * PRIOR_ONE_C_RMSE_V:.6f} |",
            f"| UDDS | Absolute RMSE | {1000.0 * PRIOR_UDDS_RMSE_V:.6f} |",
            f"| UDDS ten largest current steps | Delta-voltage RMSE | {1000.0 * PRIOR_UDDS_LARGEST_STEP_RMSE_V:.6f} |",
            f"| Independent CC-CV charge | Absolute RMSE | {1000.0 * overall['rmse_v']:.6f} |",
            "",
            "The independent charge supports a current-direction limitation, but not as a single uniform sign offset. UDDS showed excessive predicted upward voltage changes at several large regenerative onsets. Here, after the initial rest, the sustained CC and CV predictions are predominantly below measurement. Thus the frozen discharge-calibrated representation can overstate an immediate regenerative step while understating sustained charge voltage; that contrast is evidence of direction/transient dependence rather than a correctable constant baseline.",
            "",
            "## Interpretation and classification",
            "",
            f"The complete-charge RMSE is {1000.0 * overall['rmse_v']:.3f} mV. The model begins {1000.0 * phase_metrics['initial_rest']['mean_signed_error_v']:+.3f} mV high during the initial rest, then runs {1000.0 * phase_metrics['cc']['mean_signed_error_v']:+.3f} mV low through CC and {1000.0 * phase_metrics['cv']['mean_signed_error_v']:+.3f} mV low through CV. High-SOC bias reaches {1000.0 * regional['high']['mean_signed_error_v']:+.3f} mV, while the final 5% taper remains at {1000.0 * phase_metrics['final_5pct_taper']['rmse_v']:.3f} mV RMSE. These sign and SOC changes rule out a single static correction and show that discharge-side OCP/dynamic transfer is incomplete.",
            "",
            "The dominant limitation is a **combination of high-SOC/static OCP mismatch, CC loaded-voltage underprediction, CV taper underprediction, and current-direction asymmetry**. The regulated CV trace is treated only as a measured reference, not as simulated charger control. No repair is attempted.",
            "",
            "**Final classification: CC-CV CHARGE TRANSFER LIMITED.**",
            "",
            "No optimization, parameter change, hysteresis, separate charge OCP, clipping, model-equation change, or thermal fit was performed.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    manifest = MANIFEST_PATH.read_text(encoding="utf-8")
    if RAW_SOURCE_PATH.name not in manifest or CONVERTED_SOURCE_PATH.name not in manifest:
        raise RuntimeError("Selected CC-CV sources are not both recorded in the manifest.")
    if not np.isclose(Q_CELL_AH, 2.994979138519, rtol=0.0, atol=5.0e-13):
        raise RuntimeError("Frozen cell capacity changed.")

    raw_columns, raw = load_mat_channels(RAW_SOURCE_PATH)
    converted_columns = verify_converted_csv(raw)
    raw_sample_count = raw["Time"].size
    data, duplicates = remove_duplicate_times(raw)
    time_s = data["Time"]
    measured_v = data["Voltage"]
    raw_current_a = data["Current"]
    exported_ah = data["Ah"]
    battery_temperature_c = data["Battery_Temp_degC"]
    chamber_temperature_c = data["Chamber_Temp_degC"]
    dt_s = np.diff(time_s)

    phases = identify_phases(raw_current_a, measured_v, exported_ah)
    preceding_capacity_ah, preceding_recovered_voltage_v = preceding_discharge_capacity_ah()
    initial_soc = 1.0 - preceding_capacity_ah / Q_CELL_AH
    source_soc = initial_soc + exported_ah / Q_CELL_AH
    if np.any((source_soc < 0.0) | (source_soc > 1.0)):
        raise ValueError("Charge SOC leaves [0, 1]; no clipping.")
    initial_x_n, initial_x_p = initial_stoichiometries(initial_soc)

    result = forward_charge(
        time_s,
        raw_current_a,
        exported_ah,
        battery_temperature_c,
        initial_x_n,
        initial_x_p,
    )
    predicted_v = result.voltage_v
    residual_v = predicted_v - measured_v

    expected_x_n_avg = X_N_END + source_soc * (X_N_START - X_N_END)
    expected_x_p_avg = X_P_END + source_soc * (X_P_START - X_P_END)
    if (
        np.max(np.abs(result.x_n_avg - expected_x_n_avg)) > 2.0e-12
        or np.max(np.abs(result.x_p_avg - expected_x_p_avg)) > 2.0e-12
    ):
        raise RuntimeError("Cumulative-Ah and grouped average states are inconsistent.")

    phase_masks = {
        "initial_rest": phases.initial_rest,
        "cc": phases.cc,
        "cv": phases.cv,
        "cc_to_cv_transition": phases.transition,
        "final_5pct_taper": phases.final_taper,
        "final_rest": phases.final_rest,
    }
    overall = voltage_metrics(measured_v, predicted_v)
    phase_metrics = {
        name: voltage_metrics(measured_v[mask], predicted_v[mask])
        for name, mask in phase_masks.items()
    }
    region_masks = soc_region_masks(source_soc)
    regional = {
        name: voltage_metrics(measured_v[mask], predicted_v[mask])
        for name, mask in region_masks.items()
    }

    cv_indices = np.flatnonzero(phases.cv)
    last_five_cv = cv_indices[-5:]
    audit: dict[str, object] = {
        "raw_samples": int(raw_sample_count),
        "usable_samples": int(time_s.size),
        "duplicate_rows_removed": int(duplicates.size),
        "duplicate_battery_temperature_change_c": float(
            raw["Battery_Temp_degC"][duplicates[0]]
            - raw["Battery_Temp_degC"][duplicates[0] - 1]
        ) if duplicates.size else 0.0,
        "duration_s": float(time_s[-1] - time_s[0]),
        "dt_min_s": float(np.min(dt_s)),
        "dt_median_s": float(np.median(dt_s)),
        "dt_mean_s": float(np.mean(dt_s)),
        "dt_max_s": float(np.max(dt_s)),
        "raw_current_min_a": float(np.min(raw_current_a)),
        "raw_current_max_a": float(np.max(raw_current_a)),
        "voltage_min_v": float(np.min(measured_v)),
        "voltage_max_v": float(np.max(measured_v)),
        "exported_ah_min": float(np.min(exported_ah)),
        "exported_ah_max": float(np.max(exported_ah)),
        "exported_charge_ah": float(exported_ah[phases.last_charge_index] - exported_ah[0]),
        "battery_temperature_min_c": float(np.min(battery_temperature_c)),
        "battery_temperature_max_c": float(np.max(battery_temperature_c)),
        "chamber_temperature_min_c": float(np.min(chamber_temperature_c)),
        "chamber_temperature_max_c": float(np.max(chamber_temperature_c)),
        "preceding_discharge_capacity_ah": preceding_capacity_ah,
        "preceding_recovered_voltage_v": preceding_recovered_voltage_v,
        "initial_measured_voltage_v": float(measured_v[0]),
        "initial_soc": float(initial_soc),
        "final_soc": float(source_soc[-1]),
        "initial_x_n": initial_x_n,
        "initial_x_p": initial_x_p,
        "cc_current_min_a": float(np.min(raw_current_a[phases.cc])),
        "cc_current_mean_a": float(np.mean(raw_current_a[phases.cc])),
        "cc_current_max_a": float(np.max(raw_current_a[phases.cc])),
        "cv_start_voltage_v": float(measured_v[phases.cv_start_index]),
        "cv_voltage_min_v": float(np.min(measured_v[phases.cv])),
        "cv_voltage_max_v": float(np.max(measured_v[phases.cv])),
        "cv_current_start_a": float(raw_current_a[phases.cv_start_index]),
        "cv_current_end_a": float(raw_current_a[phases.last_charge_index]),
        "cv_terminal_residual_v": float(residual_v[phases.last_charge_index]),
        "cv_last_five_bias_v": float(np.mean(residual_v[last_five_cv])),
        "residual_vs_soc_spearman_rho": float(spearmanr(residual_v, source_soc).statistic),
        "x_n_avg_min": float(np.min(result.x_n_avg)),
        "x_n_avg_max": float(np.max(result.x_n_avg)),
        "x_p_avg_min": float(np.min(result.x_p_avg)),
        "x_p_avg_max": float(np.max(result.x_p_avg)),
        "x_n_surface_min": float(np.min(result.x_n_surface)),
        "x_n_surface_max": float(np.max(result.x_n_surface)),
        "x_p_surface_min": float(np.min(result.x_p_surface)),
        "x_p_surface_max": float(np.max(result.x_p_surface)),
        "ocp_domain_violation_count": 0,
        "clipping_count": 0,
    }

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    write_metrics(
        OUTPUT_DIRECTORY / "panasonic_cccv_charge_validation_metrics.csv",
        overall,
        phase_metrics,
        phases,
        regional,
        region_masks,
        audit,
    )
    save_voltage_plot(
        OUTPUT_DIRECTORY / "panasonic_cccv_charge_validation_voltage.png",
        time_s,
        measured_v,
        predicted_v,
        raw_current_a,
        phases,
    )
    save_residual_soc_plot(
        OUTPUT_DIRECTORY / "panasonic_cccv_charge_validation_residual_soc.png",
        source_soc,
        residual_v,
        phases,
    )
    save_cv_taper_plot(
        OUTPUT_DIRECTORY / "panasonic_cccv_charge_validation_cv_taper.png",
        time_s,
        raw_current_a,
        residual_v,
        phases,
    )
    write_report(
        OUTPUT_DIRECTORY / "panasonic_cccv_charge_validation_report.md",
        raw_columns,
        converted_columns,
        phases,
        overall,
        phase_metrics,
        regional,
        region_masks,
        time_s,
        raw_current_a,
        exported_ah,
        residual_v,
        audit,
    )

    print("Panasonic 25 C CC-CV charge forward-only validation complete.")
    print(
        f"Samples={time_s.size}; RMSE={1000.0 * overall['rmse_v']:.6f} mV; "
        f"bias={1000.0 * overall['mean_signed_error_v']:.6f} mV; R2={overall['r2']:.9f}"
    )
    print(
        f"CC/CV/final-taper RMSE={1000.0 * phase_metrics['cc']['rmse_v']:.6f}/"
        f"{1000.0 * phase_metrics['cv']['rmse_v']:.6f}/"
        f"{1000.0 * phase_metrics['final_5pct_taper']['rmse_v']:.6f} mV"
    )
    print("OCP-domain violations=0; clipping=0")
    print("Classification: CC-CV CHARGE TRANSFER LIMITED")


if __name__ == "__main__":
    main()
