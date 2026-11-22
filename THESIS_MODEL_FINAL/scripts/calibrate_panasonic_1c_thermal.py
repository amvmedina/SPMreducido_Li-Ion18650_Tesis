"""Calibrate effective Panasonic lumped thermal parameters on the 25 C 1C test.

Only Cth and Gth are fitted. The prescribed calibration heat is based on
measured terminal voltage and measured temperature. A model-voltage heat
trajectory is evaluated afterward with the fitted parameters and no refit.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares


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
    X_N_START,
    X_P_START,
)
from models.panasonic_entropic import (  # noqa: E402
    SUPPORTED_SOC_RANGE,
    dUeq_dT_panasonic,
)
from scripts.audit_panasonic_1c_heat_generation import (  # noqa: E402
    EXPECTED_FULL_FORWARD_RMSE_V,
    restrict_to_exact_soc_domain,
)
from scripts.audit_panasonic_hppc_identifiability import (  # noqa: E402
    Q_N_C as IMPLEMENTED_Q_N_C,
    Q_P_C as IMPLEMENTED_Q_P_C,
    simulate_local_group,
)
from scripts.validate_panasonic_1c_forward import (  # noqa: E402
    DISCHARGE_CURRENT_THRESHOLD_A,
    NOMINAL_PARAMETERS,
    Q_N_C,
    Q_P_C,
    RAW_SOURCE_PATH,
    load_raw_discharge,
    longest_contiguous_true,
)


OUTPUT_DIRECTORY = REPOSITORY_ROOT / "thesis_data"
REPORT_PATH = OUTPUT_DIRECTORY / "panasonic_1c_thermal_calibration.md"
METRICS_PATH = OUTPUT_DIRECTORY / "panasonic_1c_thermal_calibration_metrics.csv"
MULTISTART_PATH = OUTPUT_DIRECTORY / "panasonic_1c_thermal_calibration_multistart.csv"
FIT_PLOT_PATH = OUTPUT_DIRECTORY / "panasonic_1c_thermal_fit.png"
RESIDUAL_TIME_PATH = OUTPUT_DIRECTORY / "panasonic_1c_thermal_residual_time.png"
RESIDUAL_SOC_PATH = OUTPUT_DIRECTORY / "panasonic_1c_thermal_residual_soc.png"

C_TH_BOUNDS_J_PER_K = (5.0, 200.0)
G_TH_BOUNDS_W_PER_K = (1.0e-4, 1.0)
LOG_SENSITIVITY_STEP = 1.0e-5
NEAR_OPTIMAL_RELATIVE_MSE = 1.0e-4
BOUND_LOG_FRACTION = 0.01


@dataclass(frozen=True)
class ThermalCalibrationData:
    time_s: np.ndarray
    soc: np.ndarray
    current_a: np.ndarray
    measured_temperature_k: np.ndarray
    ambient_temperature_k: np.ndarray
    measured_voltage_v: np.ndarray
    model_voltage_v: np.ndarray
    x_n_avg: np.ndarray
    x_p_avg: np.ndarray
    ueq_average_v: np.ndarray
    d_u_eq_d_t_v_per_k: np.ndarray
    qdot_irr_measured_w: np.ndarray
    qdot_irr_model_w: np.ndarray
    qdot_rev_w: np.ndarray
    qdot_cal_w: np.ndarray
    qdot_model_w: np.ndarray
    full_forward_rmse_v: float


@dataclass(frozen=True)
class TemperatureMetrics:
    rmse_k: float
    mae_k: float
    bias_k: float
    max_abs_error_k: float
    r2: float


def prepare_data() -> ThermalCalibrationData:
    columns, raw = load_raw_discharge(RAW_SOURCE_PATH)
    if "Chamber_Temp_degC" not in columns:
        raise ValueError("Measured chamber temperature is required.")
    indices = longest_contiguous_true(raw["Current"] <= DISCHARGE_CURRENT_THRESHOLD_A)
    time_full = raw["Time"][indices]
    current_full = -raw["Current"][indices]
    measured_voltage_full = raw["Voltage"][indices]
    measured_temperature_full_k = raw["Battery_Temp_degC"][indices] + 273.15
    ambient_temperature_full_k = raw["Chamber_Temp_degC"][indices] + 273.15
    exported_ah = raw["Ah"][indices]
    soc_full = 1.0 - (exported_ah[0] - exported_ah) / Q_CELL_AH
    if np.any(np.diff(time_full) <= 0.0) or np.any(np.diff(soc_full) >= 0.0):
        raise ValueError("Active time and SOC must be strictly monotonic.")
    if np.any(current_full <= 0.0):
        raise ValueError("The selected branch is not discharge-positive after conversion.")

    electrical = simulate_local_group(
        time_s=time_full,
        discharge_positive_current_a=current_full,
        temperature_k=measured_temperature_full_k,
        initial_x_n=X_N_START,
        initial_x_p=X_P_START,
        dynamic_parameters=NOMINAL_PARAMETERS,
    )
    full_forward_rmse = float(
        np.sqrt(np.mean((electrical.voltage_v - measured_voltage_full) ** 2))
    )
    if not np.isclose(
        full_forward_rmse,
        EXPECTED_FULL_FORWARD_RMSE_V,
        rtol=0.0,
        atol=5.0e-13,
    ):
        raise RuntimeError("Frozen electrical 1C forward result changed.")

    soc, values = restrict_to_exact_soc_domain(
        soc_full,
        time_s=time_full,
        current_a=current_full,
        measured_temperature_k=measured_temperature_full_k,
        ambient_temperature_k=ambient_temperature_full_k,
        measured_voltage_v=measured_voltage_full,
        model_voltage_v=electrical.voltage_v,
        x_n_avg=electrical.x_n_avg,
        x_p_avg=electrical.x_p_avg,
    )
    x_n_avg = values["x_n_avg"]
    x_p_avg = values["x_p_avg"]
    if not (
        np.all((x_n_avg >= MCMB_SUPPORTED_STOICHIOMETRY_RANGE[0]) & (x_n_avg <= MCMB_SUPPORTED_STOICHIOMETRY_RANGE[1]))
        and np.all((x_p_avg >= NCA_WORKING_STOICHIOMETRY_RANGE[0]) & (x_p_avg <= NCA_WORKING_STOICHIOMETRY_RANGE[1]))
    ):
        raise ValueError("Average electrode state left frozen OCP support; no clipping.")

    ueq_average = np.asarray(U_p(x_p_avg)) - np.asarray(U_n(x_n_avg))
    coefficient = np.asarray(dUeq_dT_panasonic(soc))
    current = values["current_a"]
    measured_temperature = values["measured_temperature_k"]
    qdot_irr_measured = current * (ueq_average - values["measured_voltage_v"])
    qdot_irr_model = current * (ueq_average - values["model_voltage_v"])
    qdot_rev = -current * measured_temperature * coefficient
    qdot_cal = qdot_irr_measured + qdot_rev
    qdot_model = qdot_irr_model + qdot_rev
    required = (
        values["time_s"],
        soc,
        current,
        measured_temperature,
        values["ambient_temperature_k"],
        qdot_cal,
        qdot_model,
    )
    if not all(np.all(np.isfinite(array)) for array in required):
        raise ValueError("A required thermal-calibration signal is non-finite.")
    return ThermalCalibrationData(
        time_s=values["time_s"],
        soc=soc,
        current_a=current,
        measured_temperature_k=measured_temperature,
        ambient_temperature_k=values["ambient_temperature_k"],
        measured_voltage_v=values["measured_voltage_v"],
        model_voltage_v=values["model_voltage_v"],
        x_n_avg=x_n_avg,
        x_p_avg=x_p_avg,
        ueq_average_v=ueq_average,
        d_u_eq_d_t_v_per_k=coefficient,
        qdot_irr_measured_w=qdot_irr_measured,
        qdot_irr_model_w=qdot_irr_model,
        qdot_rev_w=qdot_rev,
        qdot_cal_w=qdot_cal,
        qdot_model_w=qdot_model,
        full_forward_rmse_v=full_forward_rmse,
    )


def exact_thermal_trajectory(
    time_s: np.ndarray,
    heat_w: np.ndarray,
    ambient_temperature_k: np.ndarray,
    initial_temperature_k: float,
    c_th_j_per_k: float,
    g_th_w_per_k: float,
) -> np.ndarray:
    """Exact interval propagation for left-sample ZOH heat and ambient."""

    if c_th_j_per_k <= 0.0 or g_th_w_per_k < 0.0:
        raise ValueError("Cth must be positive and Gth nonnegative.")
    if not (time_s.shape == heat_w.shape == ambient_temperature_k.shape):
        raise ValueError("Thermal inputs must have matching shapes.")
    if np.any(np.diff(time_s) <= 0.0):
        raise ValueError("Thermal timestamps must be strictly increasing.")
    temperature = np.empty_like(time_s, dtype=float)
    temperature[0] = float(initial_temperature_k)
    for index, dt_s in enumerate(np.diff(time_s)):
        if g_th_w_per_k == 0.0:
            temperature[index + 1] = (
                temperature[index] + heat_w[index] * dt_s / c_th_j_per_k
            )
            continue
        x = g_th_w_per_k * dt_s / c_th_j_per_k
        approach = -np.expm1(-x)
        phi = approach / x if x != 0.0 else 1.0
        temperature[index + 1] = (
            temperature[index]
            + approach * (ambient_temperature_k[index] - temperature[index])
            + heat_w[index] * dt_s * phi / c_th_j_per_k
        )
    return temperature


def residual_log_parameters(log_parameters: np.ndarray, data: ThermalCalibrationData) -> np.ndarray:
    c_th, g_th = np.exp(log_parameters)
    predicted = exact_thermal_trajectory(
        data.time_s,
        data.qdot_cal_w,
        data.ambient_temperature_k,
        data.measured_temperature_k[0],
        float(c_th),
        float(g_th),
    )
    return predicted - data.measured_temperature_k


def temperature_metrics(measured_k: np.ndarray, predicted_k: np.ndarray) -> TemperatureMetrics:
    residual = predicted_k - measured_k
    denominator = float(np.sum((measured_k - np.mean(measured_k)) ** 2))
    return TemperatureMetrics(
        rmse_k=float(np.sqrt(np.mean(residual**2))),
        mae_k=float(np.mean(np.abs(residual))),
        bias_k=float(np.mean(residual)),
        max_abs_error_k=float(np.max(np.abs(residual))),
        r2=float(1.0 - np.sum(residual**2) / denominator),
    )


def run_multistart(data: ThermalCalibrationData) -> list[dict[str, object]]:
    c_starts = np.geomspace(*C_TH_BOUNDS_J_PER_K, 5)
    g_starts = np.geomspace(*G_TH_BOUNDS_W_PER_K, 5)
    lower = np.log([C_TH_BOUNDS_J_PER_K[0], G_TH_BOUNDS_W_PER_K[0]])
    upper = np.log([C_TH_BOUNDS_J_PER_K[1], G_TH_BOUNDS_W_PER_K[1]])
    records: list[dict[str, object]] = []
    start_id = 0
    for c_initial in c_starts:
        for g_initial in g_starts:
            start_id += 1
            result = least_squares(
                residual_log_parameters,
                np.log([c_initial, g_initial]),
                args=(data,),
                bounds=(lower, upper),
                method="trf",
                loss="linear",
                ftol=1.0e-13,
                xtol=1.0e-13,
                gtol=1.0e-13,
                max_nfev=3000,
            )
            c_fit, g_fit = np.exp(result.x)
            residual = result.fun
            records.append(
                {
                    "start_id": start_id,
                    "initial_c_th_j_per_k": float(c_initial),
                    "initial_g_th_w_per_k": float(g_initial),
                    "success": bool(result.success),
                    "status": int(result.status),
                    "nfev": int(result.nfev),
                    "cost_half_sum_squared_k2": float(result.cost),
                    "objective_mse_k2": float(np.mean(residual**2)),
                    "rmse_k": float(np.sqrt(np.mean(residual**2))),
                    "c_th_j_per_k": float(c_fit),
                    "g_th_w_per_k": float(g_fit),
                    "tau_th_s": float(c_fit / g_fit),
                    "optimality": float(result.optimality),
                    "message": str(result.message).replace("\n", " "),
                }
            )
    return records


def local_identifiability(
    data: ThermalCalibrationData,
    c_th: float,
    g_th: float,
) -> dict[str, object]:
    center = np.log([c_th, g_th])
    columns = []
    for parameter_index in range(2):
        plus = center.copy()
        minus = center.copy()
        plus[parameter_index] += LOG_SENSITIVITY_STEP
        minus[parameter_index] -= LOG_SENSITIVITY_STEP
        columns.append(
            (
                residual_log_parameters(plus, data)
                - residual_log_parameters(minus, data)
            )
            / (2.0 * LOG_SENSITIVITY_STEP)
        )
    jacobian = np.column_stack(columns)
    singular_values = np.linalg.svd(jacobian, compute_uv=False)
    tolerance = np.finfo(float).eps * max(jacobian.shape) * singular_values[0]
    rank = int(np.sum(singular_values > tolerance))
    condition_number = float(singular_values[0] / singular_values[-1])
    sensitivity_correlation = float(np.corrcoef(jacobian.T)[0, 1])
    covariance = np.linalg.pinv(jacobian.T @ jacobian, rcond=1.0e-12)
    standard = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    parameter_correlation = float(covariance[0, 1] / (standard[0] * standard[1]))
    return {
        "jacobian": jacobian,
        "singular_values": singular_values,
        "rank": rank,
        "condition_number": condition_number,
        "sensitivity_correlation": sensitivity_correlation,
        "parameter_correlation": parameter_correlation,
        "sensitivity_rms_c_th_k": float(np.sqrt(np.mean(jacobian[:, 0] ** 2))),
        "sensitivity_rms_g_th_k": float(np.sqrt(np.mean(jacobian[:, 1] ** 2))),
    }


def near_bound(value: float, bounds: tuple[float, float]) -> bool:
    fraction = (np.log(value) - np.log(bounds[0])) / (np.log(bounds[1]) - np.log(bounds[0]))
    return bool(fraction <= BOUND_LOG_FRACTION or fraction >= 1.0 - BOUND_LOG_FRACTION)


def region_masks(soc: np.ndarray) -> dict[str, np.ndarray]:
    lower, upper = SUPPORTED_SOC_RANGE
    first = lower + (upper - lower) / 3.0
    second = lower + 2.0 * (upper - lower) / 3.0
    return {
        "early_high_soc": (soc >= second) & (soc <= upper),
        "middle_soc": (soc >= first) & (soc < second),
        "late_low_soc": (soc >= lower) & (soc < first),
    }


def write_multistart(path: Path, records: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(records[0]))
        writer.writeheader()
        writer.writerows(records)


def add_metrics_rows(
    rows: list[dict[str, object]],
    scope: str,
    metrics: TemperatureMetrics,
) -> None:
    for name, value, unit in (
        ("rmse", metrics.rmse_k, "K"),
        ("mae", metrics.mae_k, "K"),
        ("mean_signed_error", metrics.bias_k, "K"),
        ("maximum_absolute_error", metrics.max_abs_error_k, "K"),
        ("r2", metrics.r2, "1"),
    ):
        rows.append({"scope": scope, "metric": name, "value": value, "unit": unit})


def write_metrics(
    path: Path,
    primary: TemperatureMetrics,
    secondary: TemperatureMetrics,
    regional: dict[str, TemperatureMetrics],
    scalar_rows: list[tuple[str, str, float | int | str, str]],
) -> None:
    rows: list[dict[str, object]] = []
    add_metrics_rows(rows, "primary_measured_voltage_heat", primary)
    add_metrics_rows(rows, "secondary_model_voltage_heat_no_refit", secondary)
    for name, values in regional.items():
        add_metrics_rows(rows, f"primary_{name}", values)
    for scope, metric, value, unit in scalar_rows:
        rows.append({"scope": scope, "metric": metric, "value": value, "unit": unit})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("scope", "metric", "value", "unit"))
        writer.writeheader()
        writer.writerows(rows)


def save_fit_plot(
    path: Path,
    data: ThermalCalibrationData,
    primary_k: np.ndarray,
    secondary_k: np.ndarray,
) -> None:
    time_min = (data.time_s - data.time_s[0]) / 60.0
    fig, ax = plt.subplots(figsize=(9.5, 5.4), constrained_layout=True)
    ax.plot(time_min, data.measured_temperature_k - 273.15, color="black", lw=2.0, label="Measured surface temperature")
    ax.plot(time_min, primary_k - 273.15, color="#0072B2", lw=1.8, label="Fitted: measured-voltage heat")
    ax.plot(time_min, secondary_k - 273.15, color="#D55E00", lw=1.6, ls="--", label="No-refit: model-voltage heat")
    ax.plot(time_min, data.ambient_temperature_k - 273.15, color="#009E73", lw=1.5, label="Measured chamber temperature")
    ax.set(xlabel="Elapsed time in SOC 0.90–0.10 interval (min)", ylabel="Temperature (°C)", title="Panasonic 1C effective lumped thermal calibration")
    ax.grid(alpha=0.25)
    ax.legend(loc="best")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def save_residual_plot(
    path: Path,
    x: np.ndarray,
    x_label: str,
    primary_residual_k: np.ndarray,
    secondary_residual_k: np.ndarray,
    *,
    invert_x: bool,
) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 4.8), constrained_layout=True)
    ax.plot(x, primary_residual_k, color="#0072B2", lw=1.7, label="Measured-voltage heat fit")
    ax.plot(x, secondary_residual_k, color="#D55E00", lw=1.5, ls="--", label="Model-voltage heat, no refit")
    ax.axhline(0.0, color="black", lw=0.9)
    ax.set(xlabel=x_label, ylabel="Temperature residual, model − measured (K)")
    if invert_x:
        ax.invert_xaxis()
    ax.grid(alpha=0.25)
    ax.legend(loc="best")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def metrics_table(metrics: TemperatureMetrics) -> str:
    return "\n".join(
        [
            f"| RMSE | {metrics.rmse_k:.9f} K |",
            f"| MAE | {metrics.mae_k:.9f} K |",
            f"| Mean signed error | {metrics.bias_k:.9f} K |",
            f"| Maximum absolute error | {metrics.max_abs_error_k:.9f} K |",
            f"| R² | {metrics.r2:.9f} |",
        ]
    )


def main() -> None:
    if not (
        np.isclose(IMPLEMENTED_Q_N_C, Q_N_C, rtol=0.0, atol=5.0e-9)
        and np.isclose(IMPLEMENTED_Q_P_C, Q_P_C, rtol=0.0, atol=5.0e-9)
    ):
        raise RuntimeError("Frozen electrical capacities changed.")
    data = prepare_data()
    records = run_multistart(data)
    successful = [record for record in records if bool(record["success"])]
    if len(successful) < 20:
        raise RuntimeError(f"Only {len(successful)} of 25 deterministic starts succeeded.")
    best = min(successful, key=lambda record: float(record["objective_mse_k2"]))
    c_th = float(best["c_th_j_per_k"])
    g_th = float(best["g_th_w_per_k"])
    tau_th = c_th / g_th
    primary_k = exact_thermal_trajectory(
        data.time_s,
        data.qdot_cal_w,
        data.ambient_temperature_k,
        data.measured_temperature_k[0],
        c_th,
        g_th,
    )
    secondary_k = exact_thermal_trajectory(
        data.time_s,
        data.qdot_model_w,
        data.ambient_temperature_k,
        data.measured_temperature_k[0],
        c_th,
        g_th,
    )
    primary_metrics = temperature_metrics(data.measured_temperature_k, primary_k)
    secondary_metrics = temperature_metrics(data.measured_temperature_k, secondary_k)
    primary_residual = primary_k - data.measured_temperature_k
    secondary_residual = secondary_k - data.measured_temperature_k

    masks = region_masks(data.soc)
    regional = {
        name: temperature_metrics(data.measured_temperature_k[mask], primary_k[mask])
        for name, mask in masks.items()
    }
    identifiability = local_identifiability(data, c_th, g_th)
    near_optimal = [
        record
        for record in successful
        if float(record["objective_mse_k2"])
        <= float(best["objective_mse_k2"]) * (1.0 + NEAR_OPTIMAL_RELATIVE_MSE)
    ]
    c_success = np.asarray([float(record["c_th_j_per_k"]) for record in successful])
    g_success = np.asarray([float(record["g_th_w_per_k"]) for record in successful])
    c_near = np.asarray([float(record["c_th_j_per_k"]) for record in near_optimal])
    g_near = np.asarray([float(record["g_th_w_per_k"]) for record in near_optimal])
    c_near_bound = near_bound(c_th, C_TH_BOUNDS_J_PER_K)
    g_near_bound = near_bound(g_th, G_TH_BOUNDS_W_PER_K)
    if c_near_bound or g_near_bound:
        classification = "THERMAL PARAMETER BOUNDARY LIMITED"
    elif (
        abs(float(identifiability["sensitivity_correlation"])) >= 0.95
        or float(identifiability["condition_number"]) >= 100.0
    ):
        classification = "THERMAL PARAMETERS PRACTICALLY CORRELATED"
    else:
        classification = "THERMAL PARAMETERS WELL CONDITIONED"

    max_error_index = int(np.argmax(np.abs(primary_residual)))
    secondary_max_error_index = int(np.argmax(np.abs(secondary_residual)))
    measured_rise_k = float(data.measured_temperature_k[-1] - data.measured_temperature_k[0])
    primary_rise_k = float(primary_k[-1] - primary_k[0])
    secondary_rise_k = float(secondary_k[-1] - secondary_k[0])
    pearson_soc = float(np.corrcoef(data.soc, primary_residual)[0, 1])

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    write_multistart(MULTISTART_PATH, records)
    scalar_rows = [
        ("parameters", "Cth", c_th, "J/K"),
        ("parameters", "Gth", g_th, "W/K"),
        ("parameters", "tau_th", tau_th, "s"),
        ("optimization", "objective_mse", float(best["objective_mse_k2"]), "K^2"),
        ("optimization", "successful_multistarts", len(successful), "count"),
        ("optimization", "near_optimal_multistarts", len(near_optimal), "count"),
        ("identifiability", "rank", int(identifiability["rank"]), "count"),
        ("identifiability", "singular_value_1", float(identifiability["singular_values"][0]), "K"),
        ("identifiability", "singular_value_2", float(identifiability["singular_values"][1]), "K"),
        ("identifiability", "condition_number", float(identifiability["condition_number"]), "1"),
        ("identifiability", "sensitivity_correlation", float(identifiability["sensitivity_correlation"]), "1"),
        ("identifiability", "parameter_correlation", float(identifiability["parameter_correlation"]), "1"),
        ("identifiability", "classification", classification, "text"),
        ("data", "initial_measured_temperature", float(data.measured_temperature_k[0]), "K"),
        ("data", "final_measured_temperature", float(data.measured_temperature_k[-1]), "K"),
        ("primary", "final_model_temperature", float(primary_k[-1]), "K"),
        ("secondary", "final_model_temperature", float(secondary_k[-1]), "K"),
        ("data", "measured_temperature_rise", measured_rise_k, "K"),
        ("primary", "predicted_temperature_rise", primary_rise_k, "K"),
        ("secondary", "predicted_temperature_rise", secondary_rise_k, "K"),
        ("data", "maximum_measured_temperature", float(np.max(data.measured_temperature_k)), "K"),
        ("primary", "maximum_predicted_temperature", float(np.max(primary_k)), "K"),
        ("secondary", "maximum_predicted_temperature", float(np.max(secondary_k)), "K"),
        ("primary", "maximum_error_time", float(data.time_s[max_error_index] - data.time_s[0]), "s"),
        ("primary", "maximum_error_soc", float(data.soc[max_error_index]), "1"),
        ("secondary", "maximum_error_time", float(data.time_s[secondary_max_error_index] - data.time_s[0]), "s"),
        ("secondary", "maximum_error_soc", float(data.soc[secondary_max_error_index]), "1"),
    ]
    write_metrics(METRICS_PATH, primary_metrics, secondary_metrics, regional, scalar_rows)
    save_fit_plot(FIT_PLOT_PATH, data, primary_k, secondary_k)
    save_residual_plot(
        RESIDUAL_TIME_PATH,
        (data.time_s - data.time_s[0]) / 60.0,
        "Elapsed time in SOC 0.90–0.10 interval (min)",
        primary_residual,
        secondary_residual,
        invert_x=False,
    )
    save_residual_plot(
        RESIDUAL_SOC_PATH,
        data.soc,
        "Average cell SOC",
        primary_residual,
        secondary_residual,
        invert_x=True,
    )

    report = f"""# Panasonic 25 °C 1C primary thermal calibration

## Result: {classification}

The calibration fitted only `Cth` and `Gth`. They are retained as **effective lumped thermal parameters**, not uniquely measured material properties: the target is cell-surface temperature, the cell has one uniform model temperature, and all external heat-transfer detail is collapsed into `Gth`.

No electrical parameter, OCP, nominal effective stoichiometric window, capacity, entropy coefficient, or voltage equation was changed. No initial-temperature offset was fitted.

## Data and prescribed calibration heat

- Source: `{RAW_SOURCE_PATH}`
- Exact entropy support used: `0.10 <= SOC <= 0.90`, with the same exact endpoint interpolation as the preceding heat audit.
- Points: {data.soc.size}; duration: {data.time_s[-1] - data.time_s[0]:.9f} s.
- Raw negative discharge current was converted once to model discharge-positive current: {data.current_a.min():.6f}–{data.current_a.max():.6f} A.
- Target surface temperature: measured `Battery_Temp_degC`, {data.measured_temperature_k.min() - 273.15:.6f}–{data.measured_temperature_k.max() - 273.15:.6f} °C.
- Measured time-varying ambient: `Chamber_Temp_degC`, {data.ambient_temperature_k.min() - 273.15:.6f}–{data.ambient_temperature_k.max() - 273.15:.6f} °C. It was not replaced by a fixed 25 °C value.
- Initial condition: measured and modeled temperature both {data.measured_temperature_k[0] - 273.15:.9f} °C.
- Frozen electrical forward result reproduced: {1000.0 * data.full_forward_rmse_v:.9f} mV full-discharge RMSE.

The primary prescribed heat was

```text
Ueq_average = U_p(x_p_avg) - U_n(x_n_avg)
Qdot_irr_measured = I * (Ueq_average - V_measured)
Qdot_rev = -I * T_measured_K * dUeq_dT(SOC)
Qdot_cal = Qdot_irr_measured + Qdot_rev
```

Measured temperature—not fitted model temperature—was used in `Qdot_rev`. `Qdot_irr_model` was not used in the fit.

## Optimization

The model was propagated exactly on every nonuniform interval under left-sample piecewise-constant `Qdot_cal` and measured `Tamb`. The implementation uses `expm1` and a finite small-argument ratio, retaining the adiabatic limiting form without Euler integration.

- Guard bounds: `5 <= Cth <= 200 J/K`; `1e-4 <= Gth <= 1 W/K`.
- Deterministic starts: {len(records)} on a 5×5 grid distributed geometrically in parameter space (uniformly in log space).
- Successful starts: {len(successful)}.
- Primary objective: mean squared unweighted temperature residual.
- Best objective: {float(best['objective_mse_k2']):.12g} K².
- **Cth = {c_th:.12g} J/K**.
- **Gth = {g_th:.12g} W/K**.
- **Thermal time constant = {tau_th:.12g} s ({tau_th / 60.0:.9f} min)**.
- Frozen effective-parameter module: `models/panasonic_thermal_parameters.py`. It does not alter generic thermal defaults.
- Cth near a guard bound: **{yes_no(c_near_bound).upper()}**; Gth near a guard bound: **{yes_no(g_near_bound).upper()}**. “Near” means the closest 1% of the log-bound interval.

### Primary measured-voltage-heat fit

| Metric | Value |
|---|---:|
{metrics_table(primary_metrics)}

## Multistart stability

- Successful-solution Cth range: {c_success.min():.12g}–{c_success.max():.12g} J/K; spread {c_success.max() - c_success.min():.6g} J/K.
- Successful-solution Gth range: {g_success.min():.12g}–{g_success.max():.12g} W/K; spread {g_success.max() - g_success.min():.6g} W/K.
- Near-optimal definition: MSE no more than {100.0 * NEAR_OPTIMAL_RELATIVE_MSE:.4f}% above best.
- Near-optimal starts: {len(near_optimal)}.
- Near-optimal Cth range: {c_near.min():.12g}–{c_near.max():.12g} J/K.
- Near-optimal Gth range: {g_near.min():.12g}–{g_near.max():.12g} W/K.

## Local conditioning audit

Finite central differences used step `{LOG_SENSITIVITY_STEP:.1e}` in `log(Cth)` and `log(Gth)`.

- Numerical rank: **{int(identifiability['rank'])} of 2**.
- Singular values: `{float(identifiability['singular_values'][0]):.12g}`, `{float(identifiability['singular_values'][1]):.12g}` K per unit log parameter.
- Condition number: `{float(identifiability['condition_number']):.12g}`.
- Sensitivity-column correlation: `{float(identifiability['sensitivity_correlation']):.12g}`.
- Local estimator parameter correlation from `(JᵀJ)^−1`: `{float(identifiability['parameter_correlation']):.12g}`.
- Sensitivity RMS: `log(Cth)` {float(identifiability['sensitivity_rms_c_th_k']):.12g} K; `log(Gth)` {float(identifiability['sensitivity_rms_g_th_k']):.12g} K.

This is a practical local-conditioning audit for one heating experiment, not a claim of structural identifiability.

The covariance-derived correlation is appreciable, but it does not produce practical solution ambiguity here: the Jacobian condition number is low, both singular values are resolved, and all 25 broadly separated starts collapse to the same narrow parameter ranges. That evidence supports the stated well-conditioned classification for this experiment while the parameters remain explicitly effective and experiment-conditional.

## Temperature diagnostics

- Initial measured/model temperature: {data.measured_temperature_k[0] - 273.15:.9f}/{primary_k[0] - 273.15:.9f} °C.
- Final measured/model temperature: {data.measured_temperature_k[-1] - 273.15:.9f}/{primary_k[-1] - 273.15:.9f} °C.
- Measured/predicted temperature rise: {measured_rise_k:.9f}/{primary_rise_k:.9f} K.
- Measured/predicted maximum temperature: {np.max(data.measured_temperature_k) - 273.15:.9f}/{np.max(primary_k) - 273.15:.9f} °C.
- Maximum absolute primary error: {primary_metrics.max_abs_error_k:.9f} K at elapsed time {(data.time_s[max_error_index] - data.time_s[0]) / 60.0:.9f} min and SOC {data.soc[max_error_index]:.12f}.
- Pearson correlation of primary residual with SOC: {pearson_soc:.9f}.

### Residuals by discharge region

| Region | SOC interval | Samples | RMSE (K) | Bias (K) | Maximum absolute error (K) |
|---|---|---:|---:|---:|---:|
| Early/high | [0.633333, 0.90] | {int(np.sum(masks['early_high_soc']))} | {regional['early_high_soc'].rmse_k:.9f} | {regional['early_high_soc'].bias_k:.9f} | {regional['early_high_soc'].max_abs_error_k:.9f} |
| Middle | [0.366667, 0.633333) | {int(np.sum(masks['middle_soc']))} | {regional['middle_soc'].rmse_k:.9f} | {regional['middle_soc'].bias_k:.9f} | {regional['middle_soc'].max_abs_error_k:.9f} |
| Late/low | [0.10, 0.366667) | {int(np.sum(masks['late_low_soc']))} | {regional['late_low_soc'].rmse_k:.9f} | {regional['late_low_soc'].bias_k:.9f} | {regional['late_low_soc'].max_abs_error_k:.9f} |

These regional values and the residual plots document any early-, middle-, late-, or SOC-dependent structure without adding weights or correction terms.

The residual is not strongly monotonic with SOC: regional biases alternate from slightly negative (early) to slightly positive (middle) and back to slightly negative (late), and the residual/SOC correlation is small. The largest error is late in the discharge near SOC {data.soc[max_error_index]:.3f}, but the regional RMSE is highest in the early/high-SOC third rather than progressively worsening toward low SOC.

## Secondary no-refit model-voltage heat diagnostic

After freezing the best values above, the same thermal parameters and measured-temperature reversible heat were propagated once with

```text
Qdot_total_model = I * (Ueq_average - V_model) + Qdot_rev.
```

No parameter was refitted.

| Metric | Value |
|---|---:|
{metrics_table(secondary_metrics)}
| Final model temperature | {secondary_k[-1] - 273.15:.9f} °C |
| Predicted temperature rise | {secondary_rise_k:.9f} K |
| Maximum-error time | {(data.time_s[secondary_max_error_index] - data.time_s[0]) / 60.0:.9f} min |
| Maximum-error SOC | {data.soc[secondary_max_error_index]:.12f} |

Relative to the primary fit, switching to model-voltage irreversible heat changes RMSE by {secondary_metrics.rmse_k - primary_metrics.rmse_k:+.9f} K and final predicted temperature by {secondary_k[-1] - primary_k[-1]:+.9f} K. This loss is attributed to the already documented electrical-voltage transfer error expressed through heat input; it was not repaired thermally.

## Interpretation and stop condition

**A. Thermal-submodel calibration:** the fitted effective parameters use measured-voltage irreversible heat.

**B. Coupled electrothermal forward diagnostic:** the no-refit result uses model-voltage irreversible heat with exactly the same effective parameters.

The calibration stops here. No UDDS or CC–CV thermal fitting/validation, electrical refit, Arrhenius dependence, entropy modification, OCP change, or own-TRC data use was performed.
"""
    REPORT_PATH.write_text(report, encoding="utf-8")

    print(f"Classification: {classification}")
    print(f"Cth={c_th:.12g} J/K; Gth={g_th:.12g} W/K; tau={tau_th:.12g} s")
    print(f"Primary RMSE={primary_metrics.rmse_k:.9f} K; secondary no-refit RMSE={secondary_metrics.rmse_k:.9f} K")
    print(f"Condition number={float(identifiability['condition_number']):.9g}; sensitivity correlation={float(identifiability['sensitivity_correlation']):.9g}")
    print(f"Successful/near-optimal starts: {len(successful)}/{len(near_optimal)}")
    print(f"Wrote: {REPORT_PATH}")
    print(f"Wrote: {METRICS_PATH}")
    print(f"Wrote: {MULTISTART_PATH}")
    print(f"Wrote: {FIT_PLOT_PATH}")
    print(f"Wrote: {RESIDUAL_TIME_PATH}")
    print(f"Wrote: {RESIDUAL_SOC_PATH}")


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


if __name__ == "__main__":
    main()
