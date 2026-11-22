"""Calibrate and validate the clean Panasonic two-node thermal model.

Four positive legacy-equivalent physical parameters are fitted only on the
closed 25 C 1C entropy-supported interval. The independent UDDS evaluation is
strictly forward-only. Heat is measured-voltage irreversible plus verified
reversible heat; no empirical heat or temperature correction exists here.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from models.panasonic_thermal_parameters import C_TH_J_PER_K, G_TH_W_PER_K  # noqa: E402
from models.two_node_thermal_physical import (  # noqa: E402
    TwoNodeThermalParameters,
    simulate_two_node_physical,
    thermal_eigenvalues,
    thermal_time_constants_s,
)
from scripts.calibrate_panasonic_1c_thermal import (  # noqa: E402
    RAW_SOURCE_PATH as ONE_C_RAW_SOURCE_PATH,
    exact_thermal_trajectory as one_node_trajectory,
    prepare_data as prepare_one_c_data,
)
from scripts.validate_panasonic_udds_thermal import (  # noqa: E402
    RAW_SOURCE_PATH as UDDS_RAW_SOURCE_PATH,
    prepare_data as prepare_udds_data,
)


OUTPUT_DIRECTORY = REPOSITORY_ROOT / "thesis_data"
REPORT_PATH = OUTPUT_DIRECTORY / "panasonic_two_node_physical_calibration.md"
METRICS_PATH = OUTPUT_DIRECTORY / "panasonic_two_node_physical_metrics.csv"
MULTISTART_PATH = OUTPUT_DIRECTORY / "panasonic_two_node_physical_multistart.csv"
FIT_PLOT_PATH = OUTPUT_DIRECTORY / "panasonic_two_node_1c_fit.png"
UDDS_PLOT_PATH = OUTPUT_DIRECTORY / "panasonic_two_node_udds_validation.png"
LEGACY_AUDIT_PATH = OUTPUT_DIRECTORY / "two_node_thermal_legacy_audit.md"

PARAMETER_NAMES = ("c_c_j_per_k", "c_s_j_per_k", "g_cs_w_per_k", "hA_w_per_k")
LOWER_BOUNDS = np.asarray([1.0e-3, 1.0e-4, 1.0e-5, 1.0e-5])
UPPER_BOUNDS = np.asarray([5000.0, 2000.0, 100.0, 10.0])
N_MULTISTARTS = 40
RANDOM_SEED = 20260913
LOG_SENSITIVITY_STEP = 1.0e-5
NEAR_OPTIMAL_RELATIVE_RMSE = 1.0e-3

ONE_NODE_REFERENCE = {
    "rmse_k": 0.143622065,
    "mae_k": 0.119569322,
    "mean_signed_error_k": -0.008530653,
    "r2": 0.964702003,
}
ONE_NODE_UDDS_REFERENCE = {
    "rmse_k": 0.526697652,
    "mae_k": 0.512960900,
    "mean_signed_error_k": -0.512714512,
    "max_abs_error_k": 0.875239690,
    "r2": -4.119410369,
    "measured_final_temperature_c": 26.469807973,
    "predicted_final_temperature_c": 26.112601243,
    "measured_temperature_rise_k": 0.636968211,
    "predicted_temperature_rise_k": 0.279761481,
}


def vector_to_parameters(values: np.ndarray) -> TwoNodeThermalParameters:
    return TwoNodeThermalParameters(**dict(zip(PARAMETER_NAMES, map(float, values), strict=True)))


def simulate_surface(
    log_parameters: np.ndarray,
    time_s: np.ndarray,
    heat_w: np.ndarray,
    ambient_k: np.ndarray,
    measured_surface_k: np.ndarray,
) -> np.ndarray:
    parameters = vector_to_parameters(np.exp(log_parameters))
    result = simulate_two_node_physical(
        time_s=time_s,
        q_gen_w=heat_w,
        t_ambient_k=ambient_k,
        parameters=parameters,
        initial_core_temperature_k=float(measured_surface_k[0]),
        initial_surface_temperature_k=float(measured_surface_k[0]),
    )
    return result.t_surface_k


def residual_log_parameters(
    log_parameters: np.ndarray,
    time_s: np.ndarray,
    heat_w: np.ndarray,
    ambient_k: np.ndarray,
    measured_surface_k: np.ndarray,
) -> np.ndarray:
    return simulate_surface(log_parameters, time_s, heat_w, ambient_k, measured_surface_k) - measured_surface_k


def metrics(measured: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    residual = predicted - measured
    denominator = float(np.sum((measured - np.mean(measured)) ** 2))
    return {
        "rmse_k": float(np.sqrt(np.mean(residual**2))),
        "mae_k": float(np.mean(np.abs(residual))),
        "mean_signed_error_k": float(np.mean(residual)),
        "max_abs_error_k": float(np.max(np.abs(residual))),
        "r2": float(1.0 - np.sum(residual**2) / denominator),
        "measured_final_temperature_c": float(measured[-1] - 273.15),
        "predicted_final_temperature_c": float(predicted[-1] - 273.15),
        "measured_temperature_rise_k": float(measured[-1] - measured[0]),
        "predicted_temperature_rise_k": float(predicted[-1] - predicted[0]),
    }


def deterministic_starts() -> np.ndarray:
    named = np.asarray(
        [
            [79.266263, 21.462628, 0.422606, 0.105224],
            [78.0436825, 23.59375, 0.466711436, 0.105648286],
            [90.0, 30.0, 0.40, 0.10],
            [60.0, 60.0, 0.20, 0.10],
            [110.0, 10.0, 1.00, 0.10],
            [50.0, 70.0, 0.05, 0.20],
        ]
    )
    rng = np.random.default_rng(RANDOM_SEED)
    random_log = rng.uniform(
        np.log(LOWER_BOUNDS), np.log(UPPER_BOUNDS), size=(N_MULTISTARTS - named.shape[0], 4)
    )
    return np.vstack((named, np.exp(random_log)))


def run_multistart(data: object) -> tuple[list[dict[str, object]], dict[str, object]]:
    rows: list[dict[str, object]] = []
    best: dict[str, object] | None = None
    arguments = (
        data.time_s,
        data.qdot_cal_w,
        data.ambient_temperature_k,
        data.measured_temperature_k,
    )
    for start_id, start in enumerate(deterministic_starts(), start=1):
        result = least_squares(
            residual_log_parameters,
            np.log(start),
            args=arguments,
            bounds=(np.log(LOWER_BOUNDS), np.log(UPPER_BOUNDS)),
            method="trf",
            jac="3-point",
            x_scale="jac",
            ftol=1.0e-12,
            xtol=1.0e-12,
            gtol=1.0e-12,
            max_nfev=3000,
        )
        fitted = np.exp(result.x)
        residual = np.asarray(result.fun)
        row: dict[str, object] = {
            "start_id": start_id,
            "success": bool(result.success),
            "status": int(result.status),
            "message": str(result.message).replace("\n", " "),
            "nfev": int(result.nfev),
            "cost": float(result.cost),
            "rmse_k": float(np.sqrt(np.mean(residual**2))),
            "mae_k": float(np.mean(np.abs(residual))),
            "optimality": float(result.optimality),
        }
        for index, name in enumerate(PARAMETER_NAMES):
            row[f"start_{name}"] = float(start[index])
            row[f"fit_{name}"] = float(fitted[index])
        rows.append(row)
        if result.success and (best is None or float(row["rmse_k"]) < float(best["row"]["rmse_k"])):
            best = {"row": row, "result": result}
    if best is None:
        raise RuntimeError("No two-node multistart optimization succeeded.")
    return rows, best


def numerical_sensitivity(
    best_log: np.ndarray, data: object
) -> tuple[np.ndarray, np.ndarray, int, float, np.ndarray, np.ndarray]:
    arguments = (
        data.time_s,
        data.qdot_cal_w,
        data.ambient_temperature_k,
        data.measured_temperature_k,
    )
    jacobian = np.empty((data.time_s.size, 4))
    for column in range(4):
        plus = best_log.copy()
        minus = best_log.copy()
        plus[column] += LOG_SENSITIVITY_STEP
        minus[column] -= LOG_SENSITIVITY_STEP
        jacobian[:, column] = (
            residual_log_parameters(plus, *arguments)
            - residual_log_parameters(minus, *arguments)
        ) / (2.0 * LOG_SENSITIVITY_STEP)
    singular_values = np.linalg.svd(jacobian, full_matrices=False, compute_uv=False)
    threshold = np.finfo(float).eps * max(jacobian.shape) * singular_values[0]
    rank = int(np.count_nonzero(singular_values > threshold))
    condition = float(singular_values[0] / singular_values[-1])
    sensitivity_correlation = np.corrcoef(jacobian, rowvar=False)
    covariance_shape = np.linalg.pinv(jacobian.T @ jacobian)
    scale = np.sqrt(np.diag(covariance_shape))
    parameter_correlation = covariance_shape / np.outer(scale, scale)
    return (
        jacobian,
        singular_values,
        rank,
        condition,
        sensitivity_correlation,
        parameter_correlation,
    )


def multistart_summary(rows: list[dict[str, object]]) -> dict[str, object]:
    successful = [row for row in rows if bool(row["success"])]
    best_rmse = min(float(row["rmse_k"]) for row in successful)
    tolerance = max(1.0e-8, NEAR_OPTIMAL_RELATIVE_RMSE * best_rmse)
    near = [row for row in successful if float(row["rmse_k"]) <= best_rmse + tolerance]
    summary: dict[str, object] = {
        "successful_count": len(successful),
        "near_optimal_count": len(near),
        "best_rmse_k": best_rmse,
        "near_optimal_rmse_tolerance_k": tolerance,
    }
    for name in PARAMETER_NAMES:
        all_values = np.asarray([float(row[f"fit_{name}"]) for row in successful])
        near_values = np.asarray([float(row[f"fit_{name}"]) for row in near])
        summary[f"successful_{name}_min"] = float(np.min(all_values))
        summary[f"successful_{name}_max"] = float(np.max(all_values))
        summary[f"near_optimal_{name}_min"] = float(np.min(near_values))
        summary[f"near_optimal_{name}_max"] = float(np.max(near_values))
        summary[f"near_optimal_{name}_relative_range"] = float(
            (np.max(near_values) - np.min(near_values)) / np.median(near_values)
        )
    outside = [row for row in successful if row not in near]
    if outside:
        alternative = min(outside, key=lambda row: float(row["rmse_k"]))
        summary["alternative_basin_rmse_k"] = float(alternative["rmse_k"])
        summary["alternative_basin_rmse_gap_k"] = float(alternative["rmse_k"]) - best_rmse
        for name in PARAMETER_NAMES:
            summary[f"alternative_basin_{name}"] = float(alternative[f"fit_{name}"])
    return summary


def bound_proximity(parameters: np.ndarray) -> dict[str, float]:
    log_values = np.log(parameters)
    log_lower = np.log(LOWER_BOUNDS)
    log_upper = np.log(UPPER_BOUNDS)
    fraction = np.minimum(log_values - log_lower, log_upper - log_values) / (log_upper - log_lower)
    return {f"{name}_nearest_log_bound_fraction": float(value) for name, value in zip(PARAMETER_NAMES, fraction, strict=True)}


def conditioning_classification(
    rank: int,
    condition: float,
    parameter_correlation: np.ndarray,
    spread: dict[str, object],
    proximity: dict[str, float],
) -> str:
    max_correlation = float(np.max(np.abs(parameter_correlation - np.eye(4))))
    max_spread = max(float(spread[f"near_optimal_{name}_relative_range"]) for name in PARAMETER_NAMES)
    min_bound_fraction = min(proximity.values())
    boundary_collapse = min_bound_fraction <= 1.0e-6 and condition >= 1.0e4
    if rank < 4 or condition >= 1.0e8 or max_spread >= 1.0 or boundary_collapse:
        return "TWO-NODE PARAMETERIZATION NOT RECOVERABLE FROM THIS EXPERIMENT"
    if condition >= 1.0e4 or max_correlation >= 0.95 or max_spread >= 0.20 or min_bound_fraction <= 0.01:
        return "TWO-NODE PARAMETERS PRACTICALLY CORRELATED"
    return "TWO-NODE PARAMETERS WELL CONDITIONED"


def write_multistart(rows: list[dict[str, object]]) -> None:
    with MULTISTART_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def add_metric_rows(
    rows: list[dict[str, object]], scope: str, subset: str, values: dict[str, object], n: int | str
) -> None:
    for name, value in values.items():
        unit = "1"
        if name.endswith("_j_per_k"):
            unit = "J/K"
        elif name.endswith("_w_per_k"):
            unit = "W/K"
        elif name.endswith("_k_per_w"):
            unit = "K/W"
        elif name.endswith("_k"):
            unit = "K"
        elif name.endswith("_c"):
            unit = "degC"
        elif name.endswith("_per_s"):
            unit = "1/s"
        elif name.endswith("_s"):
            unit = "s"
        rows.append({"scope": scope, "subset": subset, "metric": name, "value": value, "unit": unit, "n_samples": n})


def write_metrics(
    one_c_metrics: dict[str, float],
    udds_metrics: dict[str, float],
    parameters: np.ndarray,
    time_constants: np.ndarray,
    eigenvalues: np.ndarray,
    spread: dict[str, object],
    proximity: dict[str, float],
    singular_values: np.ndarray,
    rank: int,
    condition: float,
    sensitivity_correlation: np.ndarray,
    parameter_correlation: np.ndarray,
    jacobian: np.ndarray,
    one_c_samples: int,
    udds_samples: int,
) -> None:
    rows: list[dict[str, object]] = []
    add_metric_rows(rows, "fit", "one_c_soc_0.10_to_0.90", one_c_metrics, one_c_samples)
    add_metric_rows(rows, "forward_validation", "udds_soc_0.90_to_0.10", udds_metrics, udds_samples)
    add_metric_rows(rows, "reference", "one_node_one_c", ONE_NODE_REFERENCE, one_c_samples)
    add_metric_rows(rows, "reference", "one_node_udds_validation_a", ONE_NODE_UDDS_REFERENCE, udds_samples)
    parameter_metrics = dict(zip(PARAMETER_NAMES, parameters, strict=True))
    parameter_metrics["r_cs_k_per_w"] = 1.0 / parameters[2]
    parameter_metrics["r_sa_k_per_w"] = 1.0 / parameters[3]
    add_metric_rows(rows, "parameters", "best", parameter_metrics, one_c_samples)
    add_metric_rows(rows, "dynamics", "best", {
        "fast_time_constant_s": time_constants[0], "slow_time_constant_s": time_constants[1],
        "fast_eigenvalue_per_s": eigenvalues[0], "slow_eigenvalue_per_s": eigenvalues[1],
    }, one_c_samples)
    add_metric_rows(rows, "multistart", "summary", spread, one_c_samples)
    add_metric_rows(rows, "bounds", "best", proximity, one_c_samples)
    add_metric_rows(rows, "identifiability", "local_log_parameter_jacobian", {
        "numerical_rank": rank, "condition_number": condition,
        **{f"singular_value_{index+1}": value for index, value in enumerate(singular_values)},
    }, one_c_samples)
    for row_index, row_name in enumerate(PARAMETER_NAMES):
        for column_index, column_name in enumerate(PARAMETER_NAMES):
            rows.append({"scope": "sensitivity_correlation", "subset": row_name, "metric": column_name, "value": sensitivity_correlation[row_index, column_index], "unit": "1", "n_samples": one_c_samples})
            rows.append({"scope": "parameter_correlation", "subset": row_name, "metric": column_name, "value": parameter_correlation[row_index, column_index], "unit": "1", "n_samples": one_c_samples})
    for sample_index, sensitivity_row in enumerate(jacobian):
        for column_index, name in enumerate(PARAMETER_NAMES):
            rows.append({"scope": "sensitivity_matrix", "subset": f"sample_{sample_index}", "metric": f"d_residual_d_log_{name}", "value": sensitivity_row[column_index], "unit": "K", "n_samples": one_c_samples})
    with METRICS_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def save_plots(
    one_c: object,
    one_c_two: object,
    one_c_one: np.ndarray,
    udds: object,
    udds_two: object,
    udds_one: np.ndarray,
) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(11.5, 7.4), sharex=True, constrained_layout=True)
    minutes = (one_c.time_s - one_c.time_s[0]) / 60.0
    axes[0].plot(minutes, one_c.measured_temperature_k - 273.15, color="black", lw=1.2, label="Measured surface")
    axes[0].plot(minutes, one_c_two.t_surface_k - 273.15, color="#0072B2", lw=1.2, label="Clean two-node surface")
    axes[0].plot(minutes, one_c_two.t_core_k - 273.15, color="#D55E00", lw=1.0, label="Clean two-node core")
    axes[0].plot(minutes, one_c_one - 273.15, color="#CC79A7", lw=0.9, ls="--", label="Frozen one-node")
    axes[0].plot(minutes, one_c.ambient_temperature_k - 273.15, color="#009E73", lw=0.9, label="Chamber")
    axes[0].set(ylabel="Temperature [degC]", title="Panasonic 25 C 1C clean two-node physical fit")
    axes[0].legend(ncol=2)
    axes[0].grid(alpha=0.22)
    axes[1].plot(minutes, one_c_two.t_surface_k - one_c.measured_temperature_k, color="#0072B2", lw=1.0, label="Two-node residual")
    axes[1].plot(minutes, one_c_one - one_c.measured_temperature_k, color="#CC79A7", lw=0.9, ls="--", label="One-node residual")
    axes[1].axhline(0.0, color="black", lw=0.7)
    axes[1].set(xlabel="Elapsed time [min]", ylabel="Predicted - measured [K]")
    axes[1].legend()
    axes[1].grid(alpha=0.22)
    fig.savefig(FIT_PLOT_PATH, dpi=180)
    plt.close(fig)

    shown = np.arange(0, udds.time_s.size, max(1, udds.time_s.size // 45000))
    hours = udds.time_s / 3600.0
    fig, axes = plt.subplots(2, 1, figsize=(12.0, 7.5), sharex=True, constrained_layout=True)
    axes[0].plot(hours[shown], udds.measured_temperature_k[shown] - 273.15, color="black", lw=1.1, label="Measured surface")
    axes[0].plot(hours[shown], udds_two.t_surface_k[shown] - 273.15, color="#0072B2", lw=1.0, label="Clean two-node surface")
    axes[0].plot(hours[shown], udds_one[shown] - 273.15, color="#CC79A7", lw=0.9, ls="--", label="Existing frozen one-node")
    axes[0].plot(hours[shown], udds.ambient_temperature_k[shown] - 273.15, color="#009E73", lw=0.9, label="Chamber")
    axes[0].plot(hours[shown], udds_two.t_core_k[shown] - 273.15, color="#D55E00", lw=0.8, alpha=0.8, label="Clean two-node core")
    axes[0].set(ylabel="Temperature [degC]", title="Panasonic 25 C UDDS thermal transfer: no refit, measured-voltage heat")
    axes[0].legend(ncol=2)
    axes[0].grid(alpha=0.22)
    axes[1].plot(hours[shown], udds_two.t_surface_k[shown] - udds.measured_temperature_k[shown], color="#0072B2", lw=0.9, label="Two-node residual")
    axes[1].plot(hours[shown], udds_one[shown] - udds.measured_temperature_k[shown], color="#CC79A7", lw=0.9, ls="--", label="One-node residual")
    axes[1].axhline(0.0, color="black", lw=0.7)
    axes[1].set(xlabel="Elapsed time in SOC 0.90-to-0.10 interval [h]", ylabel="Predicted - measured [K]")
    axes[1].legend()
    axes[1].grid(alpha=0.22)
    fig.savefig(UDDS_PLOT_PATH, dpi=180)
    plt.close(fig)


def matrix_markdown(matrix: np.ndarray) -> list[str]:
    lines = ["| | " + " | ".join(PARAMETER_NAMES) + " |", "|---|" + "---:|" * 4]
    for name, row in zip(PARAMETER_NAMES, matrix, strict=True):
        lines.append(f"| {name} | " + " | ".join(f"{value:+.6f}" for value in row) + " |")
    return lines


def write_report(
    one_c_metrics: dict[str, float],
    udds_metrics: dict[str, float],
    parameters: np.ndarray,
    time_constants: np.ndarray,
    eigenvalues: np.ndarray,
    spread: dict[str, object],
    proximity: dict[str, float],
    jacobian: np.ndarray,
    singular_values: np.ndarray,
    rank: int,
    condition: float,
    sensitivity_correlation: np.ndarray,
    parameter_correlation: np.ndarray,
    classification: str,
    one_c_two: object,
    udds_two: object,
    one_c_min_dt_s: float,
) -> None:
    one_c_improvement = 1.0 - one_c_metrics["rmse_k"] / ONE_NODE_REFERENCE["rmse_k"]
    udds_improvement = 1.0 - udds_metrics["rmse_k"] / ONE_NODE_UDDS_REFERENCE["rmse_k"]
    clean_outperforms_udds = udds_improvement > 0.0
    justify = clean_outperforms_udds and udds_metrics["r2"] > 0.0 and "NOT RECOVERABLE" not in classification
    correction_dependency = (
        "The archived SOC multiplier improved its legacy 1C RMSE (0.727691 to 0.655786 K), "
        "but worsened legacy UDDS RMSE (0.204643 to 0.222634 K). Therefore the old low UDDS "
        "error did not depend on the correction; the apparent late-1C improvement did."
    )
    lines = [
        "# Panasonic clean physical two-node calibration and validation",
        "",
        f"## {classification}",
        "",
        "Only the four legacy physical parameters were fitted, and only on the closed Panasonic 25 C 1C interval `0.10 <= SOC <= 0.90`. The independent UDDS interval was evaluated once after freezing the optimum. No electrical/OCP/entropy parameter, initial temperature, empirical heat term, or correction factor was fitted or modified.",
        "",
        "## Clean architecture",
        "",
        "```text",
        "C_c dT_c/dt = Q_gen - G_cs (T_c - T_s)",
        "C_s dT_s/dt = G_cs (T_c - T_s) - hA (T_s - T_amb)",
        "measured output = T_s",
        "```",
        "",
        "The clean module preserves the legacy conductance parameterization (`G_cs=1/R_cs`, `hA=1/R_sa`) and uses an exact zero-order-hold state transition at every actual nonuniform interval. `T_c(t0)=T_s(t0)=T_measured(t0)` in calibration and validation, matching the legacy initialization. There is no interface for SOC multipliers, regime gates, heat clipping, additive heat, or temperature bias.",
        "",
        "Heat used for both 1C fitting and UDDS validation is exactly `I*(Ueq_average-V_measured) - I*T_measured*dUeq/dT`, with discharge-positive internal current, frozen average electrode states/OCPs, and the verified Panasonic entropy table.",
        "",
        "## Datasets",
        "",
        f"- 1C calibration source: `{ONE_C_RAW_SOURCE_PATH}`. The exact entropy-supported discharge segment runs from SOC 0.90 to 0.10 (`0.10 <= SOC <= 0.90`), from source times {one_c_two.time_s[0]:.9f} to {one_c_two.time_s[-1]:.9f} s ({one_c_two.time_s.size} samples; duration {one_c_two.time_s[-1]-one_c_two.time_s[0]:.9f} s).",
        f"- Independent UDDS validation source: `{UDDS_RAW_SOURCE_PATH}`. The exact continuous SOC 0.90-to-0.10 interval contains {udds_two.time_s.size:,} samples over {udds_two.time_s[-1]-udds_two.time_s[0]:.9f} s.",
        "",
        "## 1C calibration result",
        "",
        f"Forty deterministic log-space least-squares starts were run with broad positivity guard bounds. {spread['successful_count']} succeeded; {spread['near_optimal_count']} are within {100*NEAR_OPTIMAL_RELATIVE_RMSE:.3f}% relative RMSE of the best.",
        "",
        "| Parameter | Fitted value | Guard bounds | Nearest log-bound fraction |",
        "|---|---:|---:|---:|",
    ]
    bounds_text = ("[1e-3, 5000] J/K", "[1e-4, 2000] J/K", "[1e-5, 100] W/K", "[1e-5, 10] W/K")
    for index, name in enumerate(PARAMETER_NAMES):
        lines.append(f"| `{name}` | {parameters[index]:.12g} | {bounds_text[index]} | {proximity[f'{name}_nearest_log_bound_fraction']:.6f} |")
    lines.extend([
        "",
        f"Modal eigenvalues are {eigenvalues[0]:.12g} and {eigenvalues[1]:.12g} 1/s, giving fast and slow time constants {time_constants[0]:.9f} and {time_constants[1]:.9f} s.",
        "",
        "| Metric | Clean two-node 1C | Frozen one-node 1C |",
        "|---|---:|---:|",
        f"| RMSE [K] | {one_c_metrics['rmse_k']:.9f} | {ONE_NODE_REFERENCE['rmse_k']:.9f} |",
        f"| MAE [K] | {one_c_metrics['mae_k']:.9f} | {ONE_NODE_REFERENCE['mae_k']:.9f} |",
        f"| Bias [K] | {one_c_metrics['mean_signed_error_k']:+.9f} | {ONE_NODE_REFERENCE['mean_signed_error_k']:+.9f} |",
        f"| Maximum absolute error [K] | {one_c_metrics['max_abs_error_k']:.9f} | not supplied in closure prompt |",
        f"| R2 | {one_c_metrics['r2']:.9f} | {ONE_NODE_REFERENCE['r2']:.9f} |",
        "",
        f"The corresponding resistances are `R_cs={1.0/parameters[2]:.9f} K/W` and `R_sa={1.0/parameters[3]:.9f} K/W`.",
        "",
        f"The clean two-node RMSE change relative to one-node is {100*one_c_improvement:+.3f}%. This is a small in-sample change, and the independent-transfer result below - not visual appearance - is used to decide whether the added state is worthwhile.",
        "",
        "## Multistart stability",
        "",
        "Near-optimal means RMSE no more than 0.1% above the best; the complete successful-start table is in the CSV output.",
        "",
        "| Parameter | All-successful min | All-successful max | Near-optimal min | Near-optimal max | Near-optimal relative range |",
        "|---|---:|---:|---:|---:|---:|",
    ])
    for name in PARAMETER_NAMES:
        lines.append(
            f"| `{name}` | {spread[f'successful_{name}_min']:.9g} | {spread[f'successful_{name}_max']:.9g} | {spread[f'near_optimal_{name}_min']:.9g} | {spread[f'near_optimal_{name}_max']:.9g} | {100*spread[f'near_optimal_{name}_relative_range']:.6f}% |"
        )
    lines.extend([
        "",
        "## Local practical conditioning",
        "",
        f"The numerical sensitivity matrix is the `{jacobian.shape[0]} x {jacobian.shape[1]}` Jacobian of surface-temperature residual [K] with respect to dimensionless log-parameters. All 1,196 matrix entries are stored under `scope=sensitivity_matrix` in `panasonic_two_node_physical_metrics.csv`. Its numerical rank is **{rank}/4**, singular values are `" + ", ".join(f"{value:.9e}" for value in singular_values) + f"`, and its condition number is **{condition:.9e}**.",
        "",
        "Sensitivity-column correlation matrix:",
        "",
        *matrix_markdown(sensitivity_correlation),
        "",
        "Local parameter correlation matrix inferred from the pseudoinverse of `J^T J`:",
        "",
        *matrix_markdown(parameter_correlation),
        "",
        "This is a local practical-conditioning audit based on one measured surface-temperature trace. It is not a structural-identifiability claim.",
        "",
        f"The fitted fast mode is {time_constants[0]:.9g} s, whereas the smallest retained 1C interval is {one_c_min_dt_s:.9g} s. The optimizer therefore pushes the surface node below the experiment's temporal resolution instead of identifying a resolved second thermal mode.",
        "",
        f"A distinct secondary multistart basin has RMSE {spread.get('alternative_basin_rmse_k', float('nan')):.9f} K (only {spread.get('alternative_basin_rmse_gap_k', float('nan')):.9f} K above the best) at `C_c={spread.get('alternative_basin_c_c_j_per_k', float('nan')):.6g}`, `C_s={spread.get('alternative_basin_c_s_j_per_k', float('nan')):.6g}`, `G_cs={spread.get('alternative_basin_g_cs_w_per_k', float('nan')):.6g}`, `hA={spread.get('alternative_basin_hA_w_per_k', float('nan')):.6g}`. Its surface trajectory differs from the best by only {spread.get('alternative_basin_surface_difference_rmse_k', float('nan')):.9f} K RMSE (maximum {spread.get('alternative_basin_surface_difference_max_k', float('nan')):.9f} K), despite radically different parameters. The four physical quantities are therefore not practically recoverable from this surface-only experiment.",
        "",
        "## Critical no-refit UDDS validation",
        "",
        "Both models use the same measured-voltage/reversible measured-temperature heat and exact SOC 0.90-to-0.10 interval. All fitted two-node parameters were frozen before this evaluation.",
        "",
        "| Metric | Clean two-node UDDS | Existing one-node UDDS A |",
        "|---|---:|---:|",
        f"| RMSE [K] | {udds_metrics['rmse_k']:.9f} | {ONE_NODE_UDDS_REFERENCE['rmse_k']:.9f} |",
        f"| MAE [K] | {udds_metrics['mae_k']:.9f} | {ONE_NODE_UDDS_REFERENCE['mae_k']:.9f} |",
        f"| Bias [K] | {udds_metrics['mean_signed_error_k']:+.9f} | {ONE_NODE_UDDS_REFERENCE['mean_signed_error_k']:+.9f} |",
        f"| Maximum absolute error [K] | {udds_metrics['max_abs_error_k']:.9f} | {ONE_NODE_UDDS_REFERENCE['max_abs_error_k']:.9f} |",
        f"| R2 | {udds_metrics['r2']:.9f} | {ONE_NODE_UDDS_REFERENCE['r2']:.9f} |",
        f"| Final measured surface [degC] | {udds_metrics['measured_final_temperature_c']:.9f} | {udds_metrics['measured_final_temperature_c']:.9f} |",
        f"| Final predicted [degC] | {udds_metrics['predicted_final_temperature_c']:.9f} | {ONE_NODE_UDDS_REFERENCE['predicted_final_temperature_c']:.9f} |",
        f"| Measured rise [K] | {udds_metrics['measured_temperature_rise_k']:+.9f} | {udds_metrics['measured_temperature_rise_k']:+.9f} |",
        f"| Predicted rise [K] | {udds_metrics['predicted_temperature_rise_k']:+.9f} | {ONE_NODE_UDDS_REFERENCE['predicted_temperature_rise_k']:+.9f} |",
        "",
        f"The clean two-node UDDS RMSE change relative to one-node is {100*udds_improvement:+.3f}%. Maximum predicted `T_c-T_s` is {np.max(udds_two.t_core_k-udds_two.t_surface_k):.9f} K (UDDS) and {np.max(one_c_two.t_core_k-one_c_two.t_surface_k):.9f} K (1C).",
        "",
        "## Final decisions",
        "",
        f"1. **Did the old apparent thermal improvement depend on the empirical SOC/regime correction?** {correction_dependency}",
        "",
        f"2. **Does the clean physical two-node model outperform the current one-node model on independent UDDS?** **{'Yes' if clean_outperforms_udds else 'No'}**; its RMSE change is {100*udds_improvement:+.3f}% with no refit or heat correction.",
        "",
        f"3. **Is the improvement large enough to justify the two-node model in the thesis?** **{'Yes' if justify else 'No'}**. The decision is based on independent UDDS transfer and recoverability: the clean model must actually improve held-out behavior and supply a defensible finite parameterization.",
        "",
        f"4. **Are the parameters practically recoverable from the available temperature data?** **{classification}.** The conclusion follows from rank, singular spectrum/condition number, correlations, multistart spread, and bound proximity; no structural-identifiability claim is made.",
        "",
        "No model-voltage heat validation, CC-CV thermal validation, own-TRC analysis, or empirical repair was performed.",
    ])
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    if not LEGACY_AUDIT_PATH.is_file():
        raise FileNotFoundError("Legacy audit must exist before calibration.")
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    one_c = prepare_one_c_data()
    rows, best = run_multistart(one_c)
    write_multistart(rows)
    best_result = best["result"]
    best_log = np.asarray(best_result.x)
    fitted_values = np.exp(best_log)
    fitted_parameters = vector_to_parameters(fitted_values)
    one_c_two = simulate_two_node_physical(
        one_c.time_s,
        one_c.qdot_cal_w,
        one_c.ambient_temperature_k,
        fitted_parameters,
        float(one_c.measured_temperature_k[0]),
        float(one_c.measured_temperature_k[0]),
    )
    one_c_one = one_node_trajectory(
        one_c.time_s,
        one_c.qdot_cal_w,
        one_c.ambient_temperature_k,
        float(one_c.measured_temperature_k[0]),
        C_TH_J_PER_K,
        G_TH_W_PER_K,
    )
    one_c_metrics = metrics(one_c.measured_temperature_k, one_c_two.t_surface_k)
    one_c_one_metrics = metrics(one_c.measured_temperature_k, one_c_one)
    for key in ("rmse_k", "mae_k", "mean_signed_error_k", "r2"):
        if not np.isclose(one_c_one_metrics[key], ONE_NODE_REFERENCE[key], atol=8.0e-9, rtol=0.0):
            raise RuntimeError(f"Frozen one-node 1C reference changed for {key}.")

    udds, _ = prepare_udds_data()
    qirr_udds = udds.current_a * (udds.ueq_average_v - udds.measured_voltage_v)
    qrev_udds = -udds.current_a * udds.measured_temperature_k * udds.dudt_v_per_k
    heat_udds = qirr_udds + qrev_udds
    udds_two = simulate_two_node_physical(
        udds.time_s,
        heat_udds,
        udds.ambient_temperature_k,
        fitted_parameters,
        float(udds.measured_temperature_k[0]),
        float(udds.measured_temperature_k[0]),
    )
    udds_one = one_node_trajectory(
        udds.time_s,
        heat_udds,
        udds.ambient_temperature_k,
        float(udds.measured_temperature_k[0]),
        C_TH_J_PER_K,
        G_TH_W_PER_K,
    )
    udds_metrics = metrics(udds.measured_temperature_k, udds_two.t_surface_k)
    udds_one_metrics = metrics(udds.measured_temperature_k, udds_one)
    for key, reference in ONE_NODE_UDDS_REFERENCE.items():
        if not np.isclose(udds_one_metrics[key], reference, atol=8.0e-9, rtol=0.0):
            raise RuntimeError(f"Frozen one-node UDDS reference changed for {key}.")

    jacobian, singular_values, rank, condition, sensitivity_correlation, parameter_correlation = numerical_sensitivity(best_log, one_c)
    spread = multistart_summary(rows)
    if "alternative_basin_rmse_k" in spread:
        alternative_values = np.asarray(
            [spread[f"alternative_basin_{name}"] for name in PARAMETER_NAMES], dtype=float
        )
        alternative = simulate_two_node_physical(
            one_c.time_s,
            one_c.qdot_cal_w,
            one_c.ambient_temperature_k,
            vector_to_parameters(alternative_values),
            float(one_c.measured_temperature_k[0]),
            float(one_c.measured_temperature_k[0]),
        )
        trajectory_difference = alternative.t_surface_k - one_c_two.t_surface_k
        spread["alternative_basin_surface_difference_rmse_k"] = float(
            np.sqrt(np.mean(trajectory_difference**2))
        )
        spread["alternative_basin_surface_difference_max_k"] = float(
            np.max(np.abs(trajectory_difference))
        )
    proximity = bound_proximity(fitted_values)
    classification = conditioning_classification(rank, condition, parameter_correlation, spread, proximity)
    eigenvalues = thermal_eigenvalues(fitted_parameters)
    time_constants = thermal_time_constants_s(fitted_parameters)

    save_plots(one_c, one_c_two, one_c_one, udds, udds_two, udds_one)
    write_metrics(
        one_c_metrics, udds_metrics, fitted_values, time_constants, eigenvalues,
        spread, proximity, singular_values, rank, condition,
        sensitivity_correlation, parameter_correlation, jacobian, one_c.time_s.size, udds.time_s.size,
    )
    write_report(
        one_c_metrics, udds_metrics, fitted_values, time_constants, eigenvalues,
        spread, proximity, jacobian, singular_values, rank, condition,
        sensitivity_correlation, parameter_correlation, classification, one_c_two, udds_two,
        float(np.min(np.diff(one_c.time_s))),
    )
    print("Clean two-node Panasonic calibration/validation complete.")
    print("Parameters: " + ", ".join(f"{name}={value:.12g}" for name, value in zip(PARAMETER_NAMES, fitted_values, strict=True)))
    print(f"1C RMSE={one_c_metrics['rmse_k']:.9f} K; one-node={ONE_NODE_REFERENCE['rmse_k']:.9f} K")
    print(f"UDDS RMSE={udds_metrics['rmse_k']:.9f} K; one-node={ONE_NODE_UDDS_REFERENCE['rmse_k']:.9f} K")
    print(f"Rank={rank}/4; condition={condition:.6e}; {classification}")


if __name__ == "__main__":
    main()
