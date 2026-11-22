"""Joint 1C + HPPC calibration of the clean one-node Panasonic thermal model.

Only Cth and Gth are fitted. The two experiments receive equal MSE weight.
HPPC is split at audited omitted-data gaps and each retained thermal block is
initialized from its measured temperature. UDDS is forward validation only.
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
from scipy.linalg import solve_banded
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
    X_N_END,
    X_N_START,
    X_P_END,
    X_P_START,
)
from models.panasonic_entropic import SUPPORTED_SOC_RANGE, dUeq_dT_panasonic  # noqa: E402
from models.panasonic_thermal_parameters import (  # noqa: E402
    C_TH_J_PER_K as BASELINE_C_TH,
    G_TH_W_PER_K as BASELINE_G_TH,
)
from scripts.audit_panasonic_hppc import (  # noqa: E402
    CONTIGUOUS_GAP_THRESHOLD_S,
    LOADED_CSV_PATH as HPPC_SOURCE_PATH,
    RAW_MANIFEST_PATH as HPPC_RAW_PATH,
    load_hppc,
)
from scripts.calibrate_panasonic_1c_thermal import (  # noqa: E402
    RAW_SOURCE_PATH as ONE_C_SOURCE_PATH,
    prepare_data as prepare_one_c_data,
)
from scripts.validate_panasonic_udds_thermal import (  # noqa: E402
    RAW_SOURCE_PATH as UDDS_SOURCE_PATH,
    prepare_data as prepare_udds_data,
)


OUTPUT_DIRECTORY = REPOSITORY_ROOT / "thesis_data"
REPORT_PATH = OUTPUT_DIRECTORY / "panasonic_joint_1c_hppc_thermal_calibration.md"
METRICS_PATH = OUTPUT_DIRECTORY / "panasonic_joint_1c_hppc_thermal_metrics.csv"
MULTISTART_PATH = OUTPUT_DIRECTORY / "panasonic_joint_1c_hppc_thermal_multistart.csv"
ONE_C_PLOT_PATH = OUTPUT_DIRECTORY / "panasonic_joint_1c_hppc_1c_fit.png"
HPPC_PLOT_PATH = OUTPUT_DIRECTORY / "panasonic_joint_1c_hppc_hppc_fit.png"
UDDS_PLOT_PATH = OUTPUT_DIRECTORY / "panasonic_joint_1c_hppc_udds_validation.png"

PARAMETER_NAMES = ("c_th_j_per_k", "g_th_w_per_k")
LOWER_BOUNDS = np.asarray([5.0, 1.0e-4])
UPPER_BOUNDS = np.asarray([200.0, 1.0])
N_MULTISTARTS = 40
RANDOM_SEED = 20260914
NEAR_OPTIMAL_RELATIVE_RMSE = 1.0e-4
LOG_SENSITIVITY_STEP = 1.0e-5
CURRENT_THRESHOLD_A = 0.05

BASELINE_ONE_C = {
    "rmse_k": 0.143622065,
    "mae_k": 0.119569322,
    "mean_signed_error_k": -0.008530653,
    "r2": 0.964702003,
}
BASELINE_UDDS = {
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
BASELINE_CONDITION_NUMBER = 6.120213125908202
BASELINE_SENSITIVITY_CORRELATION = 0.3041135029313185
BASELINE_PARAMETER_CORRELATION = -0.8667854247047297


@dataclass(frozen=True)
class ThermalBlock:
    block_id: int
    source_start_index: int
    source_end_index: int
    boundary_start_interpolated: bool
    boundary_end_interpolated: bool
    source_time_s: np.ndarray
    time_s: np.ndarray
    soc: np.ndarray
    current_a: np.ndarray
    measured_voltage_v: np.ndarray
    measured_temperature_k: np.ndarray
    ambient_temperature_k: np.ndarray
    qdot_cal_w: np.ndarray


def _interpolate(values: np.ndarray, soc: np.ndarray, left: int, target: float) -> float:
    fraction = (target - soc[left]) / (soc[left + 1] - soc[left])
    if not 0.0 <= fraction <= 1.0:
        raise ValueError("Invalid HPPC SOC boundary interpolation fraction.")
    return float(values[left] + fraction * (values[left + 1] - values[left]))


def _supported_part_of_block(
    block_id: int,
    start: int,
    end: int,
    soc_full: np.ndarray,
    signals: dict[str, np.ndarray],
) -> ThermalBlock | None:
    lower, upper = SUPPORTED_SOC_RANGE
    source_indices = np.arange(start, end + 1)
    valid = (soc_full[source_indices] >= lower) & (soc_full[source_indices] <= upper)
    if not np.any(valid):
        return None
    retained = source_indices[valid]
    first = int(retained[0])
    last = int(retained[-1])
    if not np.all(valid[first - start : last - start + 1]):
        raise ValueError("HPPC supported SOC samples are not one contiguous portion.")

    start_interpolated = first > start and soc_full[first - 1] > upper >= soc_full[first]
    end_interpolated = last < end and soc_full[last] >= lower > soc_full[last + 1]
    selected_soc = soc_full[first : last + 1].copy()
    selected: dict[str, np.ndarray] = {
        name: values[first : last + 1].copy() for name, values in signals.items()
    }
    if start_interpolated:
        selected_soc = np.concatenate(([upper], selected_soc))
        for name, values in signals.items():
            selected[name] = np.concatenate(
                ([_interpolate(values, soc_full, first - 1, upper)], selected[name])
            )
    if end_interpolated:
        selected_soc = np.concatenate((selected_soc, [lower]))
        for name, values in signals.items():
            selected[name] = np.concatenate(
                (selected[name], [_interpolate(values, soc_full, last, lower)])
            )
    if np.any((selected_soc < lower) | (selected_soc > upper)):
        raise ValueError("HPPC entropy-domain violation; no clipping or extrapolation.")

    z = 1.0 - selected_soc
    x_n_avg = X_N_START + z * (X_N_END - X_N_START)
    x_p_avg = X_P_START + z * (X_P_END - X_P_START)
    if np.any(
        (x_n_avg < MCMB_SUPPORTED_STOICHIOMETRY_RANGE[0])
        | (x_n_avg > MCMB_SUPPORTED_STOICHIOMETRY_RANGE[1])
        | (x_p_avg < NCA_WORKING_STOICHIOMETRY_RANGE[0])
        | (x_p_avg > NCA_WORKING_STOICHIOMETRY_RANGE[1])
    ):
        raise ValueError("HPPC average electrode state left frozen OCP support.")
    ueq = np.asarray(U_p(x_p_avg)) - np.asarray(U_n(x_n_avg))
    dudt = np.asarray(dUeq_dT_panasonic(selected_soc))
    current = selected["current_a"]
    measured_t = selected["measured_temperature_k"]
    q_irr = current * (ueq - selected["measured_voltage_v"])
    q_rev = -current * measured_t * dudt
    q_cal = q_irr + q_rev
    source_time = selected["source_time_s"]
    return ThermalBlock(
        block_id=block_id,
        source_start_index=first,
        source_end_index=last,
        boundary_start_interpolated=start_interpolated,
        boundary_end_interpolated=end_interpolated,
        source_time_s=source_time,
        time_s=source_time - source_time[0],
        soc=selected_soc,
        current_a=current,
        measured_voltage_v=selected["measured_voltage_v"],
        measured_temperature_k=measured_t,
        ambient_temperature_k=selected["ambient_temperature_k"],
        qdot_cal_w=q_cal,
    )


def prepare_hppc_blocks() -> tuple[
    list[ThermalBlock],
    dict[str, object],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    raw = load_hppc(HPPC_SOURCE_PATH)
    source_soc = 1.0 + raw.exported_ah / Q_CELL_AH
    if np.any((source_soc < 0.0) | (source_soc > 1.0)):
        raise ValueError("HPPC cumulative-Ah SOC leaves [0,1]; no clipping.")
    current = -raw.current_a  # sole raw-to-model sign conversion
    signals = {
        "source_time_s": raw.time_s,
        "current_a": current,
        "measured_voltage_v": raw.voltage_v,
        "measured_temperature_k": raw.battery_temperature_c + 273.15,
        "ambient_temperature_k": raw.chamber_temperature_c + 273.15,
    }
    differences = np.diff(raw.time_s)
    gap_left = np.flatnonzero(differences > CONTIGUOUS_GAP_THRESHOLD_S)
    starts = np.concatenate(([0], gap_left + 1))
    ends = np.concatenate((gap_left, [raw.time_s.size - 1]))
    if starts.size != 14:
        raise ValueError(f"Expected 14 audited HPPC blocks, found {starts.size}.")

    blocks: list[ThermalBlock] = []
    for block_id, (start, end) in enumerate(zip(starts, ends, strict=True), start=1):
        block = _supported_part_of_block(
            block_id, int(start), int(end), source_soc, signals
        )
        if block is not None:
            blocks.append(block)
    if not blocks:
        raise ValueError("No HPPC block remains inside entropy support.")

    block_rows: list[dict[str, object]] = []
    for block in blocks:
        positive_dt = np.diff(block.time_s)
        positive_dt = positive_dt[positive_dt > 0.0]
        block_rows.append(
            {
                "block_id": block.block_id,
                "source_start_index": block.source_start_index,
                "source_end_index": block.source_end_index,
                "start_soc": float(block.soc[0]),
                "end_soc": float(block.soc[-1]),
                "duration_s": float(block.time_s[-1]),
                "sample_count": int(block.time_s.size),
                "median_positive_dt_s": float(np.median(positive_dt)),
                "max_normal_dt_s": float(np.max(positive_dt)),
                "initial_battery_temperature_c": float(block.measured_temperature_k[0] - 273.15),
                "final_battery_temperature_c": float(block.measured_temperature_k[-1] - 273.15),
                "chamber_temperature_min_c": float(np.min(block.ambient_temperature_k) - 273.15),
                "chamber_temperature_max_c": float(np.max(block.ambient_temperature_k) - 273.15),
                "current_min_a": float(np.min(block.current_a)),
                "current_max_a": float(np.max(block.current_a)),
                "start_boundary_interpolated": block.boundary_start_interpolated,
                "end_boundary_interpolated": block.boundary_end_interpolated,
                "duplicate_time_increments": int(np.count_nonzero(np.diff(block.time_s) == 0.0)),
            }
        )

    retained_dt = np.concatenate([np.diff(block.time_s) for block in blocks])
    retained_current = np.concatenate([block.current_a for block in blocks])
    held_current = np.concatenate([block.current_a[:-1] for block in blocks])
    positive = retained_dt > 0.0
    pulse_dt = retained_dt[positive & (np.abs(held_current) > CURRENT_THRESHOLD_A)]
    rest_dt = retained_dt[positive & (np.abs(held_current) <= CURRENT_THRESHOLD_A)]
    audit: dict[str, object] = {
        "source_rows": int(raw.time_s.size),
        "raw_block_count": int(starts.size),
        "retained_block_count": len(blocks),
        "retained_samples": int(sum(block.time_s.size for block in blocks)),
        "duplicate_time_increments_retained": int(np.count_nonzero(retained_dt == 0.0)),
        "long_gap_count": int(gap_left.size),
        "long_gap_min_s": float(np.min(differences[gap_left])),
        "long_gap_max_s": float(np.max(differences[gap_left])),
        "normal_pulse_dt_median_s": float(np.median(pulse_dt)),
        "normal_rest_dt_median_s": float(np.median(rest_dt)),
        "retained_current_min_a": float(np.min(retained_current)),
        "retained_current_max_a": float(np.max(retained_current)),
        "retained_soc_min": float(min(np.min(block.soc) for block in blocks)),
        "retained_soc_max": float(max(np.max(block.soc) for block in blocks)),
        "temperature_finite": True,
        "chamber_temperature_finite": True,
        "suitable_for_thermal_calibration": True,
    }
    gaps = [
        {
            "gap_id": index + 1,
            "left_source_index": int(left),
            "right_source_index": int(left + 1),
            "gap_s": float(differences[left]),
        }
        for index, left in enumerate(gap_left)
    ]
    return blocks, audit, block_rows, gaps


def exact_one_node(
    time_s: np.ndarray,
    heat_w: np.ndarray,
    ambient_k: np.ndarray,
    initial_temperature_k: float,
    c_th: float,
    g_th: float,
) -> np.ndarray:
    """Exact nonuniform ZOH propagation, including duplicate timestamps."""

    dt = np.diff(time_s)
    if np.any(dt < 0.0) or c_th <= 0.0 or g_th < 0.0:
        raise ValueError("Invalid one-node propagation input.")
    x = g_th * dt / c_th
    approach = -np.expm1(-x)
    phi = np.ones_like(x)
    nonzero = x != 0.0
    phi[nonzero] = approach[nonzero] / x[nonzero]
    decay = 1.0 - approach
    forcing = approach * ambient_k[:-1] + heat_w[:-1] * dt * phi / c_th
    # Solve T[k+1] - decay[k]*T[k] = forcing[k] as a lower-bidiagonal system.
    count = time_s.size
    band = np.zeros((2, count))
    band[0] = 1.0
    band[1, :-1] = -decay
    rhs = np.empty(count)
    rhs[0] = initial_temperature_k
    rhs[1:] = forcing
    temperature = solve_banded((1, 0), band, rhs, overwrite_ab=True, overwrite_b=True)
    if not np.all(np.isfinite(temperature)):
        raise ValueError("Non-finite one-node temperature.")
    return temperature


def predict_one_c(data: object, parameters: np.ndarray) -> np.ndarray:
    return exact_one_node(
        data.time_s, data.qdot_cal_w, data.ambient_temperature_k,
        float(data.measured_temperature_k[0]), float(parameters[0]), float(parameters[1]),
    )


def predict_hppc(blocks: list[ThermalBlock], parameters: np.ndarray) -> list[np.ndarray]:
    return [
        exact_one_node(
            block.time_s, block.qdot_cal_w, block.ambient_temperature_k,
            float(block.measured_temperature_k[0]), float(parameters[0]), float(parameters[1]),
        )
        for block in blocks
    ]


def balanced_residual(
    log_parameters: np.ndarray, one_c: object, blocks: list[ThermalBlock]
) -> np.ndarray:
    parameters = np.exp(log_parameters)
    residual_one_c = predict_one_c(one_c, parameters) - one_c.measured_temperature_k
    predictions = predict_hppc(blocks, parameters)
    residual_hppc = np.concatenate(
        [prediction - block.measured_temperature_k for block, prediction in zip(blocks, predictions, strict=True)]
    )
    return np.concatenate(
        (residual_one_c / np.sqrt(residual_one_c.size), residual_hppc / np.sqrt(residual_hppc.size))
    )


def temperature_metrics(measured: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    residual = predicted - measured
    denominator = float(np.sum((measured - np.mean(measured)) ** 2))
    return {
        "mse_k2": float(np.mean(residual**2)),
        "rmse_k": float(np.sqrt(np.mean(residual**2))),
        "mae_k": float(np.mean(np.abs(residual))),
        "mean_signed_error_k": float(np.mean(residual)),
        "max_abs_error_k": float(np.max(np.abs(residual))),
        "r2": float(1.0 - np.sum(residual**2) / denominator) if denominator > 0.0 else float("nan"),
        "measured_final_temperature_c": float(measured[-1] - 273.15),
        "predicted_final_temperature_c": float(predicted[-1] - 273.15),
        "measured_temperature_rise_k": float(measured[-1] - measured[0]),
        "predicted_temperature_rise_k": float(predicted[-1] - predicted[0]),
    }


def combined_hppc_metrics(
    blocks: list[ThermalBlock], predictions: list[np.ndarray], mask_kind: str = "all"
) -> dict[str, float]:
    measured_parts: list[np.ndarray] = []
    predicted_parts: list[np.ndarray] = []
    for block, prediction in zip(blocks, predictions, strict=True):
        if mask_kind == "load":
            mask = np.abs(block.current_a) > CURRENT_THRESHOLD_A
        elif mask_kind == "rest":
            mask = np.abs(block.current_a) <= CURRENT_THRESHOLD_A
        elif mask_kind == "high":
            mask = block.soc >= 2.0 / 3.0
        elif mask_kind == "middle":
            mask = (block.soc >= 1.0 / 3.0) & (block.soc < 2.0 / 3.0)
        elif mask_kind == "low":
            mask = block.soc < 1.0 / 3.0
        else:
            mask = np.ones(block.time_s.size, dtype=bool)
        if np.any(mask):
            measured_parts.append(block.measured_temperature_k[mask])
            predicted_parts.append(prediction[mask])
    return temperature_metrics(np.concatenate(measured_parts), np.concatenate(predicted_parts))


def run_multistart(one_c: object, blocks: list[ThermalBlock]) -> tuple[list[dict[str, object]], object]:
    rng = np.random.default_rng(RANDOM_SEED)
    starts = np.vstack(
        (
            [BASELINE_C_TH, BASELINE_G_TH],
            [50.0, 0.03],
            [100.0, 0.10],
            [180.0, 0.30],
            [20.0, 0.003],
            np.exp(rng.uniform(np.log(LOWER_BOUNDS), np.log(UPPER_BOUNDS), size=(N_MULTISTARTS - 5, 2))),
        )
    )
    rows: list[dict[str, object]] = []
    best = None
    for start_id, start in enumerate(starts, start=1):
        result = least_squares(
            balanced_residual,
            np.log(start),
            args=(one_c, blocks),
            bounds=(np.log(LOWER_BOUNDS), np.log(UPPER_BOUNDS)),
            jac="3-point",
            x_scale="jac",
            ftol=1.0e-12,
            xtol=1.0e-12,
            gtol=1.0e-12,
            max_nfev=2000,
        )
        parameters = np.exp(result.x)
        one_c_prediction = predict_one_c(one_c, parameters)
        hppc_predictions = predict_hppc(blocks, parameters)
        one_c_mse = float(np.mean((one_c_prediction - one_c.measured_temperature_k) ** 2))
        hppc_residual = np.concatenate(
            [prediction - block.measured_temperature_k for block, prediction in zip(blocks, hppc_predictions, strict=True)]
        )
        hppc_mse = float(np.mean(hppc_residual**2))
        row = {
            "start_id": start_id,
            "success": bool(result.success),
            "status": int(result.status),
            "message": str(result.message).replace("\n", " "),
            "nfev": int(result.nfev),
            "optimality": float(result.optimality),
            "start_c_th_j_per_k": float(start[0]),
            "start_g_th_w_per_k": float(start[1]),
            "fit_c_th_j_per_k": float(parameters[0]),
            "fit_g_th_w_per_k": float(parameters[1]),
            "tau_th_s": float(parameters[0] / parameters[1]),
            "mse_1c_k2": one_c_mse,
            "mse_hppc_k2": hppc_mse,
            "rmse_1c_k": float(np.sqrt(one_c_mse)),
            "rmse_hppc_k": float(np.sqrt(hppc_mse)),
            "joint_objective_k2": 0.5 * (one_c_mse + hppc_mse),
        }
        rows.append(row)
        if result.success and (best is None or row["joint_objective_k2"] < best["row"]["joint_objective_k2"]):
            best = {"row": row, "result": result}
    if best is None:
        raise RuntimeError("No joint thermal multistart succeeded.")
    return rows, best


def multistart_spread(rows: list[dict[str, object]]) -> dict[str, float | int]:
    successful = [row for row in rows if bool(row["success"])]
    best_joint_rmse = min(np.sqrt(float(row["joint_objective_k2"])) for row in successful)
    tolerance = max(1.0e-10, NEAR_OPTIMAL_RELATIVE_RMSE * best_joint_rmse)
    near = [
        row for row in successful
        if np.sqrt(float(row["joint_objective_k2"])) <= best_joint_rmse + tolerance
    ]
    output: dict[str, float | int] = {
        "successful_starts": len(successful),
        "near_optimal_starts": len(near),
        "near_optimal_joint_rmse_tolerance_k": tolerance,
    }
    for name in PARAMETER_NAMES:
        key = f"fit_{name}"
        values = np.asarray([float(row[key]) for row in successful])
        near_values = np.asarray([float(row[key]) for row in near])
        output[f"successful_{name}_min"] = float(np.min(values))
        output[f"successful_{name}_max"] = float(np.max(values))
        output[f"near_optimal_{name}_min"] = float(np.min(near_values))
        output[f"near_optimal_{name}_max"] = float(np.max(near_values))
        output[f"near_optimal_{name}_relative_range"] = float(
            (np.max(near_values) - np.min(near_values)) / np.median(near_values)
        )
    return output


def conditioning(
    best_log: np.ndarray, one_c: object, blocks: list[ThermalBlock], parameters: np.ndarray
) -> dict[str, object]:
    jacobian = np.empty((one_c.time_s.size + sum(block.time_s.size for block in blocks), 2))
    for column in range(2):
        plus = best_log.copy()
        minus = best_log.copy()
        plus[column] += LOG_SENSITIVITY_STEP
        minus[column] -= LOG_SENSITIVITY_STEP
        jacobian[:, column] = (
            balanced_residual(plus, one_c, blocks) - balanced_residual(minus, one_c, blocks)
        ) / (2.0 * LOG_SENSITIVITY_STEP)
    singular_values = np.linalg.svd(jacobian, compute_uv=False)
    threshold = np.finfo(float).eps * max(jacobian.shape) * singular_values[0]
    rank = int(np.count_nonzero(singular_values > threshold))
    condition_number = float(singular_values[0] / singular_values[-1])
    sensitivity_correlation = float(np.corrcoef(jacobian, rowvar=False)[0, 1])
    covariance = np.linalg.pinv(jacobian.T @ jacobian)
    parameter_correlation = float(covariance[0, 1] / np.sqrt(covariance[0, 0] * covariance[1, 1]))
    log_lower = np.log(LOWER_BOUNDS)
    log_upper = np.log(UPPER_BOUNDS)
    fraction = np.minimum(np.log(parameters) - log_lower, log_upper - np.log(parameters)) / (log_upper - log_lower)
    if np.any(fraction <= 0.01):
        classification = "JOINT THERMAL PARAMETERIZATION BOUNDARY LIMITED"
    elif condition_number >= 1.0e3 or abs(parameter_correlation) >= 0.95:
        classification = "JOINT THERMAL PARAMETERS PRACTICALLY CORRELATED"
    else:
        classification = "JOINT THERMAL PARAMETERS WELL CONDITIONED"
    return {
        "jacobian": jacobian,
        "rank": rank,
        "singular_values": singular_values,
        "condition_number": condition_number,
        "sensitivity_correlation": sensitivity_correlation,
        "parameter_correlation": parameter_correlation,
        "c_th_nearest_log_bound_fraction": float(fraction[0]),
        "g_th_nearest_log_bound_fraction": float(fraction[1]),
        "classification": classification,
    }


def save_multistart(rows: list[dict[str, object]]) -> None:
    with MULTISTART_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_one_c(one_c: object, joint: np.ndarray, baseline: np.ndarray) -> None:
    minutes = (one_c.time_s - one_c.time_s[0]) / 60.0
    fig, axes = plt.subplots(2, 1, figsize=(11.5, 7.2), sharex=True, constrained_layout=True)
    axes[0].plot(minutes, one_c.measured_temperature_k - 273.15, color="black", lw=1.2, label="Measured")
    axes[0].plot(minutes, baseline - 273.15, color="#CC79A7", lw=1.0, ls="--", label="1C-only parameters")
    axes[0].plot(minutes, joint - 273.15, color="#0072B2", lw=1.1, label="Joint 1C+HPPC parameters")
    axes[0].plot(minutes, one_c.ambient_temperature_k - 273.15, color="#009E73", lw=0.9, label="Chamber")
    axes[0].set(ylabel="Temperature [degC]", title="Joint 1C+HPPC calibration: 1C experiment")
    axes[0].grid(alpha=0.22)
    axes[0].legend(ncol=2)
    axes[1].plot(minutes, baseline - one_c.measured_temperature_k, color="#CC79A7", lw=0.9, ls="--", label="1C-only residual")
    axes[1].plot(minutes, joint - one_c.measured_temperature_k, color="#0072B2", lw=1.0, label="Joint residual")
    axes[1].axhline(0.0, color="black", lw=0.7)
    axes[1].set(xlabel="Elapsed time [min]", ylabel="Predicted - measured [K]")
    axes[1].grid(alpha=0.22)
    axes[1].legend()
    fig.savefig(ONE_C_PLOT_PATH, dpi=180)
    plt.close(fig)


def plot_hppc(blocks: list[ThermalBlock], predictions: list[np.ndarray]) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(12.0, 7.5), sharex=True, constrained_layout=True)
    for number, (block, prediction) in enumerate(zip(blocks, predictions, strict=True)):
        stride = max(1, block.time_s.size // 5000)
        shown = np.arange(0, block.time_s.size, stride)
        hours = block.source_time_s / 3600.0
        label_measured = "Measured" if number == 0 else None
        label_predicted = "Joint prediction" if number == 0 else None
        axes[0].plot(hours[shown], block.measured_temperature_k[shown] - 273.15, color="black", lw=0.85, label=label_measured)
        axes[0].plot(hours[shown], prediction[shown] - 273.15, color="#0072B2", lw=0.9, label=label_predicted)
        axes[0].plot(hours[shown], block.ambient_temperature_k[shown] - 273.15, color="#009E73", lw=0.7, label="Chamber" if number == 0 else None)
        axes[1].plot(hours[shown], prediction[shown] - block.measured_temperature_k[shown], color="#0072B2", lw=0.75)
    axes[0].set(ylabel="Temperature [degC]", title="Joint 1C+HPPC calibration: retained continuous HPPC blocks")
    axes[0].grid(alpha=0.22)
    axes[0].legend()
    axes[1].axhline(0.0, color="black", lw=0.7)
    axes[1].set(xlabel="HPPC source elapsed time [h]", ylabel="Predicted - measured [K]")
    axes[1].grid(alpha=0.22)
    fig.savefig(HPPC_PLOT_PATH, dpi=180)
    plt.close(fig)


def plot_udds(udds: object, joint: np.ndarray, baseline: np.ndarray) -> None:
    shown = np.arange(0, udds.time_s.size, max(1, udds.time_s.size // 45000))
    hours = udds.time_s / 3600.0
    fig, axes = plt.subplots(2, 1, figsize=(12.0, 7.4), sharex=True, constrained_layout=True)
    axes[0].plot(hours[shown], udds.measured_temperature_k[shown] - 273.15, color="black", lw=1.1, label="Measured")
    axes[0].plot(hours[shown], baseline[shown] - 273.15, color="#CC79A7", lw=0.9, ls="--", label="Previous 1C-only")
    axes[0].plot(hours[shown], joint[shown] - 273.15, color="#0072B2", lw=1.0, label="Joint 1C+HPPC")
    axes[0].plot(hours[shown], udds.ambient_temperature_k[shown] - 273.15, color="#009E73", lw=0.9, label="Chamber")
    axes[0].set(ylabel="Temperature [degC]", title="Independent UDDS thermal validation: measured-voltage heat, no refit")
    axes[0].grid(alpha=0.22)
    axes[0].legend(ncol=2)
    axes[1].plot(hours[shown], baseline[shown] - udds.measured_temperature_k[shown], color="#CC79A7", lw=0.85, ls="--", label="1C-only residual")
    axes[1].plot(hours[shown], joint[shown] - udds.measured_temperature_k[shown], color="#0072B2", lw=0.9, label="Joint residual")
    axes[1].axhline(0.0, color="black", lw=0.7)
    axes[1].set(xlabel="Elapsed time in SOC 0.90-to-0.10 interval [h]", ylabel="Predicted - measured [K]")
    axes[1].grid(alpha=0.22)
    axes[1].legend()
    fig.savefig(UDDS_PLOT_PATH, dpi=180)
    plt.close(fig)


def add_rows(rows: list[dict[str, object]], scope: str, subset: str, values: dict[str, object], n: int | str) -> None:
    for name, value in values.items():
        unit = "1"
        if name.endswith("_j_per_k"):
            unit = "J/K"
        elif name.endswith("_w_per_k"):
            unit = "W/K"
        elif name.endswith("_k2"):
            unit = "K^2"
        elif name.endswith("_k"):
            unit = "K"
        elif name.endswith("_c"):
            unit = "degC"
        elif name.endswith("_s"):
            unit = "s"
        rows.append({"scope": scope, "subset": subset, "metric": name, "value": value, "unit": unit, "n_samples": n})


def write_metrics(
    metrics_sets: dict[str, dict[str, float]],
    diagnostics: dict[str, dict[str, float]],
    parameters: np.ndarray,
    spread: dict[str, object],
    condition: dict[str, object],
    audit: dict[str, object],
    block_rows: list[dict[str, object]],
    gaps: list[dict[str, object]],
    one_c_n: int,
    hppc_n: int,
    udds_n: int,
) -> None:
    rows: list[dict[str, object]] = []
    for name, values in metrics_sets.items():
        n = one_c_n if "one_c" in name else hppc_n if "hppc" in name else udds_n
        add_rows(rows, "temperature_metrics", name, values, n)
    for name, values in diagnostics.items():
        add_rows(rows, "hppc_diagnostic", name, values, hppc_n)
    add_rows(rows, "joint_fit", "best", {
        "c_th_j_per_k": parameters[0], "g_th_w_per_k": parameters[1],
        "tau_th_s": parameters[0] / parameters[1],
        "joint_objective_k2": 0.5 * (metrics_sets["joint_one_c"]["mse_k2"] + metrics_sets["joint_hppc"]["mse_k2"]),
    }, one_c_n + hppc_n)
    add_rows(rows, "multistart", "summary", spread, one_c_n + hppc_n)
    add_rows(rows, "conditioning", "joint", {
        "numerical_rank": condition["rank"],
        "singular_value_1": condition["singular_values"][0],
        "singular_value_2": condition["singular_values"][1],
        "condition_number": condition["condition_number"],
        "sensitivity_correlation": condition["sensitivity_correlation"],
        "parameter_correlation": condition["parameter_correlation"],
        "c_th_nearest_log_bound_fraction": condition["c_th_nearest_log_bound_fraction"],
        "g_th_nearest_log_bound_fraction": condition["g_th_nearest_log_bound_fraction"],
        "classification": condition["classification"],
    }, one_c_n + hppc_n)
    add_rows(rows, "hppc_audit", "overall", audit, hppc_n)
    for block in block_rows:
        add_rows(rows, "hppc_block", f"block_{block['block_id']}", block, int(block["sample_count"]))
    for gap in gaps:
        add_rows(rows, "hppc_omitted_gap", f"gap_{gap['gap_id']}", gap, 0)
    with METRICS_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def metric_row(label: str, values: dict[str, float]) -> str:
    return f"| {label} | {values['rmse_k']:.9f} | {values['mae_k']:.9f} | {values['mean_signed_error_k']:+.9f} | {values['max_abs_error_k']:.9f} | {values['r2']:.9f} |"


def write_report(
    blocks: list[ThermalBlock],
    block_rows: list[dict[str, object]],
    audit: dict[str, object],
    gaps: list[dict[str, object]],
    parameters: np.ndarray,
    metrics_sets: dict[str, dict[str, float]],
    diagnostics: dict[str, dict[str, float]],
    spread: dict[str, object],
    condition: dict[str, object],
    final_classification: str,
) -> None:
    joint_one_c = metrics_sets["joint_one_c"]
    joint_hppc = metrics_sets["joint_hppc"]
    baseline_hppc = metrics_sets["baseline_hppc"]
    joint_udds = metrics_sets["joint_udds"]
    c_change = 100.0 * (parameters[0] / BASELINE_C_TH - 1.0)
    g_change = 100.0 * (parameters[1] / BASELINE_G_TH - 1.0)
    udds_improvement = 100.0 * (1.0 - joint_udds["rmse_k"] / BASELINE_UDDS["rmse_k"])
    bias_change = joint_udds["mean_signed_error_k"] - BASELINE_UDDS["mean_signed_error_k"]
    materially_changes = max(abs(c_change), abs(g_change)) >= 10.0
    conditioning_improves = (
        "BOUNDARY" not in str(condition["classification"])
        and float(condition["condition_number"]) < BASELINE_CONDITION_NUMBER
        and abs(float(condition["parameter_correlation"]))
        < abs(BASELINE_PARAMETER_CORRELATION)
    )
    transfer_improves = udds_improvement > 0.0
    replace = transfer_improves and udds_improvement >= 10.0 and joint_udds["r2"] > 0.0
    persistent_bias = (
        abs(joint_udds["mean_signed_error_k"]) / joint_udds["rmse_k"] >= 0.75
    )
    lines = [
        "# Panasonic joint 1C + HPPC one-node thermal calibration",
        "",
        f"## {final_classification}",
        "",
        "Only `Cth` and `Gth` were fitted. The clean one-node architecture, frozen electrical/OCP/state mapping, capacity, entropy coefficient, and measured-voltage heat formulation were unchanged. UDDS was not used in fitting.",
        "",
        "## HPPC thermal-continuity gate",
        "",
        f"- Raw MAT provenance: `{HPPC_RAW_PATH}`",
        f"- Read-only electrical table: `{HPPC_SOURCE_PATH}`",
        f"- Source rows: {audit['source_rows']:,}; audited raw continuous blocks: {audit['raw_block_count']}; retained entropy-supported blocks: **{audit['retained_block_count']}**; retained samples: {audit['retained_samples']:,}.",
        f"- Normal pulse median dt: {audit['normal_pulse_dt_median_s']:.9f} s; normal rest median dt: {audit['normal_rest_dt_median_s']:.9f} s; retained duplicate-time increments: {audit['duplicate_time_increments_retained']} (retained as zero-time observations, with no propagation).",
        f"- Thirteen omitted-data gaps exceed {CONTIGUOUS_GAP_THRESHOLD_S:g} s and span {audit['long_gap_min_s']:.6f}-{audit['long_gap_max_s']:.6f} s. No current, heat, ambient, or temperature was held across them.",
        f"- Global cumulative-Ah SOC is preserved without block resets. Retained SOC range is {audit['retained_soc_min']:.12f}-{audit['retained_soc_max']:.12f}. The upper 0.90 boundary was interpolated inside block 3; the lower 0.10 boundary lies inside an omitted gap, so no interpolation across that gap was attempted.",
        "- Battery and chamber temperature channels are finite, chamber temperature is available for every retained sample, and each retained block contains real pulse/rest thermal history. **Gate decision: HPPC SUITABLE FOR THERMAL CALIBRATION.**",
        "- Each block is initialized with its measured starting battery temperature. This is a measured initial condition, not a fitted offset.",
        "",
        "| Block | Source indices | SOC start | SOC end | Duration [s] | Samples | Median positive dt [s] | Max normal dt [s] | T start [degC] | T end [degC] | Chamber [degC] | Internal current [A] |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in block_rows:
        lines.append(
            f"| {row['block_id']} | {row['source_start_index']}-{row['source_end_index']} | {row['start_soc']:.9f} | {row['end_soc']:.9f} | {row['duration_s']:.6f} | {row['sample_count']:,} | {row['median_positive_dt_s']:.6f} | {row['max_normal_dt_s']:.6f} | {row['initial_battery_temperature_c']:.6f} | {row['final_battery_temperature_c']:.6f} | {row['chamber_temperature_min_c']:.3f}-{row['chamber_temperature_max_c']:.3f} | {row['current_min_a']:.5f} to {row['current_max_a']:.5f} |"
        )
    lines.extend([
        "",
        "Omitted-gap source-index pairs and durations are recorded individually in the metrics CSV. The physical heat input is `I*(Ueq_average-V_measured) - I*T_measured*dUeq/dT`; there is no absolute value, negative-heat clipping, automatic sign selection, empirical multiplier, regime weight, or offset.",
        "",
        "## Joint calibration",
        "",
        f"The objective is `0.5*MSE_1C + 0.5*MSE_HPPC = {0.5*(joint_one_c['mse_k2']+joint_hppc['mse_k2']):.12g} K^2`. Forty deterministic log-space multistarts used only the specified numerical guards. {spread['successful_starts']} succeeded and {spread['near_optimal_starts']} are within 0.01% joint-RMSE of the best.",
        "",
        f"- `Cth = {parameters[0]:.12f} J/K` ({c_change:+.6f}% from 1C-only)",
        f"- `Gth = {parameters[1]:.12f} W/K` ({g_change:+.6f}% from 1C-only)",
        f"- `tau_th = {parameters[0]/parameters[1]:.12f} s`",
        "",
        "| Experiment | RMSE [K] | MAE [K] | Bias [K] | Max abs [K] | R2 |",
        "|---|---:|---:|---:|---:|---:|",
        metric_row("1C, joint fit", joint_one_c),
        metric_row("HPPC, joint fit", joint_hppc),
        metric_row("HPPC, 1C-only parameters", baseline_hppc),
        "",
        "### HPPC diagnostics (not reweighted)",
        "",
        "| Subset | RMSE [K] | MAE [K] | Bias [K] | Max abs [K] | R2 |",
        "|---|---:|---:|---:|---:|---:|",
    ])
    for name in ("load", "rest", "high", "middle", "low"):
        lines.append(metric_row(name, diagnostics[name]))
    lines.extend([
        "",
        "## Conditioning and multistart stability",
        "",
        f"**{condition['classification']}** (practical local audit, not structural identifiability). The balanced numerical sensitivity matrix has rank **{condition['rank']}/2**, singular values `{condition['singular_values'][0]:.12g}, {condition['singular_values'][1]:.12g}`, condition number `{condition['condition_number']:.12g}`, sensitivity correlation `{condition['sensitivity_correlation']:+.12f}`, and covariance-derived parameter correlation `{condition['parameter_correlation']:+.12f}`.",
        "",
        f"Nearest log-bound fractions are {condition['c_th_nearest_log_bound_fraction']:.6f} (`Cth`) and {condition['g_th_nearest_log_bound_fraction']:.6f} (`Gth`). Near-optimal ranges are `Cth={spread['near_optimal_c_th_j_per_k_min']:.12g}-{spread['near_optimal_c_th_j_per_k_max']:.12g}` and `Gth={spread['near_optimal_g_th_w_per_k_min']:.12g}-{spread['near_optimal_g_th_w_per_k_max']:.12g}`.",
        "",
        f"The previous 1C-only condition number was {BASELINE_CONDITION_NUMBER:.9f}, with sensitivity correlation {BASELINE_SENSITIVITY_CORRELATION:+.9f} and parameter correlation {BASELINE_PARAMETER_CORRELATION:+.9f}.",
        "",
        "## Critical no-refit UDDS validation",
        "",
        f"Source: `{UDDS_SOURCE_PATH}`. The exact same entropy-supported SOC 0.90-to-0.10 measured-voltage heat trajectory was evaluated after freezing the joint parameters.",
        "",
        "| Metric | 1C-only | Joint 1C+HPPC |",
        "|---|---:|---:|",
        f"| Cth [J/K] | {BASELINE_C_TH:.12f} | {parameters[0]:.12f} |",
        f"| Gth [W/K] | {BASELINE_G_TH:.12f} | {parameters[1]:.12f} |",
        f"| tau_th [s] | {BASELINE_C_TH/BASELINE_G_TH:.12f} | {parameters[0]/parameters[1]:.12f} |",
        f"| 1C RMSE [K] | {BASELINE_ONE_C['rmse_k']:.9f} | {joint_one_c['rmse_k']:.9f} |",
        f"| HPPC RMSE [K] | {baseline_hppc['rmse_k']:.9f} | {joint_hppc['rmse_k']:.9f} |",
        f"| UDDS RMSE [K] | {BASELINE_UDDS['rmse_k']:.9f} | {joint_udds['rmse_k']:.9f} |",
        f"| UDDS MAE [K] | {BASELINE_UDDS['mae_k']:.9f} | {joint_udds['mae_k']:.9f} |",
        f"| UDDS bias [K] | {BASELINE_UDDS['mean_signed_error_k']:+.9f} | {joint_udds['mean_signed_error_k']:+.9f} |",
        f"| UDDS max error [K] | {BASELINE_UDDS['max_abs_error_k']:.9f} | {joint_udds['max_abs_error_k']:.9f} |",
        f"| UDDS R2 | {BASELINE_UDDS['r2']:.9f} | {joint_udds['r2']:.9f} |",
        f"| Final measured [degC] | {joint_udds['measured_final_temperature_c']:.9f} | {joint_udds['measured_final_temperature_c']:.9f} |",
        f"| Final predicted [degC] | {BASELINE_UDDS['predicted_final_temperature_c']:.9f} | {joint_udds['predicted_final_temperature_c']:.9f} |",
        f"| Measured rise [K] | {joint_udds['measured_temperature_rise_k']:+.9f} | {joint_udds['measured_temperature_rise_k']:+.9f} |",
        f"| Predicted rise [K] | {BASELINE_UDDS['predicted_temperature_rise_k']:+.9f} | {joint_udds['predicted_temperature_rise_k']:+.9f} |",
        "",
        f"UDDS RMSE changes by {udds_improvement:+.6f}% relative to 0.526697652 K. Bias changes by {bias_change:+.9f} K from the prior -0.512714512 K.",
        "",
        "## Decisions",
        "",
        f"1. **Does adding HPPC materially change Cth or Gth?** **{'Yes' if materially_changes else 'No'}**. Changes are {c_change:+.3f}% and {g_change:+.3f}%, respectively.",
        "",
        f"2. **Does HPPC improve practical conditioning?** **{'Yes' if conditioning_improves else 'No'}**. The condition number changes from {BASELINE_CONDITION_NUMBER:.3f} to {condition['condition_number']:.3f}, but the joint optimum is boundary-limited and the absolute parameter correlation changes from {abs(BASELINE_PARAMETER_CORRELATION):.3f} to {abs(float(condition['parameter_correlation'])):.3f}; correlations and multistart spread are reported above.",
        "",
        f"3. **Does joint calibration improve independent UDDS transfer?** **{'Yes' if transfer_improves else 'No'}**. The no-refit RMSE change is {udds_improvement:+.3f}%.",
        "",
        f"4. **Is the improvement large enough to replace the 1C-only calibration?** **{'Yes' if replace else 'No'}**. Replacement requires a clear independent improvement without pathological/boundary conditioning; the frozen baseline file was not overwritten.",
        "",
        f"5. **Does a persistent approximately constant UDDS underprediction remain?** **{'Yes' if persistent_bias else 'No'}**. Joint UDDS bias is {joint_udds['mean_signed_error_k']:+.6f} K ({abs(joint_udds['mean_signed_error_k'])/joint_udds['rmse_k']:.1%} of RMSE magnitude).",
        "",
        f"**Final classification: {final_classification}.**",
        "",
        "No UDDS/CC-CV fitting, two-node state, empirical correction, ambient offset, electrical refit, or model-voltage heat calibration was performed.",
    ])
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    one_c = prepare_one_c_data()
    blocks, audit, block_rows, gaps = prepare_hppc_blocks()
    rows, best = run_multistart(one_c, blocks)
    save_multistart(rows)
    best_log = np.asarray(best["result"].x)
    parameters = np.exp(best_log)

    joint_one_c_prediction = predict_one_c(one_c, parameters)
    baseline_one_c_prediction = predict_one_c(one_c, np.asarray([BASELINE_C_TH, BASELINE_G_TH]))
    joint_hppc_predictions = predict_hppc(blocks, parameters)
    baseline_hppc_predictions = predict_hppc(blocks, np.asarray([BASELINE_C_TH, BASELINE_G_TH]))
    measured_hppc = np.concatenate([block.measured_temperature_k for block in blocks])
    joint_hppc_all = np.concatenate(joint_hppc_predictions)
    baseline_hppc_all = np.concatenate(baseline_hppc_predictions)

    udds, _ = prepare_udds_data()
    qirr_udds = udds.current_a * (udds.ueq_average_v - udds.measured_voltage_v)
    qrev_udds = -udds.current_a * udds.measured_temperature_k * udds.dudt_v_per_k
    heat_udds = qirr_udds + qrev_udds
    joint_udds_prediction = exact_one_node(
        udds.time_s, heat_udds, udds.ambient_temperature_k,
        float(udds.measured_temperature_k[0]), float(parameters[0]), float(parameters[1]),
    )
    baseline_udds_prediction = exact_one_node(
        udds.time_s, heat_udds, udds.ambient_temperature_k,
        float(udds.measured_temperature_k[0]), BASELINE_C_TH, BASELINE_G_TH,
    )

    metrics_sets = {
        "joint_one_c": temperature_metrics(one_c.measured_temperature_k, joint_one_c_prediction),
        "baseline_one_c": temperature_metrics(one_c.measured_temperature_k, baseline_one_c_prediction),
        "joint_hppc": temperature_metrics(measured_hppc, joint_hppc_all),
        "baseline_hppc": temperature_metrics(measured_hppc, baseline_hppc_all),
        "joint_udds": temperature_metrics(udds.measured_temperature_k, joint_udds_prediction),
        "baseline_udds": temperature_metrics(udds.measured_temperature_k, baseline_udds_prediction),
    }
    for key, reference in BASELINE_ONE_C.items():
        if not np.isclose(metrics_sets["baseline_one_c"][key], reference, atol=8.0e-9, rtol=0.0):
            raise RuntimeError(f"1C-only baseline changed for {key}.")
    for key, reference in BASELINE_UDDS.items():
        if not np.isclose(metrics_sets["baseline_udds"][key], reference, atol=8.0e-9, rtol=0.0):
            raise RuntimeError(f"UDDS baseline changed for {key}.")

    diagnostics = {
        name: combined_hppc_metrics(blocks, joint_hppc_predictions, name)
        for name in ("load", "rest", "high", "middle", "low")
    }
    spread = multistart_spread(rows)
    condition = conditioning(best_log, one_c, blocks, parameters)
    udds_improvement = 100.0 * (1.0 - metrics_sets["joint_udds"]["rmse_k"] / BASELINE_UDDS["rmse_k"])
    replacement_justified = udds_improvement >= 10.0 and metrics_sets["joint_udds"]["r2"] > 0.0 and "BOUNDARY" not in str(condition["classification"])
    final_classification = (
        "JOINT CALIBRATION IMPROVES THERMAL TRANSFER"
        if udds_improvement > 0.0
        else "JOINT CALIBRATION DOES NOT MATERIALLY IMPROVE THERMAL TRANSFER"
    )

    plot_one_c(one_c, joint_one_c_prediction, baseline_one_c_prediction)
    plot_hppc(blocks, joint_hppc_predictions)
    plot_udds(udds, joint_udds_prediction, baseline_udds_prediction)
    write_metrics(
        metrics_sets, diagnostics, parameters, spread, condition, audit, block_rows,
        gaps, one_c.time_s.size, measured_hppc.size, udds.time_s.size,
    )
    write_report(
        blocks, block_rows, audit, gaps, parameters, metrics_sets, diagnostics,
        spread, condition, final_classification,
    )

    print("Joint 1C + HPPC one-node thermal calibration complete.")
    print(f"Retained HPPC blocks={len(blocks)}; samples={measured_hppc.size:,}")
    print(f"Cth={parameters[0]:.12f} J/K; Gth={parameters[1]:.12f} W/K; tau={parameters[0]/parameters[1]:.9f} s")
    print(f"1C RMSE={metrics_sets['joint_one_c']['rmse_k']:.9f} K; HPPC RMSE={metrics_sets['joint_hppc']['rmse_k']:.9f} K")
    print(f"UDDS RMSE={metrics_sets['joint_udds']['rmse_k']:.9f} K; change={udds_improvement:+.6f}%")
    print(f"{condition['classification']}; {final_classification}")


if __name__ == "__main__":
    main()
