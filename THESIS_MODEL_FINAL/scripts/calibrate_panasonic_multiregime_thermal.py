"""Common-parameter Panasonic one-node thermal calibration.

Fit only Cth and Gth to the accepted 1C, HPPC, and CC-CV preprocessing,
with equal experiment-level MSE weights. UDDS is strictly withheld and is
evaluated only after selecting the calibration optimum. This script does not
write any frozen model or parameter module.
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


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from models.ocp_panasonic_qi_frozen import U_n, U_p  # noqa: E402
from models.panasonic_entropic import dUeq_dT_panasonic  # noqa: E402
from models.panasonic_thermal_parameters import (  # noqa: E402
    C_TH_J_PER_K as BASELINE_C,
    G_TH_W_PER_K as BASELINE_G,
)
from scripts.calibrate_panasonic_1c_thermal import (  # noqa: E402
    prepare_data as prepare_one_c,
)
from scripts.calibrate_panasonic_joint_1c_hppc_thermal import (  # noqa: E402
    combined_hppc_metrics,
    exact_one_node,
    predict_hppc,
    predict_one_c,
    prepare_hppc_blocks,
    temperature_metrics,
)
from scripts.validate_panasonic_cccv_thermal import (  # noqa: E402
    prepare_data as prepare_cccv,
)
from scripts.validate_panasonic_udds_thermal import (  # noqa: E402
    prepare_data as prepare_udds,
    terminal_voltage_components,
)


DATA_DIR = ROOT / "thesis_data"
THESIS_FIG_DIR = ROOT / "thesis_figures" / "final"
PRESENTATION_FIG_DIR = ROOT / "presentation_assets" / "final"
REPORT_PATH = DATA_DIR / "panasonic_multiregime_thermal_calibration.md"
METRICS_PATH = DATA_DIR / "panasonic_multiregime_thermal_metrics.csv"
MULTISTART_PATH = DATA_DIR / "panasonic_multiregime_thermal_multistart.csv"

LOWER = np.asarray([5.0, 1.0e-4], dtype=float)
UPPER = np.asarray([1000.0, 1.0], dtype=float)
LOG_LOWER = np.log(LOWER)
LOG_UPPER = np.log(UPPER)
N_GRID_PER_PARAMETER = 8  # 8 x 8 = 64 deterministic starts
FD_STEP = 1.0e-5
NEAR_OPTIMAL_RELATIVE_OBJECTIVE = 1.0e-6

PREVIOUS_JOINT_C = 200.0
PREVIOUS_JOINT_G = 0.07845155487126856

ACCEPTED_BASELINE = {
    "one_c_rmse_k": 0.143622065213,
    "hppc_rmse_k": 0.617273179233,
    "udds_rmse_k": 0.526697652104,
    "udds_bias_k": -0.512714512143,
}
ACCEPTED_PREVIOUS_JOINT = {
    "one_c_rmse_k": 0.254140708884,
    "hppc_rmse_k": 0.512740397251,
    "udds_rmse_k": 0.378519120183,
    "udds_bias_k": -0.353251679028,
}


def primary_cccv_heat(data: object) -> np.ndarray:
    ueq = np.asarray(U_p(data.x_p_avg)) - np.asarray(U_n(data.x_n_avg))
    dudt = np.asarray(dUeq_dT_panasonic(data.soc))
    qirr = data.model_current_a * (ueq - data.measured_voltage_v)
    qrev = -data.model_current_a * data.measured_temperature_k * dudt
    heat = qirr + qrev
    if not np.all(np.isfinite(heat)):
        raise ValueError("Non-finite CC-CV primary heat.")
    return heat


def primary_udds_heat(data: object) -> np.ndarray:
    qirr = data.current_a * (data.ueq_average_v - data.measured_voltage_v)
    qrev = -data.current_a * data.measured_temperature_k * data.dudt_v_per_k
    heat = qirr + qrev
    if not np.all(np.isfinite(heat)):
        raise ValueError("Non-finite UDDS Validation-A heat.")
    return heat


def predict_cccv(data: object, heat: np.ndarray, parameters: np.ndarray) -> np.ndarray:
    return exact_one_node(
        data.time_s,
        heat,
        data.ambient_temperature_k,
        float(data.measured_temperature_k[0]),
        float(parameters[0]),
        float(parameters[1]),
    )


def predict_udds_a(data: object, heat: np.ndarray, parameters: np.ndarray) -> np.ndarray:
    return exact_one_node(
        data.time_s,
        heat,
        data.ambient_temperature_k,
        float(data.measured_temperature_k[0]),
        float(parameters[0]),
        float(parameters[1]),
    )


def predict_udds_fully_predictive(
    data: object, parameters: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Couple frozen voltage, predicted temperature, and heat with supplied C/G."""

    c_th, g_th = map(float, parameters)
    temperature = np.empty_like(data.time_s)
    voltage = np.empty_like(data.time_s)
    heat = np.empty_like(data.time_s)
    temperature[0] = data.measured_temperature_k[0]
    base_voltage, kinetic_per_k = terminal_voltage_components(
        data.current_a, data.x_n_surface, data.x_p_surface
    )
    for i in range(data.time_s.size):
        voltage[i] = base_voltage[i] + temperature[i] * kinetic_per_k[i]
        qirr = data.current_a[i] * (data.ueq_average_v[i] - voltage[i])
        qrev = -data.current_a[i] * temperature[i] * data.dudt_v_per_k[i]
        heat[i] = qirr + qrev
        if i + 1 == data.time_s.size:
            continue
        dt = data.time_s[i + 1] - data.time_s[i]
        x = g_th * dt / c_th
        approach = -np.expm1(-x)
        phi = approach / x if x != 0.0 else 1.0
        temperature[i + 1] = (
            temperature[i]
            + approach * (data.ambient_temperature_k[i] - temperature[i])
            + heat[i] * dt * phi / c_th
        )
    if not all(np.all(np.isfinite(v)) for v in (temperature, voltage, heat)):
        raise ValueError("Non-finite fully predictive UDDS result.")
    return temperature, voltage, heat


def common_residual(
    log_parameters: np.ndarray,
    one_c: object,
    hppc_blocks: list[object],
    cccv: object,
    cccv_heat: np.ndarray,
) -> np.ndarray:
    """Residual whose squared norm is exactly mean(MSE_1C, MSE_HPPC, MSE_CCCV)."""

    parameters = np.exp(log_parameters)
    r_one_c = predict_one_c(one_c, parameters) - one_c.measured_temperature_k
    hppc_predictions = predict_hppc(hppc_blocks, parameters)
    r_hppc = np.concatenate(
        [
            predicted - block.measured_temperature_k
            for block, predicted in zip(hppc_blocks, hppc_predictions, strict=True)
        ]
    )
    r_cccv = predict_cccv(cccv, cccv_heat, parameters) - cccv.measured_temperature_k
    return np.concatenate(
        (
            r_one_c / np.sqrt(3.0 * r_one_c.size),
            r_hppc / np.sqrt(3.0 * r_hppc.size),
            r_cccv / np.sqrt(3.0 * r_cccv.size),
        )
    )


def deterministic_starts() -> list[np.ndarray]:
    c_values = np.geomspace(LOWER[0] * 1.001, UPPER[0] / 1.001, N_GRID_PER_PARAMETER)
    g_values = np.geomspace(LOWER[1] * 1.001, UPPER[1] / 1.001, N_GRID_PER_PARAMETER)
    return [np.log([c_value, g_value]) for c_value in c_values for g_value in g_values]


def run_multistart(
    one_c: object,
    hppc_blocks: list[object],
    cccv: object,
    cccv_heat: np.ndarray,
) -> tuple[list[dict[str, object]], object]:
    rows: list[dict[str, object]] = []
    solutions: list[object] = []
    for start_id, initial in enumerate(deterministic_starts(), start=1):
        result = least_squares(
            common_residual,
            initial,
            bounds=(LOG_LOWER, LOG_UPPER),
            args=(one_c, hppc_blocks, cccv, cccv_heat),
            method="trf",
            jac="2-point",
            ftol=1.0e-12,
            xtol=1.0e-12,
            gtol=1.0e-12,
            max_nfev=1000,
        )
        parameters = np.exp(result.x)
        p_one_c = predict_one_c(one_c, parameters)
        p_hppc = predict_hppc(hppc_blocks, parameters)
        p_cccv = predict_cccv(cccv, cccv_heat, parameters)
        m_one_c = temperature_metrics(one_c.measured_temperature_k, p_one_c)
        m_hppc = combined_hppc_metrics(hppc_blocks, p_hppc)
        m_cccv = temperature_metrics(cccv.measured_temperature_k, p_cccv)
        objective = (m_one_c["mse_k2"] + m_hppc["mse_k2"] + m_cccv["mse_k2"]) / 3.0
        rows.append(
            {
                "start_id": start_id,
                "initial_c_th_j_per_k": float(np.exp(initial[0])),
                "initial_g_th_w_per_k": float(np.exp(initial[1])),
                "success": bool(result.success),
                "status": int(result.status),
                "nfev": int(result.nfev),
                "optimality": float(result.optimality),
                "c_th_j_per_k": float(parameters[0]),
                "g_th_w_per_k": float(parameters[1]),
                "tau_s": float(parameters[0] / parameters[1]),
                "objective_k2": float(objective),
                "one_c_rmse_k": m_one_c["rmse_k"],
                "hppc_rmse_k": m_hppc["rmse_k"],
                "cccv_rmse_k": m_cccv["rmse_k"],
                "message": str(result.message).replace("\n", " "),
            }
        )
        solutions.append(result)
    successful = [
        (row, result)
        for row, result in zip(rows, solutions, strict=True)
        if row["success"] and np.isfinite(row["objective_k2"])
    ]
    if not successful:
        raise RuntimeError("No deterministic multistart converged.")
    best_row, best_result = min(successful, key=lambda item: item[0]["objective_k2"])
    best_row["selected_best"] = True
    for row in rows:
        row.setdefault("selected_best", False)
    return rows, best_result


def conditioning(
    best_log: np.ndarray,
    one_c: object,
    hppc_blocks: list[object],
    cccv: object,
    cccv_heat: np.ndarray,
) -> dict[str, object]:
    base = common_residual(best_log, one_c, hppc_blocks, cccv, cccv_heat)
    jacobian = np.empty((base.size, 2))
    for column in range(2):
        plus = best_log.copy()
        minus = best_log.copy()
        plus[column] += FD_STEP
        minus[column] -= FD_STEP
        jacobian[:, column] = (
            common_residual(plus, one_c, hppc_blocks, cccv, cccv_heat)
            - common_residual(minus, one_c, hppc_blocks, cccv, cccv_heat)
        ) / (2.0 * FD_STEP)
    singular = np.linalg.svd(jacobian, compute_uv=False)
    tolerance = max(jacobian.shape) * np.finfo(float).eps * singular[0]
    rank = int(np.count_nonzero(singular > tolerance))
    condition_number = float(singular[0] / singular[-1]) if singular[-1] > 0.0 else float("inf")
    sensitivity_correlation = float(np.corrcoef(jacobian[:, 0], jacobian[:, 1])[0, 1])
    covariance = np.linalg.pinv(jacobian.T @ jacobian)
    covariance_correlation = float(
        covariance[0, 1] / np.sqrt(covariance[0, 0] * covariance[1, 1])
    )
    return {
        "rank": rank,
        "singular_values": singular,
        "condition_number": condition_number,
        "sensitivity_correlation": sensitivity_correlation,
        "covariance_parameter_correlation": covariance_correlation,
    }


def metric_subset(measured: np.ndarray, predicted: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    return temperature_metrics(measured[mask], predicted[mask])


def multistart_summary(rows: list[dict[str, object]]) -> dict[str, float | int]:
    successful = [row for row in rows if bool(row["success"])]
    best = min(float(row["objective_k2"]) for row in successful)
    near = [
        row
        for row in successful
        if float(row["objective_k2"]) <= best * (1.0 + NEAR_OPTIMAL_RELATIVE_OBJECTIVE)
    ]
    return {
        "attempted": len(rows),
        "successful": len(successful),
        "near_optimal": len(near),
        "all_c_min": min(float(row["c_th_j_per_k"]) for row in successful),
        "all_c_max": max(float(row["c_th_j_per_k"]) for row in successful),
        "all_g_min": min(float(row["g_th_w_per_k"]) for row in successful),
        "all_g_max": max(float(row["g_th_w_per_k"]) for row in successful),
        "near_c_min": min(float(row["c_th_j_per_k"]) for row in near),
        "near_c_max": max(float(row["c_th_j_per_k"]) for row in near),
        "near_g_min": min(float(row["g_th_w_per_k"]) for row in near),
        "near_g_max": max(float(row["g_th_w_per_k"]) for row in near),
    }


def bound_proximity(parameters: np.ndarray) -> dict[str, float | bool]:
    log_p = np.log(parameters)
    fractions = np.minimum(log_p - LOG_LOWER, LOG_UPPER - log_p) / (LOG_UPPER - LOG_LOWER)
    return {
        "c_log_fraction_to_nearest_bound": float(fractions[0]),
        "g_log_fraction_to_nearest_bound": float(fractions[1]),
        "c_within_one_percent_log_range": bool(fractions[0] <= 0.01),
        "g_within_one_percent_log_range": bool(fractions[1] <= 0.01),
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def append_metric_rows(
    rows: list[dict[str, object]], scope: str, subset: str, values: dict[str, object], n: int
) -> None:
    explicit_units = {
        "c_th_j_per_k": "J/K",
        "g_th_w_per_k": "W/K",
        "tau_s": "s",
        "objective_k2": "K2",
        "measured_final_temperature_c": "degC",
        "predicted_final_temperature_c": "degC",
    }
    for metric, value in values.items():
        if np.isscalar(value):
            unit = explicit_units.get(
                metric,
                "K" if metric.endswith("_k") else ("K2" if metric.endswith("_k2") else "1"),
            )
            rows.append(
                {
                    "scope": scope,
                    "subset": subset,
                    "metric": metric,
                    "value": value,
                    "unit": unit,
                    "n_samples": n,
                }
            )


def save_figure(fig: plt.Figure, stem: str) -> None:
    for directory, dpi in ((THESIS_FIG_DIR, 320), (PRESENTATION_FIG_DIR, 220)):
        directory.mkdir(parents=True, exist_ok=True)
        fig.savefig(directory / f"{stem}.png", dpi=dpi, bbox_inches="tight", facecolor="white")
        fig.savefig(directory / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def calibration_figure(
    one_c: object,
    p_one_c: np.ndarray,
    blocks: list[object],
    p_hppc: list[np.ndarray],
    cccv: object,
    p_cccv: np.ndarray,
) -> None:
    plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.2})
    fig, axes = plt.subplots(3, 1, figsize=(7.25, 8.0), constrained_layout=True)
    axes[0].plot((one_c.time_s - one_c.time_s[0]) / 60.0, one_c.measured_temperature_k - 273.15,
                 color="black", lw=1.25, label="Measured")
    axes[0].plot((one_c.time_s - one_c.time_s[0]) / 60.0, p_one_c - 273.15,
                 color="#0072B2", lw=1.25, ls="--", label="Common model")
    axes[0].set(title="(a) 1C discharge, SOC 0.90 to 0.10", ylabel="Battery temperature (°C)")
    axes[0].legend(frameon=False, ncol=2)

    t0 = min(float(block.source_time_s[0]) for block in blocks)
    for i, (block, prediction) in enumerate(zip(blocks, p_hppc, strict=True)):
        time_h = (block.source_time_s - t0) / 3600.0
        axes[1].plot(time_h, block.measured_temperature_k - 273.15, color="black", lw=0.65,
                     label="Measured" if i == 0 else None)
        axes[1].plot(time_h, prediction - 273.15, color="#0072B2", lw=0.75, ls="--",
                     label="Common model" if i == 0 else None)
    axes[1].set(title="(b) HPPC retained thermal blocks", ylabel="Battery temperature (°C)",
                xlabel="Source elapsed time (h)")

    axes[2].plot(cccv.time_s / 60.0, cccv.measured_temperature_k - 273.15,
                 color="black", lw=1.25, label="Measured")
    axes[2].plot(cccv.time_s / 60.0, p_cccv - 273.15,
                 color="#0072B2", lw=1.25, ls="--", label="Common model")
    axes[2].axvline(cccv.time_s[cccv.cv_start_retained_index] / 60.0,
                    color="0.55", lw=0.9, ls=":", label="CC-CV transition")
    axes[2].set(title="(c) CC-CV charge, SOC 0.10 to 0.90", ylabel="Battery temperature (°C)",
                xlabel="Elapsed time (min)")
    fig.suptitle("Common one-node thermal calibration across three operating regimes", fontsize=13)
    save_figure(fig, "multiregime_thermal_calibration_summary")


def udds_figure(
    udds: object,
    baseline_prediction: np.ndarray,
    common_prediction: np.ndarray,
) -> None:
    time_min = udds.time_s / 60.0
    measured_c = udds.measured_temperature_k - 273.15
    fig, axes = plt.subplots(2, 1, figsize=(7.25, 5.8), sharex=True,
                             gridspec_kw={"height_ratios": [2.0, 1.0]}, constrained_layout=True)
    axes[0].plot(time_min, measured_c, color="black", lw=1.0, label="Measured")
    axes[0].plot(time_min, baseline_prediction - 273.15, color="#CC79A7", lw=1.0, ls=":",
                 label="Accepted 1C-only pair")
    axes[0].plot(time_min, common_prediction - 273.15, color="#0072B2", lw=1.1, ls="--",
                 label="New common pair")
    axes[0].set_ylabel("Battery temperature (°C)")
    axes[0].legend(frameon=False, ncol=3, fontsize=8)
    axes[0].set_title("(a) Strict withheld UDDS Validation A")
    axes[1].plot(time_min, baseline_prediction - udds.measured_temperature_k,
                 color="#CC79A7", lw=0.8, ls=":", label="1C-only residual")
    axes[1].plot(time_min, common_prediction - udds.measured_temperature_k,
                 color="#0072B2", lw=0.9, label="Common-pair residual")
    axes[1].axhline(0.0, color="0.25", lw=0.7)
    axes[1].set(xlabel="Elapsed time (min)", ylabel="Predicted - measured (K)",
                title="(b) Validation residual")
    for axis in axes:
        axis.grid(alpha=0.2)
    fig.suptitle("Withheld UDDS thermal transfer comparison", fontsize=13)
    save_figure(fig, "multiregime_thermal_udds_validation")


def comparison_case(
    name: str,
    parameters: np.ndarray,
    one_c: object,
    blocks: list[object],
    cccv: object,
    cccv_heat: np.ndarray,
    udds: object,
    udds_heat: np.ndarray,
) -> tuple[dict[str, object], dict[str, np.ndarray | list[np.ndarray]]]:
    p1 = predict_one_c(one_c, parameters)
    ph = predict_hppc(blocks, parameters)
    pc = predict_cccv(cccv, cccv_heat, parameters)
    pu = predict_udds_a(udds, udds_heat, parameters)
    return (
        {
            "case": name,
            "c_th_j_per_k": float(parameters[0]),
            "g_th_w_per_k": float(parameters[1]),
            "tau_s": float(parameters[0] / parameters[1]),
            "one_c_rmse_k": temperature_metrics(one_c.measured_temperature_k, p1)["rmse_k"],
            "hppc_rmse_k": combined_hppc_metrics(blocks, ph)["rmse_k"],
            "cccv_a_rmse_k": temperature_metrics(cccv.measured_temperature_k, pc)["rmse_k"],
            "udds_a_rmse_k": temperature_metrics(udds.measured_temperature_k, pu)["rmse_k"],
            "udds_a_bias_k": temperature_metrics(udds.measured_temperature_k, pu)["mean_signed_error_k"],
        },
        {"one_c": p1, "hppc": ph, "cccv": pc, "udds": pu},
    )


def f(value: object, digits: int = 9) -> str:
    return f"{float(value):.{digits}g}"


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    one_c = prepare_one_c()
    blocks, hppc_audit, block_rows, gaps = prepare_hppc_blocks()
    cccv = prepare_cccv()
    udds, udds_audit = prepare_udds()
    cccv_heat = primary_cccv_heat(cccv)
    udds_heat = primary_udds_heat(udds)

    rows, best_result = run_multistart(one_c, blocks, cccv, cccv_heat)
    write_csv(MULTISTART_PATH, rows)
    best_parameters = np.exp(best_result.x)
    c_best, g_best = map(float, best_parameters)
    tau_best = c_best / g_best

    common_case, common_predictions = comparison_case(
        "New 1C+HPPC+CC-CV common pair", best_parameters,
        one_c, blocks, cccv, cccv_heat, udds, udds_heat,
    )
    baseline_case, baseline_predictions = comparison_case(
        "Accepted 1C-only pair", np.asarray([BASELINE_C, BASELINE_G]),
        one_c, blocks, cccv, cccv_heat, udds, udds_heat,
    )
    previous_case, _ = comparison_case(
        "Previous 1C+HPPC joint optimum", np.asarray([PREVIOUS_JOINT_C, PREVIOUS_JOINT_G]),
        one_c, blocks, cccv, cccv_heat, udds, udds_heat,
    )
    # Guard against silently changing the accepted old-case preprocessing/results.
    for calculated, reference, label in (
        (baseline_case["one_c_rmse_k"], ACCEPTED_BASELINE["one_c_rmse_k"], "baseline 1C"),
        (baseline_case["hppc_rmse_k"], ACCEPTED_BASELINE["hppc_rmse_k"], "baseline HPPC"),
        (baseline_case["udds_a_rmse_k"], ACCEPTED_BASELINE["udds_rmse_k"], "baseline UDDS"),
        (previous_case["one_c_rmse_k"], ACCEPTED_PREVIOUS_JOINT["one_c_rmse_k"], "joint 1C"),
        (previous_case["hppc_rmse_k"], ACCEPTED_PREVIOUS_JOINT["hppc_rmse_k"], "joint HPPC"),
        (previous_case["udds_a_rmse_k"], ACCEPTED_PREVIOUS_JOINT["udds_rmse_k"], "joint UDDS"),
    ):
        if not np.isclose(calculated, reference, rtol=0.0, atol=2.0e-9):
            raise RuntimeError(f"Accepted {label} result changed: {calculated} vs {reference}.")

    p_one_c = np.asarray(common_predictions["one_c"])
    p_hppc = common_predictions["hppc"]
    p_cccv = np.asarray(common_predictions["cccv"])
    p_udds = np.asarray(common_predictions["udds"])
    m_one_c = temperature_metrics(one_c.measured_temperature_k, p_one_c)
    m_hppc = combined_hppc_metrics(blocks, p_hppc)
    m_cccv = temperature_metrics(cccv.measured_temperature_k, p_cccv)
    m_cc = metric_subset(cccv.measured_temperature_k, p_cccv, cccv.cc_mask)
    m_cv = metric_subset(cccv.measured_temperature_k, p_cccv, cccv.cv_mask)
    m_udds = temperature_metrics(udds.measured_temperature_k, p_udds)
    fully_temperature, fully_voltage, fully_heat = predict_udds_fully_predictive(udds, best_parameters)
    m_udds_fully = temperature_metrics(udds.measured_temperature_k, fully_temperature)
    objective = (m_one_c["mse_k2"] + m_hppc["mse_k2"] + m_cccv["mse_k2"]) / 3.0
    condition = conditioning(best_result.x, one_c, blocks, cccv, cccv_heat)
    spread = multistart_summary(rows)
    proximity = bound_proximity(best_parameters)

    # Evidence-based classification; UDDS is evaluated, never optimized or selected upon.
    well_defined = (
        condition["rank"] == 2
        and np.isfinite(condition["condition_number"])
        and condition["condition_number"] < 100.0
        and not proximity["c_within_one_percent_log_range"]
        and not proximity["g_within_one_percent_log_range"]
        and spread["near_optimal"] == spread["successful"]
    )
    udds_improvement = 1.0 - m_udds["rmse_k"] / ACCEPTED_BASELINE["udds_rmse_k"]
    bias_improvement = 1.0 - abs(m_udds["mean_signed_error_k"]) / abs(ACCEPTED_BASELINE["udds_bias_k"])
    if not well_defined:
        classification = "COMMON THERMAL PARAMETERIZATION NOT WELL DETERMINED"
    elif udds_improvement >= 0.10 and bias_improvement > 0.0:
        classification = "MULTI-REGIME THERMAL CALIBRATION IMPROVES TRANSFER"
    else:
        classification = "MULTI-REGIME THERMAL CALIBRATION DOES NOT IMPROVE TRANSFER"

    metrics_rows: list[dict[str, object]] = []
    for subset, values, n in (
        ("1C", m_one_c, one_c.time_s.size),
        ("HPPC", m_hppc, sum(block.time_s.size for block in blocks)),
        ("CC-CV", m_cccv, cccv.time_s.size),
        ("CC", m_cc, int(np.count_nonzero(cccv.cc_mask))),
        ("CV", m_cv, int(np.count_nonzero(cccv.cv_mask))),
        ("UDDS Validation A", m_udds, udds.time_s.size),
        ("UDDS fully predictive", m_udds_fully, udds.time_s.size),
    ):
        append_metric_rows(metrics_rows, "temperature", subset, values, n)
    scalar_summary = {
        "c_th_j_per_k": c_best,
        "g_th_w_per_k": g_best,
        "tau_s": tau_best,
        "objective_k2": objective,
        "conditioning_rank": condition["rank"],
        "conditioning_singular_value_1": condition["singular_values"][0],
        "conditioning_singular_value_2": condition["singular_values"][1],
        "conditioning_condition_number": condition["condition_number"],
        "sensitivity_correlation": condition["sensitivity_correlation"],
        "covariance_parameter_correlation": condition["covariance_parameter_correlation"],
        **proximity,
        **{f"multistart_{key}": value for key, value in spread.items()},
        "udds_rmse_relative_improvement": udds_improvement,
        "udds_abs_bias_relative_improvement": bias_improvement,
    }
    append_metric_rows(metrics_rows, "calibration", "summary", scalar_summary,
                       one_c.time_s.size + sum(block.time_s.size for block in blocks) + cccv.time_s.size)
    for case in (baseline_case, previous_case, common_case):
        append_metric_rows(metrics_rows, "comparison", str(case["case"]),
                           {k: v for k, v in case.items() if k != "case"}, 0)
    write_csv(METRICS_PATH, metrics_rows)

    calibration_figure(one_c, p_one_c, blocks, p_hppc, cccv, p_cccv)
    udds_figure(udds, np.asarray(baseline_predictions["udds"]), p_udds)

    changed_c_pct = 100.0 * (c_best / BASELINE_C - 1.0)
    changed_g_pct = 100.0 * (g_best / BASELINE_G - 1.0)
    common_reasonable = (
        m_one_c["rmse_k"] < 0.30
        and m_cccv["rmse_k"] < 0.30
        and m_hppc["rmse_k"] < ACCEPTED_BASELINE["hppc_rmse_k"]
        and abs(m_hppc["mean_signed_error_k"]) < 0.50
    )
    electrical_limited = m_udds_fully["rmse_k"] > 1.5 * m_udds["rmse_k"]
    replace_evidence = well_defined and udds_improvement >= 0.10 and bias_improvement > 0.0 and common_reasonable

    def comparison_row(case: dict[str, object]) -> str:
        return (
            f"| {case['case']} | {float(case['c_th_j_per_k']):.9f} | "
            f"{float(case['g_th_w_per_k']):.9f} | {float(case['tau_s']):.6f} | "
            f"{float(case['one_c_rmse_k']):.9f} | {float(case['hppc_rmse_k']):.9f} | "
            f"{float(case['cccv_a_rmse_k']):.9f} | {float(case['udds_a_rmse_k']):.9f} | "
            f"{float(case['udds_a_bias_k']):.9f} |"
        )

    report = f"""# Panasonic common-parameter multiregime thermal calibration

## Scope and result

This was the single requested experiment. Only `Cth` and `Gth` were fitted; no model equation, OCP, entropy coefficient, electrical parameter, ambient offset, heat sign, or frozen parameter file was changed. UDDS remained fully withheld from calibration.

**Classification: {classification}**

## Objective and preprocessing

The minimized objective was exactly

`J = (MSE_1C + MSE_HPPC + MSE_CC-CV) / 3 = {objective:.12g} K^2`.

- 1C: accepted SOC 0.90 to 0.10 interval; {one_c.time_s.size} samples.
- HPPC: {len(blocks)} retained blocks and {sum(block.time_s.size for block in blocks)} retained samples; no propagation across {len(gaps)} omitted gaps; every block initialized at its measured starting temperature. The HPPC experiment contributes one retained-sample MSE, matching the preceding joint calibration convention.
- CC-CV: accepted SOC 0.10 to 0.90 interval; duplicate removal and existing CC/CV phase masks retained; {cccv.time_s.size} samples.
- Heat for all three calibration experiments: measured-voltage irreversible heat plus measured-temperature reversible heat.
- Bounds were numerical guards: `5 <= Cth <= 1000 J/K`, `1e-4 <= Gth <= 1 W/K`. Optimization was in log space with {spread['attempted']} deterministic starts.

## Best common parameters

| Quantity | Result |
|---|---:|
| Cth [J/K] | {c_best:.12g} |
| Gth [W/K] | {g_best:.12g} |
| tau = Cth/Gth [s] | {tau_best:.12g} |
| Objective J [K²] | {objective:.12g} |

Relative to the accepted 1C-only pair, Cth changed by {changed_c_pct:+.3f}% and Gth by {changed_g_pct:+.3f}%.

## Calibration performance

| Experiment/phase | RMSE [K] | MAE [K] | Bias (predicted-measured) [K] |
|---|---:|---:|---:|
| 1C | {m_one_c['rmse_k']:.9f} | {m_one_c['mae_k']:.9f} | {m_one_c['mean_signed_error_k']:.9f} |
| HPPC | {m_hppc['rmse_k']:.9f} | {m_hppc['mae_k']:.9f} | {m_hppc['mean_signed_error_k']:.9f} |
| CC-CV total | {m_cccv['rmse_k']:.9f} | {m_cccv['mae_k']:.9f} | {m_cccv['mean_signed_error_k']:.9f} |
| CC phase | {m_cc['rmse_k']:.9f} | {m_cc['mae_k']:.9f} | {m_cc['mean_signed_error_k']:.9f} |
| CV phase | {m_cv['rmse_k']:.9f} | {m_cv['mae_k']:.9f} | {m_cv['mean_signed_error_k']:.9f} |

## Bounds, multistart stability, and practical conditioning

- Log-range distance to nearest bound: Cth {proximity['c_log_fraction_to_nearest_bound']:.6f}; Gth {proximity['g_log_fraction_to_nearest_bound']:.6f}. Neither parameter is within 1% of its log-bound range: {not proximity['c_within_one_percent_log_range'] and not proximity['g_within_one_percent_log_range']}.
- Successful starts: {spread['successful']} / {spread['attempted']}.
- All successful Cth range: [{spread['all_c_min']:.12g}, {spread['all_c_max']:.12g}] J/K; Gth range: [{spread['all_g_min']:.12g}, {spread['all_g_max']:.12g}] W/K.
- Near-optimal (within {NEAR_OPTIMAL_RELATIVE_OBJECTIVE:.1e} relative objective) starts: {spread['near_optimal']}; Cth range [{spread['near_c_min']:.12g}, {spread['near_c_max']:.12g}], Gth range [{spread['near_g_min']:.12g}, {spread['near_g_max']:.12g}].
- Numerical rank: {condition['rank']}.
- Singular values of the equal-experiment residual Jacobian in log parameters: {condition['singular_values'][0]:.12g}, {condition['singular_values'][1]:.12g}.
- Condition number: {condition['condition_number']:.9g}.
- Sensitivity-column correlation: {condition['sensitivity_correlation']:.9g}.
- Covariance-derived log-parameter correlation: {condition['covariance_parameter_correlation']:.9g}.

These are practical local-conditioning diagnostics, not a structural-identifiability claim. The accepted 1C-only analysis had condition number 6.120213126, sensitivity correlation 0.304113503, and covariance-derived parameter correlation -0.866785425.

## Strict withheld UDDS validation

Validation A used the same audited SOC 0.90 to 0.10 interval, chamber temperature, initialization, measured-voltage irreversible heat, and measured-temperature reversible heat, with no refit.

| Metric | Accepted 1C-only | New common pair |
|---|---:|---:|
| RMSE [K] | {ACCEPTED_BASELINE['udds_rmse_k']:.9f} | {m_udds['rmse_k']:.9f} |
| MAE [K] | 0.512960900 | {m_udds['mae_k']:.9f} |
| Bias [K] | {ACCEPTED_BASELINE['udds_bias_k']:.9f} | {m_udds['mean_signed_error_k']:.9f} |
| Maximum absolute error [K] | 0.875239690 | {m_udds['max_abs_error_k']:.9f} |
| R² | -4.119410369 | {m_udds['r2']:.9f} |
| Final measured temperature [°C] | {m_udds['measured_final_temperature_c']:.9f} | {m_udds['measured_final_temperature_c']:.9f} |
| Final predicted temperature [°C] | 26.112601243 | {m_udds['predicted_final_temperature_c']:.9f} |
| Measured temperature rise [K] | {m_udds['measured_temperature_rise_k']:.9f} | {m_udds['measured_temperature_rise_k']:.9f} |
| Predicted temperature rise [K] | 0.279761481 | {m_udds['predicted_temperature_rise_k']:.9f} |

Relative to the accepted pair, the withheld UDDS RMSE {'improved' if udds_improvement >= 0.0 else 'worsened'} by {abs(100.0*udds_improvement):.3f}%, and the absolute bias {'improved' if bias_improvement >= 0.0 else 'worsened'} by {abs(100.0*bias_improvement):.3f}%.

## Secondary fully predictive UDDS diagnostic

| Metric | Result |
|---|---:|
| RMSE [K] | {m_udds_fully['rmse_k']:.9f} |
| Bias [K] | {m_udds_fully['mean_signed_error_k']:.9f} |
| Maximum absolute error [K] | {m_udds_fully['max_abs_error_k']:.9f} |

This path used frozen model voltage and recursively predicted temperature in reversible heat. It did not influence fitting. Its RMSE is {m_udds_fully['rmse_k']/m_udds['rmse_k']:.3f} times the primary Validation-A RMSE.

## Compact comparison

The old pairs were not refitted. Their 1C, HPPC, and UDDS entries reproduce the accepted results. The previously unreported old-pair CC-CV entries are forward-only evaluations using the already-audited CC-CV Validation-A preprocessing.

| Case | Cth [J/K] | Gth [W/K] | tau [s] | 1C RMSE [K] | HPPC RMSE [K] | CC-CV A RMSE [K] | UDDS A RMSE [K] | UDDS bias [K] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
{comparison_row(baseline_case)}
{comparison_row(previous_case)}
{comparison_row(common_case)}

## Required decisions

1. **Does multiregime calibration materially change Cth and/or Gth?** {'Yes' if abs(changed_c_pct) >= 10.0 or abs(changed_g_pct) >= 10.0 else 'No'}: changes are Cth {changed_c_pct:+.3f}% and Gth {changed_g_pct:+.3f}%.
2. **Can one common pair represent 1C, HPPC, and CC-CV reasonably?** {'Yes' if common_reasonable else 'No'}: 1C and CC-CV remain below 0.30 K RMSE, but HPPC has {m_hppc['rmse_k']:.9f} K RMSE and {m_hppc['mean_signed_error_k']:.9f} K bias, both worse than required for accepting a common representation. The individual results remain visible and are not hidden by the mean objective.
3. **Is the new pair numerically well-defined?** {'Yes' if well_defined else 'No'}: rank, condition, bounds, and multistart convergence are reported above. This is a practical statement only.
4. **Does withheld UDDS Validation A improve relative to 1C-only?** {'Yes' if udds_improvement > 0.0 else 'No'}; RMSE changes from {ACCEPTED_BASELINE['udds_rmse_k']:.9f} to {m_udds['rmse_k']:.9f} K.
5. **Does systematic negative UDDS bias decrease?** {'Yes' if abs(m_udds['mean_signed_error_k']) < abs(ACCEPTED_BASELINE['udds_bias_k']) else 'No'}; bias changes from {ACCEPTED_BASELINE['udds_bias_k']:.9f} to {m_udds['mean_signed_error_k']:.9f} K.
6. **Is the fully predictive path still primarily limited by electrical heat-generation error?** {'Yes' if electrical_limited else 'No clear evidence that it is the primary limitation'}: its RMSE is {m_udds_fully['rmse_k']/m_udds['rmse_k']:.3f} times Validation A, so {'the additional electrical-heat error dominates' if electrical_limited else 'electrical heat adds error, but the Validation-A thermal-transfer discrepancy remains dominant'}. No repair was attempted.
7. **Is there sufficient evidence to replace the frozen 1C-only pair?** {'Yes, as evidence for a later explicit model-freeze decision; this script did not replace it.' if replace_evidence else 'No.'}

The currently frozen thermal parameter file was not overwritten.
"""
    REPORT_PATH.write_text(report, encoding="utf-8")

    print(classification)
    print(f"Cth={c_best:.12g} J/K; Gth={g_best:.12g} W/K; tau={tau_best:.12g} s")
    print(f"J={objective:.12g} K^2; starts={spread['successful']}/{spread['attempted']}")
    print(f"RMSE: 1C={m_one_c['rmse_k']:.9f}, HPPC={m_hppc['rmse_k']:.9f}, CC-CV={m_cccv['rmse_k']:.9f} K")
    print(f"UDDS A RMSE={m_udds['rmse_k']:.9f} K; bias={m_udds['mean_signed_error_k']:.9f} K")
    print(f"UDDS fully predictive RMSE={m_udds_fully['rmse_k']:.9f} K")


if __name__ == "__main__":
    main()
