"""Boundary/profile audit for the Panasonic joint 1C + HPPC thermal fit.

This diagnostic changes only the numerical upper guard on Cth from 200 to
1000 J/K.  UDDS remains withheld and is evaluated only after each calibration
or fixed-Cth profile optimization has completed.
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

from scripts.calibrate_panasonic_joint_1c_hppc_thermal import (  # noqa: E402
    BASELINE_C_TH,
    BASELINE_G_TH,
    BASELINE_UDDS,
    LOG_SENSITIVITY_STEP,
    ThermalBlock,
    balanced_residual,
    exact_one_node,
    predict_hppc,
    predict_one_c,
    prepare_hppc_blocks,
    prepare_one_c_data,
    prepare_udds_data,
    temperature_metrics,
)


OUTPUT_DIRECTORY = REPOSITORY_ROOT / "thesis_data"
REPORT_PATH = OUTPUT_DIRECTORY / "panasonic_joint_boundary_audit.md"
METRICS_PATH = OUTPUT_DIRECTORY / "panasonic_joint_boundary_audit_metrics.csv"
PROFILE_PATH = OUTPUT_DIRECTORY / "panasonic_joint_Cth_profile.csv"
PROFILE_PLOT_PATH = OUTPUT_DIRECTORY / "panasonic_joint_Cth_profile.png"
PROFILE_UDDS_PLOT_PATH = OUTPUT_DIRECTORY / "panasonic_joint_Cth_profile_udds.png"

LOWER_BOUNDS = np.asarray([5.0, 1.0e-4])
UPPER_BOUNDS = np.asarray([1000.0, 1.0])
N_MULTISTARTS = 40
RANDOM_SEED = 20260915
PROFILE_C_TH = np.asarray([50.0, 75.0, 100.0, BASELINE_C_TH, 150.0, 200.0, 300.0, 500.0, 750.0, 1000.0])
PREVIOUS_BOUNDARY_PARAMETERS = np.asarray([200.0, 0.07845155487126856])
NEAR_OPTIMAL_RELATIVE_RMSE = 1.0e-4
FLAT_PROFILE_RELATIVE_OBJECTIVE = 0.01


def experiment_metrics(one_c: object, blocks: list[ThermalBlock], parameters: np.ndarray) -> tuple[dict[str, float], dict[str, float]]:
    prediction_one_c = predict_one_c(one_c, parameters)
    predictions_hppc = predict_hppc(blocks, parameters)
    measured_hppc = np.concatenate([block.measured_temperature_k for block in blocks])
    prediction_hppc = np.concatenate(predictions_hppc)
    return (
        temperature_metrics(one_c.measured_temperature_k, prediction_one_c),
        temperature_metrics(measured_hppc, prediction_hppc),
    )


def run_expanded_multistart(one_c: object, blocks: list[ThermalBlock]) -> tuple[list[dict[str, object]], dict[str, object]]:
    rng = np.random.default_rng(RANDOM_SEED)
    fixed_starts = np.asarray(
        [
            [BASELINE_C_TH, BASELINE_G_TH],
            PREVIOUS_BOUNDARY_PARAMETERS,
            [350.0, 0.06],
            [700.0, 0.04],
            [950.0, 0.02],
        ]
    )
    random_starts = np.exp(
        rng.uniform(
            np.log(LOWER_BOUNDS),
            np.log(UPPER_BOUNDS),
            size=(N_MULTISTARTS - fixed_starts.shape[0], 2),
        )
    )
    starts = np.vstack((fixed_starts, random_starts))
    rows: list[dict[str, object]] = []
    best: dict[str, object] | None = None
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
        metrics_one_c, metrics_hppc = experiment_metrics(one_c, blocks, parameters)
        objective = 0.5 * (metrics_one_c["mse_k2"] + metrics_hppc["mse_k2"])
        row = {
            "start_id": start_id,
            "success": bool(result.success),
            "status": int(result.status),
            "nfev": int(result.nfev),
            "optimality": float(result.optimality),
            "start_c_th_j_per_k": float(start[0]),
            "start_g_th_w_per_k": float(start[1]),
            "fit_c_th_j_per_k": float(parameters[0]),
            "fit_g_th_w_per_k": float(parameters[1]),
            "tau_th_s": float(parameters[0] / parameters[1]),
            "joint_objective_k2": objective,
            "rmse_1c_k": metrics_one_c["rmse_k"],
            "rmse_hppc_k": metrics_hppc["rmse_k"],
        }
        rows.append(row)
        if bool(result.success) and (best is None or objective < float(best["objective"])):
            best = {"result": result, "parameters": parameters, "objective": objective}
    if best is None:
        raise RuntimeError("No expanded-bound joint multistart succeeded.")
    return rows, best


def optimize_profile_g(one_c: object, blocks: list[ThermalBlock], c_th: float) -> tuple[np.ndarray, object]:
    log_lower = np.log(LOWER_BOUNDS[1:])
    log_upper = np.log(UPPER_BOUNDS[1:])
    starts = np.asarray([BASELINE_G_TH, PREVIOUS_BOUNDARY_PARAMETERS[1], 0.01, 0.3])
    best = None
    for start in starts:
        result = least_squares(
            lambda log_g: balanced_residual(
                np.log(np.asarray([c_th, float(np.exp(log_g[0]))])), one_c, blocks
            ),
            np.log([start]),
            bounds=(log_lower, log_upper),
            jac="3-point",
            x_scale="jac",
            ftol=1.0e-12,
            xtol=1.0e-12,
            gtol=1.0e-12,
            max_nfev=1000,
        )
        objective = float(result.cost)
        if result.success and (best is None or objective < best["objective"]):
            best = {"result": result, "objective": objective}
    if best is None:
        raise RuntimeError(f"No Gth profile optimization succeeded at Cth={c_th:g} J/K.")
    g_th = float(np.exp(best["result"].x[0]))
    return np.asarray([c_th, g_th]), best["result"]


def prepare_udds_heat() -> tuple[object, np.ndarray]:
    udds, _ = prepare_udds_data()
    q_irr = udds.current_a * (udds.ueq_average_v - udds.measured_voltage_v)
    q_rev = -udds.current_a * udds.measured_temperature_k * udds.dudt_v_per_k
    return udds, q_irr + q_rev


def evaluate_udds(udds: object, heat_w: np.ndarray, parameters: np.ndarray) -> dict[str, float]:
    prediction = exact_one_node(
        udds.time_s,
        heat_w,
        udds.ambient_temperature_k,
        float(udds.measured_temperature_k[0]),
        float(parameters[0]),
        float(parameters[1]),
    )
    return temperature_metrics(udds.measured_temperature_k, prediction)


def build_profile(one_c: object, blocks: list[ThermalBlock], udds: object, udds_heat: np.ndarray) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for c_th in PROFILE_C_TH:
        parameters, result = optimize_profile_g(one_c, blocks, float(c_th))
        one_c_metrics, hppc_metrics = experiment_metrics(one_c, blocks, parameters)
        udds_metrics = evaluate_udds(udds, udds_heat, parameters)
        rows.append(
            {
                "c_th_j_per_k": float(parameters[0]),
                "g_th_w_per_k": float(parameters[1]),
                "tau_th_s": float(parameters[0] / parameters[1]),
                "joint_objective_k2": 0.5 * (one_c_metrics["mse_k2"] + hppc_metrics["mse_k2"]),
                "rmse_1c_k": one_c_metrics["rmse_k"],
                "rmse_hppc_k": hppc_metrics["rmse_k"],
                "udds_rmse_k": udds_metrics["rmse_k"],
                "udds_bias_k": udds_metrics["mean_signed_error_k"],
                "udds_max_abs_error_k": udds_metrics["max_abs_error_k"],
                "g_fit_success": bool(result.success),
                "g_fit_nfev": int(result.nfev),
            }
        )
    return rows


def conditioning(log_parameters: np.ndarray, one_c: object, blocks: list[ThermalBlock]) -> dict[str, object]:
    row_count = one_c.time_s.size + sum(block.time_s.size for block in blocks)
    jacobian = np.empty((row_count, 2))
    for column in range(2):
        plus = log_parameters.copy()
        minus = log_parameters.copy()
        plus[column] += LOG_SENSITIVITY_STEP
        minus[column] -= LOG_SENSITIVITY_STEP
        jacobian[:, column] = (
            balanced_residual(plus, one_c, blocks)
            - balanced_residual(minus, one_c, blocks)
        ) / (2.0 * LOG_SENSITIVITY_STEP)
    singular_values = np.linalg.svd(jacobian, compute_uv=False)
    threshold = np.finfo(float).eps * max(jacobian.shape) * singular_values[0]
    covariance = np.linalg.pinv(jacobian.T @ jacobian)
    return {
        "rank": int(np.count_nonzero(singular_values > threshold)),
        "singular_value_1": float(singular_values[0]),
        "singular_value_2": float(singular_values[1]),
        "condition_number": float(singular_values[0] / singular_values[-1]),
        "sensitivity_correlation": float(np.corrcoef(jacobian, rowvar=False)[0, 1]),
        "parameter_correlation": float(
            covariance[0, 1] / np.sqrt(covariance[0, 0] * covariance[1, 1])
        ),
    }


def multistart_summary(rows: list[dict[str, object]]) -> dict[str, object]:
    successful = [row for row in rows if bool(row["success"])]
    joint_rmse = np.asarray([np.sqrt(float(row["joint_objective_k2"])) for row in successful])
    best_rmse = float(np.min(joint_rmse))
    tolerance = max(1.0e-10, NEAR_OPTIMAL_RELATIVE_RMSE * best_rmse)
    near = [row for row, rmse in zip(successful, joint_rmse, strict=True) if rmse <= best_rmse + tolerance]
    return {
        "successful_starts": len(successful),
        "near_optimal_starts": len(near),
        "near_optimal_joint_rmse_tolerance_k": tolerance,
        "successful_c_th_min": min(float(row["fit_c_th_j_per_k"]) for row in successful),
        "successful_c_th_max": max(float(row["fit_c_th_j_per_k"]) for row in successful),
        "successful_g_th_min": min(float(row["fit_g_th_w_per_k"]) for row in successful),
        "successful_g_th_max": max(float(row["fit_g_th_w_per_k"]) for row in successful),
        "near_optimal_c_th_min": min(float(row["fit_c_th_j_per_k"]) for row in near),
        "near_optimal_c_th_max": max(float(row["fit_c_th_j_per_k"]) for row in near),
        "near_optimal_g_th_min": min(float(row["fit_g_th_w_per_k"]) for row in near),
        "near_optimal_g_th_max": max(float(row["fit_g_th_w_per_k"]) for row in near),
    }


def parameter_bound_proximity(parameters: np.ndarray) -> dict[str, float]:
    log_parameters = np.log(parameters)
    log_lower = np.log(LOWER_BOUNDS)
    log_upper = np.log(UPPER_BOUNDS)
    nearest_fraction = np.minimum(log_parameters - log_lower, log_upper - log_parameters) / (log_upper - log_lower)
    return {
        "c_th_nearest_log_bound_fraction": float(nearest_fraction[0]),
        "g_th_nearest_log_bound_fraction": float(nearest_fraction[1]),
        "c_th_fraction_of_upper_bound": float(parameters[0] / UPPER_BOUNDS[0]),
        "g_th_fraction_of_upper_bound": float(parameters[1] / UPPER_BOUNDS[1]),
    }


def save_profile(rows: list[dict[str, object]]) -> None:
    with PROFILE_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def plot_profiles(
    rows: list[dict[str, object]],
    best_parameters: np.ndarray,
    best_objective: float,
    best_udds_rmse: float,
) -> None:
    c_th = np.asarray([float(row["c_th_j_per_k"]) for row in rows])
    objective = np.asarray([float(row["joint_objective_k2"]) for row in rows])
    udds_rmse = np.asarray([float(row["udds_rmse_k"]) for row in rows])

    fig, axis = plt.subplots(figsize=(9.5, 5.8), constrained_layout=True)
    axis.plot(c_th, objective, marker="o", color="#0072B2", label="Joint 1C+HPPC profile")
    axis.scatter([best_parameters[0]], [best_objective], marker="*", s=150, color="#D55E00", zorder=4, label="Expanded-bound optimum")
    axis.axvline(PREVIOUS_BOUNDARY_PARAMETERS[0], color="#CC79A7", ls="--", label="Previous 200 J/K guard")
    axis.axvline(float(best_parameters[0]), color="#D55E00", ls=":")
    axis.set_xscale("log")
    axis.set(xlabel="Fixed effective Cth [J/K]", ylabel="Joint objective [K²]", title="Joint 1C+HPPC objective profile (Gth optimized at each point)")
    axis.grid(alpha=0.25, which="both")
    axis.legend()
    fig.savefig(PROFILE_PLOT_PATH, dpi=180)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(9.5, 5.8), constrained_layout=True)
    axis.plot(c_th, udds_rmse, marker="o", color="#009E73", label="Withheld UDDS Validation A")
    axis.scatter([best_parameters[0]], [best_udds_rmse], marker="*", s=150, color="#D55E00", zorder=4, label="Expanded-bound joint optimum")
    axis.axhline(BASELINE_UDDS["rmse_k"], color="black", ls="--", label="1C-only UDDS RMSE")
    axis.axvline(PREVIOUS_BOUNDARY_PARAMETERS[0], color="#CC79A7", ls="--", label="Previous 200 J/K guard")
    axis.set_xscale("log")
    axis.set(xlabel="Fixed effective Cth [J/K]", ylabel="UDDS RMSE [K]", title="Withheld-validation performance (not used for parameter selection)")
    axis.grid(alpha=0.25, which="both")
    axis.legend()
    fig.savefig(PROFILE_UDDS_PLOT_PATH, dpi=180)
    plt.close(fig)


def append_metric(rows: list[dict[str, object]], scope: str, item: str, metric: str, value: object, unit: str = "1") -> None:
    rows.append({"scope": scope, "item": item, "metric": metric, "value": value, "unit": unit})


def save_metrics(
    multistarts: list[dict[str, object]],
    parameters: np.ndarray,
    one_c_metrics: dict[str, float],
    hppc_metrics: dict[str, float],
    udds_metrics: dict[str, float],
    summary: dict[str, object],
    bounds: dict[str, float],
    current_condition: dict[str, object],
    previous_condition: dict[str, object],
    decision: dict[str, object],
) -> None:
    rows: list[dict[str, object]] = []
    append_metric(rows, "expanded_optimum", "parameters", "c_th_j_per_k", parameters[0], "J/K")
    append_metric(rows, "expanded_optimum", "parameters", "g_th_w_per_k", parameters[1], "W/K")
    append_metric(rows, "expanded_optimum", "parameters", "tau_th_s", parameters[0] / parameters[1], "s")
    append_metric(rows, "expanded_optimum", "fit", "joint_objective_k2", 0.5 * (one_c_metrics["mse_k2"] + hppc_metrics["mse_k2"]), "K^2")
    for label, values in (("one_c", one_c_metrics), ("hppc", hppc_metrics), ("withheld_udds", udds_metrics)):
        for metric in ("rmse_k", "mae_k", "mean_signed_error_k", "max_abs_error_k", "r2"):
            append_metric(rows, "expanded_optimum", label, metric, values[metric], "K" if metric.endswith("_k") else "1")
    for metric, value in summary.items():
        append_metric(rows, "multistart_summary", "expanded_bound", metric, value)
    for metric, value in bounds.items():
        append_metric(rows, "bound_proximity", "expanded_optimum", metric, value)
    for label, values in (("expanded_optimum", current_condition), ("previous_200_boundary", previous_condition)):
        for metric, value in values.items():
            append_metric(rows, "conditioning", label, metric, value)
    for start in multistarts:
        for metric, value in start.items():
            append_metric(rows, "multistart", f"start_{start['start_id']}", metric, value)
    for metric, value in decision.items():
        append_metric(rows, "decision", "boundary_audit", metric, value)
    with METRICS_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def profile_markdown(rows: list[dict[str, object]]) -> list[str]:
    output = [
        "| Fixed Cth [J/K] | Optimal Gth [W/K] | tau [s] | Joint objective [K^2] | 1C RMSE [K] | HPPC RMSE [K] | Withheld UDDS RMSE [K] | UDDS bias [K] | UDDS max abs [K] |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        output.append(
            f"| {float(row['c_th_j_per_k']):.9f} | {float(row['g_th_w_per_k']):.12f} | {float(row['tau_th_s']):.6f} | {float(row['joint_objective_k2']):.12f} | {float(row['rmse_1c_k']):.9f} | {float(row['rmse_hppc_k']):.9f} | {float(row['udds_rmse_k']):.9f} | {float(row['udds_bias_k']):+.9f} | {float(row['udds_max_abs_error_k']):.9f} |"
        )
    return output


def save_report(
    parameters: np.ndarray,
    one_c_metrics: dict[str, float],
    hppc_metrics: dict[str, float],
    udds_metrics: dict[str, float],
    profile: list[dict[str, object]],
    summary: dict[str, object],
    bounds: dict[str, float],
    current_condition: dict[str, object],
    previous_condition: dict[str, object],
    decision: dict[str, object],
    hppc_audit: dict[str, object],
) -> None:
    previous_profile = next(row for row in profile if np.isclose(float(row["c_th_j_per_k"]), 200.0))
    high_profile = [row for row in profile if float(row["c_th_j_per_k"]) >= 300.0]
    high_objectives = np.asarray([float(row["joint_objective_k2"]) for row in high_profile])
    lines = [
        "# Panasonic joint thermal boundary/profile audit",
        "",
        f"## {decision['classification']}",
        "",
        "This is a numerical-boundary diagnostic of the unchanged clean one-node model. Only the Cth numerical guard was expanded from 200 to 1000 J/K; Gth retained its 1e-4-to-1 W/K guard. The objective, datasets, heat generation, entropy support, HPPC block handling, and withheld UDDS procedure were unchanged.",
        "",
        "## Expanded-bound joint optimum",
        "",
        f"- Effective lumped `Cth = {parameters[0]:.12f} J/K`",
        f"- `Gth = {parameters[1]:.12f} W/K`",
        f"- `tau_th = {parameters[0]/parameters[1]:.12f} s`",
        f"- Joint objective: `{0.5*(one_c_metrics['mse_k2']+hppc_metrics['mse_k2']):.12f} K^2`",
        f"- 1C RMSE: `{one_c_metrics['rmse_k']:.9f} K`; HPPC RMSE: `{hppc_metrics['rmse_k']:.9f} K`",
        f"- Cth/Gth nearest log-bound fractions: `{bounds['c_th_nearest_log_bound_fraction']:.9f}`, `{bounds['g_th_nearest_log_bound_fraction']:.9f}`",
        f"- Successful multistarts: {summary['successful_starts']}/{N_MULTISTARTS}; near-optimal: {summary['near_optimal_starts']}. Successful ranges: Cth `{summary['successful_c_th_min']:.12g}-{summary['successful_c_th_max']:.12g} J/K`, Gth `{summary['successful_g_th_min']:.12g}-{summary['successful_g_th_max']:.12g} W/K`. Near-optimal ranges: Cth `{summary['near_optimal_c_th_min']:.12g}-{summary['near_optimal_c_th_max']:.12g}`, Gth `{summary['near_optimal_g_th_min']:.12g}-{summary['near_optimal_g_th_max']:.12g}`.",
        f"- HPPC continuity handling remained fixed at {hppc_audit['retained_block_count']} blocks and {hppc_audit['retained_samples']:,} samples; no omitted gap was propagated.",
        "- **Profile outcome A:** Cth converges to an interior optimum; it neither reaches the new 1000 J/K guard nor remains practically flat over the sampled high-Cth region.",
        "",
        "Cth is an **effective lumped thermal capacitance**, not a measured material heat capacity. A very large effective value can compensate for incomplete heat-generation representation, boundary-condition mismatch, experiment-to-experiment differences, or lack of transferability of one common lumped parameter set. None of those effects was repaired here.",
        "",
        "## Fixed-Cth objective profile",
        "",
        *profile_markdown(profile),
        "",
        f"Across the sampled high-Cth region (300-1000 J/K), the joint-objective span relative to its minimum is `{(np.max(high_objectives)-np.min(high_objectives))/np.min(high_objectives):.6%}`. The stated practical-flatness diagnostic threshold is {FLAT_PROFILE_RELATIVE_OBJECTIVE:.1%}.",
        "At fixed Cth values of 500 J/K and above, the profile drives Gth to its unchanged 1e-4 W/K lower numerical guard while the joint objective rises. Those points do not constitute alternative joint optima.",
        "",
        "UDDS values in this table were calculated only after Gth had been optimized against 1C+HPPC. They did not influence any optimization or parameter selection.",
        "",
        "## Conditioning",
        "",
        "| Point | Rank | Singular values | Condition number | Sensitivity correlation | Parameter correlation |",
        "|---|---:|---:|---:|---:|---:|",
        f"| Previous 200 J/K boundary | {previous_condition['rank']}/2 | {previous_condition['singular_value_1']:.12g}, {previous_condition['singular_value_2']:.12g} | {previous_condition['condition_number']:.12g} | {previous_condition['sensitivity_correlation']:+.12f} | {previous_condition['parameter_correlation']:+.12f} |",
        f"| Expanded-bound optimum | {current_condition['rank']}/2 | {current_condition['singular_value_1']:.12g}, {current_condition['singular_value_2']:.12g} | {current_condition['condition_number']:.12g} | {current_condition['sensitivity_correlation']:+.12f} | {current_condition['parameter_correlation']:+.12f} |",
        "",
        "## Withheld UDDS at expanded-bound optimum",
        "",
        f"RMSE `{udds_metrics['rmse_k']:.9f} K`, bias `{udds_metrics['mean_signed_error_k']:+.9f} K`, maximum absolute error `{udds_metrics['max_abs_error_k']:.9f} K`. The RMSE change relative to the frozen 1C-only value `{BASELINE_UDDS['rmse_k']:.9f} K` is `{100.0*(1.0-udds_metrics['rmse_k']/BASELINE_UDDS['rmse_k']):+.6f}%`.",
        f"The improvement therefore persists relative to the frozen 1C-only baseline, but it shrinks from the previous 200 J/K profile-point RMSE `{float(previous_profile['udds_rmse_k']):.9f} K` to `{udds_metrics['rmse_k']:.9f} K`; the bias also reverses from `{float(previous_profile['udds_bias_k']):+.9f} K` to `{udds_metrics['mean_signed_error_k']:+.9f} K`.",
        "",
        "## Decisions",
        "",
        f"1. **Was Cth=200 J/K artificially truncated?** **{'Yes' if decision['previous_guard_truncated'] else 'No'}**. At fixed 200 J/K the profile objective is {float(previous_profile['joint_objective_k2']):.12f} K^2; the expanded optimum is {0.5*(one_c_metrics['mse_k2']+hppc_metrics['mse_k2']):.12f} K^2.",
        f"2. **Is there a finite interior joint optimum?** **{'Yes' if decision['finite_interior_optimum'] else 'No'}**.",
        f"3. **Is the objective effectively flat toward large Cth?** **{'Yes' if decision['practically_flat'] else 'No'}**, using the reported 300-1000 J/K profile span and the explicit 1% diagnostic threshold.",
        f"4. **Does withheld UDDS improvement persist after relaxing the bound?** **{'Yes' if decision['udds_improvement_persists'] else 'No'}**.",
        f"5. **Is the joint parameterization numerically well-defined?** **{'Yes' if decision['numerically_well_defined'] else 'No'}**.",
        "",
        f"**Final classification: {decision['classification']}.**",
        "",
        "The frozen thermal parameter module was not modified. No UDDS fitting, architecture change, empirical correction, electrical refit, OCP change, or entropy-coefficient change was performed.",
    ]
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    one_c = prepare_one_c_data()
    blocks, hppc_audit, _block_rows, _gaps = prepare_hppc_blocks()
    udds, udds_heat = prepare_udds_heat()

    multistarts, best = run_expanded_multistart(one_c, blocks)
    parameters = np.asarray(best["parameters"])
    one_c_metrics, hppc_metrics = experiment_metrics(one_c, blocks, parameters)
    udds_metrics = evaluate_udds(udds, udds_heat, parameters)
    profile = build_profile(one_c, blocks, udds, udds_heat)
    summary = multistart_summary(multistarts)
    bounds = parameter_bound_proximity(parameters)
    current_condition = conditioning(np.log(parameters), one_c, blocks)
    previous_condition = conditioning(np.log(PREVIOUS_BOUNDARY_PARAMETERS), one_c, blocks)

    best_objective = 0.5 * (one_c_metrics["mse_k2"] + hppc_metrics["mse_k2"])
    profile_200 = next(row for row in profile if np.isclose(float(row["c_th_j_per_k"]), 200.0))
    if not np.isclose(float(profile_200["g_th_w_per_k"]), PREVIOUS_BOUNDARY_PARAMETERS[1], atol=2.0e-9, rtol=0.0):
        raise RuntimeError("The 200 J/K profile point did not reproduce the previous Gth result.")
    if not np.isclose(float(profile_200["joint_objective_k2"]), 0.16374510744244164, atol=2.0e-10, rtol=0.0):
        raise RuntimeError("The 200 J/K profile point did not reproduce the previous objective.")
    high_profile = [row for row in profile if float(row["c_th_j_per_k"]) >= 300.0]
    high_objectives = np.asarray([float(row["joint_objective_k2"]) for row in high_profile])
    high_relative_span = float((np.max(high_objectives) - np.min(high_objectives)) / np.min(high_objectives))
    at_upper = bool(parameters[0] >= UPPER_BOUNDS[0] * (1.0 - 1.0e-7))
    practically_flat = bool(high_relative_span <= FLAT_PROFILE_RELATIVE_OBJECTIVE)
    finite_interior = bool(not at_upper and bounds["c_th_nearest_log_bound_fraction"] > 0.01)
    if at_upper:
        classification = "JOINT OPTIMUM REMAINS BOUNDARY LIMITED"
    elif practically_flat:
        classification = "JOINT OPTIMUM PRACTICALLY FLAT / WEAKLY DETERMINED"
    else:
        classification = "JOINT OPTIMUM RESOLVED"
    decision = {
        "previous_guard_truncated": bool(
            parameters[0] > 200.0 * (1.0 + 1.0e-6)
            and best_objective < float(profile_200["joint_objective_k2"])
        ),
        "finite_interior_optimum": finite_interior,
        "practically_flat": practically_flat,
        "udds_improvement_persists": bool(udds_metrics["rmse_k"] < BASELINE_UDDS["rmse_k"]),
        "numerically_well_defined": bool(
            finite_interior
            and not practically_flat
            and current_condition["rank"] == 2
            and float(current_condition["condition_number"]) < 1.0e3
            and abs(float(current_condition["parameter_correlation"])) < 0.95
        ),
        "high_c_th_profile_relative_objective_span": high_relative_span,
        "classification": classification,
    }

    save_profile(profile)
    plot_profiles(profile, parameters, best_objective, udds_metrics["rmse_k"])
    save_metrics(
        multistarts,
        parameters,
        one_c_metrics,
        hppc_metrics,
        udds_metrics,
        summary,
        bounds,
        current_condition,
        previous_condition,
        decision,
    )
    save_report(
        parameters,
        one_c_metrics,
        hppc_metrics,
        udds_metrics,
        profile,
        summary,
        bounds,
        current_condition,
        previous_condition,
        decision,
        hppc_audit,
    )

    print(f"{classification}")
    print(f"Cth={parameters[0]:.12f} J/K; Gth={parameters[1]:.12f} W/K; tau={parameters[0]/parameters[1]:.9f} s")
    print(f"Objective={best_objective:.12f} K^2; 1C RMSE={one_c_metrics['rmse_k']:.9f} K; HPPC RMSE={hppc_metrics['rmse_k']:.9f} K")
    print(f"UDDS RMSE={udds_metrics['rmse_k']:.9f} K; bias={udds_metrics['mean_signed_error_k']:+.9f} K")


if __name__ == "__main__":
    main()
