"""Identify a Panasonic positive OCP using Guo MCMB and the Qi-style method.

The objective uses only 0.05 <= normalized discharged capacity <= 0.95.
Fitted variables are four stoichiometric-window endpoints and thirteen fixed-
coordinate positive-electrode OCP knot voltages.  No dynamic terms are used.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, minimize


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from models.ocp_panasonic_qi import (  # noqa: E402
    MCMB_SUPPORTED_STOICHIOMETRY_RANGE,
    NCA_WORKING_STOICHIOMETRY_RANGE,
    N_UP_KNOTS,
    UP_KNOT_STOICHIOMETRIES,
    UP_VOLTAGE_BOUNDS_V,
    U_n,
    U_p,
    validate_positive_knot_voltages,
)


INPUT_PATH = Path(
    r"C:\Users\amvel\Documents\Andes\Tesis\datasets\converted_csv"
    r"\c20_ocv_25C_converted.csv"
)
OUTPUT_DIRECTORY = REPOSITORY_ROOT / "thesis_data"
OUTPUT_STEM = "panasonic_c20_qi_up"

Q_CELL_AH = 2.994979138519
DISCHARGE_START_INDEX = 6
DISCHARGE_END_INDEX_INCLUSIVE = 1246
FIT_Z_MIN = 0.05
FIT_Z_MAX = 0.95
NUMBER_OF_STARTS = 40
MULTISTART_SEED = 20260913
MIN_DIRECTION_GAP = 1.0e-8
GOOD_SOLUTION_RMSE_DELTA_V = 5.0e-4
MATERIAL_WINDOW_DIFFERENCE = 0.05
NEAR_STOICHIOMETRY_BOUND = 1.0e-3
NEAR_VOLTAGE_BOUND_V = 5.0e-3
FLAT_KNOT_DROP_TOLERANCE_V = 1.0e-4
LONG_FLAT_RUN_ADJACENCIES = 3

# Diagnostic acceptance thresholds, applied jointly rather than using RMSE
# alone.  The bin criteria reject the broad residual lobes seen with the
# Karthikeyan pair.
MAX_FIT_RMSE_V = 0.015
MAX_FIT_MAE_V = 0.012
MIN_FIT_R_SQUARED = 0.995
MAX_ABS_INTERIOR_BIN_MEAN_V = 0.020
MAX_INTERIOR_BIN_MEAN_SPAN_V = 0.035


@dataclass(frozen=True)
class DischargeData:
    source_indices: np.ndarray
    time_s: np.ndarray
    voltage_v: np.ndarray
    current_a: np.ndarray
    exported_ah: np.ndarray
    discharged_capacity_ah: np.ndarray
    z: np.ndarray
    fit_mask: np.ndarray


@dataclass(frozen=True)
class Metrics:
    rmse_v: float
    mae_v: float
    r_squared: float


def load_data(path: Path = INPUT_PATH) -> DischargeData:
    table = np.genfromtxt(path, delimiter=",", names=True, dtype=float, encoding="utf-8")
    required = {"Time", "Voltage", "Current", "Ah"}
    missing = sorted(required - set(table.dtype.names or ()))
    if missing:
        raise ValueError(f"Panasonic C/20 CSV is missing columns: {missing}")
    selection = slice(DISCHARGE_START_INDEX, DISCHARGE_END_INDEX_INCLUSIVE + 1)
    time_s = np.asarray(table["Time"][selection], dtype=float)
    voltage_v = np.asarray(table["Voltage"][selection], dtype=float)
    current_a = np.asarray(table["Current"][selection], dtype=float)
    exported_ah = np.asarray(table["Ah"][selection], dtype=float)
    if not all(
        np.all(np.isfinite(values))
        for values in (time_s, voltage_v, current_a, exported_ah)
    ):
        raise ValueError("Panasonic discharge branch contains non-finite data.")
    if np.any(np.diff(time_s) <= 0.0):
        raise ValueError("Time must increase strictly through the discharge branch.")
    if np.any(current_a >= 0.0):
        raise ValueError("Selected branch must contain only raw negative discharge current.")

    increments_ah = (
        0.5
        * (np.abs(current_a[:-1]) + np.abs(current_a[1:]))
        * np.diff(time_s)
        / 3600.0
    )
    discharged_capacity_ah = np.concatenate(([0.0], np.cumsum(increments_ah)))
    final_capacity_ah = float(discharged_capacity_ah[-1])
    if not np.isclose(final_capacity_ah, Q_CELL_AH, rtol=0.0, atol=1.0e-9):
        raise ValueError(
            f"Integrated capacity {final_capacity_ah:.12f} Ah does not match the "
            f"Panasonic capacity {Q_CELL_AH:.12f} Ah."
        )
    z = discharged_capacity_ah / final_capacity_ah
    fit_mask = (z >= FIT_Z_MIN) & (z <= FIT_Z_MAX)
    return DischargeData(
        source_indices=np.arange(
            DISCHARGE_START_INDEX, DISCHARGE_END_INDEX_INCLUSIVE + 1, dtype=int
        ),
        time_s=time_s,
        voltage_v=voltage_v,
        current_a=current_a,
        exported_ah=exported_ah,
        discharged_capacity_ah=discharged_capacity_ah,
        z=z,
        fit_mask=fit_mask,
    )


def unpack_parameters(parameters: np.ndarray) -> tuple[float, float, float, float, np.ndarray]:
    values = np.asarray(parameters, dtype=float)
    if values.shape != (4 + N_UP_KNOTS,):
        raise ValueError(f"Exactly {4 + N_UP_KNOTS} fit parameters are required.")
    return values[0], values[1], values[2], values[3], values[4:]


def is_feasible(parameters: np.ndarray, tolerance: float = 1.0e-10) -> bool:
    try:
        xn_start, xn_end, xp_start, xp_end, up_knots = unpack_parameters(parameters)
    except ValueError:
        return False
    if not np.all(np.isfinite(parameters)):
        return False
    xn_lower, xn_upper = MCMB_SUPPORTED_STOICHIOMETRY_RANGE
    xp_lower, xp_upper = NCA_WORKING_STOICHIOMETRY_RANGE
    voltage_lower, voltage_upper = UP_VOLTAGE_BOUNDS_V
    return bool(
        xn_lower - tolerance <= xn_start <= xn_upper + tolerance
        and xn_lower - tolerance <= xn_end <= xn_upper + tolerance
        and xp_lower - tolerance <= xp_start <= xp_upper + tolerance
        and xp_lower - tolerance <= xp_end <= xp_upper + tolerance
        and xn_start - xn_end >= MIN_DIRECTION_GAP - tolerance
        and xp_end - xp_start >= MIN_DIRECTION_GAP - tolerance
        and np.all(up_knots >= voltage_lower - tolerance)
        and np.all(up_knots <= voltage_upper + tolerance)
        and np.all(np.diff(up_knots) <= tolerance)
    )


def evaluate(
    parameters: np.ndarray,
    z: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if not is_feasible(parameters):
        raise ValueError("Infeasible OCP parameters were rejected; no clipping was used.")
    xn_start, xn_end, xp_start, xp_end, up_knots = unpack_parameters(parameters)
    validate_positive_knot_voltages(up_knots)
    x_n = xn_start + z * (xn_end - xn_start)
    x_p = xp_start + z * (xp_end - xp_start)
    voltage_v = U_p(x_p, up_knots) - U_n(x_n)
    return np.asarray(voltage_v), x_n, x_p


def metrics(measured_v: np.ndarray, fitted_v: np.ndarray) -> Metrics:
    residual_v = measured_v - fitted_v
    ss_res = float(np.sum(residual_v**2))
    ss_total = float(np.sum((measured_v - np.mean(measured_v)) ** 2))
    return Metrics(
        rmse_v=float(np.sqrt(np.mean(residual_v**2))),
        mae_v=float(np.mean(np.abs(residual_v))),
        r_squared=1.0 - ss_res / ss_total,
    )


def deterministic_windows() -> list[np.ndarray]:
    explicit = [
        (0.78, 0.02, 0.20, 0.94),
        (0.75, 0.05, 0.22, 0.92),
        (0.70, 0.10, 0.25, 0.90),
        (0.65, 0.05, 0.18, 0.88),
        (0.80, 0.15, 0.30, 0.96),
        (0.60, 0.10, 0.20, 0.80),
        (0.50, 0.02, 0.35, 0.95),
        (0.79, 0.30, 0.18, 0.75),
    ]
    windows = [np.asarray(values, dtype=float) for values in explicit]
    rng = np.random.default_rng(MULTISTART_SEED)
    xn_lower, xn_upper = MCMB_SUPPORTED_STOICHIOMETRY_RANGE
    xp_lower, xp_upper = NCA_WORKING_STOICHIOMETRY_RANGE
    while len(windows) < NUMBER_OF_STARTS:
        xn_end = rng.uniform(xn_lower, xn_upper - 0.06)
        xn_start = rng.uniform(xn_end + 0.05, xn_upper)
        xp_start = rng.uniform(xp_lower, xp_upper - 0.06)
        xp_end = rng.uniform(xp_start + 0.05, xp_upper)
        windows.append(np.asarray([xn_start, xn_end, xp_start, xp_end]))
    return windows


def initial_positive_knots(window: np.ndarray, data: DischargeData) -> np.ndarray:
    xn_start, xn_end, xp_start, xp_end = window
    x_n = xn_start + data.z * (xn_end - xn_start)
    x_p = xp_start + data.z * (xp_end - xp_start)
    fit_xp = x_p[data.fit_mask]
    target_up = data.voltage_v[data.fit_mask] + U_n(x_n[data.fit_mask])
    knot_values = np.interp(UP_KNOT_STOICHIOMETRIES, fit_xp, target_up)

    # Extend unused initial-knot regions with the local data slope.  This only
    # creates a feasible starting point; no extrapolation is used by U_p.
    slope_left = np.polyfit(fit_xp[:40], target_up[:40], 1)[0]
    slope_right = np.polyfit(fit_xp[-40:], target_up[-40:], 1)[0]
    left = UP_KNOT_STOICHIOMETRIES < fit_xp[0]
    right = UP_KNOT_STOICHIOMETRIES > fit_xp[-1]
    knot_values[left] = target_up[0] + slope_left * (
        UP_KNOT_STOICHIOMETRIES[left] - fit_xp[0]
    )
    knot_values[right] = target_up[-1] + slope_right * (
        UP_KNOT_STOICHIOMETRIES[right] - fit_xp[-1]
    )
    knot_values = np.minimum.accumulate(knot_values)
    lower, upper = UP_VOLTAGE_BOUNDS_V
    if np.any((knot_values <= lower) | (knot_values >= upper)):
        knot_values = np.linspace(4.35, 2.70, N_UP_KNOTS)
    return knot_values


def run_multistart(data: DischargeData) -> list[dict[str, object]]:
    xn_lower, xn_upper = MCMB_SUPPORTED_STOICHIOMETRY_RANGE
    xp_lower, xp_upper = NCA_WORKING_STOICHIOMETRY_RANGE
    voltage_lower, voltage_upper = UP_VOLTAGE_BOUNDS_V
    lower_bounds = np.asarray(
        [xn_lower, xn_lower, xp_lower, xp_lower]
        + [voltage_lower] * N_UP_KNOTS
    )
    upper_bounds = np.asarray(
        [xn_upper, xn_upper, xp_upper, xp_upper]
        + [voltage_upper] * N_UP_KNOTS
    )
    bounds = Bounds(lower_bounds, upper_bounds, keep_feasible=True)

    constraint_matrix = np.zeros((2 + N_UP_KNOTS - 1, 4 + N_UP_KNOTS))
    constraint_matrix[0, 0:2] = [1.0, -1.0]
    constraint_matrix[1, 2:4] = [-1.0, 1.0]
    for index in range(N_UP_KNOTS - 1):
        constraint_matrix[2 + index, 4 + index] = 1.0
        constraint_matrix[2 + index, 4 + index + 1] = -1.0
    constraint_lower = np.concatenate(
        ([MIN_DIRECTION_GAP, MIN_DIRECTION_GAP], np.zeros(N_UP_KNOTS - 1))
    )
    constraints = LinearConstraint(constraint_matrix, constraint_lower, np.inf)

    def objective(parameters: np.ndarray) -> float:
        xn_start, xn_end, xp_start, xp_end, up_knots = unpack_parameters(parameters)
        if not (
            xn_start > xn_end
            and xp_start < xp_end
            and np.all(np.diff(up_knots) <= 1.0e-9)
        ):
            return 1.0e6
        x_n = xn_start + data.z[data.fit_mask] * (xn_end - xn_start)
        x_p = xp_start + data.z[data.fit_mask] * (xp_end - xp_start)
        try:
            fitted_v = U_p(x_p, up_knots) - U_n(x_n)
        except ValueError:
            return 1.0e6
        residual_v = data.voltage_v[data.fit_mask] - fitted_v
        return float(np.mean(residual_v**2))

    records: list[dict[str, object]] = []
    for start_id, window in enumerate(deterministic_windows(), start=1):
        initial_knots = initial_positive_knots(window, data)
        initial = np.concatenate((window, initial_knots))
        result = minimize(
            objective,
            initial,
            method="SLSQP",
            bounds=bounds,
            constraints=[constraints],
            options={"ftol": 1.0e-14, "maxiter": 3000, "disp": False},
        )
        parameters = np.asarray(result.x, dtype=float)
        feasible = is_feasible(parameters)
        successful = bool(result.success and feasible and np.isfinite(result.fun))
        record: dict[str, object] = {
            "start_id": start_id,
            "optimizer_success": bool(result.success),
            "feasible": feasible,
            "successful_solution": successful,
            "status": int(result.status),
            "message": str(result.message),
            "iterations": int(result.nit),
            "function_evaluations": int(result.nfev),
            "objective_mse_v2": float(result.fun),
            "fit_rmse_v": np.nan,
            "fit_mae_v": np.nan,
            "fit_r_squared": np.nan,
            "full_rmse_v": np.nan,
            "full_mae_v": np.nan,
            "full_r_squared": np.nan,
            "x_n_start": parameters[0],
            "x_n_end": parameters[1],
            "x_p_start": parameters[2],
            "x_p_end": parameters[3],
            "is_good_solution": False,
        }
        for index, value in enumerate(parameters[4:]):
            record[f"u_p_knot_{index:02d}_v"] = value
        if successful:
            fitted_v, _, _ = evaluate(parameters, data.z)
            fit_metrics = metrics(
                data.voltage_v[data.fit_mask], fitted_v[data.fit_mask]
            )
            full_metrics = metrics(data.voltage_v, fitted_v)
            record.update(
                fit_rmse_v=fit_metrics.rmse_v,
                fit_mae_v=fit_metrics.mae_v,
                fit_r_squared=fit_metrics.r_squared,
                full_rmse_v=full_metrics.rmse_v,
                full_mae_v=full_metrics.mae_v,
                full_r_squared=full_metrics.r_squared,
            )
        records.append(record)

    successful_records = [record for record in records if record["successful_solution"]]
    if not successful_records:
        raise RuntimeError("No feasible Qi-style multistart solution completed.")
    best_rmse = min(float(record["fit_rmse_v"]) for record in successful_records)
    for record in successful_records:
        record["is_good_solution"] = (
            float(record["fit_rmse_v"]) <= best_rmse + GOOD_SOLUTION_RMSE_DELTA_V
        )
    return records


def interior_bin_statistics(
    data: DischargeData, residual_v: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    edges = np.linspace(FIT_Z_MIN, FIT_Z_MAX, 19)
    centers = 0.5 * (edges[:-1] + edges[1:])
    selected_z = data.z[data.fit_mask]
    selected_residual = residual_v[data.fit_mask]
    indices = np.minimum(np.digitize(selected_z, edges) - 1, 17)
    means = np.asarray(
        [np.mean(selected_residual[indices == index]) for index in range(18)]
    )
    return centers, means


def longest_flat_run(knot_voltages_v: np.ndarray) -> int:
    flat = np.abs(np.diff(knot_voltages_v)) <= FLAT_KNOT_DROP_TOLERANCE_V
    longest = current = 0
    for is_flat in flat:
        current = current + 1 if is_flat else 0
        longest = max(longest, current)
    return longest


def _parameters_from_record(record: dict[str, object]) -> np.ndarray:
    return np.asarray(
        [
            record["x_n_start"],
            record["x_n_end"],
            record["x_p_start"],
            record["x_p_end"],
        ]
        + [record[f"u_p_knot_{index:02d}_v"] for index in range(N_UP_KNOTS)],
        dtype=float,
    )


def _csv_value(value: object) -> object:
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.12g}"
    if isinstance(value, (int, np.integer)) and not isinstance(value, bool):
        return int(value)
    return value


def write_multistart(path: Path, records: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        for record in records:
            writer.writerow({name: _csv_value(value) for name, value in record.items()})


def write_fit_data(
    path: Path,
    data: DischargeData,
    fitted_v: np.ndarray,
    residual_v: np.ndarray,
    x_n: np.ndarray,
    x_p: np.ndarray,
) -> None:
    header = [
        "source_row_index",
        "time_s",
        "raw_current_a",
        "exported_capacity_ah",
        "discharged_capacity_ah",
        "normalized_discharged_capacity",
        "in_primary_fit_region",
        "measured_voltage_v",
        "fitted_equilibrium_voltage_v",
        "residual_measured_minus_fitted_v",
        "x_n",
        "x_p",
    ]
    rows = zip(
        data.source_indices,
        data.time_s,
        data.current_a,
        data.exported_ah,
        data.discharged_capacity_ah,
        data.z,
        data.fit_mask,
        data.voltage_v,
        fitted_v,
        residual_v,
        x_n,
        x_p,
        strict=True,
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for row in rows:
            writer.writerow([_csv_value(value) for value in row])


def write_knots(path: Path, knot_voltages_v: np.ndarray) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["knot_index", "x_p", "u_p_v", "drop_to_next_knot_v"])
        for index, (x_p, voltage) in enumerate(
            zip(UP_KNOT_STOICHIOMETRIES, knot_voltages_v, strict=True)
        ):
            drop = voltage - knot_voltages_v[index + 1] if index + 1 < N_UP_KNOTS else ""
            writer.writerow([index, f"{x_p:.12g}", f"{voltage:.12g}", _csv_value(drop)])


def save_figures(
    data: DischargeData,
    fitted_v: np.ndarray,
    residual_v: np.ndarray,
    x_n: np.ndarray,
    x_p: np.ndarray,
    knot_voltages_v: np.ndarray,
    bin_centers: np.ndarray,
    bin_means: np.ndarray,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fit_figure, fit_axis = plt.subplots(figsize=(10.5, 5.8), constrained_layout=True)
    fit_axis.plot(data.z, data.voltage_v, linewidth=2.0, label="Measured Panasonic C/20")
    fit_axis.plot(data.z, fitted_v, linewidth=1.8, label="$U_p(x_p)-U_n(x_n)$")
    fit_axis.axvspan(FIT_Z_MIN, FIT_Z_MAX, color="tab:green", alpha=0.09, label="Fit region")
    fit_axis.set(
        xlabel="Normalized discharged capacity, $z$",
        ylabel="Voltage [V]",
        title="Panasonic C/20 Qi-style OCP identification",
    )
    fit_axis.grid(alpha=0.25)
    fit_axis.legend()
    fit_figure.savefig(OUTPUT_DIRECTORY / f"{OUTPUT_STEM}_fit.png", dpi=180)
    plt.close(fit_figure)

    residual_figure, residual_axes = plt.subplots(
        2, 1, figsize=(10.5, 7.0), constrained_layout=True
    )
    residual_axes[0].plot(data.z, 1000 * residual_v, linewidth=1.2)
    residual_axes[0].axvspan(FIT_Z_MIN, FIT_Z_MAX, color="tab:green", alpha=0.09)
    residual_axes[0].axhline(0.0, color="black", linewidth=0.8)
    residual_axes[0].set(
        ylabel="Measured - fitted [mV]",
        title="Full-discharge residual; objective uses 5--95% only",
    )
    residual_axes[1].plot(
        bin_centers, 1000 * bin_means, marker="o", linewidth=1.2,
        label="5%-wide fit-region bin means"
    )
    residual_axes[1].axhline(0.0, color="black", linewidth=0.8)
    residual_axes[1].set(
        xlabel="Normalized discharged capacity, $z$",
        ylabel="Mean residual [mV]",
    )
    residual_axes[1].grid(alpha=0.25)
    residual_axes[1].legend()
    residual_axes[0].grid(alpha=0.25)
    residual_figure.savefig(OUTPUT_DIRECTORY / f"{OUTPUT_STEM}_residual.png", dpi=180)
    plt.close(residual_figure)

    ocp_figure, ocp_axes = plt.subplots(1, 2, figsize=(11.0, 4.8), constrained_layout=True)
    dense_xn = np.linspace(*MCMB_SUPPORTED_STOICHIOMETRY_RANGE, 2000)
    dense_xp = np.linspace(*NCA_WORKING_STOICHIOMETRY_RANGE, 2000)
    ocp_axes[0].plot(dense_xn, U_n(dense_xn), color="tab:blue", linewidth=1.7)
    ocp_axes[0].axvspan(np.min(x_n), np.max(x_n), color="tab:blue", alpha=0.08)
    ocp_axes[0].set(
        xlabel="$x_n$", ylabel="$U_n$ [V vs Li/Li$^+$]",
        title="Fixed Guo MCMB discharge OCP"
    )
    ocp_axes[1].plot(dense_xp, U_p(dense_xp, knot_voltages_v), color="tab:red", linewidth=1.7)
    ocp_axes[1].plot(UP_KNOT_STOICHIOMETRIES, knot_voltages_v, "o", color="tab:red", label="Fitted knots")
    ocp_axes[1].axvspan(np.min(x_p), np.max(x_p), color="tab:red", alpha=0.08)
    ocp_axes[1].set(
        xlabel="$x_p$", ylabel="$U_p$ [V vs Li/Li$^+$]",
        title="Identified monotonic positive OCP"
    )
    ocp_axes[1].legend()
    for axis in ocp_axes:
        axis.grid(alpha=0.25)
    ocp_figure.savefig(OUTPUT_DIRECTORY / f"{OUTPUT_STEM}_ocps.png", dpi=180)
    plt.close(ocp_figure)

    trajectory_figure, trajectory_axes = plt.subplots(
        1, 2, figsize=(11.0, 4.6), constrained_layout=True
    )
    trajectory_axes[0].plot(data.z, x_n, color="tab:blue", linewidth=1.8)
    trajectory_axes[0].set(xlabel="$z$", ylabel="$x_n$", title="MCMB trajectory")
    trajectory_axes[1].plot(data.z, x_p, color="tab:red", linewidth=1.8)
    trajectory_axes[1].set(xlabel="$z$", ylabel="$x_p$", title="Positive-electrode trajectory")
    for axis in trajectory_axes:
        axis.axvspan(FIT_Z_MIN, FIT_Z_MAX, color="tab:green", alpha=0.09)
        axis.grid(alpha=0.25)
    trajectory_figure.savefig(
        OUTPUT_DIRECTORY / f"{OUTPUT_STEM}_stoichiometries.png", dpi=180
    )
    plt.close(trajectory_figure)


def write_report(
    path: Path,
    data: DischargeData,
    best: dict[str, object],
    records: list[dict[str, object]],
    residual_v: np.ndarray,
    fit_metrics: Metrics,
    full_metrics: Metrics,
    bin_means: np.ndarray,
    decision: str,
    checks: dict[str, bool],
    diagnostics: dict[str, object],
) -> None:
    parameters = _parameters_from_record(best)
    xn_start, xn_end, xp_start, xp_end, knot_voltages = unpack_parameters(parameters)
    successful = [record for record in records if record["successful_solution"]]
    good = [record for record in records if record["is_good_solution"]]
    fit_indices = np.flatnonzero(data.fit_mask)
    lines = [
        "# Panasonic NCR18650PF Qi-style OCP identification",
        "",
        f"## Decision: **{decision}**",
        "",
        "The negative OCP is the fixed Guo--White MCMB discharge spline. The "
        "positive OCP is represented by 13 fitted voltages at fixed "
        "stoichiometric knots with monotonic linear interpolation. Only the four "
        "stoichiometric endpoints and those 13 voltages were fitted.",
        "",
        "## Data and objective",
        "",
        f"- Panasonic capacity: `{Q_CELL_AH:.12f} Ah`",
        f"- Discharge rows: {DISCHARGE_START_INDEX}--{DISCHARGE_END_INDEX_INCLUSIVE}, inclusive",
        f"- Full samples: {data.z.size}",
        f"- Primary objective samples: {fit_indices.size}, from "
        f"`z={data.z[fit_indices[0]]:.6f}` to `z={data.z[fit_indices[-1]]:.6f}`",
        "- No clipping, resistance, kinetics, diffusion, overpotential, or thermal terms",
        "",
        "## Fit-region metrics: 5--95%",
        "",
        "| Metric | Value |",
        "| --- | ---: |",
        f"| RMSE | {1000 * fit_metrics.rmse_v:.3f} mV |",
        f"| MAE | {1000 * fit_metrics.mae_v:.3f} mV |",
        f"| R² | {fit_metrics.r_squared:.9f} |",
        f"| Mean-residual bin range | {1000 * np.min(bin_means):.3f} to "
        f"{1000 * np.max(bin_means):.3f} mV |",
        f"| Mean-residual bin span | {1000 * np.ptp(bin_means):.3f} mV |",
        "",
        "## Full-discharge evaluation: 0--100%",
        "",
        "| Metric | Value |",
        "| --- | ---: |",
        f"| RMSE | {1000 * full_metrics.rmse_v:.3f} mV |",
        f"| MAE | {1000 * full_metrics.mae_v:.3f} mV |",
        f"| R² | {full_metrics.r_squared:.9f} |",
        f"| Residual at `z=0` | {1000 * residual_v[0]:.3f} mV |",
        f"| Residual at `z=1` | {1000 * residual_v[-1]:.3f} mV |",
        "",
        "## Fitted stoichiometric windows",
        "",
        "| Endpoint | Value |",
        "| --- | ---: |",
        f"| `x_n_start` | {xn_start:.9f} |",
        f"| `x_n_end` | {xn_end:.9f} |",
        f"| `x_p_start` | {xp_start:.9f} |",
        f"| `x_p_end` | {xp_end:.9f} |",
        "",
        "## Identified positive-OCP knots",
        "",
        "| Knot | `x_p` | `U_p` [V] | Drop to next [mV] |",
        "| ---: | ---: | ---: | ---: |",
    ]
    for index, (x_value, voltage) in enumerate(
        zip(UP_KNOT_STOICHIOMETRIES, knot_voltages, strict=True)
    ):
        drop = (
            1000 * (voltage - knot_voltages[index + 1])
            if index + 1 < N_UP_KNOTS
            else np.nan
        )
        drop_text = "--" if np.isnan(drop) else f"{drop:.3f}"
        lines.append(f"| {index} | {x_value:.6f} | {voltage:.9f} | {drop_text} |")
    lines.extend(
        [
            "",
            "## Constraint and pathology audit",
            "",
            "| Check | Finding |",
            "| --- | --- |",
            f"| Physical directions and domains | {'PASS' if checks['domains'] else 'FAIL'} |",
            f"| Monotonic positive OCP | {'PASS' if checks['monotonic'] else 'FAIL'} |",
            f"| Stoichiometric endpoint near a domain bound | "
            f"{'YES' if diagnostics['stoichiometry_near_bound'] else 'NO'} |",
            f"| Positive-OCP knot near 2.5 or 4.4 V | "
            f"{'YES' if diagnostics['voltage_near_bound'] else 'NO'} |",
            f"| Long flat knot run | {'YES' if diagnostics['long_flat_run'] else 'NO'} "
            f"(longest adjacent run: {diagnostics['longest_flat_run']}) |",
            f"| Broad systematic fit-region residual lobes | "
            f"{'YES' if not checks['residual_shape'] else 'NO'} |",
            "",
            "## Multistart identifiability",
            "",
            f"- Deterministic starts: {len(records)}",
            f"- Successful feasible solutions: {len(successful)}",
            f"- Good solutions within {1000 * GOOD_SOLUTION_RMSE_DELTA_V:.1f} mV "
            f"RMSE of best: {len(good)}",
            f"- Largest endpoint span across good solutions: "
            f"{diagnostics['largest_good_endpoint_span']:.9f}",
            f"- Largest fit-region voltage difference among good solutions: "
            f"{1000 * diagnostics['largest_good_curve_difference_v']:.3f} mV",
            f"- Materially different windows with nearly identical voltage error: "
            f"**{'YES' if diagnostics['materially_different_windows'] else 'NO'}**",
            "",
            "| Start | RMSE [mV] | MAE [mV] | R² | `x_n_start` | `x_n_end` | `x_p_start` | `x_p_end` |",
            "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for record in sorted(good, key=lambda item: float(item["fit_rmse_v"])):
        lines.append(
            f"| {record['start_id']} | {1000 * float(record['fit_rmse_v']):.3f} | "
            f"{1000 * float(record['fit_mae_v']):.3f} | "
            f"{float(record['fit_r_squared']):.9f} | "
            f"{float(record['x_n_start']):.9f} | {float(record['x_n_end']):.9f} | "
            f"{float(record['x_p_start']):.9f} | {float(record['x_p_end']):.9f} |"
        )
    lines.extend(
        [
            "",
            "Every successful result, its metrics, four endpoints, and 13 knot "
            "voltages is recorded in `panasonic_c20_qi_up_multistart.csv`.",
            "",
            "## Decision checks",
            "",
            "| Requirement | Result |",
            "| --- | --- |",
        ]
    )
    for name, passed in checks.items():
        lines.append(f"| `{name}` | {'PASS' if passed else 'FAIL'} |")
    lines.extend(["", "## Interpretation", ""])
    if decision == "ACCEPTABLE OCP IDENTIFICATION":
        lines.append(
            "The Guo-negative/Qi-style-positive construction reproduces the "
            "Panasonic 5--95% C/20 voltage shape without broad systematic "
            "residual mismatch and passes the OCP compatibility checks."
        )
        if diagnostics["materially_different_windows"]:
            lines.append(
                "**Voltage/OCP representation acceptable, stoichiometric "
                "parameters weakly identifiable.** Materially different "
                "near-optimal windows exist, so the stoichiometric solution "
                "is not unique."
            )
    else:
        failures = [name for name, passed in checks.items() if not passed]
        lines.append(
            "The identified OCP construction does not pass this gate. Failed "
            "checks: " + ", ".join(failures) + "."
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    data = load_data()
    records = run_multistart(data)
    successful = [record for record in records if record["successful_solution"]]
    best = min(successful, key=lambda record: float(record["fit_rmse_v"]))
    best_parameters = _parameters_from_record(best)
    fitted_v, x_n, x_p = evaluate(best_parameters, data.z)
    residual_v = data.voltage_v - fitted_v
    fit_metrics = metrics(data.voltage_v[data.fit_mask], fitted_v[data.fit_mask])
    full_metrics = metrics(data.voltage_v, fitted_v)
    bin_centers, bin_means = interior_bin_statistics(data, residual_v)
    _, _, _, _, knot_voltages = unpack_parameters(best_parameters)

    good = [record for record in records if record["is_good_solution"]]
    endpoint_matrix = np.asarray(
        [
            [
                record["x_n_start"], record["x_n_end"],
                record["x_p_start"], record["x_p_end"]
            ]
            for record in good
        ],
        dtype=float,
    )
    endpoint_spans = np.ptp(endpoint_matrix, axis=0)
    good_curves = [evaluate(_parameters_from_record(record), data.z)[0] for record in good]
    largest_curve_difference_v = 0.0
    for first in range(len(good_curves)):
        for second in range(first + 1, len(good_curves)):
            difference = np.max(
                np.abs(
                    good_curves[first][data.fit_mask]
                    - good_curves[second][data.fit_mask]
                )
            )
            largest_curve_difference_v = max(largest_curve_difference_v, float(difference))

    xn_lower, xn_upper = MCMB_SUPPORTED_STOICHIOMETRY_RANGE
    xp_lower, xp_upper = NCA_WORKING_STOICHIOMETRY_RANGE
    voltage_lower, voltage_upper = UP_VOLTAGE_BOUNDS_V
    stoichiometry_distances = np.asarray(
        [
            best_parameters[0] - xn_lower, xn_upper - best_parameters[0],
            best_parameters[1] - xn_lower, xn_upper - best_parameters[1],
            best_parameters[2] - xp_lower, xp_upper - best_parameters[2],
            best_parameters[3] - xp_lower, xp_upper - best_parameters[3],
        ]
    )
    voltage_distances = np.minimum(
        knot_voltages - voltage_lower, voltage_upper - knot_voltages
    )
    flat_run = longest_flat_run(knot_voltages)
    diagnostics: dict[str, object] = {
        "stoichiometry_near_bound": bool(
            np.any(stoichiometry_distances <= NEAR_STOICHIOMETRY_BOUND)
        ),
        "voltage_near_bound": bool(np.any(voltage_distances <= NEAR_VOLTAGE_BOUND_V)),
        "longest_flat_run": flat_run,
        "long_flat_run": flat_run >= LONG_FLAT_RUN_ADJACENCIES,
        "largest_good_endpoint_span": float(np.max(endpoint_spans)),
        "largest_good_curve_difference_v": largest_curve_difference_v,
        "materially_different_windows": bool(
            np.max(endpoint_spans) >= MATERIAL_WINDOW_DIFFERENCE
        ),
    }
    checks = {
        "fit_rmse": fit_metrics.rmse_v <= MAX_FIT_RMSE_V,
        "fit_mae": fit_metrics.mae_v <= MAX_FIT_MAE_V,
        "fit_r_squared": fit_metrics.r_squared >= MIN_FIT_R_SQUARED,
        "residual_shape": float(np.max(np.abs(bin_means)))
        <= MAX_ABS_INTERIOR_BIN_MEAN_V
        and float(np.ptp(bin_means)) <= MAX_INTERIOR_BIN_MEAN_SPAN_V,
        "domains": is_feasible(best_parameters),
        "monotonic": bool(np.all(np.diff(knot_voltages) <= 1.0e-12)),
        "no_long_flat_run": not diagnostics["long_flat_run"],
        "voltage_bounds_not_active": not diagnostics["voltage_near_bound"],
    }
    decision = (
        "ACCEPTABLE OCP IDENTIFICATION"
        if all(checks.values())
        else "NOT ACCEPTABLE OCP IDENTIFICATION"
    )

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    write_fit_data(
        OUTPUT_DIRECTORY / f"{OUTPUT_STEM}_fit.csv",
        data, fitted_v, residual_v, x_n, x_p
    )
    write_multistart(OUTPUT_DIRECTORY / f"{OUTPUT_STEM}_multistart.csv", records)
    write_knots(OUTPUT_DIRECTORY / f"{OUTPUT_STEM}_knots.csv", knot_voltages)
    save_figures(
        data, fitted_v, residual_v, x_n, x_p, knot_voltages, bin_centers, bin_means
    )
    write_report(
        OUTPUT_DIRECTORY / f"{OUTPUT_STEM}_report.md",
        data, best, records, residual_v, fit_metrics, full_metrics,
        bin_means, decision, checks, diagnostics
    )
    print(f"{decision}: fit-region RMSE={1000 * fit_metrics.rmse_v:.3f} mV")
    print(OUTPUT_DIRECTORY / f"{OUTPUT_STEM}_report.md")


if __name__ == "__main__":
    main()
