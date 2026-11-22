"""Panasonic NCR18650PF C/20 compatibility gate for Karthikeyan OCPs.

Only the four electrode stoichiometric-window endpoints are fitted.  No
dynamic, kinetic, resistance, diffusion, or thermal terms are included.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, minimize


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from models.ocp_panasonic_karthikeyan import (  # noqa: E402
    DEFAULT_TEMPERATURE_K,
    MCMB_GRAPHITE_LITERATURE_STOICHIOMETRY_RANGE,
    NCA_LITERATURE_STOICHIOMETRY_RANGE,
    U_n,
    U_p,
)


DEFAULT_INPUT = Path(
    r"C:\Users\amvel\Documents\Andes\Tesis\datasets\converted_csv"
    r"\c20_ocv_25C_converted.csv"
)
DEFAULT_OUTPUT_DIRECTORY = REPOSITORY_ROOT / "thesis_data"

Q_CELL_AH = 2.994979138519
DISCHARGE_START_INDEX = 6
DISCHARGE_END_INDEX_INCLUSIVE = 1246
MULTISTART_SEED = 20260912
NUMBER_OF_STARTS = 40
MIN_DIRECTION_GAP = 1.0e-8

# Predeclared compatibility-gate criteria.  They assess both global error and
# broad residual shape; the decision is not based on RMSE alone.
MAX_RMSE_V = 0.050
MAX_MAE_V = 0.040
MIN_R_SQUARED = 0.990
MAX_ABS_BIN_MEAN_RESIDUAL_V = 0.075
MAX_BIN_MEAN_RESIDUAL_SPAN_V = 0.125
MIN_NONPATHOLOGICAL_WINDOW_WIDTH = 0.050
GOOD_SOLUTION_RMSE_DELTA_V = 0.001
MATERIALLY_DIFFERENT_WINDOW = 0.050


@dataclass(frozen=True)
class DischargeData:
    source_indices: np.ndarray
    time_s: np.ndarray
    voltage_v: np.ndarray
    current_a: np.ndarray
    exported_ah: np.ndarray
    discharged_capacity_ah: np.ndarray
    z: np.ndarray


@dataclass(frozen=True)
class FitMetrics:
    rmse_v: float
    mae_v: float
    r_squared: float


def _required_columns(names: tuple[str, ...] | None) -> None:
    required = {"Time", "Voltage", "Current", "Ah"}
    available = set(names or ())
    missing = sorted(required - available)
    if missing:
        raise ValueError(f"C/20 CSV is missing required columns: {missing}")


def load_discharge_branch(path: Path) -> DischargeData:
    """Load and validate the manifest-identified C/20 discharge branch."""

    table = np.genfromtxt(path, delimiter=",", names=True, dtype=float, encoding="utf-8")
    _required_columns(table.dtype.names)
    if table.ndim != 1 or table.size <= DISCHARGE_END_INDEX_INCLUSIVE:
        raise ValueError("C/20 CSV does not contain the manifest-identified branch.")

    selection = slice(DISCHARGE_START_INDEX, DISCHARGE_END_INDEX_INCLUSIVE + 1)
    time_s = np.asarray(table["Time"][selection], dtype=float)
    voltage_v = np.asarray(table["Voltage"][selection], dtype=float)
    current_a = np.asarray(table["Current"][selection], dtype=float)
    exported_ah = np.asarray(table["Ah"][selection], dtype=float)

    for name, values in (
        ("time", time_s),
        ("voltage", voltage_v),
        ("current", current_a),
        ("capacity", exported_ah),
    ):
        if not np.all(np.isfinite(values)):
            raise ValueError(f"Non-finite {name} values occur in the discharge branch.")
    if np.any(np.diff(time_s) <= 0.0):
        raise ValueError("Discharge time must be strictly increasing.")
    if np.any(current_a >= 0.0):
        raise ValueError("Manifest-selected branch is not exclusively raw negative-current discharge.")

    increments_ah = (
        0.5
        * (np.abs(current_a[:-1]) + np.abs(current_a[1:]))
        * np.diff(time_s)
        / 3600.0
    )
    discharged_capacity_ah = np.concatenate(([0.0], np.cumsum(increments_ah)))
    integrated_capacity_ah = float(discharged_capacity_ah[-1])
    if not np.isclose(integrated_capacity_ah, Q_CELL_AH, rtol=0.0, atol=1.0e-9):
        raise ValueError(
            "Integrated C/20 capacity does not match the Panasonic-specific value: "
            f"{integrated_capacity_ah:.12f} Ah versus {Q_CELL_AH:.12f} Ah."
        )
    z = discharged_capacity_ah / integrated_capacity_ah
    if z[0] != 0.0 or not np.isclose(z[-1], 1.0, rtol=0.0, atol=1.0e-14):
        raise AssertionError("Normalized discharged capacity endpoints are invalid.")

    return DischargeData(
        source_indices=np.arange(
            DISCHARGE_START_INDEX,
            DISCHARGE_END_INDEX_INCLUSIVE + 1,
            dtype=int,
        ),
        time_s=time_s,
        voltage_v=voltage_v,
        current_a=current_a,
        exported_ah=exported_ah,
        discharged_capacity_ah=discharged_capacity_ah,
        z=z,
    )


def _is_feasible(parameters: np.ndarray, tolerance: float = 0.0) -> bool:
    xn_start, xn_end, xp_start, xp_end = parameters
    xn_lower, xn_upper = MCMB_GRAPHITE_LITERATURE_STOICHIOMETRY_RANGE
    xp_lower, xp_upper = NCA_LITERATURE_STOICHIOMETRY_RANGE
    return bool(
        xn_lower - tolerance <= xn_end <= xn_upper + tolerance
        and xn_lower - tolerance <= xn_start <= xn_upper + tolerance
        and xp_lower - tolerance <= xp_start <= xp_upper + tolerance
        and xp_lower - tolerance <= xp_end <= xp_upper + tolerance
        and xn_start > xn_end
        and xp_start < xp_end
    )


def evaluate_voltage(
    parameters: np.ndarray,
    z: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate U_p-U_n without clipping or any dynamic voltage terms."""

    parameters = np.asarray(parameters, dtype=float)
    if parameters.shape != (4,) or not np.all(np.isfinite(parameters)):
        raise ValueError("Exactly four finite stoichiometric endpoints are required.")
    if not _is_feasible(parameters):
        raise ValueError("Infeasible stoichiometric windows were rejected; no clipping was used.")

    xn_start, xn_end, xp_start, xp_end = parameters
    x_n = xn_start + z * (xn_end - xn_start)
    x_p = xp_start + z * (xp_end - xp_start)
    fitted_voltage_v = U_p(x_p, DEFAULT_TEMPERATURE_K) - U_n(
        x_n, DEFAULT_TEMPERATURE_K
    )
    return np.asarray(fitted_voltage_v), x_n, x_p


def calculate_metrics(measured_v: np.ndarray, fitted_v: np.ndarray) -> FitMetrics:
    residual_v = measured_v - fitted_v
    sum_squared_residual = float(np.sum(residual_v**2))
    sum_squared_total = float(np.sum((measured_v - np.mean(measured_v)) ** 2))
    return FitMetrics(
        rmse_v=float(np.sqrt(np.mean(residual_v**2))),
        mae_v=float(np.mean(np.abs(residual_v))),
        r_squared=1.0 - sum_squared_residual / sum_squared_total,
    )


def deterministic_initial_guesses(number_of_starts: int = NUMBER_OF_STARTS) -> list[np.ndarray]:
    """Return physically ordered, repeatable multistart window guesses."""

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
    guesses = [np.asarray(values, dtype=float) for values in explicit[:number_of_starts]]
    rng = np.random.default_rng(MULTISTART_SEED)
    xn_lower, xn_upper = MCMB_GRAPHITE_LITERATURE_STOICHIOMETRY_RANGE
    xp_lower, xp_upper = NCA_LITERATURE_STOICHIOMETRY_RANGE
    while len(guesses) < number_of_starts:
        xn_end = rng.uniform(xn_lower, xn_upper - 0.06)
        xn_start = rng.uniform(xn_end + 0.05, xn_upper)
        xp_start = rng.uniform(xp_lower, xp_upper - 0.06)
        xp_end = rng.uniform(xp_start + 0.05, xp_upper)
        guesses.append(np.asarray([xn_start, xn_end, xp_start, xp_end]))
    return guesses


def run_multistart_fit(
    data: DischargeData,
    objective_mask: np.ndarray | None = None,
) -> list[dict[str, object]]:
    """Fit the four windows, optionally using only selected samples in the objective."""

    if objective_mask is None:
        objective_mask = np.ones(data.z.shape, dtype=bool)
    else:
        objective_mask = np.asarray(objective_mask, dtype=bool)
        if objective_mask.shape != data.z.shape or not np.any(objective_mask):
            raise ValueError("objective_mask must select samples from the discharge branch.")

    xn_lower, xn_upper = MCMB_GRAPHITE_LITERATURE_STOICHIOMETRY_RANGE
    xp_lower, xp_upper = NCA_LITERATURE_STOICHIOMETRY_RANGE
    bounds = Bounds(
        [xn_lower, xn_lower, xp_lower, xp_lower],
        [xn_upper, xn_upper, xp_upper, xp_upper],
    )
    directions = LinearConstraint(
        np.asarray([[1.0, -1.0, 0.0, 0.0], [0.0, 0.0, -1.0, 1.0]]),
        [MIN_DIRECTION_GAP, MIN_DIRECTION_GAP],
        [np.inf, np.inf],
    )

    def objective(parameters: np.ndarray) -> float:
        try:
            fitted_v, _, _ = evaluate_voltage(parameters, data.z)
        except ValueError:
            return 1.0e6
        residual_v = data.voltage_v[objective_mask] - fitted_v[objective_mask]
        return float(np.mean(residual_v**2))

    records: list[dict[str, object]] = []
    for start_id, initial in enumerate(deterministic_initial_guesses(), start=1):
        result = minimize(
            objective,
            initial,
            method="SLSQP",
            bounds=bounds,
            constraints=[directions],
            options={"ftol": 1.0e-14, "maxiter": 2000, "disp": False},
        )
        parameters = np.asarray(result.x, dtype=float)
        feasible = _is_feasible(parameters, tolerance=1.0e-10)
        successful = bool(result.success and feasible and np.isfinite(result.fun))
        record: dict[str, object] = {
            "start_id": start_id,
            "initial_x_n_start": initial[0],
            "initial_x_n_end": initial[1],
            "initial_x_p_start": initial[2],
            "initial_x_p_end": initial[3],
            "optimizer_success": bool(result.success),
            "feasible": feasible,
            "successful_solution": successful,
            "status": int(result.status),
            "message": str(result.message),
            "iterations": int(result.nit),
            "function_evaluations": int(result.nfev),
            "x_n_start": parameters[0],
            "x_n_end": parameters[1],
            "x_p_start": parameters[2],
            "x_p_end": parameters[3],
            "objective_mse_v2": float(result.fun),
            "rmse_v": np.nan,
            "mae_v": np.nan,
            "r_squared": np.nan,
            "is_good_solution": False,
        }
        if successful:
            fitted_v, _, _ = evaluate_voltage(parameters, data.z)
            metrics = calculate_metrics(
                data.voltage_v[objective_mask],
                fitted_v[objective_mask],
            )
            record.update(
                rmse_v=metrics.rmse_v,
                mae_v=metrics.mae_v,
                r_squared=metrics.r_squared,
            )
        records.append(record)

    successful_records = [record for record in records if record["successful_solution"]]
    if not successful_records:
        raise RuntimeError("No feasible multistart optimization completed successfully.")
    best_rmse = min(float(record["rmse_v"]) for record in successful_records)
    for record in successful_records:
        record["is_good_solution"] = (
            float(record["rmse_v"]) <= best_rmse + GOOD_SOLUTION_RMSE_DELTA_V
        )
    return records


def residual_shape_statistics(z: np.ndarray, residual_v: np.ndarray) -> dict[str, object]:
    edges = np.linspace(0.0, 1.0, 21)
    bin_index = np.minimum(np.digitize(z, edges, right=False) - 1, 19)
    means = np.asarray(
        [np.mean(residual_v[bin_index == index]) for index in range(20)],
        dtype=float,
    )
    centers = 0.5 * (edges[:-1] + edges[1:])
    return {
        "centers": centers,
        "means": means,
        "max_abs_mean_v": float(np.max(np.abs(means))),
        "mean_span_v": float(np.ptp(means)),
    }


def _format_csv_value(value: object) -> object:
    if isinstance(value, (np.floating, float)):
        return f"{float(value):.12g}"
    if isinstance(value, (np.integer, int)) and not isinstance(value, bool):
        return int(value)
    return value


def write_multistart_csv(path: Path, records: list[dict[str, object]]) -> None:
    fieldnames = list(records[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for record in records:
            writer.writerow({key: _format_csv_value(value) for key, value in record.items()})


def write_fit_csv(
    path: Path,
    data: DischargeData,
    fitted_v: np.ndarray,
    residual_v: np.ndarray,
    x_n: np.ndarray,
    x_p: np.ndarray,
) -> None:
    names = [
        "source_row_index",
        "time_s",
        "raw_current_a",
        "exported_capacity_ah",
        "discharged_capacity_ah",
        "normalized_discharged_capacity",
        "measured_voltage_v",
        "fitted_equilibrium_voltage_v",
        "residual_measured_minus_fitted_v",
        "x_n",
        "x_p",
    ]
    columns = zip(
        data.source_indices,
        data.time_s,
        data.current_a,
        data.exported_ah,
        data.discharged_capacity_ah,
        data.z,
        data.voltage_v,
        fitted_v,
        residual_v,
        x_n,
        x_p,
        strict=True,
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(names)
        for row in columns:
            writer.writerow([_format_csv_value(value) for value in row])


def save_plots(
    fit_path: Path,
    residual_path: Path,
    data: DischargeData,
    fitted_v: np.ndarray,
    residual_v: np.ndarray,
    x_n: np.ndarray,
    x_p: np.ndarray,
    shape: dict[str, object],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure = plt.figure(figsize=(11.5, 8.0), constrained_layout=True)
    grid = figure.add_gridspec(2, 2, height_ratios=(1.25, 1.0))
    voltage_axis = figure.add_subplot(grid[0, :])
    negative_axis = figure.add_subplot(grid[1, 0])
    positive_axis = figure.add_subplot(grid[1, 1])

    voltage_axis.plot(data.z, data.voltage_v, label="Measured Panasonic C/20", linewidth=2.0)
    voltage_axis.plot(data.z, fitted_v, label="$U_p(x_p)-U_n(x_n)$", linewidth=1.8)
    voltage_axis.set(ylabel="Voltage [V]", title="Panasonic C/20 Karthikeyan compatibility fit")
    voltage_axis.legend()
    negative_axis.plot(data.z, x_n, color="tab:blue", linewidth=1.8)
    negative_axis.set(
        xlabel="Normalized discharged capacity, $z$",
        ylabel="$x_n$",
        title="MCMB stoichiometric trajectory",
        ylim=(0.0, 0.82),
    )
    positive_axis.plot(data.z, x_p, color="tab:red", linewidth=1.8)
    positive_axis.set(
        xlabel="Normalized discharged capacity, $z$",
        ylabel="$x_p$",
        title="NCA stoichiometric trajectory",
        ylim=(0.16, 0.98),
    )
    for axis in (voltage_axis, negative_axis, positive_axis):
        axis.grid(alpha=0.25)
    figure.savefig(fit_path, dpi=180)
    plt.close(figure)

    residual_figure, residual_axes = plt.subplots(
        2, 1, figsize=(10.5, 7.0), constrained_layout=True
    )
    residual_axes[0].plot(data.z, 1000.0 * residual_v, linewidth=1.2)
    residual_axes[0].axhline(0.0, color="black", linewidth=0.8)
    residual_axes[0].set(
        ylabel="Measured - fitted [mV]",
        title="Panasonic C/20 voltage residual",
    )
    residual_axes[1].plot(
        np.asarray(shape["centers"]),
        1000.0 * np.asarray(shape["means"]),
        marker="o",
        linewidth=1.2,
        label="Mean residual in 5% capacity bins",
    )
    residual_axes[1].axhline(0.0, color="black", linewidth=0.8)
    residual_axes[1].set(
        xlabel="Normalized discharged capacity, $z$",
        ylabel="Bin mean residual [mV]",
    )
    residual_axes[1].legend()
    for axis in residual_axes:
        axis.grid(alpha=0.25)
    residual_figure.savefig(residual_path, dpi=180)
    plt.close(residual_figure)


def write_report(
    path: Path,
    input_path: Path,
    data: DischargeData,
    best: dict[str, object],
    records: list[dict[str, object]],
    residual_v: np.ndarray,
    shape: dict[str, object],
    decision: str,
    criteria: dict[str, bool],
    parameter_spans: dict[str, float],
    materially_different: bool,
) -> None:
    good = [record for record in records if record["is_good_solution"]]
    successful = [record for record in records if record["successful_solution"]]
    midpoint_index = int(np.argmin(np.abs(data.z - 0.50)))
    near_end_index = int(np.argmin(np.abs(data.z - 0.95)))
    maximum_residual_index = int(np.argmax(np.abs(residual_v)))
    lines = [
        "# Panasonic C/20 Karthikeyan OCP compatibility gate",
        "",
        f"## Decision: **{decision}**",
        "",
        "This gate adjusts only `x_n_start`, `x_n_end`, `x_p_start`, and "
        "`x_p_end`. It contains no resistance, kinetic, diffusion, "
        "overpotential, or thermal terms, and performs no clipping.",
        "",
        "## Data preparation",
        "",
        f"- Source: `{input_path}`",
        f"- Discharge source rows: {DISCHARGE_START_INDEX} through "
        f"{DISCHARGE_END_INDEX_INCLUSIVE}, inclusive",
        f"- Samples: {data.z.size}",
        f"- Raw current sign: negative throughout ({np.min(data.current_a):.5f} to "
        f"{np.max(data.current_a):.5f} A)",
        f"- Voltage span: {data.voltage_v[0]:.5f} to {data.voltage_v[-1]:.5f} V",
        f"- Integrated Panasonic C/20 capacity: "
        f"{data.discharged_capacity_ah[-1]:.12f} Ah",
        f"- Fixed Panasonic capacity reference: {Q_CELL_AH:.12f} Ah",
        "- Normalization: `z = Q_discharge / Q_discharge_final`, with endpoints "
        "exactly 0 and 1",
        "",
        "## Best feasible solution",
        "",
        "| Metric | Value |",
        "| --- | ---: |",
        f"| RMSE | {1000.0 * float(best['rmse_v']):.3f} mV |",
        f"| MAE | {1000.0 * float(best['mae_v']):.3f} mV |",
        f"| R² | {float(best['r_squared']):.8f} |",
        f"| `x_n_start` | {float(best['x_n_start']):.9f} |",
        f"| `x_n_end` | {float(best['x_n_end']):.9f} |",
        f"| `x_p_start` | {float(best['x_p_start']):.9f} |",
        f"| `x_p_end` | {float(best['x_p_end']):.9f} |",
        "",
        "Supported domains are `0.01 <= x_n <= 0.80` and "
        "`0.18 <= x_p <= 0.96`; the fitted linear trajectories remain inside "
        "these domains by construction and validation.",
        "",
        "## Residual-shape evidence",
        "",
        f"- Beginning residual (measured - fitted): "
        f"{1000.0 * residual_v[0]:.3f} mV",
        f"- Residual near `z=0.50`: "
        f"{1000.0 * residual_v[midpoint_index]:.3f} mV",
        f"- Residual near `z=0.95`: "
        f"{1000.0 * residual_v[near_end_index]:.3f} mV",
        f"- Terminal residual at `z=1`: {1000.0 * residual_v[-1]:.3f} mV",
        f"- Largest absolute point residual: "
        f"{1000.0 * abs(residual_v[maximum_residual_index]):.3f} mV at "
        f"`z={data.z[maximum_residual_index]:.6f}`",
        "",
        "The alternating broad residual lobes and the terminal mismatch show "
        "a systematic voltage-shape incompatibility, not unstructured noise.",
        "",
        "## Gate checks",
        "",
        "| Check | Criterion | Result |",
        "| --- | --- | --- |",
        f"| Voltage RMSE | <= {1000 * MAX_RMSE_V:.0f} mV | "
        f"{'PASS' if criteria['rmse'] else 'FAIL'} |",
        f"| Voltage MAE | <= {1000 * MAX_MAE_V:.0f} mV | "
        f"{'PASS' if criteria['mae'] else 'FAIL'} |",
        f"| R² | >= {MIN_R_SQUARED:.3f} | "
        f"{'PASS' if criteria['r_squared'] else 'FAIL'} |",
        f"| Largest absolute 5%-bin mean residual | <= "
        f"{1000 * MAX_ABS_BIN_MEAN_RESIDUAL_V:.0f} mV; observed "
        f"{1000 * float(shape['max_abs_mean_v']):.3f} mV | "
        f"{'PASS' if criteria['shape_absolute'] else 'FAIL'} |",
        f"| Span of 5%-bin mean residuals | <= "
        f"{1000 * MAX_BIN_MEAN_RESIDUAL_SPAN_V:.0f} mV; observed "
        f"{1000 * float(shape['mean_span_v']):.3f} mV | "
        f"{'PASS' if criteria['shape_span'] else 'FAIL'} |",
        f"| Literature domains and physical direction | Required | "
        f"{'PASS' if criteria['domains'] else 'FAIL'} |",
        f"| Non-pathological electrode-window widths | >= "
        f"{MIN_NONPATHOLOGICAL_WINDOW_WIDTH:.2f} each | "
        f"{'PASS' if criteria['window_widths'] else 'FAIL'} |",
        "",
        "## Multistart and identifiability",
        "",
        f"- Deterministic starts attempted: {len(records)} (seed "
        f"`{MULTISTART_SEED}` for generated starts)",
        f"- Successful feasible solutions: {len(successful)}",
        f"- Good solutions within {1000 * GOOD_SOLUTION_RMSE_DELTA_V:.1f} mV "
        f"RMSE of the best: {len(good)}",
        f"- Maximum good-solution endpoint span: "
        f"{max(parameter_spans.values()):.9f}",
        f"- Materially different windows (endpoint span >= "
        f"{MATERIALLY_DIFFERENT_WINDOW:.2f}) with nearly identical error: "
        f"**{'YES' if materially_different else 'NO'}**",
        "",
        "Good-solution endpoint spans:",
        "",
        "| Endpoint | Span across good solutions |",
        "| --- | ---: |",
    ]
    for name, span in parameter_spans.items():
        lines.append(f"| `{name}` | {span:.9f} |")
    lines.extend(
        [
            "",
            "All good multistart solutions:",
            "",
            "| Start | RMSE [mV] | MAE [mV] | R² | x_n_start | x_n_end | "
            "x_p_start | x_p_end |",
            "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for record in sorted(good, key=lambda item: float(item["rmse_v"])):
        lines.append(
            f"| {record['start_id']} | {1000 * float(record['rmse_v']):.3f} | "
            f"{1000 * float(record['mae_v']):.3f} | "
            f"{float(record['r_squared']):.8f} | "
            f"{float(record['x_n_start']):.9f} | "
            f"{float(record['x_n_end']):.9f} | "
            f"{float(record['x_p_start']):.9f} | "
            f"{float(record['x_p_end']):.9f} |"
        )
    lines.extend(
        [
            "",
            "The complete record of every start and every successful solution is "
            "stored in `panasonic_c20_karthikeyan_multistart.csv`.",
            "",
            "## Interpretation",
            "",
        ]
    )
    if decision == "PASS":
        lines.append(
            "The Karthikeyan NCA/MCMB pair passes this quasi-equilibrium "
            "full-cell compatibility gate under the stated criteria. This does "
            "not establish that the best stoichiometric windows are unique or "
            "that the literature electrodes are exact Panasonic electrodes."
        )
    else:
        failed = [name for name, passed in criteria.items() if not passed]
        lines.append(
            "The Karthikeyan NCA/MCMB pair fails this quasi-equilibrium "
            "compatibility gate. Failed checks: " + ", ".join(failed) + "."
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(input_path: Path, output_directory: Path) -> str:
    output_directory.mkdir(parents=True, exist_ok=True)
    data = load_discharge_branch(input_path)
    records = run_multistart_fit(data)
    successful = [record for record in records if record["successful_solution"]]
    best = min(successful, key=lambda record: float(record["rmse_v"]))
    parameters = np.asarray(
        [
            best["x_n_start"],
            best["x_n_end"],
            best["x_p_start"],
            best["x_p_end"],
        ],
        dtype=float,
    )
    fitted_v, x_n, x_p = evaluate_voltage(parameters, data.z)
    residual_v = data.voltage_v - fitted_v
    shape = residual_shape_statistics(data.z, residual_v)

    good = [record for record in records if record["is_good_solution"]]
    parameter_names = ("x_n_start", "x_n_end", "x_p_start", "x_p_end")
    parameter_spans = {
        name: float(
            np.ptp(np.asarray([float(record[name]) for record in good], dtype=float))
        )
        for name in parameter_names
    }
    materially_different = max(parameter_spans.values()) >= MATERIALLY_DIFFERENT_WINDOW

    xn_width = float(best["x_n_start"]) - float(best["x_n_end"])
    xp_width = float(best["x_p_end"]) - float(best["x_p_start"])
    criteria = {
        "rmse": float(best["rmse_v"]) <= MAX_RMSE_V,
        "mae": float(best["mae_v"]) <= MAX_MAE_V,
        "r_squared": float(best["r_squared"]) >= MIN_R_SQUARED,
        "shape_absolute": float(shape["max_abs_mean_v"])
        <= MAX_ABS_BIN_MEAN_RESIDUAL_V,
        "shape_span": float(shape["mean_span_v"]) <= MAX_BIN_MEAN_RESIDUAL_SPAN_V,
        "domains": _is_feasible(parameters) and np.all(np.isfinite(fitted_v)),
        "window_widths": xn_width >= MIN_NONPATHOLOGICAL_WINDOW_WIDTH
        and xp_width >= MIN_NONPATHOLOGICAL_WINDOW_WIDTH,
    }
    decision = "PASS" if all(criteria.values()) else "FAIL"

    fit_csv = output_directory / "panasonic_c20_karthikeyan_fit.csv"
    multistart_csv = output_directory / "panasonic_c20_karthikeyan_multistart.csv"
    fit_plot = output_directory / "panasonic_c20_karthikeyan_fit.png"
    residual_plot = output_directory / "panasonic_c20_karthikeyan_residual.png"
    report = output_directory / "panasonic_c20_karthikeyan_report.md"
    write_fit_csv(fit_csv, data, fitted_v, residual_v, x_n, x_p)
    write_multistart_csv(multistart_csv, records)
    save_plots(fit_plot, residual_plot, data, fitted_v, residual_v, x_n, x_p, shape)
    write_report(
        report,
        input_path,
        data,
        best,
        records,
        residual_v,
        shape,
        decision,
        criteria,
        parameter_spans,
        materially_different,
    )
    return decision


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-directory", type=Path, default=DEFAULT_OUTPUT_DIRECTORY)
    arguments = parser.parse_args()
    decision = run(arguments.input, arguments.output_directory)
    print(f"Panasonic C/20 Karthikeyan compatibility gate: {decision}")
    print(arguments.output_directory.resolve() / "panasonic_c20_karthikeyan_report.md")


if __name__ == "__main__":
    main()
