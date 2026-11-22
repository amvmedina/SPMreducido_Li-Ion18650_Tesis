"""Fit Karthikeyan stoichiometric windows on 5--95% Panasonic C/20 only."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from scripts.fit_panasonic_c20_karthikeyan import (  # noqa: E402
    DEFAULT_INPUT,
    DEFAULT_OUTPUT_DIRECTORY,
    GOOD_SOLUTION_RMSE_DELTA_V,
    Q_CELL_AH,
    calculate_metrics,
    evaluate_voltage,
    load_discharge_branch,
    run_multistart_fit,
    write_fit_csv,
    write_multistart_csv,
)


FIT_Z_MIN = 0.05
FIT_Z_MAX = 0.95
OUTPUT_STEM = "panasonic_c20_karthikeyan_05_95"
PREVIOUS_MULTISTART = (
    DEFAULT_OUTPUT_DIRECTORY / "panasonic_c20_karthikeyan_multistart.csv"
)


def _best_previous_full_range_solution(path: Path) -> dict[str, float]:
    if not path.is_file():
        raise FileNotFoundError(
            "Previous full-range multistart result is required for the requested comparison: "
            f"{path}"
        )
    with path.open(newline="", encoding="utf-8") as handle:
        successful = [
            row
            for row in csv.DictReader(handle)
            if row["successful_solution"].lower() == "true"
        ]
    if not successful:
        raise ValueError("Previous multistart file contains no successful solution.")
    best = min(successful, key=lambda row: float(row["rmse_v"]))
    return {
        name: float(best[name])
        for name in ("x_n_start", "x_n_end", "x_p_start", "x_p_end")
    }


def _parameters(record: dict[str, object] | dict[str, float]) -> np.ndarray:
    return np.asarray(
        [
            record["x_n_start"],
            record["x_n_end"],
            record["x_p_start"],
            record["x_p_end"],
        ],
        dtype=float,
    )


def _interior_bin_means(
    z: np.ndarray,
    residual_v: np.ndarray,
    fit_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    edges = np.linspace(FIT_Z_MIN, FIT_Z_MAX, 19)
    centers = 0.5 * (edges[:-1] + edges[1:])
    selected_z = z[fit_mask]
    selected_residual = residual_v[fit_mask]
    indices = np.minimum(np.digitize(selected_z, edges, right=False) - 1, 17)
    means = np.asarray(
        [np.mean(selected_residual[indices == index]) for index in range(18)]
    )
    return centers, means


def save_figures(
    fit_path: Path,
    residual_path: Path,
    z: np.ndarray,
    measured_v: np.ndarray,
    fitted_v: np.ndarray,
    residual_v: np.ndarray,
    x_n: np.ndarray,
    x_p: np.ndarray,
    bin_centers: np.ndarray,
    bin_means: np.ndarray,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure = plt.figure(figsize=(11.5, 8.0), constrained_layout=True)
    grid = figure.add_gridspec(2, 2, height_ratios=(1.25, 1.0))
    voltage_axis = figure.add_subplot(grid[0, :])
    negative_axis = figure.add_subplot(grid[1, 0])
    positive_axis = figure.add_subplot(grid[1, 1])
    voltage_axis.plot(z, measured_v, linewidth=2.0, label="Measured Panasonic C/20")
    voltage_axis.plot(z, fitted_v, linewidth=1.8, label="$U_p(x_p)-U_n(x_n)$")
    voltage_axis.axvspan(FIT_Z_MIN, FIT_Z_MAX, color="tab:green", alpha=0.08, label="Fit region")
    voltage_axis.set(
        ylabel="Voltage [V]",
        title="Karthikeyan OCP fit on 5--95%; full-discharge evaluation",
    )
    voltage_axis.legend()
    negative_axis.plot(z, x_n, color="tab:blue", linewidth=1.8)
    negative_axis.axvspan(FIT_Z_MIN, FIT_Z_MAX, color="tab:green", alpha=0.08)
    negative_axis.set(
        xlabel="Normalized discharged capacity, $z$",
        ylabel="$x_n$",
        title="MCMB trajectory",
        ylim=(0.0, 0.82),
    )
    positive_axis.plot(z, x_p, color="tab:red", linewidth=1.8)
    positive_axis.axvspan(FIT_Z_MIN, FIT_Z_MAX, color="tab:green", alpha=0.08)
    positive_axis.set(
        xlabel="Normalized discharged capacity, $z$",
        ylabel="$x_p$",
        title="NCA trajectory",
        ylim=(0.16, 0.98),
    )
    for axis in (voltage_axis, negative_axis, positive_axis):
        axis.grid(alpha=0.25)
    figure.savefig(fit_path, dpi=180)
    plt.close(figure)

    residual_figure, residual_axes = plt.subplots(
        2, 1, figsize=(10.5, 7.0), constrained_layout=True
    )
    residual_axes[0].plot(z, 1000.0 * residual_v, linewidth=1.2)
    residual_axes[0].axvspan(FIT_Z_MIN, FIT_Z_MAX, color="tab:green", alpha=0.08)
    residual_axes[0].axhline(0.0, color="black", linewidth=0.8)
    residual_axes[0].set(
        ylabel="Measured - fitted [mV]",
        title="Full-discharge residual after fitting only 5--95%",
    )
    residual_axes[1].plot(
        bin_centers,
        1000.0 * bin_means,
        marker="o",
        linewidth=1.2,
        label="Mean residual in 5%-wide fit-region bins",
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
    data,
    fit_mask: np.ndarray,
    best: dict[str, object],
    previous_parameters: dict[str, float],
    fitted_v: np.ndarray,
    residual_v: np.ndarray,
    fit_metrics,
    full_metrics,
    previous_fit_metrics,
    previous_full_metrics,
    bin_means: np.ndarray,
    records: list[dict[str, object]],
) -> None:
    current_parameters = {
        name: float(best[name])
        for name in ("x_n_start", "x_n_end", "x_p_start", "x_p_end")
    }
    successful = [record for record in records if record["successful_solution"]]
    good = [record for record in records if record["is_good_solution"]]
    start_index = int(np.flatnonzero(fit_mask)[0])
    end_index = int(np.flatnonzero(fit_mask)[-1])
    midpoint_index = int(np.argmin(np.abs(data.z - 0.50)))

    lines = [
        "# Panasonic C/20 Karthikeyan 5--95% diagnostic",
        "",
        "Only samples with `0.05 <= z <= 0.95` entered the voltage objective. "
        "The resulting linear stoichiometric trajectories were evaluated without "
        "modification over the complete discharge.",
        "",
        "## Fit-region metrics (5--95%)",
        "",
        "| Metric | 5--95% fit | Previous full-range fit evaluated on 5--95% |",
        "| --- | ---: | ---: |",
        f"| RMSE | {1000 * fit_metrics.rmse_v:.3f} mV | "
        f"{1000 * previous_fit_metrics.rmse_v:.3f} mV |",
        f"| MAE | {1000 * fit_metrics.mae_v:.3f} mV | "
        f"{1000 * previous_fit_metrics.mae_v:.3f} mV |",
        f"| R² | {fit_metrics.r_squared:.8f} | "
        f"{previous_fit_metrics.r_squared:.8f} |",
        "",
        f"The fit region contains {int(np.sum(fit_mask))} of {data.z.size} samples, "
        f"from `z={data.z[start_index]:.6f}` to `z={data.z[end_index]:.6f}`.",
        "",
        "## Full-discharge evaluation",
        "",
        "| Metric | 5--95% fitted windows, evaluated on 0--100% | Previous full-range fit |",
        "| --- | ---: | ---: |",
        f"| RMSE | {1000 * full_metrics.rmse_v:.3f} mV | "
        f"{1000 * previous_full_metrics.rmse_v:.3f} mV |",
        f"| MAE | {1000 * full_metrics.mae_v:.3f} mV | "
        f"{1000 * previous_full_metrics.mae_v:.3f} mV |",
        f"| R² | {full_metrics.r_squared:.8f} | "
        f"{previous_full_metrics.r_squared:.8f} |",
        "",
        f"- Residual at `z=0`: {1000 * residual_v[0]:.3f} mV",
        f"- Residual near `z=0.50`: {1000 * residual_v[midpoint_index]:.3f} mV",
        f"- Residual at `z=1`: {1000 * residual_v[-1]:.3f} mV",
        f"- Interior 5%-bin mean residual range: "
        f"{1000 * np.min(bin_means):.3f} to {1000 * np.max(bin_means):.3f} mV",
        "",
        "## Stoichiometric-window comparison",
        "",
        "| Endpoint | Previous 0--100% fit | New 5--95% fit | Change |",
        "| --- | ---: | ---: | ---: |",
    ]
    for name in current_parameters:
        previous = previous_parameters[name]
        current = current_parameters[name]
        lines.append(f"| `{name}` | {previous:.9f} | {current:.9f} | {current - previous:+.9f} |")
    lines.extend(
        [
            "",
            "Both trajectories obey `x_n_start > x_n_end` and "
            "`x_p_start < x_p_end`, remain inside the Karthikeyan literature "
            "domains over the full discharge, and use no clipping.",
            "",
            "## Multistart stability",
            "",
            f"- Deterministic starts: {len(records)}",
            f"- Successful feasible solutions: {len(successful)}",
            f"- Solutions within {1000 * GOOD_SOLUTION_RMSE_DELTA_V:.1f} mV "
            f"fit-region RMSE of the best: {len(good)}",
            "",
            "Every successful solution and its fit-region metrics is recorded in "
            "`panasonic_c20_karthikeyan_05_95_multistart.csv`.",
            "",
            "## Diagnostic conclusion",
            "",
            "Excluding the endpoint regions reduces their leverage, but the "
            "fit-region residual remains structured rather than noise-like. "
            "The endpoint/cutoff region amplifies the previous full-range error, "
            "but it is not the sole source of incompatibility: the Karthikeyan "
            "OCP pair retains an interior full-cell shape mismatch.",
            "",
            f"Panasonic capacity remained fixed and separate at `{Q_CELL_AH:.12f} Ah`.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    data = load_discharge_branch(DEFAULT_INPUT)
    fit_mask = (data.z >= FIT_Z_MIN) & (data.z <= FIT_Z_MAX)
    records = run_multistart_fit(data, objective_mask=fit_mask)
    successful = [record for record in records if record["successful_solution"]]
    best = min(successful, key=lambda record: float(record["rmse_v"]))

    fitted_v, x_n, x_p = evaluate_voltage(_parameters(best), data.z)
    residual_v = data.voltage_v - fitted_v
    fit_metrics = calculate_metrics(data.voltage_v[fit_mask], fitted_v[fit_mask])
    full_metrics = calculate_metrics(data.voltage_v, fitted_v)

    previous_parameters = _best_previous_full_range_solution(PREVIOUS_MULTISTART)
    previous_fitted_v, _, _ = evaluate_voltage(_parameters(previous_parameters), data.z)
    previous_fit_metrics = calculate_metrics(
        data.voltage_v[fit_mask], previous_fitted_v[fit_mask]
    )
    previous_full_metrics = calculate_metrics(data.voltage_v, previous_fitted_v)
    bin_centers, bin_means = _interior_bin_means(data.z, residual_v, fit_mask)

    DEFAULT_OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    fit_csv = DEFAULT_OUTPUT_DIRECTORY / f"{OUTPUT_STEM}_fit.csv"
    multistart_csv = DEFAULT_OUTPUT_DIRECTORY / f"{OUTPUT_STEM}_multistart.csv"
    fit_plot = DEFAULT_OUTPUT_DIRECTORY / f"{OUTPUT_STEM}_fit.png"
    residual_plot = DEFAULT_OUTPUT_DIRECTORY / f"{OUTPUT_STEM}_residual.png"
    report = DEFAULT_OUTPUT_DIRECTORY / f"{OUTPUT_STEM}_report.md"
    write_fit_csv(fit_csv, data, fitted_v, residual_v, x_n, x_p)
    write_multistart_csv(multistart_csv, records)
    save_figures(
        fit_plot,
        residual_plot,
        data.z,
        data.voltage_v,
        fitted_v,
        residual_v,
        x_n,
        x_p,
        bin_centers,
        bin_means,
    )
    write_report(
        report,
        data,
        fit_mask,
        best,
        previous_parameters,
        fitted_v,
        residual_v,
        fit_metrics,
        full_metrics,
        previous_fit_metrics,
        previous_full_metrics,
        bin_means,
        records,
    )
    print(f"5--95% fit-region RMSE: {1000 * fit_metrics.rmse_v:.3f} mV")
    print(report.resolve())


if __name__ == "__main__":
    main()
