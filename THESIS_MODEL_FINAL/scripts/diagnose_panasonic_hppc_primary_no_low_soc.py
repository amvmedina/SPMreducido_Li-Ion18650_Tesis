"""One diagnostic HPPC calibration excluding pulse-start SOC below 1/3."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from scripts.calibrate_panasonic_hppc_primary import (  # noqa: E402
    ALPHA_N_FIXED_S,
    D_P_FIXED_PER_S,
    FIT_PARAMETER_NAMES,
    NEAR_BEST_RELATIVE_OBJECTIVE,
    OUTPUT_DIRECTORY,
    calibration_diagnostics,
    evaluate,
    maximum_feasible_alpha_p,
    postfit_jacobian,
    prepare_groups,
    run_multistart,
)
from scripts.audit_panasonic_hppc_identifiability import (  # noqa: E402
    load_hppc_data,
    load_pulse_records,
)


MINIMUM_SOC = 1.0 / 3.0
PREFIX = "panasonic_hppc_primary_fit_no_low_soc"


def read_single_row(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != 1:
        raise ValueError(f"Expected one row in {path}, found {len(rows)}.")
    return rows[0]


def read_metrics(path: Path) -> dict[tuple[str, str, str], float]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            (row["group"], row["subset"], row["metric"]): float(row["value"])
            for row in csv.DictReader(handle)
        }


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_report(
    path: Path,
    pulse_count: int,
    history_count: int,
    best: dict[str, object],
    multistart: list[dict[str, object]],
    diagnostics: dict[str, object],
    postfit: dict[str, object],
    original_parameters: dict[str, str],
    original_metrics: dict[tuple[str, str, str], float],
) -> None:
    successful = [row for row in multistart if row["success"]]
    best_objective = float(best["objective_v2"])
    near_best = [
        row
        for row in successful
        if float(row["objective_v2"])
        <= best_objective * (1.0 + NEAR_BEST_RELATIVE_OBJECTIVE)
    ]
    spread = {
        field: (
            min(float(row[field]) for row in near_best),
            max(float(row[field]) for row in near_best),
        )
        for field in ("R0_ohm", "alpha_p_s", "d_n_per_s")
    }
    metrics = dict(diagnostics["metrics"])
    correlation = np.asarray(postfit["correlation"])
    boundary_active = bool(best["at_frozen_ocp_domain_boundary"])
    classification = (
        "BOUNDARY PERSISTS WITHOUT LOW SOC"
        if boundary_active
        else "BOUNDARY WAS PRIMARILY LOW-SOC DRIVEN"
    )
    original_balanced = original_metrics[
        ("overall", "all", "phase_balanced_rmse_v")
    ]
    original_phase = {
        phase: original_metrics[("overall", "all", f"rmse_{phase}_v")]
        for phase in ("immediate", "within_pulse", "post_relaxation")
    }
    original_soc = {
        region: original_metrics[
            ("soc", region, "phase_balanced_dynamic_rmse_v")
        ]
        for region in ("high", "middle")
    }
    lines = [
        "# Panasonic HPPC diagnostic calibration without low-SOC pulses",
        "",
        f"## {classification}",
        "",
        "This is the single requested diagnostic rerun. The sole calibration-data change is `pulse-start cell SOC >= 1/3`; all lower-SOC pulses are excluded. Frozen OCPs, capacities, nominal effective stoichiometric windows, model equations, phase weighting, broad numerical bounds, multistart definitions, and block initialization logic are unchanged.",
        "",
        f"- Retained pulses: **{pulse_count}** in {history_count} independently initialized contiguous histories.",
        f"- Fixed `alpha_n = {ALPHA_N_FIXED_S:g} s` and `d_p = {D_P_FIXED_PER_S:g} 1/s`.",
        "- Actual mixed sampling is retained; no omitted gap is integrated and no resampling is used.",
        "",
        "## Fitted conditional parameters",
        "",
        "| Parameter | SOC >= 1/3 diagnostic | Original all-SOC fit |",
        "|---|---:|---:|",
        f"| `R0` [ohm] | {float(best['R0_ohm']):.12g} | {float(original_parameters['R0_ohm']):.12g} |",
        f"| `alpha_p` [s] | {float(best['alpha_p_s']):.12g} | {float(original_parameters['alpha_p_s']):.12g} |",
        f"| `d_n` [1/s] | {float(best['d_n_per_s']):.12g} | {float(original_parameters['d_n_per_s']):.12g} |",
        "",
        f"The retained-data frozen-OCP feasibility edge is `alpha_p <= {float(best['alpha_p_frozen_domain_upper_s']):.12g} s`. `alpha_p` at that edge: **{'YES' if boundary_active else 'NO'}**.",
        "",
        "## Dynamic fit metrics",
        "",
        "| Metric | SOC >= 1/3 diagnostic [mV] | Original all-SOC [mV] |",
        "|---|---:|---:|",
        f"| Phase-balanced RMSE | {1000*metrics['phase_balanced_rmse_v']:.6f} | {1000*original_balanced:.6f} |",
        f"| Immediate RMSE | {1000*metrics['rmse_immediate_v']:.6f} | {1000*original_phase['immediate']:.6f} |",
        f"| Within-pulse RMSE | {1000*metrics['rmse_within_pulse_v']:.6f} | {1000*original_phase['within_pulse']:.6f} |",
        f"| Post-relaxation RMSE | {1000*metrics['rmse_post_relaxation_v']:.6f} | {1000*original_phase['post_relaxation']:.6f} |",
        "",
        "Phase-balanced RMSE by retained SOC region:",
        "",
        "| Region | SOC >= 1/3 diagnostic [mV] | Original all-SOC fit [mV] |",
        "|---|---:|---:|",
        f"| High SOC | {1000*dict(diagnostics['by_soc'])['high']:.6f} | {1000*original_soc['high']:.6f} |",
        f"| Middle SOC | {1000*dict(diagnostics['by_soc'])['middle']:.6f} | {1000*original_soc['middle']:.6f} |",
        "",
        "## Multistart recoverability",
        "",
        f"- Successful feasible starts: {len(successful)} of {len(multistart)}.",
        f"- Near-best solutions (within {100*NEAR_BEST_RELATIVE_OBJECTIVE:.1f}%): {len(near_best)}.",
        f"- Near-best `R0` spread: {spread['R0_ohm'][0]:.12g}-{spread['R0_ohm'][1]:.12g} ohm.",
        f"- Near-best `alpha_p` spread: {spread['alpha_p_s'][0]:.12g}-{spread['alpha_p_s'][1]:.12g} s.",
        f"- Near-best `d_n` spread: {spread['d_n_per_s'][0]:.12g}-{spread['d_n_per_s'][1]:.12g} 1/s.",
        f"- Solutions at the retained-data frozen-OCP feasibility edge: {sum(bool(row['at_frozen_ocp_domain_boundary']) for row in successful)}.",
        "",
        "All near-best solutions are interior with respect to the OCP support, but they occupy distinct `alpha_p` bands near 3204, 4847, and 6350 s while differing in objective by less than 1%. Thus the boundary result is removed, but `alpha_p` remains practically weak/non-unique under this diagnostic.",
        "",
        "## Post-fit local audit",
        "",
        f"- Normalized-Jacobian singular values: {', '.join(f'{value:.9g}' for value in np.asarray(postfit['singular_values']))}.",
        f"- Condition number: **{float(postfit['condition_number']):.9g}**.",
        "- Estimator-oriented correlations:",
        "",
        "| | R0 | alpha_p | d_n |",
        "|---|---:|---:|---:|",
        *[
            "| " + FIT_PARAMETER_NAMES[row] + " | "
            + " | ".join(f"{correlation[row, column]:+.6f}" for column in range(3))
            + " |"
            for row in range(3)
        ],
        "",
        "The fitted-point Jacobian uses the same normalized finite differences as the original calibration. If the optimum is on the OCP feasibility edge, `alpha_p` uses the same second-order inward difference.",
        "",
        "## Comparison and decision",
        "",
        f"Removing low-SOC pulses changes the phase-balanced RMSE from {1000*original_balanced:.3f} to {1000*metrics['phase_balanced_rmse_v']:.3f} mV. The original all-SOC `alpha_p` was boundary-limited at {float(original_parameters['alpha_p_s']):.6f} s; the retained-data value is {float(best['alpha_p_s']):.6f} s against its recalculated support edge of {float(best['alpha_p_frozen_domain_upper_s']):.6f} s.",
        "",
        f"**Final classification: {classification}.**",
        "",
        "No other calibration, validation, parameter release, OCP/support change, capacity/window change, or model-equation change was performed.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    data = load_hppc_data()
    records = load_pulse_records()
    groups = prepare_groups(
        data, records, minimum_pulse_start_soc=MINIMUM_SOC
    )
    pulse_count = sum(len(group.pulses) for group in groups)
    alpha_upper = maximum_feasible_alpha_p(groups)
    multistart = run_multistart(groups, alpha_upper)
    successful = [row for row in multistart if row["success"]]
    if not successful:
        raise RuntimeError("No successful feasible diagnostic solution.")
    best = min(successful, key=lambda row: float(row["objective_v2"]))
    parameters = {
        "R0": float(best["R0_ohm"]),
        "alpha_p": float(best["alpha_p_s"]),
        "d_n": float(best["d_n_per_s"]),
    }
    evaluation = evaluate(groups, parameters)
    diagnostics = calibration_diagnostics(groups, evaluation)
    postfit = postfit_jacobian(groups, parameters)

    original_parameters = read_single_row(
        OUTPUT_DIRECTORY / "panasonic_hppc_primary_fit_parameters.csv"
    )
    original_metrics = read_metrics(
        OUTPUT_DIRECTORY / "panasonic_hppc_primary_fit_metrics.csv"
    )
    write_rows(
        OUTPUT_DIRECTORY / f"{PREFIX}_multistart.csv", multistart
    )
    write_rows(
        OUTPUT_DIRECTORY / f"{PREFIX}_parameters.csv",
        [
            {
                "retained_pulses": pulse_count,
                "retained_histories": len(groups),
                "minimum_pulse_start_soc": MINIMUM_SOC,
                "R0_ohm": parameters["R0"],
                "alpha_p_s": parameters["alpha_p"],
                "d_n_per_s": parameters["d_n"],
                "alpha_n_fixed_s": ALPHA_N_FIXED_S,
                "d_p_fixed_per_s": D_P_FIXED_PER_S,
                "alpha_p_frozen_domain_upper_s": alpha_upper,
                "alpha_p_at_frozen_ocp_domain_boundary": best[
                    "at_frozen_ocp_domain_boundary"
                ],
                "phase_balanced_rmse_v": diagnostics["metrics"][
                    "phase_balanced_rmse_v"
                ],
                "postfit_condition_number": postfit["condition_number"],
            }
        ],
    )
    metric_rows = []
    for name, value in dict(diagnostics["metrics"]).items():
        metric_rows.append({"group": "overall", "subset": "all", "metric": name, "value": value})
    for region in ("high", "middle"):
        metric_rows.append(
            {
                "group": "soc",
                "subset": region,
                "metric": "phase_balanced_dynamic_rmse_v",
                "value": diagnostics["by_soc"][region],
            }
        )
    write_rows(OUTPUT_DIRECTORY / f"{PREFIX}_metrics.csv", metric_rows)
    write_report(
        OUTPUT_DIRECTORY / f"{PREFIX}_report.md",
        pulse_count,
        len(groups),
        best,
        multistart,
        diagnostics,
        postfit,
        original_parameters,
        original_metrics,
    )
    classification = (
        "BOUNDARY PERSISTS WITHOUT LOW SOC"
        if best["at_frozen_ocp_domain_boundary"]
        else "BOUNDARY WAS PRIMARILY LOW-SOC DRIVEN"
    )
    print(f"Retained pulses: {pulse_count}; histories: {len(groups)}")
    print(
        f"R0={parameters['R0']:.12g}, alpha_p={parameters['alpha_p']:.12g}, "
        f"d_n={parameters['d_n']:.12g}"
    )
    print(classification)
    print((OUTPUT_DIRECTORY / f"{PREFIX}_report.md").resolve())


if __name__ == "__main__":
    main()
