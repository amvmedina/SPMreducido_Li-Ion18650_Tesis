"""Build the own-data case-study reference in capacity/SOC coordinates only.

This script does not run an electrochemical or thermal model, fit parameters,
correct measured voltage, or interpret the loaded 1C curve as OCV.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from analysis.experimental_data_audit import audit_experiment, extract_rest_to_load_transitions  # noqa: E402
from analysis.own_case_study_reference import (  # noqa: E402
    SOC_GRID,
    charge_curve,
    discharge_curve,
    interpolate_on_soc,
    regional_voltage_statistics,
    relative_series_path_detector,
    segment_complete_cycle,
    segment_report_row,
    voltage_statistics,
)
from scripts.audit_own_trc18650_data import classify_and_augment, load_all  # noqa: E402


REFERENCE_ID = "cycle_cell1_ambient"
AFFECTED_CYCLE_IDS = (
    "cycle_cell1_30c",
    "cycle_cell1_35c",
    "cycle_cell2_undeclared",
    "cycle_cell2_30c",
)
PULSE_IDS = ("pulse_cell1_part1", "pulse_cell1_part2")
OUTPUT_DIRECTORY = REPOSITORY_ROOT / "thesis_data"
FIGURE_DIRECTORY = OUTPUT_DIRECTORY / "own_case_study_figures"
REFERENCE_CSV = OUTPUT_DIRECTORY / "own_case_study_ambient_discharge_reference.csv"
DIRECTIONAL_CSV = OUTPUT_DIRECTORY / "own_case_study_directional_voltage.csv"
METRICS_CSV = OUTPUT_DIRECTORY / "own_case_study_cross_experiment_metrics.csv"
REPORT_PATH = OUTPUT_DIRECTORY / "own_case_study_nonideality_report.md"


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def source_data_line_offset(path: Path) -> int:
    """Return the one-based source line of raw data row zero."""

    with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
        for index, line in enumerate(handle, start=1):
            if "Voltage(V)" in line and "Current(A)" in line:
                return index + 1
    raise ValueError(f"Could not locate Chroma header in {path}")


def reference_rows(experiment, curve, interpolated) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    first_data_line = source_data_line_offset(experiment.metadata.path)
    temp_names = [name for name in experiment.numeric if name.startswith("Temp")]
    for local, raw_row in enumerate(np.asarray(curve["rows"], dtype=int)):
        row: dict[str, object] = {
            "record_type": "raw_source_point",
            "is_interpolated": False,
            "reference_soc": float(np.asarray(curve["soc"])[local]),
            "normalized_discharged_throughput": float(np.asarray(curve["z"])[local]),
            "discharged_capacity_ah": float(np.asarray(curve["throughput_ah"])[local]),
            "voltage_v": float(np.asarray(curve["voltage_v"])[local]),
            "source_file": str(experiment.metadata.path.resolve()),
            "source_row_index_zero_based": int(raw_row),
            "source_csv_line_one_based": int(first_data_line + raw_row),
            "source_row_left_zero_based": int(raw_row),
            "source_row_right_zero_based": int(raw_row),
            "soc_label": "coulomb-counted reference SOC",
            "curve_label": "clean ambient 1C discharge reference curve",
        }
        for name in temp_names:
            row[f"{name}_c"] = float(np.asarray(curve[name])[local])
        rows.append(row)
    for local, soc in enumerate(np.asarray(interpolated["soc"])):
        left = int(np.asarray(interpolated["source_row_left"])[local])
        right = int(np.asarray(interpolated["source_row_right"])[local])
        row = {
            "record_type": "uniform_soc_interpolation",
            "is_interpolated": True,
            "reference_soc": float(soc),
            "normalized_discharged_throughput": float(1.0 - soc),
            "discharged_capacity_ah": float((1.0 - soc) * float(curve["capacity_ah"])),
            "voltage_v": float(np.asarray(interpolated["voltage_v"])[local]),
            "source_file": str(experiment.metadata.path.resolve()),
            "source_row_index_zero_based": "",
            "source_csv_line_one_based": "",
            "source_row_left_zero_based": left,
            "source_row_right_zero_based": right,
            "soc_label": "coulomb-counted reference SOC",
            "curve_label": "clean ambient 1C discharge reference curve",
        }
        for name in temp_names:
            row[f"{name}_c"] = float(np.asarray(interpolated[name])[local])
        rows.append(row)
    return rows


def region_label(soc: float) -> str:
    if soc >= 2.0 / 3.0:
        return "high"
    if soc >= 1.0 / 3.0:
        return "middle"
    return "low"


def fmt(value: object, digits: int = 6) -> str:
    if isinstance(value, (float, np.floating)):
        return "not applicable" if not np.isfinite(float(value)) else f"{float(value):.{digits}f}"
    return str(value)


def interpretation(dataset_id: str) -> str:
    if dataset_id == REFERENCE_ID:
        return "CLEAN REFERENCE"
    if dataset_id in PULSE_IDS:
        return "ELEVATED SERIES-PATH RESPONSE"
    if dataset_id == "cycle_cell2_undeclared":
        return "ELEVATED RESPONSE + CELL-TO-CELL DIFFERENCE"
    if dataset_id == "cycle_cell2_30c":
        return "ELEVATED RESPONSE + TEMPERATURE DIFFERENCE + CELL-TO-CELL DIFFERENCE"
    if dataset_id in {"cycle_cell1_30c", "cycle_cell1_35c"}:
        return "ELEVATED RESPONSE + TEMPERATURE DIFFERENCE"
    return "AMBIGUOUS"


def thermal_status(experiment) -> str:
    roles = list(experiment.metadata.temperature_roles.values())
    if any("surface" in role for role in roles) and any(
        "ambient" in role or "chamber" in role for role in roles
    ):
        return "two identified surface sensors plus chamber/ambient"
    return "temperature channels present; placement/ambient undocumented"


def make_figures(experiments, curves, interps, transitions, directional) -> None:
    colors = {
        "cycle_cell1_ambient": "#000000",
        "cycle_cell1_30c": "#D55E00",
        "cycle_cell1_35c": "#CC79A7",
        "cycle_cell2_undeclared": "#0072B2",
        "cycle_cell2_30c": "#009E73",
    }
    reference = experiments[REFERENCE_ID]
    reference_curve = curves[REFERENCE_ID]
    reference_interp = interps[REFERENCE_ID]
    charge = curves["ambient_charge"]
    charge_interp = interps["ambient_charge"]

    figure, axis = plt.subplots(figsize=(8.5, 5.3), constrained_layout=True)
    axis.plot(reference_curve["soc"], reference_curve["voltage_v"], color="#0072B2", lw=1.2, label="1C discharge")
    axis.plot(charge["soc"], charge["voltage_v"], color="#D55E00", lw=1.2, label="CC-CV charge")
    axis.set(xlabel="Coulomb-counted reference SOC", ylabel="Measured voltage [V]", title="Clean ambient cell-1 cycle (loaded directional curves)")
    axis.grid(alpha=0.25)
    axis.legend()
    figure.savefig(FIGURE_DIRECTORY / "ambient_clean_cycle_voltage_vs_soc.png", dpi=180)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(9.2, 5.5), constrained_layout=True)
    for dataset_id in (REFERENCE_ID, *AFFECTED_CYCLE_IDS):
        curve = curves[dataset_id]
        axis.plot(curve["soc"], curve["voltage_v"], lw=1.1, color=colors[dataset_id], label=dataset_id)
    axis.set(xlabel="Coulomb-counted reference SOC", ylabel="Measured discharge voltage [V]", title="Own complete-cycle 1C discharge curves")
    axis.grid(alpha=0.25)
    axis.legend(fontsize=8)
    figure.savefig(FIGURE_DIRECTORY / "own_cycles_discharge_voltage_vs_soc.png", dpi=180)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(9.2, 5.5), constrained_layout=True)
    for dataset_id in AFFECTED_CYCLE_IDS:
        delta = np.asarray(interps[dataset_id]["voltage_v"]) - np.asarray(reference_interp["voltage_v"])
        axis.plot(SOC_GRID, 1000.0 * delta, lw=1.1, color=colors[dataset_id], label=dataset_id)
    axis.axhline(0.0, color="black", lw=0.8)
    axis.set(xlabel="Coulomb-counted reference SOC", ylabel="Affected - clean reference [mV]", title="Voltage deviation from clean ambient cell-1 discharge")
    axis.grid(alpha=0.25)
    axis.legend(fontsize=8)
    figure.savefig(FIGURE_DIRECTORY / "own_cycles_voltage_deviation_vs_soc.png", dpi=180)
    plt.close(figure)

    ids = [experiment.metadata.dataset_id for experiment in experiments.values()]
    values = [[float(row["r_app_ohm"]) for row in transitions if row["dataset_id"] == dataset_id] for dataset_id in ids]
    figure, axis = plt.subplots(figsize=(11.8, 6.3), constrained_layout=True)
    boxes = axis.boxplot(values, tick_labels=ids, patch_artist=True, showmeans=True)
    for box, dataset_id in zip(boxes["boxes"], ids, strict=True):
        box.set_facecolor("#56B4E9" if dataset_id == REFERENCE_ID else "#E69F00")
        box.set_alpha(0.68)
    for x, dataset_id in enumerate(ids, start=1):
        for row in (item for item in transitions if item["dataset_id"] == dataset_id):
            axis.scatter(x, row["r_app_ohm"], marker="v" if row["direction"] == "discharge" else "^", color="black", s=26, zorder=4)
    axis.set(ylabel="Apparent rest-to-load series resistance [ohm]", title="Relative series-path diagnostic (not intrinsic cell R0)")
    axis.tick_params(axis="x", rotation=28)
    axis.grid(alpha=0.25, axis="y")
    figure.savefig(FIGURE_DIRECTORY / "own_apparent_resistance_comparison.png", dpi=180)
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(12.0, 5.0), sharey=True, constrained_layout=True)
    for name in (column for column in reference.numeric if column.startswith("Temp")):
        role = reference.metadata.temperature_roles.get(name, "unknown")
        axes[0].plot(reference_curve["soc"], reference_curve[name], lw=1.0, label=f"{name}: {role}")
        axes[1].plot(charge["soc"], charge[name], lw=1.0, label=f"{name}: {role}")
    axes[0].set_title("1C discharge")
    axes[1].set_title("CC-CV recharge")
    for axis in axes:
        axis.set_xlabel("Coulomb-counted reference SOC")
        axis.grid(alpha=0.25)
        axis.legend(fontsize=7)
    axes[0].set_ylabel("Measured temperature [degC]")
    figure.suptitle("Clean ambient-cycle temperature in capacity coordinates")
    figure.savefig(FIGURE_DIRECTORY / "ambient_clean_cycle_temperature_vs_soc.png", dpi=180)
    plt.close(figure)


def build_report(
    experiments,
    audits,
    segment_rows,
    curves,
    directional_stats,
    directional_regions,
    metrics,
    pulse_rows,
) -> None:
    reference = experiments[REFERENCE_ID]
    discharge = curves[REFERENCE_ID]
    charge = curves["ambient_charge"]
    discharge_segment = next(segment for segment in segment_complete_cycle(reference) if segment.name == "main_1c_discharge")
    post_rest = next(segment for segment in segment_complete_cycle(reference) if segment.name == "post_discharge_rest")
    temp_names = [name for name in reference.numeric if name.startswith("Temp")]
    surface_names = [name for name in temp_names if "surface" in reference.metadata.temperature_roles.get(name, "")]
    ambient_names = [name for name in temp_names if "ambient" in reference.metadata.temperature_roles.get(name, "") or "chamber" in reference.metadata.temperature_roles.get(name, "")]
    first_surface = surface_names[0]
    second_surface = surface_names[1]
    first_values = np.asarray(discharge[first_surface])
    second_values = np.asarray(discharge[second_surface])
    audit_by_id = {str(row["dataset_id"]): row for row in audits}
    metric_by_id = {str(row["dataset"]): row for row in metrics}

    lines = [
        "# Own TRC18650 capacity-coordinate case-study reference",
        "",
        "This reference is constructed directly from measured Chroma signals. No SPM or thermal model was run, no parameter was fitted, no voltage was corrected, and the loaded 1C curve is **not OCV**.",
        "",
        "## Primary reference and coordinate",
        "",
        f"`{REFERENCE_ID}` is used as the clean relative reference because the completed audit placed it alone in the lowest apparent-resistance group. Its measured complete-discharge capacity is `{float(discharge['capacity_ah']):.9f} Ah`; this is used as `Q_ref_Ah`. The coordinate is explicitly **coulomb-counted reference SOC**, with `SOC_ref_dis = 1 - discharged_Ah/Q_ref_Ah`.",
        f"The post-discharge recharge capacity is `{float(charge['capacity_ah']):.9f} Ah`. Raw discharge points and the comparison-only uniform `0.001` SOC interpolation are stored separately by `record_type`; interpolation does not replace or alter the raw observations.",
        "",
        f"Its median apparent rest-to-load series resistance is `{float(audit_by_id[REFERENCE_ID]['r_app_median_ohm']):.9f} ohm`. This is a relative setup diagnostic, not intrinsic cell resistance or model R0.",
        "",
        "## Signal-derived ambient-cycle segmentation",
        "",
        "The file does not start with a rest: it begins with a partial conditioning charge. CC/CV boundaries use the upper-voltage band together with declining positive current; discharge/rest identity uses step number and measured current sign.",
        "",
        "| Segment | Step | Rows (zero-based, inclusive) | Samples | Current [A] | Voltage [V] | Step capacity [Ah] | Signed total-capacity change [Ah] | Inferred/provisional duration [s] | Temperature ranges [degC] |",
        "|---|---:|---|---:|---|---|---|---:|---:|---|",
    ]
    for row in segment_rows:
        temperatures = "; ".join(
            f"{name} {row[name+'_min_c']:.3f}-{row[name+'_max_c']:.3f}" for name in temp_names
        )
        lines.append(
            f"| {row['segment']} | {row['step_number']} | {row['row_start_zero_based']}-{row['row_end_zero_based']} | {row['sample_count']} | {row['current_min_a']:.4f} to {row['current_max_a']:.4f} | {row['voltage_min_v']:.4f} to {row['voltage_max_v']:.4f} | {row['step_capacity_min_ah']:.4f} to {row['step_capacity_max_ah']:.4f} | {row['signed_total_capacity_change_ah']:.4f} | {row['inferred_provisional_duration_s']:.2f} | {temperatures} |"
        )
    lines.extend(
        [
            "",
            "All durations above are inferred/provisional. No raw time or timestamp exists, so no duration is presented as measured.",
            "",
            "## Charge-discharge directional voltage separation",
            "",
            f"The normalized charge and discharge curves share the interval `SOC_ref = {SOC_GRID[0]:.3f} to {SOC_GRID[-1]:.3f}`. `DeltaV_dir = V_charge - V_discharge` has mean `{1000*directional_stats['mean_v']:.3f} mV`, median `{1000*directional_stats['median_v']:.3f} mV`, RMS `{1000*directional_stats['rms_v']:.3f} mV`, minimum `{1000*directional_stats['min_v']:.3f} mV`, maximum `{1000*directional_stats['max_v']:.3f} mV`, P5 `{1000*directional_stats['p5_v']:.3f} mV`, and P95 `{1000*directional_stats['p95_v']:.3f} mV`.",
            "",
            "| SOC region | Mean [mV] | Median [mV] | RMS [mV] | P5/P95 [mV] |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for region in ("high", "middle", "low"):
        stats = directional_regions[region]
        lines.append(f"| {region} | {1000*stats['mean_v']:.3f} | {1000*stats['median_v']:.3f} | {1000*stats['rms_v']:.3f} | {1000*stats['p5_v']:.3f} / {1000*stats['p95_v']:.3f} |")
    lines.extend(
        [
            "",
            "This is **charge-discharge directional voltage separation**, not thermodynamic hysteresis. It can include ohmic, kinetic, diffusion, mapping, and true hysteretic contributions.",
            "",
            "## Descriptive temperature evidence",
            "",
            f"During the clean discharge, {first_surface} starts at `{first_values[0]:.3f} degC`, reaches `{np.nanmax(first_values):.3f} degC`, and ends at `{first_values[-1]:.3f} degC`; its maximum rise above the discharge-start value is `{np.nanmax(first_values)-first_values[0]:.3f} degC`. {second_surface} starts at `{second_values[0]:.3f} degC`, reaches `{np.nanmax(second_values):.3f} degC`, and ends at `{second_values[-1]:.3f} degC`.",
            f"At the end of the post-discharge rest, the two surface sensors read `{reference.numeric[first_surface][post_rest.end]:.3f}` and `{reference.numeric[second_surface][post_rest.end]:.3f} degC`. During discharge their mean signed difference ({first_surface} - {second_surface}) is `{np.nanmean(first_values-second_values):.3f} degC`, and maximum absolute difference is `{np.nanmax(np.abs(first_values-second_values)):.3f} degC`.",
            f"The chamber/ambient range over the complete file is " + ", ".join(f"`{name}: {np.nanmin(reference.numeric[name]):.3f}-{np.nanmax(reference.numeric[name]):.3f} degC`" for name in ambient_names) + ". Ten nonfinite ambient samples identified by the prior audit remain unmodified.",
            "",
            "These are capacity-coordinate descriptions only. Thermal time constants, heat-transfer coefficients, and dynamic thermal errors cannot be inferred without measured timestamps.",
            "",
            "## Complete-cycle comparison with clean reference",
            "",
            "`DeltaV_reference = V_affected - V_clean_reference` at matched normalized coulomb-counted reference SOC. The deviation combines temperature, cell, setup, degradation/state, and electrochemical differences; it is not assigned solely to parasitic resistance.",
            "",
            "| Dataset | Capacity [Ah] | Median R_app [ohm] | R_app ratio | Mean / median deviation [mV] | RMS / max abs [mV] | P5 / P95 [mV] | High/middle/low mean [mV] | Interpretation |",
            "|---|---:|---:|---:|---:|---:|---:|---|---|",
        ]
    )
    for dataset_id in (REFERENCE_ID, *AFFECTED_CYCLE_IDS):
        row = metric_by_id[dataset_id]
        lines.append(
            f"| {dataset_id} | {fmt(row['measured_discharge_capacity_ah'], 4)} | {fmt(row['median_r_app_ohm'], 6)} | {fmt(row['r_app_ratio_to_clean_reference'], 3)}x | {1000*float(row['mean_voltage_deviation_v']):.2f} / {1000*float(row['median_voltage_deviation_v']):.2f} | {1000*float(row['rms_voltage_deviation_v']):.2f} / {1000*float(row['max_abs_voltage_deviation_v']):.2f} | {1000*float(row['p5_voltage_deviation_v']):.2f} / {1000*float(row['p95_voltage_deviation_v']):.2f} | {1000*float(row['high_soc_mean_deviation_v']):.2f} / {1000*float(row['middle_soc_mean_deviation_v']):.2f} / {1000*float(row['low_soc_mean_deviation_v']):.2f} | {row['interpretation']} |"
        )
    lines.extend(
        [
            "",
            "## Repeated-transition pulse evidence",
            "",
            "| Dataset | N | Median [ohm] | IQR [ohm] | Range [ohm] | Ratio to clean | Detector |",
            "|---|---:|---:|---:|---|---:|---|",
        ]
    )
    for row in pulse_rows:
        lines.append(f"| {row['dataset']} | {row['transition_count']} | {row['median_r_app_ohm']:.9f} | {row['iqr_r_app_ohm']:.9f} | {row['min_r_app_ohm']:.9f}-{row['max_r_app_ohm']:.9f} | {row['r_app_ratio_to_clean_reference']:.3f}x | {row['detector_flag']} |")
    part_difference = abs(pulse_rows[0]["median_r_app_ohm"] - pulse_rows[1]["median_r_app_ohm"])
    lines.extend(
        [
            "",
            f"All ten pulse transitions lie in the elevated group. Part-1 and part-2 medians differ by only `{1000*part_difference:.3f} milliohm`, so the elevation is **persistent across repeated transitions**, not isolated. Its physical source is not uniquely identified.",
            "",
            "## Nonideality detector provenance",
            "",
            "The detector reuses the audit's largest-adjacent-log-gap grouping of experiment-level R_app medians. It does not impose a universal resistance threshold. Each output row retains transition count, reference/experiment median and IQR, ratio, reference identity, and an evidence statement in the metrics CSV.",
            "",
            "## Explicit answers",
            "",
            "1. **Coherent clean reference?** Yes, as a relative experimental reference in capacity/SOC coordinates: it has a complete 2.8999 Ah discharge/recharge and the lowest transition-response group. It is not asserted to be ideal or OCV.",
            f"2. **Directional separation magnitude?** Mean `{1000*directional_stats['mean_v']:.3f} mV`, RMS `{1000*directional_stats['rms_v']:.3f} mV`, with the full distribution and SOC regions reported above.",
            "3. **Affected-cycle deviation strength?** All four affected cycles show systematic capacity-matched voltage deviations; their numerical mean/RMS/maximum and regional values are tabulated above and exported without correction.",
            "4. **Persistent pulse elevation?** Yes. Both pulse parts and every retained transition occupy the elevated relative-response group.",
            "5. **Useful nonideality case studies?** Yes. They remain valid measurements of the complete experimental chain and are useful precisely because the audit detects a reproducible nonideality; they are not treated as calibration ground truth or discarded as failed measurements.",
            "6. **Quantitative versus qualitative conclusions?** Capacity, normalized SOC coordinates, measured voltage/temperature values, voltage deviations, and adjacent-sample R_app are quantitative. Any time-dependent electrical or thermal inference, thermal constants, dynamic RMSE, or unique attribution of the excess resistance remains qualitative/not identifiable because measured time is absent.",
            "",
            "## Stop point",
            "",
            "No electrochemical or thermal model was run; no OCP, R0, diffusion, or thermal parameter was fitted; and no voltage correction was applied.",
        ]
    )
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    FIGURE_DIRECTORY.mkdir(parents=True, exist_ok=True)
    experiment_list = load_all()
    experiments = {experiment.metadata.dataset_id: experiment for experiment in experiment_list}
    transitions_by_id = {
        dataset_id: extract_rest_to_load_transitions(experiment)
        for dataset_id, experiment in experiments.items()
    }
    transitions = [row for dataset_id in experiments for row in transitions_by_id[dataset_id]]
    audits = [audit_experiment(experiment, transitions_by_id[experiment.metadata.dataset_id]) for experiment in experiment_list]
    classify_and_augment(audits, transitions)
    audit_by_id = {str(row["dataset_id"]): row for row in audits}

    cycle_ids = (REFERENCE_ID, *AFFECTED_CYCLE_IDS)
    segments = {dataset_id: segment_complete_cycle(experiments[dataset_id]) for dataset_id in cycle_ids}
    segment_rows = [segment_report_row(experiments[REFERENCE_ID], segment) for segment in segments[REFERENCE_ID]]
    curves: dict[str, dict[str, np.ndarray | float]] = {}
    interps: dict[str, dict[str, np.ndarray]] = {}
    for dataset_id in cycle_ids:
        discharge_segment = next(segment for segment in segments[dataset_id] if segment.name == "main_1c_discharge")
        curves[dataset_id] = discharge_curve(experiments[dataset_id], discharge_segment)
        interps[dataset_id] = interpolate_on_soc(curves[dataset_id])
    curves["ambient_charge"] = charge_curve(experiments[REFERENCE_ID], segments[REFERENCE_ID])
    interps["ambient_charge"] = interpolate_on_soc(curves["ambient_charge"])

    write_csv(REFERENCE_CSV, reference_rows(experiments[REFERENCE_ID], curves[REFERENCE_ID], interps[REFERENCE_ID]))
    directional = np.asarray(interps["ambient_charge"]["voltage_v"]) - np.asarray(interps[REFERENCE_ID]["voltage_v"])
    directional_stats = voltage_statistics(directional)
    directional_regions = regional_voltage_statistics(SOC_GRID, directional)
    directional_rows = [
        {
            "reference_soc": float(soc),
            "soc_region": region_label(float(soc)),
            "discharge_voltage_v": float(interps[REFERENCE_ID]["voltage_v"][index]),
            "charge_voltage_v": float(interps["ambient_charge"]["voltage_v"][index]),
            "directional_voltage_separation_v": float(directional[index]),
            "definition": "V_charge - V_discharge; not thermodynamic hysteresis",
            "soc_label": "coulomb-counted reference SOC",
        }
        for index, soc in enumerate(SOC_GRID)
    ]
    write_csv(DIRECTIONAL_CSV, directional_rows)

    metrics: list[dict[str, object]] = []
    for experiment in experiment_list:
        dataset_id = experiment.metadata.dataset_id
        audit = audit_by_id[dataset_id]
        detector = relative_series_path_detector(audit, REFERENCE_ID)
        if dataset_id in cycle_ids:
            delta = np.asarray(interps[dataset_id]["voltage_v"]) - np.asarray(interps[REFERENCE_ID]["voltage_v"])
            stats = voltage_statistics(delta)
            regions = regional_voltage_statistics(SOC_GRID, delta)
            measured_capacity = float(curves[dataset_id]["capacity_ah"])
        else:
            stats = {key: float("nan") for key in ("mean_v", "median_v", "rms_v", "min_v", "max_v", "max_abs_v", "p5_v", "p95_v")}
            regions = {name: {"mean_v": float("nan"), "rms_v": float("nan")} for name in ("low", "middle", "high")}
            measured_capacity = float("nan")
        metrics.append(
            {
                "dataset": dataset_id,
                "cell": experiment.metadata.cell_id,
                "nominal_temperature": experiment.metadata.nominal_temperature,
                "measured_discharge_capacity_ah": measured_capacity,
                "capacity_status": "complete measured discharge" if dataset_id in cycle_ids else "not applicable: incomplete pulse sequence",
                "median_r_app_ohm": detector["experiment_median_r_app_ohm"],
                "iqr_r_app_ohm": detector["experiment_iqr_r_app_ohm"],
                "delta_r_app_vs_clean_reference_ohm": detector["r_app_excess_vs_reference_ohm"],
                "r_app_ratio_to_clean_reference": detector["r_app_ratio_to_reference"],
                "rms_voltage_deviation_v": stats["rms_v"],
                "mean_voltage_deviation_v": stats["mean_v"],
                "median_voltage_deviation_v": stats["median_v"],
                "max_abs_voltage_deviation_v": stats["max_abs_v"],
                "p5_voltage_deviation_v": stats["p5_v"],
                "p95_voltage_deviation_v": stats["p95_v"],
                "high_soc_mean_deviation_v": regions["high"]["mean_v"],
                "high_soc_rms_deviation_v": regions["high"]["rms_v"],
                "middle_soc_mean_deviation_v": regions["middle"]["mean_v"],
                "middle_soc_rms_deviation_v": regions["middle"]["rms_v"],
                "low_soc_mean_deviation_v": regions["low"]["mean_v"],
                "low_soc_rms_deviation_v": regions["low"]["rms_v"],
                "thermal_channel_status": thermal_status(experiment),
                "interpretation": interpretation(dataset_id),
                **detector,
            }
        )
    write_csv(METRICS_CSV, metrics)

    pulse_rows = []
    for dataset_id in PULSE_IDS:
        audit = audit_by_id[dataset_id]
        detector = relative_series_path_detector(audit, REFERENCE_ID)
        pulse_rows.append(
            {
                "dataset": dataset_id,
                "transition_count": int(audit["r_app_transition_count"]),
                "median_r_app_ohm": float(audit["r_app_median_ohm"]),
                "iqr_r_app_ohm": float(audit["r_app_iqr_ohm"]),
                "min_r_app_ohm": float(audit["r_app_min_ohm"]),
                "max_r_app_ohm": float(audit["r_app_max_ohm"]),
                "r_app_ratio_to_clean_reference": float(audit["r_app_ratio_vs_reference"]),
                "detector_flag": detector["detector_flag"],
            }
        )
    make_figures(experiments, curves, interps, transitions, directional)
    build_report(
        experiments,
        audits,
        segment_rows,
        curves,
        directional_stats,
        directional_regions,
        metrics,
        pulse_rows,
    )
    print(f"Built capacity-coordinate own-data case study with Q_ref={curves[REFERENCE_ID]['capacity_ah']:.9f} Ah")
    print(f"Directional separation: mean={1000*directional_stats['mean_v']:.3f} mV, RMS={1000*directional_stats['rms_v']:.3f} mV")
    print(f"Outputs written under {OUTPUT_DIRECTORY}")


if __name__ == "__main__":
    main()
