"""Capacity-coordinate tools for the own-data TRC18650 case study.

These functions operate on :class:`StandardExperiment` objects from the
experimental-data audit.  They do not run or fit a battery model, alter raw
voltage, or treat the loaded 1C curve as OCV.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from analysis.experimental_data_audit import StandardExperiment


REST_CURRENT_A = 0.02
ACTIVE_CURRENT_A = 0.5
SOC_GRID = np.linspace(0.0, 1.0, 1001)


@dataclass(frozen=True)
class PhysicalSegment:
    """Inclusive raw-row bounds and an evidence-based physical label."""

    name: str
    start: int
    end: int
    step_number: int
    basis: str


def _step_runs(experiment: StandardExperiment) -> list[tuple[int, int, int]]:
    step = experiment.numeric["Step No."]
    changes = np.flatnonzero(np.diff(step) != 0) + 1
    starts = np.r_[0, changes]
    ends = np.r_[changes - 1, step.size - 1]
    return [(int(a), int(b), int(step[a])) for a, b in zip(starts, ends, strict=True)]


def _split_charge_cc_cv(experiment: StandardExperiment, start: int, end: int) -> int | None:
    """Return first CV row using simultaneous voltage and current evidence."""

    current = experiment.raw_current_a[start : end + 1]
    voltage = experiment.voltage_v[start : end + 1]
    active_current = current[current > ACTIVE_CURRENT_A]
    if active_current.size == 0:
        # A positive low-current step already at the upper voltage limit is a
        # CV-only top-up rather than a CC segment.
        return start if float(np.nanmedian(voltage)) >= 4.15 else None
    nominal_current = float(np.quantile(active_current, 0.9))
    voltage_threshold = float(np.nanmax(voltage) - 0.010)
    candidate = (voltage >= voltage_threshold) & (current < 0.98 * nominal_current)
    # Require the voltage to remain in the near-limit band for the next records;
    # this rejects a single noisy crossing without using reconstructed time.
    for local in np.flatnonzero(candidate):
        right = min(local + 5, candidate.size)
        if right - local >= 3 and np.count_nonzero(voltage[local:right] >= voltage_threshold) >= 3:
            return start + int(local)
    return None


def segment_complete_cycle(experiment: StandardExperiment) -> list[PhysicalSegment]:
    """Segment a cycle from step identity, raw current, voltage, and capacity signs."""

    runs = _step_runs(experiment)
    current = experiment.raw_current_a
    discharge_runs = [
        (a, b, step)
        for a, b, step in runs
        if float(np.nanmedian(current[a : b + 1])) < -ACTIVE_CURRENT_A
    ]
    if not discharge_runs:
        raise ValueError(f"No active discharge found in {experiment.metadata.dataset_id}")
    main_discharge = max(discharge_runs, key=lambda item: item[1] - item[0])
    discharge_start, discharge_end, _ = main_discharge
    output: list[PhysicalSegment] = []
    for start, end, step_number in runs:
        values = current[start : end + 1]
        median = float(np.nanmedian(values))
        max_abs = float(np.nanmax(np.abs(values)))
        if max_abs <= REST_CURRENT_A:
            if end < discharge_start:
                label = "pre_discharge_rest" if start > 0 else "initial_rest"
            elif start > discharge_end and not any(
                np.nanmedian(current[a : b + 1]) > ACTIVE_CURRENT_A
                for a, b, _ in runs
                if discharge_end < a < start
            ):
                label = "post_discharge_rest"
            else:
                label = "final_rest"
            output.append(PhysicalSegment(label, start, end, step_number, "zero-current step"))
        elif median < -ACTIVE_CURRENT_A:
            label = "main_1c_discharge" if (start, end, step_number) == main_discharge else "other_discharge"
            output.append(PhysicalSegment(label, start, end, step_number, "negative near-1C current step"))
        elif median > REST_CURRENT_A:
            prefix = "initial_conditioning_charge" if end < discharge_start else "recharge"
            cv_start = _split_charge_cc_cv(experiment, start, end)
            if cv_start is None:
                output.append(PhysicalSegment(prefix, start, end, step_number, "positive-current charge step"))
            elif cv_start <= start:
                output.append(
                    PhysicalSegment(
                        f"{prefix}_cv",
                        start,
                        end,
                        step_number,
                        "upper-voltage low-current charge step",
                    )
                )
            else:
                output.append(
                    PhysicalSegment(
                        f"{prefix}_cc",
                        start,
                        cv_start - 1,
                        step_number,
                        "positive near-constant current below upper-voltage band",
                    )
                )
                output.append(
                    PhysicalSegment(
                        f"{prefix}_cv",
                        cv_start,
                        end,
                        step_number,
                        "upper-voltage band with declining positive current",
                    )
                )
        else:
            output.append(PhysicalSegment("ambiguous_step", start, end, step_number, "mixed/low-current step"))
    return output


def segment_report_row(experiment: StandardExperiment, segment: PhysicalSegment) -> dict[str, object]:
    sl = slice(segment.start, segment.end + 1)
    numeric = experiment.numeric
    total_capacity = numeric["Total Capacity(Ah)"]
    result: dict[str, object] = {
        "segment": segment.name,
        "step_number": segment.step_number,
        "row_start_zero_based": segment.start,
        "row_end_zero_based": segment.end,
        "sample_count": segment.end - segment.start + 1,
        "current_min_a": float(np.nanmin(experiment.raw_current_a[sl])),
        "current_max_a": float(np.nanmax(experiment.raw_current_a[sl])),
        "voltage_min_v": float(np.nanmin(experiment.voltage_v[sl])),
        "voltage_max_v": float(np.nanmax(experiment.voltage_v[sl])),
        "step_capacity_min_ah": float(np.nanmin(numeric["Capacity(Ah)"][sl])),
        "step_capacity_max_ah": float(np.nanmax(numeric["Capacity(Ah)"][sl])),
        "signed_total_capacity_change_ah": float(total_capacity[segment.end] - total_capacity[segment.start]),
        "inferred_provisional_duration_s": float((segment.end - segment.start) * experiment.estimated_dt_s),
        "segmentation_basis": segment.basis,
    }
    for name in (column for column in numeric if column.startswith("Temp")):
        result[f"{name}_min_c"] = float(np.nanmin(numeric[name][sl]))
        result[f"{name}_max_c"] = float(np.nanmax(numeric[name][sl]))
    return result


def _monotone_capacity_coordinate(capacity: np.ndarray) -> np.ndarray:
    accumulated = np.abs(capacity - capacity[0])
    if np.any(np.diff(accumulated) < -1.0e-12):
        raise ValueError("Capacity counter is non-monotonic inside the selected segment; raw data were not repaired.")
    final = float(accumulated[-1])
    if final <= 0.0:
        raise ValueError("Segment has no capacity throughput.")
    return accumulated / final


def discharge_curve(experiment: StandardExperiment, segment: PhysicalSegment) -> dict[str, np.ndarray | float]:
    rows = np.arange(segment.start, segment.end + 1)
    capacity = experiment.numeric["Capacity(Ah)"][rows]
    z_dis = _monotone_capacity_coordinate(capacity)
    result: dict[str, np.ndarray | float] = {
        "rows": rows,
        "throughput_ah": z_dis * abs(float(capacity[-1] - capacity[0])),
        "z": z_dis,
        "soc": 1.0 - z_dis,
        "voltage_v": experiment.voltage_v[rows],
        "capacity_ah": abs(float(capacity[-1] - capacity[0])),
    }
    for name in (column for column in experiment.numeric if column.startswith("Temp")):
        result[name] = experiment.numeric[name][rows]
    return result


def charge_curve(experiment: StandardExperiment, segments: list[PhysicalSegment]) -> dict[str, np.ndarray | float]:
    charge_segments = [segment for segment in segments if segment.name.startswith("recharge")]
    if not charge_segments:
        raise ValueError(f"No post-discharge recharge found in {experiment.metadata.dataset_id}")
    start = min(segment.start for segment in charge_segments)
    end = max(segment.end for segment in charge_segments)
    rows = np.arange(start, end + 1)
    capacity = experiment.numeric["Capacity(Ah)"][rows]
    z_charge = _monotone_capacity_coordinate(capacity)
    result: dict[str, np.ndarray | float] = {
        "rows": rows,
        "throughput_ah": z_charge * abs(float(capacity[-1] - capacity[0])),
        "z": z_charge,
        "soc": z_charge,
        "voltage_v": experiment.voltage_v[rows],
        "capacity_ah": abs(float(capacity[-1] - capacity[0])),
    }
    for name in (column for column in experiment.numeric if column.startswith("Temp")):
        result[name] = experiment.numeric[name][rows]
    return result


def interpolate_on_soc(curve: dict[str, np.ndarray | float], grid: np.ndarray = SOC_GRID) -> dict[str, np.ndarray]:
    """Interpolate measured curve values after stable sorting and duplicate collapse."""

    soc = np.asarray(curve["soc"], dtype=float)
    order = np.argsort(soc, kind="stable")
    sorted_soc = soc[order]
    unique_soc, first, counts = np.unique(sorted_soc, return_index=True, return_counts=True)
    result: dict[str, np.ndarray] = {"soc": grid.copy()}
    for name, values in curve.items():
        if name in {"capacity_ah"}:
            continue
        array = np.asarray(values)
        if array.ndim != 1 or array.size != soc.size or name == "rows":
            continue
        sorted_values = array[order].astype(float)
        collapsed = np.asarray(
            [np.nanmean(sorted_values[index : index + count]) for index, count in zip(first, counts, strict=True)]
        )
        result[name] = np.interp(grid, unique_soc, collapsed)
    sorted_rows = np.asarray(curve["rows"], dtype=int)[order]
    right = np.searchsorted(unique_soc, grid, side="left")
    right = np.clip(right, 0, unique_soc.size - 1)
    left = np.maximum(right - 1, 0)
    representative_rows = sorted_rows[first]
    result["source_row_left"] = representative_rows[left]
    result["source_row_right"] = representative_rows[right]
    return result


def voltage_statistics(delta_v: np.ndarray) -> dict[str, float]:
    finite = np.asarray(delta_v, dtype=float)
    finite = finite[np.isfinite(finite)]
    return {
        "mean_v": float(np.mean(finite)),
        "median_v": float(np.median(finite)),
        "rms_v": float(np.sqrt(np.mean(finite**2))),
        "min_v": float(np.min(finite)),
        "max_v": float(np.max(finite)),
        "max_abs_v": float(np.max(np.abs(finite))),
        "p5_v": float(np.quantile(finite, 0.05)),
        "p95_v": float(np.quantile(finite, 0.95)),
    }


def regional_voltage_statistics(soc: np.ndarray, delta_v: np.ndarray) -> dict[str, dict[str, float]]:
    masks = {
        "low": soc < 1.0 / 3.0,
        "middle": (soc >= 1.0 / 3.0) & (soc < 2.0 / 3.0),
        "high": soc >= 2.0 / 3.0,
    }
    return {name: voltage_statistics(delta_v[mask]) for name, mask in masks.items()}


def relative_series_path_detector(
    audit_row: dict[str, object],
    reference_id: str,
) -> dict[str, object]:
    """Map the audit's data-driven grouping to a reusable evidence flag."""

    count = int(audit_row["r_app_transition_count"])
    dataset_id = str(audit_row["dataset_id"])
    if count < 2:
        flag = "INSUFFICIENT TRANSITIONS"
        evidence = "Fewer than two comparable adjacent rest-to-load transitions."
    elif dataset_id == reference_id:
        flag = "NO RELATIVE ELEVATION DETECTED"
        evidence = "Selected lowest-group experimental reference; this is not intrinsic-R0 certification."
    elif "elevated" in str(audit_row["resistance_group"]):
        flag = "EVIDENCE OF ELEVATED SERIES-PATH RESPONSE"
        evidence = "Experiment median lies in the audit's data-driven elevated group relative to the clean reference."
    else:
        flag = "NO RELATIVE ELEVATION DETECTED"
        evidence = "Experiment remains in the audit's lowest internally consistent group."
    return {
        "detector_flag": flag,
        "detector_evidence": evidence,
        "reference_dataset": reference_id,
        "reference_group_method": "largest adjacent gap in log experiment-median R_app from own-data audit",
        "transition_count": count,
        "experiment_median_r_app_ohm": float(audit_row["r_app_median_ohm"]),
        "experiment_iqr_r_app_ohm": float(audit_row["r_app_iqr_ohm"]),
        "reference_median_r_app_ohm": float(audit_row["reference_group_median_ohm"]),
        "reference_iqr_r_app_ohm": float(audit_row["reference_group_iqr_ohm"]),
        "r_app_ratio_to_reference": float(audit_row["r_app_ratio_vs_reference"]),
        "r_app_excess_vs_reference_ohm": float(audit_row["r_app_excess_vs_reference_ohm"]),
    }
