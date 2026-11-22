"""Reusable experimental battery-data audit utilities.

The module standardizes signals without changing raw files.  It deliberately
distinguishes a measured time channel from a provisional timebase inferred
from capacity/current sampling.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


REST_CURRENT_A = 0.02
ACTIVE_CURRENT_A = 0.5
COMPARABLE_1C_CURRENT_A = 2.5
NOMINAL_REFERENCE_AH = 2.9


@dataclass(frozen=True)
class ExperimentMetadata:
    dataset_id: str
    path: Path
    cell_id: str
    nominal_temperature: str
    encoded_date: str
    temperature_roles: dict[str, str]


@dataclass(frozen=True)
class StandardExperiment:
    metadata: ExperimentMetadata
    columns: tuple[str, ...]
    header_present: bool
    numeric: dict[str, np.ndarray]
    descriptions: np.ndarray
    estimated_dt_s: float
    estimated_time_s: np.ndarray
    time_basis: str
    test_type: str
    raw_current_sign: str
    current_sign_confidence: str

    @property
    def voltage_v(self) -> np.ndarray:
        """Canonical voltage view; the source column remains unchanged."""

        return self.numeric["Voltage(V)"]

    @property
    def raw_current_a(self) -> np.ndarray:
        """Canonical raw-current view using the acquisition sign convention."""

        return self.numeric["Current(A)"]

    @property
    def model_current_a(self) -> np.ndarray:
        """Discharge-positive view for downstream loaders; never written to source."""

        if self.raw_current_sign != "discharge_negative":
            raise ValueError("Current sign is not established; model-current conversion is unsafe.")
        return -self.raw_current_a


def _detect_header_row(path: Path) -> int | None:
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        for index, row in enumerate(csv.reader(handle)):
            if "Voltage(V)" in row and "Current(A)" in row:
                return index
            if index >= 20:
                break
    return None


def _read_rows(path: Path, fallback_columns: tuple[str, ...] | None) -> tuple[tuple[str, ...], list[list[str]], bool]:
    header_row = _detect_header_row(path)
    with path.open("r", encoding="utf-8-sig", errors="replace", newline="") as handle:
        rows = list(csv.reader(handle))
    if header_row is None:
        if fallback_columns is None:
            raise ValueError(f"No Chroma header and no fallback schema: {path}")
        columns = fallback_columns
        data_rows = rows
        header_present = False
    else:
        columns = tuple(value.strip() for value in rows[header_row])
        data_rows = rows[header_row + 1 :]
        header_present = True
    while columns and columns[-1] == "":
        columns = columns[:-1]
    cleaned_rows: list[list[str]] = []
    for row in data_rows:
        row = row[: len(columns)] + [""] * max(0, len(columns) - len(row))
        if any(value.strip() for value in row):
            cleaned_rows.append(row)
    return columns, cleaned_rows, header_present


def _float_column(rows: list[list[str]], index: int) -> np.ndarray:
    values = np.empty(len(rows), dtype=float)
    for row_index, row in enumerate(rows):
        try:
            values[row_index] = float(row[index].strip())
        except (ValueError, IndexError):
            values[row_index] = np.nan
    return values


def contiguous_runs(mask: np.ndarray) -> list[tuple[int, int]]:
    indices = np.flatnonzero(mask)
    if indices.size == 0:
        return []
    breaks = np.flatnonzero(np.diff(indices) > 1)
    starts = np.concatenate(([indices[0]], indices[breaks + 1]))
    ends = np.concatenate((indices[breaks], [indices[-1]]))
    return [(int(start), int(end)) for start, end in zip(starts, ends, strict=True)]


def estimate_sample_interval_s(
    current_a: np.ndarray, capacity_ah: np.ndarray, step_number: np.ndarray
) -> float:
    """Infer nominal sampling from active step capacity; never call it measured time."""

    estimates: list[float] = []
    weights: list[int] = []
    active = np.isfinite(current_a) & (np.abs(current_a) >= ACTIVE_CURRENT_A)
    for start, end in contiguous_runs(active):
        step_breaks = np.flatnonzero(np.diff(step_number[start : end + 1]) != 0) + start
        sub_starts = np.concatenate(([start], step_breaks + 1))
        sub_ends = np.concatenate((step_breaks, [end]))
        for left, right in zip(sub_starts, sub_ends, strict=True):
            intervals = int(right - left)
            if intervals < 5:
                continue
            delta_capacity = abs(float(capacity_ah[right] - capacity_ah[left]))
            mean_current = float(np.nanmean(np.abs(current_a[left : right + 1])))
            if delta_capacity > 0.0 and mean_current > ACTIVE_CURRENT_A:
                estimate = delta_capacity * 3600.0 / mean_current / intervals
                if 1.0e-3 <= estimate <= 120.0:
                    estimates.append(estimate)
                    weights.append(intervals)
    if not estimates:
        raise ValueError("Cannot infer a provisional sample interval from capacity/current.")
    order = np.argsort(estimates)
    sorted_values = np.asarray(estimates)[order]
    sorted_weights = np.asarray(weights)[order]
    cumulative = np.cumsum(sorted_weights)
    return float(sorted_values[np.searchsorted(cumulative, cumulative[-1] / 2.0)])


def infer_test_type(current_a: np.ndarray) -> str:
    active_positive = current_a >= ACTIVE_CURRENT_A
    active_negative = current_a <= -ACTIVE_CURRENT_A
    rest = np.abs(current_a) <= REST_CURRENT_A
    transitions = int(np.count_nonzero(rest[:-1] & (np.abs(current_a[1:]) >= ACTIVE_CURRENT_A)))
    if np.any(active_positive) and np.any(active_negative):
        return "combined cycle"
    if np.any(active_negative) and transitions >= 2:
        return "pulse discharge/rest"
    if np.any(active_negative):
        return "discharge"
    if np.any(active_positive):
        return "charge"
    if np.all(rest):
        return "rest"
    return "unknown"


def infer_raw_current_sign(current_a: np.ndarray, voltage_v: np.ndarray) -> tuple[str, str]:
    rest = np.abs(current_a[:-1]) <= REST_CURRENT_A
    loaded = np.abs(current_a[1:]) >= ACTIVE_CURRENT_A
    indices = np.flatnonzero(rest & loaded)
    evidence: list[bool] = []
    for index in indices:
        delta_i = current_a[index + 1] - current_a[index]
        delta_v = voltage_v[index + 1] - voltage_v[index]
        if delta_i < 0.0:
            evidence.append(delta_v < 0.0)
        elif delta_i > 0.0:
            evidence.append(delta_v > 0.0)
    if evidence and np.mean(evidence) >= 0.8:
        return "discharge_negative", "high"
    return "ambiguous", "low"


def load_chroma_experiment(
    metadata: ExperimentMetadata,
    fallback_columns: tuple[str, ...] | None = None,
) -> StandardExperiment:
    columns, rows, header_present = _read_rows(metadata.path, fallback_columns)
    required = ("Step No.", "Voltage(V)", "Current(A)", "Capacity(Ah)", "Mode", "Total Capacity(Ah)")
    missing = [name for name in required if name not in columns]
    if missing:
        raise ValueError(f"Missing required Chroma columns {missing}: {metadata.path}")
    numeric: dict[str, np.ndarray] = {}
    for index, name in enumerate(columns):
        if name != "Description":
            numeric[name] = _float_column(rows, index)
    description_index = columns.index("Description") if "Description" in columns else None
    descriptions = np.asarray(
        [row[description_index].strip() if description_index is not None else "" for row in rows],
        dtype=object,
    )
    dt = estimate_sample_interval_s(
        numeric["Current(A)"], numeric["Capacity(Ah)"], numeric["Step No."]
    )
    estimated_time = np.arange(len(rows), dtype=float) * dt
    test_type = infer_test_type(numeric["Current(A)"])
    sign, confidence = infer_raw_current_sign(numeric["Current(A)"], numeric["Voltage(V)"])
    return StandardExperiment(
        metadata=metadata,
        columns=columns,
        header_present=header_present,
        numeric=numeric,
        descriptions=descriptions,
        estimated_dt_s=dt,
        estimated_time_s=estimated_time,
        time_basis="provisional_uniform_from_step_capacity_and_current",
        test_type=test_type,
        raw_current_sign=sign,
        current_sign_confidence=confidence,
    )


def extract_rest_to_load_transitions(experiment: StandardExperiment) -> list[dict[str, object]]:
    current = experiment.numeric["Current(A)"]
    voltage = experiment.numeric["Voltage(V)"]
    candidates = np.flatnonzero(
        (np.abs(current[:-1]) <= REST_CURRENT_A)
        & (np.abs(current[1:]) >= COMPARABLE_1C_CURRENT_A)
        & (np.abs(current[1:] - current[:-1]) >= COMPARABLE_1C_CURRENT_A)
    )
    rows: list[dict[str, object]] = []
    for number, index in enumerate(candidates, start=1):
        delta_i = float(current[index + 1] - current[index])
        delta_v = float(voltage[index + 1] - voltage[index])
        direction = "discharge" if delta_i < 0.0 else "charge"
        expected_sign = (direction == "discharge" and delta_v < 0.0) or (
            direction == "charge" and delta_v > 0.0
        )
        rows.append(
            {
                "dataset_id": experiment.metadata.dataset_id,
                "transition_id": number,
                "pre_index": int(index),
                "loaded_index": int(index + 1),
                "estimated_time_s": float(experiment.estimated_time_s[index + 1]),
                "direction": direction,
                "pre_current_a": float(current[index]),
                "loaded_current_a": float(current[index + 1]),
                "delta_current_a": delta_i,
                "pre_voltage_v": float(voltage[index]),
                "loaded_voltage_v": float(voltage[index + 1]),
                "delta_voltage_v": delta_v,
                "r_app_ohm": abs(delta_v) / abs(delta_i),
                "voltage_response_sign_consistent": expected_sign,
                "diagnostic_label": "apparent rest-to-load series resistance; not intrinsic cell R0",
            }
        )
    return rows


def robust_summary(values: Iterable[float]) -> dict[str, float | int]:
    array = np.asarray(list(values), dtype=float)
    if array.size == 0:
        return {
            "count": 0,
            "median": float("nan"),
            "mean": float("nan"),
            "std": float("nan"),
            "min": float("nan"),
            "max": float("nan"),
            "iqr": float("nan"),
        }
    return {
        "count": int(array.size),
        "median": float(np.median(array)),
        "mean": float(np.mean(array)),
        "std": float(np.std(array, ddof=1)) if array.size > 1 else 0.0,
        "min": float(np.min(array)),
        "max": float(np.max(array)),
        "iqr": float(np.quantile(array, 0.75) - np.quantile(array, 0.25)),
    }


def largest_log_gap_reference_group(experiment_medians: dict[str, float]) -> tuple[set[str], float]:
    """Split the lowest internally coherent group using the data's largest log gap."""

    ordered = sorted(experiment_medians.items(), key=lambda item: item[1])
    if len(ordered) < 2:
        return {ordered[0][0]} if ordered else set(), float("nan")
    log_values = np.log([value for _, value in ordered])
    gaps = np.diff(log_values)
    split = int(np.argmax(gaps)) + 1
    threshold = float(np.sqrt(ordered[split - 1][1] * ordered[split][1]))
    return {name for name, _ in ordered[:split]}, threshold


def _temperature_columns(experiment: StandardExperiment) -> list[str]:
    return [name for name in experiment.numeric if re.fullmatch(r"Temp\d+", name)]


def _quantization(values: np.ndarray) -> float:
    finite_pairs = np.isfinite(values[:-1]) & np.isfinite(values[1:])
    differences = np.round(np.abs(np.diff(values)[finite_pairs]), 4)
    differences = differences[differences > 0.0]
    if differences.size == 0:
        return 0.0
    unique, counts = np.unique(differences, return_counts=True)
    return float(unique[np.argmax(counts)])


def _temperature_jump_screen(values: np.ndarray) -> tuple[float, int]:
    """Robustly flag abrupt changes; this is a screening diagnostic, not rejection."""

    finite_pairs = np.isfinite(values[:-1]) & np.isfinite(values[1:])
    differences = np.abs(np.diff(values)[finite_pairs])
    if differences.size == 0:
        return float("nan"), 0
    median = float(np.median(differences))
    mad = float(np.median(np.abs(differences - median)))
    threshold = max(0.5, median + 10.0 * 1.4826 * mad)
    return threshold, int(np.count_nonzero(differences > threshold))


def _temperature_pair_diagnostic(experiment: StandardExperiment, columns: list[str]) -> dict[str, object]:
    if len(columns) < 2:
        return {
            "temperature_pair_mean_offset_c": float("nan"),
            "temperature_pair_correlation": float("nan"),
            "temperature_pair_lag_samples": "not_assessable",
        }
    first = experiment.numeric[columns[0]]
    second = experiment.numeric[columns[1]]
    finite = np.isfinite(first) & np.isfinite(second)
    first = first[finite]
    second = second[finite]
    correlation = float(np.corrcoef(first, second)[0, 1])
    centered_a = first - np.mean(first)
    centered_b = second - np.mean(second)
    best_lag = 0
    best_corr = -np.inf
    maximum_lag = min(100, first.size // 20)
    for lag in range(-maximum_lag, maximum_lag + 1):
        if lag < 0:
            a, b = centered_a[-lag:], centered_b[:lag]
        elif lag > 0:
            a, b = centered_a[:-lag], centered_b[lag:]
        else:
            a, b = centered_a, centered_b
        if a.size > 10 and np.std(a) > 0.0 and np.std(b) > 0.0:
            value = float(np.corrcoef(a, b)[0, 1])
            if value > best_corr:
                best_corr, best_lag = value, lag
    return {
        "temperature_pair_mean_offset_c": float(np.mean(first - second)),
        "temperature_pair_correlation": correlation,
        "temperature_pair_lag_samples": best_lag,
    }


def segment_capacities(experiment: StandardExperiment) -> list[dict[str, object]]:
    current = experiment.numeric["Current(A)"]
    voltage = experiment.numeric["Voltage(V)"]
    capacity = experiment.numeric["Capacity(Ah)"]
    step = experiment.numeric["Step No."]
    # Preserve the low-current CV tail; the rest threshold separates it from rest.
    active = np.abs(current) > REST_CURRENT_A
    segments: list[dict[str, object]] = []
    for start, end in contiguous_runs(active):
        step_breaks = np.flatnonzero(np.diff(step[start : end + 1]) != 0) + start
        sub_starts = np.concatenate(([start], step_breaks + 1))
        sub_ends = np.concatenate((step_breaks, [end]))
        for left, right in zip(sub_starts, sub_ends, strict=True):
            direction = "charge" if np.nanmedian(current[left : right + 1]) > 0.0 else "discharge"
            counter_capacity = abs(float(capacity[right] - capacity[left]))
            if counter_capacity == 0.0:
                counter_capacity = float(np.nanmax(capacity[left : right + 1]) - np.nanmin(capacity[left : right + 1]))
            integrated = float(
                abs(np.trapezoid(current[left : right + 1], dx=experiment.estimated_dt_s))
                / 3600.0
            )
            preceded_rest = left > 0 and abs(current[left - 1]) <= REST_CURRENT_A
            followed_rest = right + 1 < current.size and abs(current[right + 1]) <= REST_CURRENT_A
            complete = bool(
                preceded_rest
                and followed_rest
                and (
                    (direction == "discharge" and voltage[right] <= 3.0)
                    or (direction == "charge" and voltage[right] >= 4.15)
                )
            )
            segments.append(
                {
                    "direction": direction,
                    "start_index": int(left),
                    "end_index": int(right),
                    "counter_capacity_ah": counter_capacity,
                    "integrated_capacity_ah": integrated,
                    "start_voltage_v": float(voltage[left]),
                    "end_voltage_v": float(voltage[right]),
                    "preceded_by_rest": preceded_rest,
                    "followed_by_rest": followed_rest,
                    "appears_complete": complete,
                }
            )
    return segments


def audit_experiment(experiment: StandardExperiment, transitions: list[dict[str, object]]) -> dict[str, object]:
    numeric = experiment.numeric
    current = numeric["Current(A)"]
    voltage = numeric["Voltage(V)"]
    total_capacity = numeric["Total Capacity(Ah)"]
    finite_values = np.concatenate([values for values in numeric.values()])
    finite_count = int(np.count_nonzero(np.isfinite(finite_values)))
    nonfinite_count = int(finite_values.size - finite_count)
    integrated_signed = float(np.trapezoid(current, dx=experiment.estimated_dt_s) / 3600.0)
    integrated_absolute = float(np.trapezoid(np.abs(current), dx=experiment.estimated_dt_s) / 3600.0)
    counter_signed = float(total_capacity[-1] - total_capacity[0])
    counter_absolute = float(np.sum(np.abs(np.diff(total_capacity))))
    transition_summary = robust_summary(row["r_app_ohm"] for row in transitions)
    discharge_transitions = [row for row in transitions if row["direction"] == "discharge"]
    discharge_summary = robust_summary(row["r_app_ohm"] for row in discharge_transitions)
    delta_v = np.diff(voltage)
    delta_i = np.diff(current)
    quiet = np.abs(delta_i) < 0.05
    max_unexplained_jump = float(np.max(np.abs(delta_v[quiet]))) if np.any(quiet) else float("nan")
    sign_inconsistent = int(sum(not bool(row["voltage_response_sign_consistent"]) for row in transitions))
    repeated_extrema = int(
        max(
            np.count_nonzero(np.isclose(voltage, np.min(voltage), atol=5.0e-5)),
            np.count_nonzero(np.isclose(voltage, np.max(voltage), atol=5.0e-5)),
        )
    )
    segments = segment_capacities(experiment)
    discharge_segments = [row for row in segments if row["direction"] == "discharge"]
    charge_segments = [row for row in segments if row["direction"] == "charge"]
    main_discharge = max(discharge_segments, key=lambda row: row["counter_capacity_ah"], default=None)
    main_charge = max(charge_segments, key=lambda row: row["counter_capacity_ah"], default=None)
    temp_columns = _temperature_columns(experiment)
    temp_values = np.concatenate([numeric[name] for name in temp_columns]) if temp_columns else np.asarray([])
    temp_jumps = {}
    for name in temp_columns:
        differences = np.abs(np.diff(numeric[name]))
        finite_differences = differences[np.isfinite(differences)]
        temp_jumps[name] = float(np.max(finite_differences)) if finite_differences.size else float("nan")
    temperature_diagnostics: dict[str, object] = {}
    for name in temp_columns:
        values = numeric[name]
        role = experiment.metadata.temperature_roles.get(name, "unknown")
        jump_threshold, jump_count = _temperature_jump_screen(values)
        temperature_diagnostics.update(
            {
                f"{name}_role": role,
                f"{name}_initial_c": float(values[np.flatnonzero(np.isfinite(values))[0]]),
                f"{name}_final_c": float(values[np.flatnonzero(np.isfinite(values))[-1]]),
                f"{name}_min_c": float(np.nanmin(values)),
                f"{name}_max_c": float(np.nanmax(values)),
                f"{name}_excursion_c": float(np.nanmax(values) - np.nanmin(values)),
                f"{name}_quantization_c": _quantization(values),
                f"{name}_max_adjacent_jump_c": temp_jumps[name],
                f"{name}_jump_screen_threshold_c": jump_threshold,
                f"{name}_jump_flag_count": jump_count,
                f"{name}_nonfinite_count": int(np.count_nonzero(~np.isfinite(values))),
            }
        )
    temperature_diagnostics.update(_temperature_pair_diagnostic(experiment, temp_columns))
    ambient_columns = [
        name
        for name in temp_columns
        if "ambient" in experiment.metadata.temperature_roles.get(name, "")
        or "chamber" in experiment.metadata.temperature_roles.get(name, "")
    ]
    surface_columns = [name for name in temp_columns if "surface" in experiment.metadata.temperature_roles.get(name, "")]
    if ambient_columns and surface_columns:
        ambient = numeric[ambient_columns[0]]
        surface = numeric[surface_columns[0]]
        finite = np.isfinite(surface) & np.isfinite(ambient)
        valid = np.flatnonzero(finite)
        temperature_diagnostics["surface_minus_ambient_initial_c"] = float(surface[valid[0]] - ambient[valid[0]])
        temperature_diagnostics["surface_minus_ambient_final_c"] = float(surface[valid[-1]] - ambient[valid[-1]])
    else:
        temperature_diagnostics["surface_minus_ambient_initial_c"] = float("nan")
        temperature_diagnostics["surface_minus_ambient_final_c"] = float("nan")
    active_runs = contiguous_runs(np.abs(current) >= ACTIVE_CURRENT_A)
    rest_runs = contiguous_runs(np.abs(current) <= REST_CURRENT_A)
    primary_temp = numeric[surface_columns[0] if surface_columns else temp_columns[0]] if temp_columns else np.asarray([])
    heating_runs = sum(primary_temp[end] > primary_temp[start] for start, end in active_runs if end > start) if primary_temp.size else 0
    cooling_runs = sum(primary_temp[end] < primary_temp[start] for start, end in rest_runs if end > start) if primary_temp.size else 0
    capacity_complete = bool(
        experiment.test_type != "pulse discharge/rest"
        and (
            (main_discharge and main_discharge["appears_complete"])
            or (main_charge and main_charge["appears_complete"])
        )
    )
    measured_capacity = max(
        [float(row["counter_capacity_ah"]) for row in (main_discharge, main_charge) if row is not None],
        default=float("nan"),
    )
    output: dict[str, object] = {
        "dataset_id": experiment.metadata.dataset_id,
        "full_path": str(experiment.metadata.path.resolve()),
        "filename": experiment.metadata.path.name,
        "cell_id": experiment.metadata.cell_id,
        "nominal_test_temperature": experiment.metadata.nominal_temperature,
        "encoded_date": experiment.metadata.encoded_date,
        "test_type_signal_verified": experiment.test_type,
        "sample_count": int(current.size),
        "raw_time_available": False,
        "time_basis": experiment.time_basis,
        "estimated_dt_s": experiment.estimated_dt_s,
        "estimated_duration_s": float(experiment.estimated_time_s[-1]),
        "raw_duplicate_timestamps": "not_assessable_no_time_channel",
        "raw_nonmonotonic_timestamps": "not_assessable_no_time_channel",
        "raw_long_gaps": "not_assessable_no_time_channel",
        "min_positive_dt_s": experiment.estimated_dt_s,
        "median_dt_s": experiment.estimated_dt_s,
        "max_normal_dt_s": experiment.estimated_dt_s,
        "available_channels": ";".join(experiment.columns),
        "raw_current_sign": experiment.raw_current_sign,
        "current_sign_confidence": experiment.current_sign_confidence,
        "model_loader_mapping": "I_model=-I_raw",
        "standard_channel_map": "voltage_v=Voltage(V);raw_current_a=Current(A);step_capacity_ah=Capacity(Ah);signed_total_capacity_ah=Total Capacity(Ah);temperature_c=Temp*",
        "finite_value_count": finite_count,
        "nonfinite_value_count": nonfinite_count,
        "voltage_min_v": float(np.min(voltage)),
        "voltage_max_v": float(np.max(voltage)),
        "current_min_a": float(np.min(current)),
        "current_max_a": float(np.max(current)),
        "integrated_signed_charge_ah_provisional": integrated_signed,
        "integrated_absolute_charge_ah_provisional": integrated_absolute,
        "counter_signed_change_ah": counter_signed,
        "counter_absolute_throughput_ah": counter_absolute,
        "counter_discharge_throughput_ah": float(np.sum(-np.diff(total_capacity)[np.diff(total_capacity) < 0.0])),
        "counter_charge_throughput_ah": float(np.sum(np.diff(total_capacity)[np.diff(total_capacity) > 0.0])),
        "signed_integration_minus_counter_ah": integrated_signed - counter_signed,
        "absolute_integration_minus_counter_ah": integrated_absolute - counter_absolute,
        "main_discharge_counter_capacity_ah": float(main_discharge["counter_capacity_ah"]) if main_discharge else float("nan"),
        "main_discharge_integrated_capacity_ah": float(main_discharge["integrated_capacity_ah"]) if main_discharge else float("nan"),
        "main_charge_counter_capacity_ah": float(main_charge["counter_capacity_ah"]) if main_charge else float("nan"),
        "main_charge_integrated_capacity_ah": float(main_charge["integrated_capacity_ah"]) if main_charge else float("nan"),
        "measured_capacity_ah": measured_capacity,
        "nominal_reference_difference_ah": measured_capacity - NOMINAL_REFERENCE_AH if np.isfinite(measured_capacity) else float("nan"),
        "capacity_sequence_appears_complete": capacity_complete,
        "coulomb_counted_reference_soc_available": capacity_complete,
        "soc_label": "coulomb-counted reference SOC" if capacity_complete else "not justified",
        "main_discharge_start_voltage_v": float(main_discharge["start_voltage_v"]) if main_discharge else float("nan"),
        "main_discharge_end_voltage_v": float(main_discharge["end_voltage_v"]) if main_discharge else float("nan"),
        "main_discharge_preceded_by_rest": bool(main_discharge["preceded_by_rest"]) if main_discharge else False,
        "main_discharge_followed_by_rest": bool(main_discharge["followed_by_rest"]) if main_discharge else False,
        "r_app_transition_count": transition_summary["count"],
        "r_app_median_ohm": transition_summary["median"],
        "r_app_mean_ohm": transition_summary["mean"],
        "r_app_std_ohm": transition_summary["std"],
        "r_app_min_ohm": transition_summary["min"],
        "r_app_max_ohm": transition_summary["max"],
        "r_app_iqr_ohm": transition_summary["iqr"],
        "discharge_r_app_count": discharge_summary["count"],
        "discharge_r_app_median_ohm": discharge_summary["median"],
        "discharge_r_app_iqr_ohm": discharge_summary["iqr"],
        "transition_current_min_abs_a": min((abs(float(row["delta_current_a"])) for row in transitions), default=float("nan")),
        "transition_current_max_abs_a": max((abs(float(row["delta_current_a"])) for row in transitions), default=float("nan")),
        "transition_sign_inconsistent_count": sign_inconsistent,
        "max_voltage_jump_v": float(np.max(np.abs(delta_v))),
        "max_voltage_jump_at_small_current_change_v": max_unexplained_jump,
        "voltage_extreme_near_repeat_count": repeated_extrema,
        "temperature_channels": ";".join(temp_columns),
        "temperature_min_c": float(np.nanmin(temp_values)) if temp_values.size else float("nan"),
        "temperature_max_c": float(np.nanmax(temp_values)) if temp_values.size else float("nan"),
        "temperature_active_runs_showing_net_heating": heating_runs,
        "temperature_rest_runs_showing_net_cooling": cooling_runs,
    }
    output.update(temperature_diagnostics)
    return output
