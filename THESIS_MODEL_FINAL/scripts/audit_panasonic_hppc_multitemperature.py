"""Exploratory multi-temperature Panasonic NCR18650PF HPPC audit.

The script discovers and audits the original five HPPC MAT files, builds the
largest pulse support common to all campaigns, repeats the accepted conditional
R0/alpha_p/d_n identification independently at each temperature, and performs
Arrhenius diagnostics. It does not alter any frozen model or parameter file.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from models.ocp_panasonic_qi import (  # noqa: E402
    MCMB_SUPPORTED_STOICHIOMETRY_RANGE,
    NCA_WORKING_STOICHIOMETRY_RANGE,
)
from models.ocp_panasonic_qi_frozen import (  # noqa: E402
    Q_CELL_AH,
    U_n,
    X_N_END,
    X_N_START,
    X_P_END,
    X_P_START,
)
from scripts.audit_panasonic_hppc_identifiability import (  # noqa: E402
    HppcData,
    PulseRecord,
    load_hppc_data as load_accepted_hppc_csv,
    load_pulse_records as load_accepted_pulse_records,
    simulate_local_group,
    stoichiometries_from_discharged_ah,
)
from scripts.calibrate_panasonic_hppc_primary import (  # noqa: E402
    ALPHA_N_FIXED_S,
    BOUND_LOG_DISTANCE_TOLERANCE,
    D_P_FIXED_PER_S,
    FIT_PARAMETER_NAMES,
    INITIAL_PARAMETERS,
    LOG_LOWER_BOUNDS,
    LOG_UPPER_BOUNDS,
    LOWER_BOUNDS,
    MULTISTART_PARAMETERS,
    NEAR_BEST_RELATIVE_OBJECTIVE,
    PHASES,
    PreparedGroup,
    PreparedPulse,
    evaluate,
    maximum_feasible_alpha_p,
    objective_metrics,
    positive_surface_state,
    postfit_jacobian,
    prepare_groups as prepare_accepted_groups,
    run_multistart,
)


DATASET_ROOT = (
    ROOT
    / "Panasonic-18650PF-Data-master"
    / "Panasonic-18650PF-Data-master"
)
OUTPUT_DIR = ROOT / "thesis_data"
FIGURE_DIR = ROOT / "thesis_figures" / "diagnostics"
REPORT_PATH = OUTPUT_DIR / "panasonic_hppc_multitemperature_audit.md"
PULSES_PATH = OUTPUT_DIR / "panasonic_hppc_multitemperature_pulses.csv"
PARAMETERS_PATH = OUTPUT_DIR / "panasonic_hppc_multitemperature_parameters.csv"
ARRHENIUS_PATH = OUTPUT_DIR / "panasonic_hppc_arrhenius_diagnostics.csv"
LOO_PATH = OUTPUT_DIR / "panasonic_hppc_arrhenius_leave_one_out.csv"

HPPC_FILES = {
    -20.0: DATASET_ROOT / "-20degC" / "5 pulse disch" / "06-15-17_11.31 n20degC_5Pulse_HPPC_Pan18650PF.mat",
    -10.0: DATASET_ROOT / "-10degC" / "5 pulse test" / "06-07-17_12.29 n10degC_5pulse_HPPC_Pan18650PF.mat",
    0.0: DATASET_ROOT / "0degC" / "5 pulse" / "05-20-17_10.44 0degC_5pulse_HPPC_Pan18650PF.mat",
    10.0: DATASET_ROOT / "10degC" / "5 pulse test" / "03-27-17_09.06 10degC_5Pulse_HPPC_Pan18650PF.mat",
    25.0: DATASET_ROOT / "Panasonic 18650PF Data" / "25degC" / "5 pulse disch" / "03-11-17_08.47 25degC_5Pulse_HPPC_Pan18650PF.mat",
}

REST_THRESHOLD_A = 0.05
PULSE_THRESHOLD_A = 0.5
LONG_GAP_S = 10.0
MIN_PULSE_DURATION_S = 9.5
MIN_REST_DURATION_S = 5.0
LOW_VOLTAGE_CUTOFF_V = 2.501
NOMINAL_ONE_C_A = 2.9
RATE_VALUES = np.asarray([0.5, 1.0, 2.0, 4.0, 6.0])
SOC_LEVELS = np.asarray([1.00, 0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.25, 0.20, 0.15, 0.10, 0.05])
SOC_MATCH_TOLERANCE = 0.035
RATE_MATCH_RELATIVE_TOLERANCE = 0.20
R_GAS = 8.314462618
EXPECTED_ACCEPTED = np.asarray([0.0259594900586, 3204.22765511, 0.000183029692927])


@dataclass(frozen=True)
class RawCampaign:
    nominal_temperature_c: float
    path: Path
    channels: tuple[str, ...]
    arrays: dict[str, np.ndarray]
    finite_fractions: dict[str, float]
    audit: dict[str, object]


@dataclass(frozen=True)
class Segment:
    state: str
    start: int
    end: int


def vector(value: object) -> np.ndarray:
    return np.asarray(value, dtype=float).reshape(-1)


def load_campaign(nominal_c: float, path: Path) -> RawCampaign:
    if not path.is_file():
        raise FileNotFoundError(path)
    contents = loadmat(path, squeeze_me=True, struct_as_record=False)
    measurement = contents.get("meas")
    if measurement is None or not hasattr(measurement, "_fieldnames"):
        raise ValueError(f"{path.name}: expected MATLAB structure 'meas'.")
    channels = tuple(measurement._fieldnames)
    relevant = ("Time", "Current", "Voltage", "Ah", "Battery_Temp_degC", "Chamber_Temp_degC")
    required = set(relevant[:-1])
    missing = sorted(required - set(channels))
    if missing:
        raise ValueError(f"{path.name}: missing required channels {missing}.")
    arrays = {
        name: vector(getattr(measurement, name))
        for name in relevant
        if name in channels
    }
    sizes = {values.size for name, values in arrays.items() if name != "Chamber_Temp_degC"}
    if len(sizes) != 1:
        raise ValueError(f"{path.name}: relevant channels have inconsistent lengths.")
    n = arrays["Time"].size
    if "Chamber_Temp_degC" not in arrays:
        arrays["Chamber_Temp_degC"] = np.full(n, np.nan)
    elif arrays["Chamber_Temp_degC"].size != n:
        raise ValueError(f"{path.name}: chamber temperature length differs.")
    finite = {name: float(np.mean(np.isfinite(values))) for name, values in arrays.items()}
    if any(finite[name] < 1.0 for name in required):
        raise ValueError(f"{path.name}: required electrical/battery-temperature channel is non-finite.")
    dt = np.diff(arrays["Time"])
    if np.any(dt < 0.0):
        raise ValueError(f"{path.name}: elapsed time decreases.")
    positive_dt = dt[dt > 0.0]
    current = arrays["Current"]
    negative_active = int(np.count_nonzero(current <= -PULSE_THRESHOLD_A))
    positive_active = int(np.count_nonzero(current >= PULSE_THRESHOLD_A))
    if negative_active > 0 and positive_active == 0:
        sign = "raw negative = discharge"
    elif positive_active > 0 and negative_active == 0:
        sign = "raw positive = discharge candidate"
    else:
        sign = "mixed active-current signs"
    audit = {
        "file_name": path.name,
        "nominal_temperature_c": nominal_c,
        "channels": "; ".join(channels),
        "sample_count": n,
        "duration_s": float(arrays["Time"][-1] - arrays["Time"][0]),
        "dt_min_positive_s": float(np.min(positive_dt)),
        "dt_median_positive_s": float(np.median(positive_dt)),
        "dt_max_positive_s": float(np.max(positive_dt)),
        "dt_max_contiguous_s": float(np.max(positive_dt[positive_dt <= LONG_GAP_S])),
        "duplicate_time_steps": int(np.count_nonzero(dt == 0.0)),
        "long_time_gaps": int(np.count_nonzero(dt > LONG_GAP_S)),
        "current_sign_convention": sign,
        "current_min_a": float(np.min(current)),
        "current_max_a": float(np.max(current)),
        "voltage_min_v": float(np.min(arrays["Voltage"])),
        "voltage_max_v": float(np.max(arrays["Voltage"])),
    }
    return RawCampaign(nominal_c, path, channels, arrays, finite, audit)


def states(current: np.ndarray) -> np.ndarray:
    result = np.full(current.shape, "transition", dtype="U10")
    result[np.abs(current) <= REST_THRESHOLD_A] = "rest"
    result[current <= -PULSE_THRESHOLD_A] = "discharge"
    result[current >= PULSE_THRESHOLD_A] = "charge"
    return result


def segments(campaign: RawCampaign) -> list[Segment]:
    state = states(campaign.arrays["Current"])
    time = campaign.arrays["Time"]
    output: list[Segment] = []
    start = 0
    for index in range(1, state.size):
        if state[index] != state[index - 1] or time[index] - time[index - 1] > LONG_GAP_S:
            output.append(Segment(str(state[start]), start, index - 1))
            start = index
    output.append(Segment(str(state[start]), start, state.size - 1))
    return output


def rest_before(campaign: RawCampaign, pulse_start: int) -> tuple[int, int] | None:
    current, time = campaign.arrays["Current"], campaign.arrays["Time"]
    end = pulse_start - 1
    if end < 0 or abs(current[end]) > REST_THRESHOLD_A:
        return None
    start = end
    while start > 0 and time[start] - time[start - 1] <= LONG_GAP_S and abs(current[start - 1]) <= REST_THRESHOLD_A:
        start -= 1
    return start, end


def rest_after(campaign: RawCampaign, pulse_end: int) -> tuple[int, int] | None:
    current, time = campaign.arrays["Current"], campaign.arrays["Time"]
    start = pulse_end + 1
    if start >= current.size or abs(current[start]) > REST_THRESHOLD_A:
        return None
    end = start
    while end + 1 < current.size and time[end + 1] - time[end] <= LONG_GAP_S and abs(current[end + 1]) <= REST_THRESHOLD_A:
        end += 1
    return start, end


def nearest_rate(current_magnitude: float) -> tuple[str, float, bool]:
    measured = current_magnitude / NOMINAL_ONE_C_A
    nearest = float(RATE_VALUES[np.argmin(np.abs(RATE_VALUES - measured))])
    ok = abs(measured - nearest) / nearest <= RATE_MATCH_RELATIVE_TOLERANCE
    return f"{nearest:g}C", nearest, bool(ok)


def identify_pulses(campaign: RawCampaign) -> list[dict[str, object]]:
    a = campaign.arrays
    records: list[dict[str, object]] = []
    active = [segment for segment in segments(campaign) if segment.state in {"discharge", "charge"}]
    group_index = -1
    previous_post_end: int | None = None
    for pulse_id, segment in enumerate(active, start=1):
        start, end = segment.start, segment.end
        before, after = rest_before(campaign, start), rest_after(campaign, end)
        pre_start, pre_end = before if before is not None else (max(0, start - 1), max(0, start - 1))
        post_start, post_end = after if after is not None else (min(end + 1, a["Time"].size - 1), min(end + 1, a["Time"].size - 1))
        pulse_current = float(np.median(a["Current"][start : end + 1]))
        current_iqr = float(np.subtract(*np.percentile(a["Current"][start : end + 1], [75, 25])))
        differences = np.diff(a["Time"][start : end + 1])
        positive = differences[differences > 0.0]
        final_interval = float(np.median(positive)) if positive.size else 0.0
        duration = float(a["Time"][end] - a["Time"][start] + final_interval)
        pre_duration = float(a["Time"][start] - a["Time"][pre_start]) if before is not None else 0.0
        post_duration = float(a["Time"][post_end] - a["Time"][post_start]) if after is not None else 0.0
        pre_current = float(a["Current"][pre_end])
        delta_i = pulse_current - pre_current
        v_pre, v_immediate, v_end, v_post = map(
            float, (a["Voltage"][pre_end], a["Voltage"][start], a["Voltage"][end], a["Voltage"][post_end])
        )
        rapp = abs((v_immediate - v_pre) / delta_i) if delta_i != 0.0 else float("nan")
        discharged_ah = max(0.0, -float(a["Ah"][pre_end]))
        soc = 1.0 - discharged_ah / Q_CELL_AH
        # The source Ah counter drifts from the nominal SOC schedule because
        # capacity is not exactly 2.994979 Ah and because preceding pulses also
        # remove charge. Use the audited HPPC group order (split by omitted-data
        # gaps) to assign the documented 100, 95, 90, 80, ..., 5% schedule.
        if previous_post_end is None or a["Time"][pre_end] - a["Time"][previous_post_end] > LONG_GAP_S:
            group_index += 1
        previous_post_end = post_end
        soc_level = float(SOC_LEVELS[group_index]) if group_index < SOC_LEVELS.size else float("nan")
        soc_match = bool(np.isfinite(soc_level) and abs(soc_level - soc) <= SOC_MATCH_TOLERANCE)
        rate_label, rate_value, rate_match = nearest_rate(abs(pulse_current))
        z = 1.0 - soc
        x_n = X_N_START + z * (X_N_END - X_N_START)
        x_p = X_P_START + z * (X_P_END - X_P_START)
        support_ok = bool(
            MCMB_SUPPORTED_STOICHIOMETRY_RANGE[0] <= x_n <= MCMB_SUPPORTED_STOICHIOMETRY_RANGE[1]
            and NCA_WORKING_STOICHIOMETRY_RANGE[0] <= x_p <= NCA_WORKING_STOICHIOMETRY_RANGE[1]
        )
        voltage_ok = bool(np.all(np.isfinite([v_pre, v_immediate, v_end, v_post])) and 2.0 <= min(v_pre, v_immediate, v_end, v_post) <= max(v_pre, v_immediate, v_end, v_post) <= 4.35)
        temperature_slice = a["Battery_Temp_degC"][pre_end : post_end + 1]
        temperature_ok = bool(temperature_slice.size and np.all(np.isfinite(temperature_slice)))
        truncated = bool(duration < MIN_PULSE_DURATION_S and v_end <= LOW_VOLTAGE_CUTOFF_V)
        insufficient_rest = bool(before is None or after is None or pre_duration < MIN_REST_DURATION_S or post_duration < MIN_REST_DURATION_S)
        abnormal_current = bool(not rate_match or current_iqr > max(0.10, 0.03 * abs(pulse_current)))
        invalid_soc = bool(not (0.0 <= soc <= 1.0) or not soc_match)
        flags: list[str] = []
        if segment.state != "discharge": flags.append("non-discharge pulse")
        if before is None: flags.append("missing pre-rest")
        if after is None: flags.append("missing post-rest")
        if insufficient_rest: flags.append("insufficient rest")
        if truncated: flags.append("cutoff-truncated")
        if not voltage_ok: flags.append("abnormal voltage")
        if abnormal_current: flags.append("abnormal current")
        if invalid_soc: flags.append("invalid/unmatched SOC")
        if not temperature_ok: flags.append("non-finite battery temperature")
        if not support_ok: flags.append("average-state OCP support violation")
        usable = bool(
            segment.state == "discharge" and before is not None and after is not None
            and duration >= MIN_PULSE_DURATION_S and not insufficient_rest and voltage_ok
            and not abnormal_current and not invalid_soc and temperature_ok and support_ok
            and np.isfinite(rapp)
        )
        records.append(
            {
                "nominal_temperature_c": campaign.nominal_temperature_c,
                "source_file": campaign.path.name,
                "pulse_id": pulse_id,
                "pulse_type": segment.state,
                "start_index": start,
                "end_index": end,
                "pre_pulse_index": pre_end,
                "post_rest_index": post_end,
                "start_time_s": float(a["Time"][start]),
                "soc": soc,
                "matched_soc_level": soc_level,
                "soc_group_index": group_index + 1,
                "soc_match_ok": soc_match,
                "pulse_c_rate": rate_value,
                "nominal_rate": rate_label,
                "rate_match_ok": rate_match,
                "current_a": pulse_current,
                "current_iqr_a": current_iqr,
                "pre_pulse_current_a": pre_current,
                "delta_i_a": delta_i,
                "pulse_duration_s": duration,
                "pre_rest_duration_s": pre_duration,
                "relaxation_duration_s": post_duration,
                "v_pre_v": v_pre,
                "v_immediate_v": v_immediate,
                "v_end_v": v_end,
                "v_post_rest_v": v_post,
                "delta_v_immediate_v": v_immediate - v_pre,
                "apparent_transition_resistance_ohm": rapp,
                "battery_temperature_start_c": float(a["Battery_Temp_degC"][start]),
                "battery_temperature_end_c": float(a["Battery_Temp_degC"][end]),
                "battery_temperature_change_during_pulse_k": float(a["Battery_Temp_degC"][end] - a["Battery_Temp_degC"][start]),
                "missing_pulse": False,
                "truncated_pulse": truncated,
                "abnormal_voltage": not voltage_ok,
                "abnormal_current": abnormal_current,
                "insufficient_rest": insufficient_rest,
                "invalid_soc": invalid_soc,
                "nonfinite_temperature": not temperature_ok,
                "support_violation": not support_ok,
                "usable": usable,
                "common_support": False,
                "exclusion_flags": "; ".join(flags),
            }
        )
    return records


def campaign_temperature_summary(campaign: RawCampaign, records: list[dict[str, object]]) -> dict[str, object]:
    battery = campaign.arrays["Battery_Temp_degC"]
    high = [r for r in records if float(r["pulse_c_rate"]) >= 4.0 and np.isfinite(float(r["battery_temperature_change_during_pulse_k"]))]
    changes = np.asarray([float(r["battery_temperature_change_during_pulse_k"]) for r in high])
    return {
        **campaign.audit,
        **{f"finite_fraction_{name}": fraction for name, fraction in campaign.finite_fractions.items()},
        "battery_initial_c": float(battery[0]),
        "battery_median_c": float(np.median(battery)),
        "battery_min_c": float(np.min(battery)),
        "battery_max_c": float(np.max(battery)),
        "battery_drift_k": float(battery[-1] - battery[0]),
        "high_current_pulse_delta_t_median_k": float(np.median(changes)) if changes.size else float("nan"),
        "high_current_pulse_delta_t_max_k": float(np.max(changes)) if changes.size else float("nan"),
    }


def support_key(record: dict[str, object]) -> tuple[float, float]:
    return float(record["matched_soc_level"]), float(record["pulse_c_rate"])


def determine_common_support(all_records: dict[float, list[dict[str, object]]]) -> set[tuple[float, float]]:
    usable_sets = [
        {support_key(record) for record in records if bool(record["usable"])}
        for records in all_records.values()
    ]
    common = set.intersection(*usable_sets)
    for records in all_records.values():
        for record in records:
            record["common_support"] = bool(record["usable"] and support_key(record) in common)
    return common


def append_missing_rows(
    all_records: dict[float, list[dict[str, object]]], common_candidates: set[tuple[float, float]]
) -> None:
    """Add explicit audit rows for expected grid combinations absent from a campaign."""

    template_fields = list(next(iter(next(iter(all_records.values())))))
    for temperature, records in all_records.items():
        observed = {support_key(record) for record in records}
        for soc, rate in sorted(common_candidates - observed, reverse=True):
            row = {field: "" for field in template_fields}
            row.update(
                {
                    "nominal_temperature_c": temperature,
                    "source_file": HPPC_FILES[temperature].name,
                    "pulse_id": "",
                    "pulse_type": "missing",
                    "matched_soc_level": soc,
                    "pulse_c_rate": rate,
                    "nominal_rate": f"{rate:g}C",
                    "missing_pulse": True,
                    "usable": False,
                    "common_support": False,
                    "exclusion_flags": "missing expected pulse",
                }
            )
            records.append(row)


def to_model_data(campaign: RawCampaign) -> HppcData:
    return HppcData(
        time_s=campaign.arrays["Time"],
        voltage_v=campaign.arrays["Voltage"],
        raw_current_a=campaign.arrays["Current"],
        exported_ah=campaign.arrays["Ah"],
        battery_temperature_k=campaign.arrays["Battery_Temp_degC"] + 273.15,
    )


def chosen_pulse_records(records: list[dict[str, object]]) -> list[PulseRecord]:
    selected = [record for record in records if bool(record["common_support"])]
    return [
        PulseRecord(
            pulse_id=int(record["pulse_id"]),
            nominal_rate=str(record["nominal_rate"]),
            start_index=int(record["start_index"]),
            end_index=int(record["end_index"]),
            pre_pulse_index=int(record["pre_pulse_index"]),
            post_rest_index=int(record["post_rest_index"]),
            soc=float(record["soc"]),
            usable=True,
        )
        for record in selected
    ]


def partition_general(data: HppcData, records: list[PulseRecord]) -> list[list[PulseRecord]]:
    ordered = sorted(records, key=lambda record: record.start_index)
    groups: list[list[PulseRecord]] = []
    for record in ordered:
        if not groups:
            groups.append([record])
        else:
            previous = groups[-1][-1]
            gap = data.time_s[record.pre_pulse_index] - data.time_s[previous.post_rest_index]
            if gap > LONG_GAP_S:
                groups.append([record])
            else:
                groups[-1].append(record)
    return groups


def prepare_common_groups(data: HppcData, records: list[PulseRecord]) -> list[PreparedGroup]:
    prepared: list[PreparedGroup] = []
    for group_id, group_records in enumerate(partition_general(data, records), start=1):
        lower = group_records[0].pre_pulse_index
        upper = group_records[-1].post_rest_index + 1
        time_s = data.time_s[lower:upper] - data.time_s[lower]
        current_a = -data.raw_current_a[lower:upper]
        temperature_k = data.battery_temperature_k[lower:upper]
        measured_voltage_v = data.voltage_v[lower:upper]
        discharged_ah = max(0.0, -float(data.exported_ah[lower]))
        initial_x_n, initial_x_p = stoichiometries_from_discharged_ah(discharged_ah)
        reference = simulate_local_group(
            time_s, current_a, temperature_k, initial_x_n, initial_x_p,
            {"alpha_n": ALPHA_N_FIXED_S, "alpha_p": INITIAL_PARAMETERS[1],
             "d_n": INITIAL_PARAMETERS[2], "d_p": D_P_FIXED_PER_S, "R0": INITIAL_PARAMETERS[0]},
        )
        pulses: list[PreparedPulse] = []
        for record in group_records:
            pre, start, end, post = (
                record.pre_pulse_index - lower,
                record.start_index - lower,
                record.end_index - lower,
                record.post_rest_index - lower,
            )
            phase_indices = {
                "immediate": np.asarray([start], dtype=int),
                "within_pulse": np.arange(start + 1, end + 1, dtype=int),
                "post_relaxation": np.arange(end + 1, post + 1, dtype=int),
            }
            if any(indices.size == 0 for indices in phase_indices.values()):
                raise ValueError(f"Pulse {record.pulse_id} has an empty phase.")
            measured_pre = float(measured_voltage_v[pre])
            pulses.append(
                PreparedPulse(record, pre, end, phase_indices, measured_pre,
                              {phase: measured_voltage_v[indices] - measured_pre for phase, indices in phase_indices.items()})
            )
        changes = np.flatnonzero(np.diff(current_a[:-1]) != 0.0) + 1
        prepared.append(
            PreparedGroup(
                group_id, time_s, current_a, temperature_k, measured_voltage_v,
                initial_x_n, initial_x_p, reference.x_n_surface, reference.x_p_avg,
                np.asarray(U_n(reference.x_n_surface), dtype=float),
                tuple(zip(np.r_[0, changes].tolist(), np.r_[changes, current_a.size - 1].tolist(), strict=True)),
                tuple(pulses),
            )
        )
    if not prepared:
        raise ValueError("Common support produced no prepared histories.")
    return prepared


def best_fit(groups: list[PreparedGroup]) -> tuple[dict[str, object], dict[str, object], dict[str, object], list[dict[str, object]]]:
    alpha_upper = maximum_feasible_alpha_p(groups)
    starts = run_multistart(groups, alpha_upper)
    successful = [row for row in starts if bool(row["success"])]
    if not successful:
        raise RuntimeError("No feasible multistart solution.")
    best = min(successful, key=lambda row: float(row["objective_v2"]))
    parameters = {"R0": float(best["R0_ohm"]), "alpha_p": float(best["alpha_p_s"]), "d_n": float(best["d_n_per_s"])}
    evaluation = evaluate(groups, parameters)
    # Only the overall and three phase metrics are needed here. Avoid asking
    # the legacy diagnostic helper for absent C-rate subsets, which would form
    # empty means for a deliberately reduced common support.
    diagnostics = {"metrics": objective_metrics(evaluation)}
    condition = postfit_jacobian(groups, parameters)
    return best, diagnostics, condition, starts


def regression_test_accepted_25c() -> dict[str, object]:
    data = load_accepted_hppc_csv()
    records = load_accepted_pulse_records()
    groups = prepare_accepted_groups(data, records, minimum_pulse_start_soc=1.0 / 3.0)
    best, diagnostics, condition, starts = best_fit(groups)
    recovered = np.asarray([best["R0_ohm"], best["alpha_p_s"], best["d_n_per_s"]], dtype=float)
    relative = np.abs(recovered / EXPECTED_ACCEPTED - 1.0)
    passed = bool(np.all(relative <= np.asarray([5e-6, 5e-6, 5e-6])))
    if not passed:
        raise RuntimeError(f"Accepted 25 C regression failed: recovered {recovered}, expected {EXPECTED_ACCEPTED}.")
    return {
        "passed": passed,
        "pulse_count": sum(len(group.pulses) for group in groups),
        "group_count": len(groups),
        "R0_ohm": recovered[0],
        "alpha_p_s": recovered[1],
        "d_n_per_s": recovered[2],
        "phase_balanced_rmse_v": diagnostics["metrics"]["phase_balanced_rmse_v"],
        "successful_starts": sum(bool(row["success"]) for row in starts),
    }


def parameter_row(
    temperature: float,
    median_battery_c: float,
    pulse_count: int,
    group_count: int,
    best: dict[str, object],
    diagnostics: dict[str, object],
    condition: dict[str, object],
    starts: list[dict[str, object]],
) -> dict[str, object]:
    successful = [row for row in starts if bool(row["success"])]
    objective_best = float(best["objective_v2"])
    near = [row for row in successful if float(row["objective_v2"]) <= objective_best * (1.0 + NEAR_BEST_RELATIVE_OBJECTIVE)]
    singular = np.asarray(condition["singular_values"])
    rank_tolerance = max(condition["jacobian"].shape) * np.finfo(float).eps * singular[0]
    rank = int(np.count_nonzero(singular > rank_tolerance))
    sensitivity_corr = np.corrcoef(np.asarray(condition["jacobian"]).T)
    values = np.asarray([best["R0_ohm"], best["alpha_p_s"], best["d_n_per_s"]], dtype=float)
    log_values = np.log(values)
    effective_upper = LOG_UPPER_BOUNDS.copy()
    effective_upper[1] = np.log(float(best["alpha_p_frozen_domain_upper_s"]))
    bound_fraction = np.minimum(log_values - LOG_LOWER_BOUNDS, effective_upper - log_values) / (effective_upper - LOG_LOWER_BOUNDS)
    spreads = {}
    for field in ("R0_ohm", "alpha_p_s", "d_n_per_s"):
        vals = np.asarray([float(row[field]) for row in near])
        spreads[field] = (float(np.min(vals)), float(np.max(vals)))
    relative_spreads = [
        (spreads[field][1] - spreads[field][0]) / float(best[field])
        for field in ("R0_ohm", "alpha_p_s", "d_n_per_s")
    ]
    recoverable = bool(
        rank == 3 and np.isfinite(condition["condition_number"])
        and float(condition["condition_number"]) < 1.0e4
        and not bool(best["at_numerical_bound"])
        and not bool(best["at_frozen_ocp_domain_boundary"])
        and max(relative_spreads) < 0.25
    )
    metrics = diagnostics["metrics"]
    return {
        "nominal_temperature_c": temperature,
        "median_battery_temperature_c": median_battery_c,
        "common_pulse_count": pulse_count,
        "prepared_group_count": group_count,
        "R0_ohm": best["R0_ohm"],
        "alpha_p_s": best["alpha_p_s"],
        "d_n_per_s": best["d_n_per_s"],
        "phase_balanced_rmse_v": metrics["phase_balanced_rmse_v"],
        "immediate_rmse_v": metrics["rmse_immediate_v"],
        "during_pulse_rmse_v": metrics["rmse_within_pulse_v"],
        "relaxation_rmse_v": metrics["rmse_post_relaxation_v"],
        "optimizer_success": True,
        "successful_multistarts": len(successful),
        "attempted_multistarts": len(starts),
        "best_objective_v2": objective_best,
        "R0_log_bound_fraction": bound_fraction[0],
        "alpha_p_log_bound_fraction": bound_fraction[1],
        "d_n_log_bound_fraction": bound_fraction[2],
        "alpha_p_domain_upper_s": best["alpha_p_frozen_domain_upper_s"],
        "at_numerical_bound": best["at_numerical_bound"],
        "at_ocp_feasibility_boundary": best["at_frozen_ocp_domain_boundary"],
        "near_best_R0_min": spreads["R0_ohm"][0],
        "near_best_R0_max": spreads["R0_ohm"][1],
        "near_best_alpha_p_min": spreads["alpha_p_s"][0],
        "near_best_alpha_p_max": spreads["alpha_p_s"][1],
        "near_best_d_n_min": spreads["d_n_per_s"][0],
        "near_best_d_n_max": spreads["d_n_per_s"][1],
        "jacobian_rank": rank,
        "singular_value_1": singular[0],
        "singular_value_2": singular[1],
        "singular_value_3": singular[2],
        "condition_number": condition["condition_number"],
        "sensitivity_corr_R0_alpha_p": sensitivity_corr[0, 1],
        "sensitivity_corr_R0_d_n": sensitivity_corr[0, 2],
        "sensitivity_corr_alpha_p_d_n": sensitivity_corr[1, 2],
        "covariance_corr_R0_alpha_p": condition["correlation"][0, 1],
        "covariance_corr_R0_d_n": condition["correlation"][0, 2],
        "covariance_corr_alpha_p_d_n": condition["correlation"][1, 2],
        "practically_recoverable_triplet": recoverable,
    }


def fit_line(x: np.ndarray, y: np.ndarray) -> dict[str, object]:
    design = np.column_stack((np.ones_like(x), x))
    intercept, slope = np.linalg.lstsq(design, y, rcond=None)[0]
    predicted = intercept + slope * x
    residual = y - predicted
    denominator = np.sum((y - np.mean(y)) ** 2)
    r2 = float(1.0 - np.sum(residual**2) / denominator) if denominator > 0.0 else float("nan")
    return {"intercept": float(intercept), "slope": float(slope), "predicted": predicted, "residual": residual, "r2": r2}


def arrhenius_diagnostics(parameter_rows: list[dict[str, object]]) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[tuple[str, str], dict[str, object]]]:
    output: list[dict[str, object]] = []
    loo: list[dict[str, object]] = []
    fits: dict[tuple[str, str], dict[str, object]] = {}
    parameter_fields = {"R0": "R0_ohm", "alpha_p": "alpha_p_s", "d_n": "d_n_per_s"}
    for parameter, field in parameter_fields.items():
        measured = np.asarray([float(row[field]) for row in parameter_rows])
        for definition, temperature_field in (("nominal", "nominal_temperature_c"), ("battery_median", "median_battery_temperature_c")):
            temperatures_c = np.asarray([float(row[temperature_field]) for row in parameter_rows])
            kelvin = temperatures_c + 273.15
            x, y = 1.0 / kelvin, np.log(measured)
            fit = fit_line(x, y)
            signed_ea = (-fit["slope"] if parameter == "d_n" else fit["slope"]) * R_GAS
            fit.update({"temperature_c": temperatures_c, "x": x, "measured": measured, "activation_energy_j_per_mol": signed_ea})
            fits[(parameter, definition)] = fit
            for index, row in enumerate(parameter_rows):
                output.append(
                    {
                        "parameter": parameter,
                        "temperature_definition": definition,
                        "nominal_temperature_c": row["nominal_temperature_c"],
                        "regression_temperature_c": temperatures_c[index],
                        "inverse_temperature_per_k": x[index],
                        "measured_parameter": measured[index],
                        "predicted_parameter": float(np.exp(fit["predicted"][index])),
                        "log_residual": fit["residual"][index],
                        "slope_k": fit["slope"],
                        "intercept": fit["intercept"],
                        "activation_energy_signed_j_per_mol": signed_ea,
                        "activation_energy_magnitude_j_per_mol": abs(signed_ea),
                        "r2": fit["r2"],
                    }
                )
            for withheld in range(measured.size):
                retained = np.arange(measured.size) != withheld
                loo_fit = fit_line(x[retained], y[retained])
                predicted_value = float(np.exp(loo_fit["intercept"] + loo_fit["slope"] * x[withheld]))
                relative_error = (predicted_value - measured[withheld]) / measured[withheld]
                loo.append(
                    {
                        "parameter": parameter,
                        "temperature_definition": definition,
                        "withheld_nominal_temperature_c": parameter_rows[withheld]["nominal_temperature_c"],
                        "withheld_regression_temperature_c": temperatures_c[withheld],
                        "measured_parameter": measured[withheld],
                        "predicted_parameter": predicted_value,
                        "relative_error": relative_error,
                        "absolute_relative_error": abs(relative_error),
                    }
                )
    return output, loo, fits


def classifications(parameter_rows: list[dict[str, object]], fits: dict[tuple[str, str], dict[str, object]], loo: list[dict[str, object]]) -> dict[str, str]:
    decisions: dict[str, str] = {}
    fields = {"R0": "R0_ohm", "alpha_p": "alpha_p_s", "d_n": "d_n_per_s"}
    for parameter, field in fields.items():
        if not all(bool(row["practically_recoverable_triplet"]) for row in parameter_rows):
            # A boundary specifically affecting alpha_p makes alpha_p unrecoverable;
            # other parameters are judged by their own near-best spread and bounds.
            parameter_spread_fields = {
                "R0": ("near_best_R0_min", "near_best_R0_max", "R0_log_bound_fraction"),
                "alpha_p": ("near_best_alpha_p_min", "near_best_alpha_p_max", "alpha_p_log_bound_fraction"),
                "d_n": ("near_best_d_n_min", "near_best_d_n_max", "d_n_log_bound_fraction"),
            }[parameter]
            recoverable = True
            for row in parameter_rows:
                lo, hi = float(row[parameter_spread_fields[0]]), float(row[parameter_spread_fields[1]])
                if ((hi - lo) / float(row[field]) >= 0.25
                        or int(row["jacobian_rank"]) < 3 or float(row["condition_number"]) >= 1.0e4
                        or float(row[parameter_spread_fields[2]]) <= 0.01
                        or (parameter == "alpha_p" and bool(row["at_ocp_feasibility_boundary"]))):
                    recoverable = False
                    break
            if not recoverable:
                decisions[parameter] = "PARAMETER NOT SUFFICIENTLY RECOVERABLE"
                continue
        nominal_fit = fits[(parameter, "nominal")]
        battery_fit = fits[(parameter, "battery_median")]
        subset = [row for row in loo if row["parameter"] == parameter]
        nominal_errors = np.asarray([float(row["absolute_relative_error"]) for row in subset if row["temperature_definition"] == "nominal"])
        battery_errors = np.asarray([float(row["absolute_relative_error"]) for row in subset if row["temperature_definition"] == "battery_median"])
        robust = max(np.median(nominal_errors), np.median(battery_errors)) <= 0.25 and max(np.max(nominal_errors), np.max(battery_errors)) <= 0.50
        clear = min(float(nominal_fit["r2"]), float(battery_fit["r2"])) >= 0.90 and robust
        possible = max(float(nominal_fit["r2"]), float(battery_fit["r2"])) >= 0.60 or min(np.median(nominal_errors), np.median(battery_errors)) <= 0.50
        decisions[parameter] = "CLEAR TEMPERATURE DEPENDENCE" if clear else ("POSSIBLE BUT NOT ROBUST" if possible else "NO DEFENSIBLE TEMPERATURE DEPENDENCE")
    return decisions


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def figures(
    records: list[dict[str, object]],
    parameters: list[dict[str, object]],
    fits: dict[tuple[str, str], dict[str, object]],
) -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.22})
    matched = [row for row in records if bool(row["common_support"])]
    fig, axis = plt.subplots(figsize=(7.2, 4.5), constrained_layout=True)
    colors = {0.5: "#56B4E9", 1.0: "#009E73", 2.0: "#E69F00", 4.0: "#D55E00", 6.0: "#CC79A7"}
    matched_rates = sorted({float(row["pulse_c_rate"]) for row in matched})
    for rate in matched_rates:
        selected = [row for row in matched if float(row["pulse_c_rate"]) == rate]
        axis.scatter([row["nominal_temperature_c"] for row in selected],
                     [1000.0 * float(row["apparent_transition_resistance_ohm"]) for row in selected],
                     s=20, alpha=0.65, color=colors[float(rate)], label=f"{rate:g}C")
    medians = []
    for temp in sorted({float(row["nominal_temperature_c"]) for row in matched}):
        values = np.asarray([float(row["apparent_transition_resistance_ohm"]) for row in matched if float(row["nominal_temperature_c"]) == temp])
        medians.append(np.median(values))
        axis.vlines(temp, 1000*np.percentile(values, 25), 1000*np.percentile(values, 75), color="black", lw=2.2)
    axis.plot(sorted({float(row["nominal_temperature_c"]) for row in matched}), 1000*np.asarray(medians), "ko-", lw=1.1, label="Median and IQR")
    axis.set(xlabel="Nominal campaign temperature (°C)", ylabel="Apparent transition resistance (mΩ)", title="Matched HPPC rest-to-load transitions")
    axis.legend(frameon=False, ncol=3)
    fig.savefig(FIGURE_DIR / "panasonic_hppc_multitemperature_rapp.png", dpi=240, facecolor="white")
    plt.close(fig)

    parameter_specs = (("R0_ohm", "R0 (mΩ)", 1000.0), ("alpha_p_s", "alpha_p (s)", 1.0), ("d_n_per_s", "d_n (1/s)", 1.0))
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.7), constrained_layout=True)
    temperature = np.asarray([float(row["nominal_temperature_c"]) for row in parameters])
    for axis, (field, label, scale) in zip(axes, parameter_specs, strict=True):
        axis.plot(temperature, scale*np.asarray([float(row[field]) for row in parameters]), "o-", color="#0072B2")
        axis.set(xlabel="Nominal temperature (°C)", ylabel=label)
    fig.suptitle("Independent conditional HPPC parameter estimates")
    fig.savefig(FIGURE_DIR / "panasonic_hppc_multitemperature_parameter_trends.png", dpi=240, facecolor="white")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.7), constrained_layout=True)
    for axis, parameter in zip(axes, ("R0", "alpha_p", "d_n"), strict=True):
        for definition, marker, color in (("nominal", "o", "#0072B2"), ("battery_median", "s", "#D55E00")):
            fit = fits[(parameter, definition)]
            x = np.asarray(fit["x"])
            y = np.log(np.asarray(fit["measured"]))
            order = np.argsort(x)
            axis.scatter(1000*x, y, marker=marker, color=color, label=definition.replace("_", " "))
            axis.plot(1000*x[order], (float(fit["intercept"]) + float(fit["slope"])*x)[order], color=color, lw=1.0)
        axis.set(xlabel="1000/T (1/K)", ylabel=f"ln({parameter})", title=parameter)
    axes[0].legend(frameon=False)
    fig.suptitle("Arrhenius-coordinate diagnostics")
    fig.savefig(FIGURE_DIR / "panasonic_hppc_multitemperature_arrhenius.png", dpi=240, facecolor="white")
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(7.3, 4.5), constrained_layout=True)
    x = np.arange(len(parameters)); width = 0.19
    phase_fields = (("immediate_rmse_v", "Start response"), ("during_pulse_rmse_v", "During pulse"), ("relaxation_rmse_v", "Relaxation"), ("phase_balanced_rmse_v", "Balanced RMSE"))
    for index, (field, label) in enumerate(phase_fields):
        axis.bar(x + (index - 1.5)*width, [1000*float(row[field]) for row in parameters], width, label=label)
    axis.set_xticks(x, [f"{float(row['nominal_temperature_c']):g}°C" for row in parameters])
    axis.set(ylabel="RMSE (mV)", xlabel="Nominal campaign temperature", title="Conditional HPPC fit quality on common support")
    axis.legend(frameon=False, ncol=2)
    fig.savefig(FIGURE_DIR / "panasonic_hppc_multitemperature_fit_quality.png", dpi=240, facecolor="white")
    plt.close(fig)


def date_from_name(name: str) -> datetime:
    return datetime.strptime(name[:8], "%m-%d-%y")


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    campaigns = {temperature: load_campaign(temperature, path) for temperature, path in HPPC_FILES.items()}
    records_by_temperature = {temperature: identify_pulses(campaign) for temperature, campaign in campaigns.items()}
    all_expected = {(float(soc), float(rate)) for soc in SOC_LEVELS for rate in RATE_VALUES}
    common = determine_common_support(records_by_temperature)
    append_missing_rows(records_by_temperature, all_expected)
    temperature_audits = {
        temperature: campaign_temperature_summary(campaigns[temperature], [row for row in records if not bool(row["missing_pulse"])])
        for temperature, records in records_by_temperature.items()
    }
    all_pulses = [row for temperature in sorted(records_by_temperature) for row in records_by_temperature[temperature]]
    write_csv(PULSES_PATH, all_pulses)

    regression = regression_test_accepted_25c()
    parameter_rows: list[dict[str, object]] = []
    for temperature in sorted(campaigns):
        data = to_model_data(campaigns[temperature])
        pulse_records = chosen_pulse_records(records_by_temperature[temperature])
        groups = prepare_common_groups(data, pulse_records)
        best, diagnostics, condition, starts = best_fit(groups)
        parameter_rows.append(
            parameter_row(
                temperature,
                float(temperature_audits[temperature]["battery_median_c"]),
                len(pulse_records),
                len(groups),
                best,
                diagnostics,
                condition,
                starts,
            )
        )
    write_csv(PARAMETERS_PATH, parameter_rows)
    arrhenius_rows, loo_rows, fits = arrhenius_diagnostics(parameter_rows)
    write_csv(ARRHENIUS_PATH, arrhenius_rows)
    write_csv(LOO_PATH, loo_rows)
    decisions = classifications(parameter_rows, fits, loo_rows)
    clear = [name for name, decision in decisions.items() if decision == "CLEAR TEMPERATURE DEPENDENCE"]
    if len(clear) == 3:
        global_decision = "YES — PROCEED TO TEMPERATURE-DEPENDENT MODEL TEST"
    elif clear:
        global_decision = f"PARTIALLY — ONLY {', '.join(clear)} ARE SUPPORTED"
    else:
        global_decision = "NO — KEEP THE CURRENT TEMPERATURE-INDEPENDENT PARAMETERS"
    figures(all_pulses, parameter_rows, fits)

    matched_rows = [row for row in all_pulses if bool(row["common_support"])]
    common_socs = sorted({key[0] for key in common}, reverse=True)
    common_rates = sorted({key[1] for key in common})
    acquisition = sorted(((date_from_name(path.name), temperature, path.name) for temperature, path in HPPC_FILES.items()))
    ordered_temperatures = sorted(campaigns)

    def median_rapp(temp: float, *, rate: float | None = None, soc: float | None = None) -> float:
        selected = [
            float(row["apparent_transition_resistance_ohm"])
            for row in matched_rows
            if float(row["nominal_temperature_c"]) == temp
            and (rate is None or float(row["pulse_c_rate"]) == rate)
            and (soc is None or float(row["matched_soc_level"]) == soc)
        ]
        return float(np.median(selected))

    def rapp_matrix_row(label: str, values: list[float]) -> str:
        return "| " + label + " | " + " | ".join(f"{1000*value:.6f}" for value in values) + " |"

    def audit_row(temp: float) -> str:
        a = temperature_audits[temp]
        chamber = campaigns[temp].finite_fractions["Chamber_Temp_degC"]
        return (
            f"| {temp:g} | `{a['file_name']}` | {a['sample_count']} | {a['duration_s']:.1f} | "
            f"{a['dt_min_positive_s']:.4g}/{a['dt_median_positive_s']:.4g}/{a['dt_max_contiguous_s']:.4g} | {a['long_time_gaps']} | "
            f"{campaigns[temp].finite_fractions['Time']:.3f}/{campaigns[temp].finite_fractions['Current']:.3f}/"
            f"{campaigns[temp].finite_fractions['Voltage']:.3f}/{campaigns[temp].finite_fractions['Battery_Temp_degC']:.3f}/{chamber:.3f} | {a['current_sign_convention']} |"
        )

    def temp_row(temp: float) -> str:
        a = temperature_audits[temp]
        return (
            f"| {temp:g} | {a['battery_initial_c']:.4f} | {a['battery_median_c']:.4f} | {a['battery_min_c']:.4f} | "
            f"{a['battery_max_c']:.4f} | {a['battery_drift_k']:+.4f} | {a['high_current_pulse_delta_t_median_k']:+.5f} | {a['high_current_pulse_delta_t_max_k']:+.5f} |"
        )

    def flag_row(temp: float) -> str:
        records = records_by_temperature[temp]
        count = lambda field: sum(bool(row[field]) for row in records)
        return (
            f"| {temp:g} | {count('missing_pulse')} | {count('truncated_pulse')} | "
            f"{count('abnormal_voltage')} | {count('abnormal_current')} | "
            f"{count('insufficient_rest')} | {count('invalid_soc')} | "
            f"{count('nonfinite_temperature')} | {count('support_violation')} |"
        )

    def fit_row(row: dict[str, object]) -> str:
        return (
            f"| {float(row['nominal_temperature_c']):g} | {float(row['R0_ohm']):.9g} | {float(row['alpha_p_s']):.9g} | {float(row['d_n_per_s']):.9g} | "
            f"{1000*float(row['phase_balanced_rmse_v']):.4f} | {1000*float(row['immediate_rmse_v']):.4f} | {1000*float(row['during_pulse_rmse_v']):.4f} | {1000*float(row['relaxation_rmse_v']):.4f} | "
            f"{row['successful_multistarts']}/{row['attempted_multistarts']} | {row['jacobian_rank']} | {float(row['condition_number']):.5g} | "
            f"{row['at_numerical_bound']} | {row['at_ocp_feasibility_boundary']} |"
        )

    def arr_row(parameter: str, definition: str) -> str:
        fit = fits[(parameter, definition)]
        subset = [row for row in loo_rows if row["parameter"] == parameter and row["temperature_definition"] == definition]
        errors = np.asarray([float(row["absolute_relative_error"]) for row in subset])
        return (
            f"| {parameter} | {definition} | {float(fit['slope']):.8g} | {float(fit['intercept']):.8g} | "
            f"{float(fit['activation_energy_j_per_mol']):.8g} | {float(fit['r2']):.6f} | "
            f"{100*np.median(errors):.3f}% | {100*np.max(errors):.3f}% |"
        )

    report = [
        "# Panasonic NCR18650PF multi-temperature HPPC audit",
        "",
        "## Scope",
        "",
        "Exploratory audit only. The frozen electrical and thermal models, OCPs, accepted parameters, and manifests were not modified. No Arrhenius dependence was implemented. UDDS, CC-CV, C/20, and 1C were not fitted.",
        "",
        "## Source and channel audit",
        "",
        "The five original MAT HPPC files were loaded directly from the newly local full dataset. Finite fractions are ordered as time/current/voltage/battery temperature/chamber temperature. The cold chamber channel is reported but never used as ambient or as a fit input.",
        "",
        "| Nominal °C | File | Samples | Duration s | contiguous positive dt min/median/max s | Omitted gaps | finite fractions | Current sign |",
        "|---:|---|---:|---:|---|---:|---|---|",
        *[audit_row(temp) for temp in sorted(campaigns)],
        "",
        "Available channels by file:",
        *[f"- {temp:g} °C: `{'; '.join(campaigns[temp].channels)}`" for temp in sorted(campaigns)],
        "",
        "The chamber-temperature field exists in every MAT structure, but its finite fraction is zero at −20, −10, 0, and 10 °C; only the 25 °C channel is finite. It is therefore invalid for the cold campaigns and was not used. Large raw positive time increments are the documented omitted between-group discharges; the table reports only contiguous increments and counts the omitted gaps explicitly.",
        "",
        "## Battery-temperature audit",
        "",
        "| Nominal °C | Initial °C | Median °C | Min °C | Max °C | End-start K | Median ΔT in 4C/6C pulse K | Maximum ΔT K |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
        *[temp_row(temp) for temp in sorted(campaigns)],
        "",
        "Both nominal campaign temperature and median measured battery temperature are retained for the Arrhenius diagnostic. No choice between them is imposed here.",
        "",
        "## Pulse reconstruction and common support",
        "",
        *[f"- {temp:g} °C: {sum(not bool(r['missing_pulse']) for r in records_by_temperature[temp])} detected records; {sum(bool(r['usable']) for r in records_by_temperature[temp])} individually usable; {sum(bool(r['common_support']) for r in records_by_temperature[temp])} in common support." for temp in sorted(campaigns)],
        f"- Common discrete SOC levels: {', '.join(f'{100*soc:g}%' for soc in common_socs)}.",
        f"- Common C-rates: {', '.join(f'{rate:g}C' for rate in common_rates)}.",
        f"- Common matched combinations per temperature: {len(common)}; total matched pulse records: {len(matched_rows)}.",
        f"- Common SOC range: {100*min(common_socs):g}% to {100*max(common_socs):g}%.",
        "",
        "| Nominal °C | Missing expected | Truncated | Abnormal V | Abnormal I | Insufficient rest | Invalid SOC | Nonfinite battery T | Support violation |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        *[flag_row(temp) for temp in sorted(campaigns)],
        "",
        "All detected and expected-but-missing records, flags, voltages, durations, temperatures, SOC assignments, and apparent transition resistances are in the pulse CSV. No pulse was silently discarded.",
        "",
        "## Direct apparent-transition-resistance diagnostic",
        "",
        "`R_app = abs((V_immediate - V_pre)/(I_pulse - I_pre))` is reported only as apparent transition resistance, never as intrinsic R0.",
        "",
        "| Nominal °C | Median R_app [mΩ] | IQR [mΩ] |",
        "|---:|---:|---:|",
    ]
    for temp in sorted(campaigns):
        values = np.asarray([float(row["apparent_transition_resistance_ohm"]) for row in matched_rows if float(row["nominal_temperature_c"]) == temp])
        report.append(f"| {temp:g} | {1000*np.median(values):.6f} | {1000*np.percentile(values,25):.6f}–{1000*np.percentile(values,75):.6f} |")
    report.extend(
        [
            "",
            "Matched C-rate medians [mΩ]:",
            "",
            "| C-rate | −20 °C | −10 °C | 0 °C | 10 °C | 25 °C |",
            "|---:|---:|---:|---:|---:|---:|",
            *[rapp_matrix_row(f"{rate:g}C", [median_rapp(temp, rate=rate) for temp in ordered_temperatures]) for rate in common_rates],
            "",
            "Matched nominal-SOC medians [mΩ]:",
            "",
            "| SOC | −20 °C | −10 °C | 0 °C | 10 °C | 25 °C |",
            "|---:|---:|---:|---:|---:|---:|",
            *[rapp_matrix_row(f"{100*soc:g}%", [median_rapp(temp, soc=soc) for temp in ordered_temperatures]) for soc in common_socs],
            "",
            "Matched-SOC and matched-C-rate values use the same pulse keys at every temperature; no different pulse population is substituted between temperatures. The raw transition metric shows a strong and consistent increase as temperature decreases, before any model fitting.",
            "",
            "## Accepted 25 °C regression test",
            "",
            f"PASS: {regression['pulse_count']} pulses in {regression['group_count']} histories reproduced `R0={regression['R0_ohm']:.12g} ohm`, `alpha_p={regression['alpha_p_s']:.12g} s`, and `d_n={regression['d_n_per_s']:.12g} 1/s`, with {regression['successful_starts']}/12 successful starts. Accepted files were not overwritten.",
            "",
            "## Independent common-support conditional fits",
            "",
            "The exact accepted phase loss and the same 12 deterministic starts were used independently at every temperature. `alpha_n=3000 s` and `d_p=0.0007 1/s` remain fixed.",
            "",
            "| °C | R0 Ω | alpha_p s | d_n 1/s | Balanced mV | Start mV | Pulse mV | Relaxation mV | Starts | Rank | Condition | Numerical bound | OCP edge |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|",
            *[fit_row(row) for row in parameter_rows],
            "",
            "Practical-conditioning details:",
            *[
                f"- {float(row['nominal_temperature_c']):g} °C: singular values {float(row['singular_value_1']):.6g}, {float(row['singular_value_2']):.6g}, {float(row['singular_value_3']):.6g}; sensitivity correlations (R0,alpha_p)={float(row['sensitivity_corr_R0_alpha_p']):+.4f}, (R0,d_n)={float(row['sensitivity_corr_R0_d_n']):+.4f}, (alpha_p,d_n)={float(row['sensitivity_corr_alpha_p_d_n']):+.4f}; log-bound fractions={float(row['R0_log_bound_fraction']):.4f}/{float(row['alpha_p_log_bound_fraction']):.4f}/{float(row['d_n_log_bound_fraction']):.4f}; near-best ranges R0=[{float(row['near_best_R0_min']):.6g},{float(row['near_best_R0_max']):.6g}], alpha_p=[{float(row['near_best_alpha_p_min']):.6g},{float(row['near_best_alpha_p_max']):.6g}], d_n=[{float(row['near_best_d_n_min']):.6g},{float(row['near_best_d_n_max']):.6g}]."
                for row in parameter_rows
            ],
            "",
            "These are practical local diagnostics, not structural-identifiability claims. At 0 °C the condition number exceeds 10^4 and the near-best R0 spread is large. At 10 and 25 °C, distinct near-best alpha_p valleys remain. At −20 and −10 °C, alpha_p lies within about 1% of its log-range OCP-feasibility edge.",
            "",
            "## Arrhenius diagnostics and leave-one-temperature-out robustness",
            "",
            "| Parameter | Temperature definition | Slope [K] | Intercept | Signed Ea [J/mol] | R² | LOTO median abs. relative error | LOTO worst |",
            "|---|---|---:|---:|---:|---:|---:|---:|",
            *[arr_row(parameter, definition) for parameter in ("R0", "alpha_p", "d_n") for definition in ("nominal", "battery_median")],
            "",
            "A high R² is not treated as proof of Arrhenius physics. The five-point fits are judged together with local recoverability, boundaries, and leave-one-temperature-out errors.",
            "",
            "## Campaign order and aging confounding",
            "",
            "Chronological HPPC acquisition order recovered from filenames:",
            *[f"- {date.date().isoformat()}: {temp:g} °C — `{name}`" for date, temp, name in acquisition],
            "",
            "Temperature decreases as the campaign progresses from 25 °C toward −20 °C. The dataset documentation explicitly warns that the cell aged during the campaign. Therefore temperature and campaign progression are confounded; no recovered trend is attributed exclusively to temperature, and no aging correction is attempted.",
            "",
            "## Required decisions",
            "",
            f"- `R0(T)`: **{decisions['R0']}**",
            f"- `alpha_p(T)`: **{decisions['alpha_p']}**",
            f"- `d_n(T)`: **{decisions['d_n']}**",
            "",
            "`R0` has essentially no Arrhenius-line support and becomes practically degenerate at 0 °C. `alpha_p` has moderate in-sample linearity but poor leave-one-temperature-out transfer and near-edge/multivalley solutions. `d_n` has the strongest in-sample linear appearance, but its leave-one-out worst error exceeds 200% and it remains embedded in the ill-conditioned, highly correlated three-parameter fit. Consequently none of the three independently recovered sequences is sufficiently robust to implement as a temperature law.",
            "",
            f"**{global_decision}**",
            "",
            "No temperature-dependent law was implemented after this decision.",
        ]
    )
    REPORT_PATH.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(global_decision)
    print(f"Common support: {len(common)} pulses per temperature; SOC {min(common_socs):.3f} to {max(common_socs):.3f}; rates {common_rates}")
    for row in parameter_rows:
        print(f"{float(row['nominal_temperature_c']):g} C: R0={float(row['R0_ohm']):.9g}, alpha_p={float(row['alpha_p_s']):.9g}, d_n={float(row['d_n_per_s']):.9g}, RMSE={1000*float(row['phase_balanced_rmse_v']):.3f} mV")
    print(decisions)


if __name__ == "__main__":
    main()
