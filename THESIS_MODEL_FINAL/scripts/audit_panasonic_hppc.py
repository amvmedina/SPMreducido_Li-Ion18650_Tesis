"""Read-only audit and pulse preparation for the Panasonic 25 degC HPPC test.

This script performs segmentation and diagnostic calculations only.  It does
not optimize dynamic parameters or modify any model/OCP parameterization.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from models.ocp_panasonic_qi import (  # noqa: E402
    MCMB_SUPPORTED_STOICHIOMETRY_RANGE,
    NCA_WORKING_STOICHIOMETRY_RANGE,
)
from models.ocp_panasonic_qi_frozen import (  # noqa: E402
    Q_CELL_AH,
    X_N_END,
    X_N_START,
    X_P_END,
    X_P_START,
)


RAW_MANIFEST_PATH = Path(
    r"C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1"
    r"\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC"
    r"\5 pulse disch\03-11-17_08.47 25degC_5Pulse_HPPC_Pan18650PF.mat"
)
LOADED_CSV_PATH = Path(
    r"C:\Users\amvel\Documents\Andes\Tesis\datasets\converted_csv"
    r"\hppc_25C_converted.csv"
)
OUTPUT_DIRECTORY = REPOSITORY_ROOT / "thesis_data"

REST_CURRENT_THRESHOLD_A = 0.05
PULSE_CURRENT_THRESHOLD_A = 0.5
CONTIGUOUS_GAP_THRESHOLD_S = 10.0
NOMINAL_PULSE_DURATION_S = 10.0
MIN_USABLE_PULSE_DURATION_S = 9.5
LOW_VOLTAGE_CUTOFF_V = 2.5
NOMINAL_ONE_C_CURRENT_A = 2.9


@dataclass(frozen=True)
class HppcData:
    columns: tuple[str, ...]
    time_s: np.ndarray
    voltage_v: np.ndarray
    current_a: np.ndarray
    exported_ah: np.ndarray
    battery_temperature_c: np.ndarray
    chamber_temperature_c: np.ndarray


@dataclass(frozen=True)
class Segment:
    state: str
    start_index: int
    end_index: int


def load_hppc(path: Path = LOADED_CSV_PATH) -> HppcData:
    table = np.genfromtxt(path, delimiter=",", names=True, dtype=float, encoding="utf-8")
    columns = tuple(table.dtype.names or ())
    required = {
        "Time",
        "Voltage",
        "Current",
        "Ah",
        "Battery_Temp_degC",
        "Chamber_Temp_degC",
    }
    missing = sorted(required - set(columns))
    if missing:
        raise ValueError(f"HPPC CSV is missing required columns: {missing}")
    arrays = {
        "time_s": np.asarray(table["Time"], dtype=float),
        "voltage_v": np.asarray(table["Voltage"], dtype=float),
        "current_a": np.asarray(table["Current"], dtype=float),
        "exported_ah": np.asarray(table["Ah"], dtype=float),
        "battery_temperature_c": np.asarray(table["Battery_Temp_degC"], dtype=float),
        "chamber_temperature_c": np.asarray(table["Chamber_Temp_degC"], dtype=float),
    }
    if not all(np.all(np.isfinite(values)) for values in arrays.values()):
        raise ValueError("HPPC CSV contains non-finite values.")
    if np.any(np.diff(arrays["time_s"]) < 0.0):
        raise ValueError("HPPC elapsed time decreases.")
    return HppcData(columns=columns, **arrays)


def current_state(current_a: np.ndarray) -> np.ndarray:
    states = np.full(current_a.shape, "transition", dtype="U10")
    states[np.abs(current_a) <= REST_CURRENT_THRESHOLD_A] = "rest"
    states[current_a <= -PULSE_CURRENT_THRESHOLD_A] = "discharge"
    states[current_a >= PULSE_CURRENT_THRESHOLD_A] = "charge"
    return states


def segment_states(data: HppcData) -> list[Segment]:
    states = current_state(data.current_a)
    segments: list[Segment] = []
    start = 0
    for index in range(1, states.size):
        time_gap = data.time_s[index] - data.time_s[index - 1]
        if states[index] != states[index - 1] or time_gap > CONTIGUOUS_GAP_THRESHOLD_S:
            segments.append(Segment(states[start], start, index - 1))
            start = index
    segments.append(Segment(states[start], start, states.size - 1))
    return segments


def _rest_bounds_before(data: HppcData, pulse_start: int) -> tuple[int, int] | None:
    end = pulse_start - 1
    if end < 0 or abs(data.current_a[end]) > REST_CURRENT_THRESHOLD_A:
        return None
    start = end
    while start > 0:
        gap = data.time_s[start] - data.time_s[start - 1]
        if gap > CONTIGUOUS_GAP_THRESHOLD_S or abs(data.current_a[start - 1]) > REST_CURRENT_THRESHOLD_A:
            break
        start -= 1
    return start, end


def _rest_bounds_after(data: HppcData, pulse_end: int) -> tuple[int, int] | None:
    start = pulse_end + 1
    if start >= data.time_s.size or abs(data.current_a[start]) > REST_CURRENT_THRESHOLD_A:
        return None
    end = start
    while end + 1 < data.time_s.size:
        gap = data.time_s[end + 1] - data.time_s[end]
        if gap > CONTIGUOUS_GAP_THRESHOLD_S or abs(data.current_a[end + 1]) > REST_CURRENT_THRESHOLD_A:
            break
        end += 1
    return start, end


def _nominal_rate_label(current_magnitude_a: float) -> str:
    c_rate = current_magnitude_a / NOMINAL_ONE_C_CURRENT_A
    candidates = np.asarray([0.5, 1.0, 2.0, 4.0, 6.0])
    nearest = float(candidates[np.argmin(np.abs(candidates - c_rate))])
    return f"{nearest:g}C"


def identify_pulses(data: HppcData, segments: list[Segment]) -> list[dict[str, object]]:
    pulse_segments = [segment for segment in segments if segment.state in {"discharge", "charge"}]
    records: list[dict[str, object]] = []
    for pulse_id, segment in enumerate(pulse_segments, start=1):
        start = segment.start_index
        end = segment.end_index
        pre_bounds = _rest_bounds_before(data, start)
        post_bounds = _rest_bounds_after(data, end)
        if pre_bounds is None or post_bounds is None:
            raise ValueError(f"Pulse {pulse_id} lacks a measurable pre/post rest segment.")
        pre_start, pre_end = pre_bounds
        post_start, post_end = post_bounds

        pulse_current = float(np.median(data.current_a[start : end + 1]))
        current_magnitude = abs(pulse_current)
        active_differences = np.diff(data.time_s[start : end + 1])
        positive_active_differences = active_differences[active_differences > 0.0]
        final_active_interval_s = (
            float(np.median(positive_active_differences))
            if positive_active_differences.size
            else 0.0
        )
        pulse_duration_s = float(
            data.time_s[end] - data.time_s[start] + final_active_interval_s
        )
        pre_rest_duration_s = float(data.time_s[start] - data.time_s[pre_start])
        post_rest_duration_s = float(data.time_s[post_end] - data.time_s[post_start])
        pre_current = float(data.current_a[pre_end])
        delta_current = pulse_current - pre_current
        v_pre = float(data.voltage_v[pre_end])
        v_immediate = float(data.voltage_v[start])
        v_end = float(data.voltage_v[end])
        v_post_rest = float(data.voltage_v[post_end])
        r_app = abs((v_immediate - v_pre) / delta_current)

        # The source Ah counter is negative during discharge and retains the
        # cumulative removed capacity across the omitted inter-group records.
        discharged_capacity_ah = max(0.0, -float(data.exported_ah[pre_end]))
        z = discharged_capacity_ah / Q_CELL_AH
        soc = 1.0 - z
        x_n = X_N_START + z * (X_N_END - X_N_START)
        x_p = X_P_START + z * (X_P_END - X_P_START)
        domain_ok = bool(
            MCMB_SUPPORTED_STOICHIOMETRY_RANGE[0]
            <= x_n
            <= MCMB_SUPPORTED_STOICHIOMETRY_RANGE[1]
            and NCA_WORKING_STOICHIOMETRY_RANGE[0]
            <= x_p
            <= NCA_WORKING_STOICHIOMETRY_RANGE[1]
        )
        truncated_by_cutoff = bool(
            pulse_duration_s < MIN_USABLE_PULSE_DURATION_S
            and v_end <= LOW_VOLTAGE_CUTOFF_V + 0.001
        )
        usable = bool(
            segment.state == "discharge"
            and pulse_duration_s >= MIN_USABLE_PULSE_DURATION_S
            and domain_ok
            and np.isfinite(r_app)
        )
        records.append(
            {
                "pulse_id": pulse_id,
                "pulse_type": segment.state,
                "nominal_rate": _nominal_rate_label(current_magnitude),
                "start_index": start,
                "end_index": end,
                "pre_pulse_index": pre_end,
                "post_rest_index": post_end,
                "start_time_s": float(data.time_s[start]),
                "end_time_s": float(data.time_s[end]),
                "discharged_capacity_ah": discharged_capacity_ah,
                "soc": soc,
                "x_n": x_n,
                "x_p": x_p,
                "current_a": pulse_current,
                "current_magnitude_a": current_magnitude,
                "pulse_duration_s": pulse_duration_s,
                "pre_rest_duration_s": pre_rest_duration_s,
                "post_rest_duration_s": post_rest_duration_s,
                "v_pre_v": v_pre,
                "v_immediate_v": v_immediate,
                "v_end_v": v_end,
                "v_post_rest_v": v_post_rest,
                "delta_v_immediate_v": v_immediate - v_pre,
                "delta_i_a": delta_current,
                "r_app_ohm": r_app,
                "supported_stoichiometric_domains": domain_ok,
                "truncated_by_low_voltage_cutoff": truncated_by_cutoff,
                "usable_for_dynamic_calibration": usable,
            }
        )
    return records


def _csv_value(value: object) -> object:
    if isinstance(value, (float, np.floating)):
        return f"{float(value):.12g}"
    if isinstance(value, (int, np.integer)) and not isinstance(value, bool):
        return int(value)
    return value


def write_summary(path: Path, records: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        for record in records:
            writer.writerow({name: _csv_value(value) for name, value in record.items()})


def save_current_voltage_plot(data: HppcData, records: list[dict[str, object]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    time_h = data.time_s / 3600.0
    figure, axes = plt.subplots(3, 1, figsize=(12.0, 9.0), constrained_layout=True)
    axes[0].plot(time_h, data.current_a, linewidth=0.65)
    axes[0].set(ylabel="Raw current [A]", title="Panasonic 25 °C HPPC current")
    axes[1].plot(time_h, data.voltage_v, linewidth=0.65, color="tab:orange")
    axes[1].set(ylabel="Voltage [V]", title="Panasonic 25 °C HPPC voltage")

    starts_h = np.asarray([float(record["start_time_s"]) / 3600.0 for record in records])
    soc_percent = 100.0 * np.asarray([float(record["soc"]) for record in records])
    magnitudes = np.asarray([float(record["current_magnitude_a"]) for record in records])
    scatter = axes[2].scatter(starts_h, soc_percent, c=magnitudes, cmap="viridis", s=34)
    axes[2].set(
        xlabel="Elapsed time [h]",
        ylabel="Pulse-start SOC [%]",
        title="Detected pulse locations and SOC",
    )
    colorbar = figure.colorbar(scatter, ax=axes[2], pad=0.01)
    colorbar.set_label("Current magnitude [A]")
    for axis in axes:
        axis.grid(alpha=0.25)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def save_rapp_plot(records: list[dict[str, object]], path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(10.0, 6.0), constrained_layout=True)
    labels = ["0.5C", "1C", "2C", "4C", "6C"]
    for label in labels:
        selected = [record for record in records if record["nominal_rate"] == label]
        axis.scatter(
            [100.0 * float(record["soc"]) for record in selected],
            [1000.0 * float(record["r_app_ohm"]) for record in selected],
            s=34,
            label=label,
        )
    truncated = [record for record in records if record["truncated_by_low_voltage_cutoff"]]
    axis.scatter(
        [100.0 * float(record["soc"]) for record in truncated],
        [1000.0 * float(record["r_app_ohm"]) for record in truncated],
        marker="x",
        s=80,
        linewidth=1.8,
        color="black",
        label="Cutoff-truncated",
    )
    axis.set(
        xlabel="Pulse-start SOC [%]",
        ylabel="$R_{app}$ [mΩ]",
        title="Apparent instantaneous resistance diagnostic",
    )
    axis.grid(alpha=0.25)
    axis.legend(ncol=3)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def save_representative_pulses(
    data: HppcData, records: list[dict[str, object]], path: Path
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    one_c = [
        record
        for record in records
        if record["nominal_rate"] == "1C" and record["usable_for_dynamic_calibration"]
    ]
    targets = np.asarray([1.0, 0.8, 0.6, 0.4, 0.2, 0.08])
    selected: list[dict[str, object]] = []
    for target in targets:
        closest = min(one_c, key=lambda record: abs(float(record["soc"]) - target))
        if closest not in selected:
            selected.append(closest)

    figure, axes = plt.subplots(2, 3, figsize=(13.0, 7.5), constrained_layout=True)
    for axis, record in zip(axes.flat, selected, strict=False):
        start = int(record["start_index"])
        start_time = float(record["start_time_s"])
        lower = int(np.searchsorted(data.time_s, start_time - 5.0, side="left"))
        upper = int(np.searchsorted(data.time_s, start_time + 70.0, side="right"))
        relative_time = data.time_s[lower:upper] - start_time
        axis.plot(relative_time, data.voltage_v[lower:upper], linewidth=1.3)
        axis.axvspan(0.0, float(record["pulse_duration_s"]), color="tab:red", alpha=0.10)
        axis.set(
            title=f"Pulse {record['pulse_id']}: {100 * float(record['soc']):.1f}% SOC",
            xlabel="Time from pulse start [s]",
            ylabel="Voltage [V]",
        )
        axis.grid(alpha=0.25)
    for axis in axes.flat[len(selected) :]:
        axis.set_visible(False)
    figure.suptitle("Representative 1C HPPC voltage responses")
    figure.savefig(path, dpi=180)
    plt.close(figure)


def _duration_statistics(data: HppcData, segments: list[Segment], state: str) -> np.ndarray:
    durations = []
    for segment in segments:
        if segment.state != state:
            continue
        durations.append(data.time_s[segment.end_index] - data.time_s[segment.start_index])
    return np.asarray(durations, dtype=float)


def write_audit(
    path: Path,
    data: HppcData,
    segments: list[Segment],
    records: list[dict[str, object]],
) -> None:
    time_differences = np.diff(data.time_s)
    positive_differences = time_differences[time_differences > 0.0]
    active_pairs = (
        (np.abs(data.current_a[:-1]) > PULSE_CURRENT_THRESHOLD_A)
        & (np.abs(data.current_a[1:]) > PULSE_CURRENT_THRESHOLD_A)
        & (time_differences > 0.0)
        & (time_differences <= CONTIGUOUS_GAP_THRESHOLD_S)
    )
    rest_pairs = (
        (np.abs(data.current_a[:-1]) <= REST_CURRENT_THRESHOLD_A)
        & (np.abs(data.current_a[1:]) <= REST_CURRENT_THRESHOLD_A)
        & (time_differences > 0.0)
        & (time_differences <= CONTIGUOUS_GAP_THRESHOLD_S)
    )
    long_gaps = time_differences[time_differences > CONTIGUOUS_GAP_THRESHOLD_S]
    rest_durations = _duration_statistics(data, segments, "rest")
    discharge_records = [record for record in records if record["pulse_type"] == "discharge"]
    charge_records = [record for record in records if record["pulse_type"] == "charge"]
    usable = [record for record in records if record["usable_for_dynamic_calibration"]]
    truncated = [record for record in records if record["truncated_by_low_voltage_cutoff"]]
    domain_failures = [record for record in records if not record["supported_stoichiometric_domains"]]
    usable_soc = np.asarray([float(record["soc"]) for record in usable])
    usable_rapp = np.asarray([float(record["r_app_ohm"]) for record in usable])
    pulse_durations = np.asarray([float(record["pulse_duration_s"]) for record in records])
    post_rest_durations = np.asarray([float(record["post_rest_duration_s"]) for record in records])

    rate_counts = {
        label: sum(record["nominal_rate"] == label for record in discharge_records)
        for label in ("0.5C", "1C", "2C", "4C", "6C")
    }
    lines = [
        "# Panasonic 25 °C HPPC data audit and preparation",
        "",
        "## Scope",
        "",
        "Read-only segmentation and diagnostics only. `R_app` is an apparent "
        "instantaneous step-response quantity and is **not** automatically the "
        "model parameter `R0`. No dynamic parameters were optimized.",
        "",
        "## Sources and signals",
        "",
        f"- Raw MAT path recorded in the manifest: `{RAW_MANIFEST_PATH}`",
        f"- Read-only converted CSV loaded: `{LOADED_CSV_PATH}`",
        f"- Columns: `{', '.join(data.columns)}`",
        f"- Rows: {data.time_s.size}",
        f"- Elapsed-time range: {data.time_s[0]:.6f} to {data.time_s[-1]:.6f} s",
        f"- Total logged duration: {data.time_s[-1] - data.time_s[0]:.6f} s "
        f"({(data.time_s[-1] - data.time_s[0]) / 3600.0:.6f} h)",
        f"- Voltage range: {np.min(data.voltage_v):.5f} to {np.max(data.voltage_v):.5f} V",
        f"- Current range: {np.min(data.current_a):.5f} to {np.max(data.current_a):.5f} A",
        "- Current convention: raw negative current is discharge; raw positive "
        "current would be charge; zero is rest",
        f"- Battery case temperature: {np.min(data.battery_temperature_c):.5f} to "
        f"{np.max(data.battery_temperature_c):.5f} °C",
        f"- Chamber temperature: {np.min(data.chamber_temperature_c):.2f} to "
        f"{np.max(data.chamber_temperature_c):.2f} °C",
        "",
        "## Sampling and gaps",
        "",
        f"- Median positive sample interval overall: {np.median(positive_differences):.6f} s",
        f"- Median within-pulse sample interval: {np.median(time_differences[active_pairs]):.6f} s",
        f"- Median contiguous-rest sample interval: {np.median(time_differences[rest_pairs]):.6f} s",
        f"- Duplicate elapsed-time increments: {int(np.sum(time_differences == 0.0))}",
        f"- Long inter-group gaps greater than {CONTIGUOUS_GAP_THRESHOLD_S:g} s: "
        f"{long_gaps.size}, spanning {np.min(long_gaps):.3f} to {np.max(long_gaps):.3f} s",
        "- The long gaps separate pulse groups. The cumulative `Ah` counter "
        "jumps consistently across them and is therefore used for SOC; direct "
        "integration of this pulse-only CSV would miss the omitted between-group discharge.",
        "",
        "## Reproducible segmentation",
        "",
        f"- Rest: `abs(I) <= {REST_CURRENT_THRESHOLD_A:.2f} A`",
        f"- Discharge pulse: `I <= -{PULSE_CURRENT_THRESHOLD_A:.2f} A`",
        f"- Charge pulse: `I >= +{PULSE_CURRENT_THRESHOLD_A:.2f} A`",
        f"- A time jump greater than {CONTIGUOUS_GAP_THRESHOLD_S:g} s starts a new segment",
        f"- Rest segments: {sum(segment.state == 'rest' for segment in segments)}; "
        f"duration range {np.min(rest_durations):.3f} to {np.max(rest_durations):.3f} s",
        f"- Discharge pulse segments: {len(discharge_records)}",
        f"- Charge/recovery pulse segments: {len(charge_records)} (recovery occurs by zero-current rest)",
        f"- Pulse counts by nominal magnitude: "
        + ", ".join(f"{label}: {count}" for label, count in rate_counts.items()),
        f"- Pulse-duration range: {np.min(pulse_durations):.3f} to "
        f"{np.max(pulse_durations):.3f} s",
        f"- Post-pulse rest-duration range: {np.min(post_rest_durations):.3f} to "
        f"{np.max(post_rest_durations):.3f} s",
        "",
        "Pulse indices in `panasonic_hppc_pulse_summary.csv` are zero-based CSV "
        "data-row indices (the header is not counted). `V_pre` is the last rest "
        "sample before the pulse; `V_immediate` is the first loaded sample; "
        "`V_end` is the last loaded sample; `V_post_rest` is the final contiguous "
        "rest sample before the next pulse or long inter-group gap. `R_app` uses "
        "the median pulse current relative to the pre-pulse rest current.",
        "",
        "## SOC and frozen nominal electrode mapping",
        "",
        f"- Fixed Panasonic capacity: `{Q_CELL_AH:.12f} Ah`",
        "- `Q_discharge = -Ah` and `SOC = 1 - Q_discharge/Q_cell`",
        f"- Frozen nominal effective negative window: {X_N_START:.12f} to {X_N_END:.12f}",
        f"- Frozen nominal effective positive window: {X_P_START:.12f} to {X_P_END:.12f}",
        f"- Pulses outside frozen OCP supported domains: {len(domain_failures)}",
        "- These mappings use nominal effective windows, not uniquely identified "
        "physical electrode stoichiometries.",
        "",
        "## Usability assessment",
        "",
        f"- Detected discharge pulses: {len(discharge_records)}",
        f"- Usable complete pulses: **{len(usable)}**",
        f"- Low-voltage-cutoff-truncated pulses: {len(truncated)} "
        f"(IDs {', '.join(str(record['pulse_id']) for record in truncated)})",
        f"- Usable pulse-start SOC range: {100 * np.min(usable_soc):.3f}% to "
        f"{100 * np.max(usable_soc):.3f}%",
        f"- Usable `R_app` range: {1000 * np.min(usable_rapp):.3f} to "
        f"{1000 * np.max(usable_rapp):.3f} mΩ",
        "",
        "**Calibration suitability:** The 64 complete pulses provide repeated "
        "0.5C--6C excitation over a broad SOC range and appear suitable for a "
        "subsequent dynamic-parameter calibration after the three cutoff-truncated "
        "pulses are excluded or treated separately. The dataset also contains "
        "long omitted-data gaps, variable 0.1/1 s sampling, duplicate timestamps, "
        "and progressively missing high-current levels near low SOC; loaders must "
        "segment groups and must not integrate current across missing intervals.",
        "",
        "No HPPC parameter fitting, UDDS, 1C validation, or model-equation change "
        "was performed in this audit.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    data = load_hppc()
    segments = segment_states(data)
    records = identify_pulses(data, segments)
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    write_summary(OUTPUT_DIRECTORY / "panasonic_hppc_pulse_summary.csv", records)
    save_current_voltage_plot(
        data, records, OUTPUT_DIRECTORY / "panasonic_hppc_current_voltage.png"
    )
    save_rapp_plot(records, OUTPUT_DIRECTORY / "panasonic_hppc_rapp_vs_soc.png")
    save_representative_pulses(
        data, records, OUTPUT_DIRECTORY / "panasonic_hppc_representative_pulses.png"
    )
    write_audit(
        OUTPUT_DIRECTORY / "panasonic_hppc_audit.md", data, segments, records
    )
    usable = [record for record in records if record["usable_for_dynamic_calibration"]]
    print(f"Detected pulses: {len(records)}; usable complete pulses: {len(usable)}")
    print((OUTPUT_DIRECTORY / "panasonic_hppc_audit.md").resolve())


if __name__ == "__main__":
    main()
