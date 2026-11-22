"""Compare Panasonic C/20 charge and discharge as a data-only diagnostic.

The primary comparison uses the same frozen cell capacity for both directions
and is restricted to their common coulomb-counted SOC support.  No model or OCP
code is imported, no parameters are fitted, and no extrapolation is permitted.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.io import loadmat


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SOURCE_PATH = Path(
    r"C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1"
    r"\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC"
    r"\C20 OCV and 1C discharge tests_start_of_tests"
    r"\05-08-17_13.26 C20 OCV Test_C20_25dC.mat"
)
OUTPUT_DIRECTORY = REPOSITORY_ROOT / "thesis_data"

Q_CELL_AH = 2.994979138519
ACTIVE_CURRENT_THRESHOLD_A = 0.05
COMMON_GRID_POINTS = 2001
CCCV_CC_BIAS_MV = -50.2
CCCV_CV_BIAS_MV = -68.0


@dataclass(frozen=True)
class Segment:
    kind: str
    start: int
    end: int

    @property
    def indices(self) -> np.ndarray:
        return np.arange(self.start, self.end + 1)


def _vector(value: object, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float).reshape(-1)
    if array.size == 0 or not np.all(np.isfinite(array)):
        raise ValueError(f"MAT field {name!r} is empty or non-finite.")
    return array


def load_source(path: Path) -> tuple[tuple[str, ...], dict[str, np.ndarray]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    measurement = loadmat(path, squeeze_me=True, struct_as_record=False).get("meas")
    if measurement is None or not hasattr(measurement, "_fieldnames"):
        raise ValueError("Expected a MATLAB structure named 'meas'.")
    columns = tuple(measurement._fieldnames)
    required = (
        "Time",
        "Voltage",
        "Current",
        "Ah",
        "Battery_Temp_degC",
        "Chamber_Temp_degC",
    )
    missing = sorted(set(required) - set(columns))
    if missing:
        raise ValueError(f"Missing required fields: {missing}")
    arrays = {name: _vector(getattr(measurement, name), name) for name in required}
    if len({values.size for values in arrays.values()}) != 1:
        raise ValueError("MAT fields have inconsistent lengths.")
    if np.any(np.diff(arrays["Time"]) < 0.0):
        raise ValueError("Elapsed time decreases in the source data.")
    return columns, arrays


def identify_segments(current_a: np.ndarray) -> list[Segment]:
    """Run-length segment the measured current; signs are data classifications."""

    state = np.zeros(current_a.size, dtype=int)
    state[current_a < -ACTIVE_CURRENT_THRESHOLD_A] = -1
    state[current_a > ACTIVE_CURRENT_THRESHOLD_A] = 1
    edges = np.r_[0, np.flatnonzero(np.diff(state) != 0) + 1, state.size]
    names = {-1: "C/20 discharge", 0: "rest", 1: "C/20 charge"}
    segments = [
        Segment(names[int(state[start])], int(start), int(end - 1))
        for start, end in zip(edges[:-1], edges[1:])
    ]
    expected = ["rest", "C/20 discharge", "rest", "C/20 charge", "rest"]
    if [segment.kind for segment in segments] != expected:
        raise ValueError(
            "Signal-based segmentation did not yield the expected five-part structure: "
            + ", ".join(segment.kind for segment in segments)
        )
    return segments


def cumulative_capacity_ah(time_s: np.ndarray, current_a: np.ndarray) -> np.ndarray:
    if np.any(np.diff(time_s) <= 0.0):
        raise ValueError("Active-branch times must be strictly increasing.")
    return np.r_[0.0, cumulative_trapezoid(current_a, time_s) / 3600.0]


def gap_statistics(values_v: np.ndarray) -> dict[str, float]:
    return {
        "mean_v": float(np.mean(values_v)),
        "median_v": float(np.median(values_v)),
        "rms_v": float(np.sqrt(np.mean(values_v**2))),
        "min_v": float(np.min(values_v)),
        "max_v": float(np.max(values_v)),
        "p5_v": float(np.percentile(values_v, 5.0)),
        "p95_v": float(np.percentile(values_v, 95.0)),
    }


def region_name(soc: float) -> str:
    if soc < 1.0 / 3.0:
        return "low [0, 1/3)"
    if soc < 2.0 / 3.0:
        return "middle [1/3, 2/3)"
    return "high [2/3, overlap max]"


def fmt_mv(value_v: float) -> str:
    return f"{1000.0 * value_v:.3f}"


def write_csv(
    path: Path,
    soc: np.ndarray,
    discharge_v: np.ndarray,
    charge_v: np.ndarray,
    delta_v: np.ndarray,
    shape_coordinate: np.ndarray,
    shape_discharge_v: np.ndarray,
    shape_charge_v: np.ndarray,
    shape_delta_v: np.ndarray,
) -> None:
    fields = (
        "soc_common_capacity",
        "voltage_c20_discharge_v",
        "voltage_c20_charge_v",
        "delta_v_qocv_v",
        "soc_region",
        "shape_only_normalized_coordinate",
        "shape_only_discharge_v",
        "shape_only_charge_v",
        "shape_only_delta_v_qocv_v",
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for i in range(soc.size):
            writer.writerow(
                {
                    "soc_common_capacity": f"{soc[i]:.12g}",
                    "voltage_c20_discharge_v": f"{discharge_v[i]:.12g}",
                    "voltage_c20_charge_v": f"{charge_v[i]:.12g}",
                    "delta_v_qocv_v": f"{delta_v[i]:.12g}",
                    "soc_region": region_name(float(soc[i])),
                    "shape_only_normalized_coordinate": f"{shape_coordinate[i]:.12g}",
                    "shape_only_discharge_v": f"{shape_discharge_v[i]:.12g}",
                    "shape_only_charge_v": f"{shape_charge_v[i]:.12g}",
                    "shape_only_delta_v_qocv_v": f"{shape_delta_v[i]:.12g}",
                }
            )


def make_voltage_figure(
    path: Path,
    soc: np.ndarray,
    discharge_v: np.ndarray,
    charge_v: np.ndarray,
    shape_coordinate: np.ndarray,
    shape_discharge_v: np.ndarray,
    shape_charge_v: np.ndarray,
) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(9.0, 8.2), constrained_layout=True)
    axes[0].plot(soc, discharge_v, lw=2.0, label="C/20 discharge")
    axes[0].plot(soc, charge_v, lw=2.0, label="C/20 charge")
    axes[0].set(
        xlabel="SOC from common frozen capacity",
        ylabel="Measured voltage (V)",
        title="Primary common-capacity comparison (overlap only)",
        xlim=(0.0, 1.0),
    )
    axes[0].axvspan(float(soc[-1]), 1.0, color="0.9", label="No charge coverage")
    axes[0].grid(alpha=0.25)
    axes[0].legend(loc="best")

    axes[1].plot(shape_coordinate, shape_discharge_v, lw=2.0, label="C/20 discharge")
    axes[1].plot(shape_coordinate, shape_charge_v, lw=2.0, label="C/20 charge")
    axes[1].set(
        xlabel="Branch-normalized capacity coordinate",
        ylabel="Measured voltage (V)",
        title="SHAPE-ONLY NORMALIZED COMPARISON (secondary diagnostic)",
        xlim=(0.0, 1.0),
    )
    axes[1].grid(alpha=0.25)
    axes[1].legend(loc="best")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def make_gap_figure(
    path: Path,
    soc: np.ndarray,
    delta_v: np.ndarray,
    regional: dict[str, dict[str, float]],
) -> None:
    fig, ax = plt.subplots(figsize=(9.0, 5.2), constrained_layout=True)
    ax.plot(soc, 1000.0 * delta_v, color="tab:purple", lw=2.0)
    ax.axvline(1.0 / 3.0, color="0.5", ls="--", lw=1.0)
    ax.axvline(2.0 / 3.0, color="0.5", ls="--", lw=1.0)
    colors = {"low": "tab:blue", "middle": "tab:orange", "high": "tab:green"}
    spans = {
        "low": (0.0, min(1.0 / 3.0, float(soc[-1]))),
        "middle": (1.0 / 3.0, min(2.0 / 3.0, float(soc[-1]))),
        "high": (2.0 / 3.0, float(soc[-1])),
    }
    for name, (left, right) in spans.items():
        if right > left:
            ax.hlines(
                1000.0 * regional[name]["mean_v"],
                left,
                right,
                color=colors[name],
                lw=2.5,
                label=f"{name.capitalize()}-SOC mean",
            )
    ax.set(
        xlabel="SOC from common frozen capacity",
        ylabel=r"$V_{charge}-V_{discharge}$ (mV)",
        title="Low-rate charge/discharge quasi-OCV gap",
        xlim=(0.0, float(soc[-1])),
    )
    ax.grid(alpha=0.25)
    ax.legend(loc="best", ncol=2)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    columns, data = load_source(SOURCE_PATH)
    segments = identify_segments(data["Current"])
    initial_rest, discharge, intermediate_rest, charge, final_rest = segments
    d = discharge.indices
    c = charge.indices

    q_discharge_ah = cumulative_capacity_ah(data["Time"][d], -data["Current"][d])
    q_charge_ah = cumulative_capacity_ah(data["Time"][c], data["Current"][c])
    soc_discharge = 1.0 - q_discharge_ah / Q_CELL_AH
    soc_charge = q_charge_ah / Q_CELL_AH
    if np.any(np.diff(soc_discharge) >= 0.0) or np.any(np.diff(soc_charge) <= 0.0):
        raise ValueError("Coulomb-counted branch SOC coordinates are not strictly monotonic.")

    common_soc_min = max(float(np.min(soc_discharge)), float(np.min(soc_charge)))
    common_soc_max = min(float(np.max(soc_discharge)), float(np.max(soc_charge)))
    if common_soc_min >= common_soc_max:
        raise ValueError("The branches have no common SOC support.")
    common_soc = np.linspace(common_soc_min, common_soc_max, COMMON_GRID_POINTS)
    # np.interp requires increasing abscissas. Reversing the discharge does not
    # normalize or otherwise alter its common-capacity SOC definition.
    discharge_v = np.interp(common_soc, soc_discharge[::-1], data["Voltage"][d][::-1])
    charge_v = np.interp(common_soc, soc_charge, data["Voltage"][c])
    if common_soc[0] < soc_discharge.min() or common_soc[-1] > soc_discharge.max():
        raise ValueError("Discharge interpolation would extrapolate.")
    if common_soc[0] < soc_charge.min() or common_soc[-1] > soc_charge.max():
        raise ValueError("Charge interpolation would extrapolate.")
    delta_v = charge_v - discharge_v

    masks = {
        "low": common_soc < 1.0 / 3.0,
        "middle": (common_soc >= 1.0 / 3.0) & (common_soc < 2.0 / 3.0),
        "high": common_soc >= 2.0 / 3.0,
    }
    overall = gap_statistics(delta_v)
    regional = {
        name: {**gap_statistics(delta_v[mask]), "n_grid": int(np.sum(mask))}
        for name, mask in masks.items()
        if np.any(mask)
    }

    # Explicitly secondary: independently normalized transferred-capacity axes.
    shape_coordinate = np.linspace(0.0, 1.0, COMMON_GRID_POINTS)
    discharge_shape_soc = 1.0 - q_discharge_ah / q_discharge_ah[-1]
    charge_shape_soc = q_charge_ah / q_charge_ah[-1]
    shape_discharge_v = np.interp(
        shape_coordinate, discharge_shape_soc[::-1], data["Voltage"][d][::-1]
    )
    shape_charge_v = np.interp(shape_coordinate, charge_shape_soc, data["Voltage"][c])
    shape_delta_v = shape_charge_v - shape_discharge_v
    shape_stats = gap_statistics(shape_delta_v)

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    csv_path = OUTPUT_DIRECTORY / "panasonic_c20_charge_discharge_comparison.csv"
    voltage_plot_path = OUTPUT_DIRECTORY / "panasonic_c20_charge_discharge_voltage.png"
    gap_plot_path = OUTPUT_DIRECTORY / "panasonic_c20_charge_discharge_gap.png"
    report_path = OUTPUT_DIRECTORY / "panasonic_c20_charge_discharge_comparison_report.md"
    write_csv(
        csv_path,
        common_soc,
        discharge_v,
        charge_v,
        delta_v,
        shape_coordinate,
        shape_discharge_v,
        shape_charge_v,
        shape_delta_v,
    )
    make_voltage_figure(
        voltage_plot_path,
        common_soc,
        discharge_v,
        charge_v,
        shape_coordinate,
        shape_discharge_v,
        shape_charge_v,
    )
    make_gap_figure(gap_plot_path, common_soc, delta_v, regional)

    qd_final = float(q_discharge_ah[-1])
    qc_final = float(q_charge_ah[-1])
    final_charge_soc = qc_final / Q_CELL_AH
    ratio = qc_final / qd_final
    segment_rows = []
    for segment in segments:
        idx = segment.indices
        segment_rows.append(
            f"| {segment.kind} | {segment.start}–{segment.end} | {idx.size} | "
            f"{data['Time'][idx[0]]:.6f}–{data['Time'][idx[-1]]:.6f} | "
            f"{data['Current'][idx].min():.5f} to {data['Current'][idx].max():.5f} | "
            f"{data['Voltage'][idx].min():.5f} to {data['Voltage'][idx].max():.5f} |"
        )

    stat_order = (
        ("Mean", "mean_v"),
        ("Median", "median_v"),
        ("RMS", "rms_v"),
        ("Minimum", "min_v"),
        ("Maximum", "max_v"),
        ("P5", "p5_v"),
        ("P95", "p95_v"),
    )
    overall_rows = "\n".join(
        f"| {label} | {fmt_mv(overall[key])} |" for label, key in stat_order
    )
    region_rows = "\n".join(
        "| " + name.capitalize() + " | "
        + ("[0, 1/3)" if name == "low" else "[1/3, 2/3)" if name == "middle" else f"[2/3, {common_soc_max:.9f}]")
        + f" | {values['n_grid']} | {fmt_mv(values['mean_v'])} | "
        + f"{fmt_mv(values['median_v'])} | {fmt_mv(values['rms_v'])} | "
        + f"{fmt_mv(values['min_v'])} | {fmt_mv(values['max_v'])} | "
        + f"{fmt_mv(values['p5_v'])} | {fmt_mv(values['p95_v'])} |"
        for name, values in regional.items()
    )

    # The first rest samples after each active branch are also shown so that the
    # long-rest endpoint is not mistaken for the immediate current interruption.
    report = f"""# Panasonic C/20 charge/discharge quasi-OCV comparison

## Scope and source

This is a data-only comparison. No electrochemical model or OCP function was run, no parameter was fitted, and no hysteresis term or branch correction was introduced.

- Source: `{SOURCE_PATH}`
- MAT fields: `{', '.join(columns)}`
- Samples: {data['Time'].size}
- Frozen common capacity: {Q_CELL_AH:.12f} Ah
- Sign observed from the signals: negative current is discharge; positive current is charge.
- Active/rest threshold used only for segmentation: |I| > {ACTIVE_CURRENT_THRESHOLD_A:.3f} A.
- Source duration: {data['Time'][-1] - data['Time'][0]:.6f} s ({(data['Time'][-1] - data['Time'][0]) / 3600.0:.6f} h).
- Battery-temperature range: {data['Battery_Temp_degC'].min():.5f} to {data['Battery_Temp_degC'].max():.5f} °C.
- Chamber-temperature range: {data['Chamber_Temp_degC'].min():.5f} to {data['Chamber_Temp_degC'].max():.5f} °C.

## Reproducible signal segmentation

The script classifies every sample from measured current and run-length encodes the result. It verifies the five-part order rather than using the expected indices as inputs. Indices below are zero-based and inclusive.

| Segment | Indices | Samples | Time span (s) | Current range (A) | Voltage range (V) |
|---|---:|---:|---:|---:|---:|
{chr(10).join(segment_rows)}

The detected active branches are exactly discharge indices {discharge.start}–{discharge.end} and charge indices {charge.start}–{charge.end}.

## Common-capacity SOC comparison

The primary axes are

`SOC_dis = 1 - Q_discharge_since_start / Q_cell_Ah`

and

`SOC_chg = Q_charge_since_charge_start / Q_cell_Ah`.

The charge starts at nominal SOC = 0 because it follows the complete low-rate discharge. Its endpoint was **not** forced to SOC = 1. Both measured curves were linearly interpolated on {COMMON_GRID_POINTS} points only over the exact common interval **[{common_soc_min:.12g}, {common_soc_max:.12g}]**. No extrapolation was used.

For `deltaV_qocv = V_C20_charge - V_C20_discharge`:

| Statistic | Gap (mV) |
|---|---:|
{overall_rows}

### Gap by available SOC region

| Region | SOC interval | Grid points | Mean (mV) | Median (mV) | RMS (mV) | Min (mV) | Max (mV) | P5 (mV) | P95 (mV) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
{region_rows}

The gap is positive at every comparison point and is strongly SOC-dependent. It has a pronounced endpoint excursion near nominal SOC = 0; that point compares two loaded measurements separated by the intermediate rest and must not be interpreted as an equilibrium material property.

## Endpoint and capacity diagnostic

| Quantity | Value |
|---|---:|
| Integrated discharge capacity | {qd_final:.12f} Ah |
| Integrated charge capacity | {qc_final:.12f} Ah |
| Charge/discharge capacity ratio | {ratio:.9f} ({100.0 * ratio:.6f}%) |
| Final charge SOC from frozen capacity | {final_charge_soc:.12f} ({100.0 * final_charge_soc:.6f}%) |
| Discharge final loaded voltage | {data['Voltage'][discharge.end]:.5f} V |
| Discharge first post-current rest voltage | {data['Voltage'][intermediate_rest.start]:.5f} V |
| Discharge recovered voltage before charge | {data['Voltage'][intermediate_rest.end]:.5f} V |
| Charge initial loaded voltage | {data['Voltage'][charge.start]:.5f} V |
| Charge final loaded voltage | {data['Voltage'][charge.end]:.5f} V |
| Charge first post-current rest voltage | {data['Voltage'][final_rest.start]:.5f} V |
| Charge final recorded rest voltage | {data['Voltage'][final_rest.end]:.5f} V |

The 4.2 V loaded endpoint does not make the charge branch complete: coulomb counting reaches only {100.0 * final_charge_soc:.3f}% of the frozen capacity. The integrated discharge capacity differs from the frozen value by {1.0e9 * (qd_final - Q_CELL_AH):.3f} nAh.

## SHAPE-ONLY NORMALIZED COMPARISON

As a secondary diagnostic only, each branch was mapped by its own transferred capacity. This is not the primary SOC comparison and is not used to redefine capacity. Its gap statistics are: mean {fmt_mv(shape_stats['mean_v'])} mV, median {fmt_mv(shape_stats['median_v'])} mV, RMS {fmt_mv(shape_stats['rms_v'])} mV, minimum {fmt_mv(shape_stats['min_v'])} mV, maximum {fmt_mv(shape_stats['max_v'])} mV, P5 {fmt_mv(shape_stats['p5_v'])} mV, and P95 {fmt_mv(shape_stats['p95_v'])} mV. The lower central gap after this remapping shows that unequal effective capacity mapping contributes to the primary difference, but it does not remove the directional voltage separation.

## Relation to the existing CC–CV charge validation

1. **Systematic gap:** Yes. The measured low-rate charge voltage exceeds the discharge voltage over the entire shared SOC interval (minimum {fmt_mv(overall['min_v'])} mV; mean {fmt_mv(overall['mean_v'])} mV).
2. **Comparable magnitude:** Yes. The mean measured gap is the same order as, and larger than, the prior CC bias magnitude ({abs(CCCV_CC_BIAS_MV):.1f} mV) and CV bias magnitude ({abs(CCCV_CV_BIAS_MV):.1f} mV: mean-gap ratios {1000.0 * overall['mean_v'] / abs(CCCV_CC_BIAS_MV):.2f} and {1000.0 * overall['mean_v'] / abs(CCCV_CV_BIAS_MV):.2f}, respectively).
3. **SOC dependence:** Strong. Regional mean gaps are {fmt_mv(regional['low']['mean_v'])} mV (low), {fmt_mv(regional['middle']['mean_v'])} mV (middle), and {fmt_mv(regional['high']['mean_v'])} mV (high), with the largest endpoint difference near nominal SOC = 0.
4. **OCP implication:** The positive, sustained low-rate charge/discharge quasi-OCV gap is consistent with a discharge-derived OCP producing a substantial part of the observed negative charge-direction voltage bias. It **could explain** a substantial part, but this diagnostic does not prove causation or uniquely allocate the error: both branches remain under approximately C/20 load, the endpoint protocols include rests, and the charge transfers less capacity.

This result is therefore described only as a **low-rate charge/discharge quasi-OCV gap**, not thermodynamic hysteresis.
"""
    report_path.write_text(report, encoding="utf-8")

    print(f"Source: {SOURCE_PATH}")
    print(f"Segments: {[(s.kind, s.start, s.end) for s in segments]}")
    print(f"Common SOC: [{common_soc_min:.12g}, {common_soc_max:.12g}]")
    print(f"Discharge capacity: {qd_final:.12f} Ah")
    print(f"Charge capacity: {qc_final:.12f} Ah")
    print(f"Final charge SOC: {final_charge_soc:.12f}")
    print(f"Mean/median/RMS gap: {fmt_mv(overall['mean_v'])}/"
          f"{fmt_mv(overall['median_v'])}/{fmt_mv(overall['rms_v'])} mV")
    print(f"Wrote: {csv_path}")
    print(f"Wrote: {report_path}")
    print(f"Wrote: {voltage_plot_path}")
    print(f"Wrote: {gap_plot_path}")


if __name__ == "__main__":
    main()
