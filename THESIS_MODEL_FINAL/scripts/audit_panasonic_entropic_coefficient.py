"""Audit whether Panasonic dUeq/dT(SOC) is identifiable from local data.

This script is intentionally data- and source-audit only. It scans every MAT
profile in the Panasonic tree from measured signals, but it does not import or
run the electrochemical model, fit parameters, or modify an OCP.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.io import loadmat


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DATASET_ROOT = Path(
    r"C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1"
    r"\Panasonic 18650PF Data"
)
REFERENCE_PATH = DATASET_ROOT / (
    r"Panasonic 18650PF Data\25degC"
    r"\C20 OCV and 1C discharge tests_start_of_tests"
    r"\05-08-17_13.26 C20 OCV Test_C20_25dC.mat"
)
README_PATH = DATASET_ROOT / "Readme file - desc of tests performed.txt"
REPORT_PATH = REPOSITORY_ROOT / "thesis_data" / "panasonic_entropic_data_audit.md"

Q_CELL_AH = 2.994979138519
ACTIVE_THRESHOLD_A = 0.02
CANDIDATE_MAX_MEDIAN_CURRENT_A = 0.35
CANDIDATE_MIN_DURATION_H = 1.0
CANDIDATE_MIN_CAPACITY_AH = 0.10
MOST_SOC_FRACTION = 0.75
CONSTANT_CURRENT_MAX_CV = 0.05
CONSTANT_CURRENT_MAX_P95_P5_RATIO = 1.10
LONG_REST_MIN_H = 0.90


@dataclass(frozen=True)
class MatProfile:
    path: Path
    time_s: np.ndarray
    current_a: np.ndarray
    voltage_v: np.ndarray
    ah: np.ndarray
    battery_temp_c: np.ndarray
    chamber_temp_c: np.ndarray


@dataclass(frozen=True)
class CandidateBranch:
    path: Path
    start: int
    end: int
    direction: str
    nominal_temperature: str
    duration_h: float
    transferred_ah: float
    ah_counter_change: float
    ah_relative_difference: float
    current_min_a: float
    current_median_a: float
    current_mean_a: float
    current_max_a: float
    current_p5_a: float
    current_p95_a: float
    current_cv: float
    voltage_min_v: float
    voltage_max_v: float
    battery_temp_min_c: float
    battery_temp_max_c: float
    battery_temp_mean_c: float
    chamber_temp_min_c: float
    chamber_temp_max_c: float
    chamber_temp_mean_c: float
    chamber_finite_fraction: float
    max_rest_h: float
    spans_most_soc: bool
    has_long_rest: bool
    ah_counter_usable: bool
    constant_current: bool
    controlled_temperature: bool
    usable_low_rate: bool


def finite_stats(values: np.ndarray) -> tuple[float, float, float, float]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return np.nan, np.nan, np.nan, 0.0
    return (
        float(np.min(finite)),
        float(np.max(finite)),
        float(np.mean(finite)),
        float(finite.size / values.size),
    )


def nominal_temperature(path: Path) -> str:
    relative = str(path.relative_to(DATASET_ROOT)).lower()
    if "-20degc" in relative:
        label = "-20 °C"
    elif "-10degc" in relative:
        label = "-10 °C"
    elif "10degc" in relative:
        label = "10 °C"
    elif "0degc" in relative:
        label = "0 °C"
    elif "25degc" in relative:
        label = "25 °C"
    else:
        label = "not encoded"
    if "trise" in relative:
        label += " (temperature-rise test)"
    return label


def load_profile(path: Path) -> MatProfile:
    source = loadmat(
        path,
        squeeze_me=True,
        struct_as_record=False,
        variable_names=["meas"],
    )
    measurement = source.get("meas")
    if measurement is None or not hasattr(measurement, "_fieldnames"):
        raise ValueError(f"No 'meas' structure in {path}")
    required = (
        "Time",
        "Current",
        "Voltage",
        "Ah",
        "Battery_Temp_degC",
        "Chamber_Temp_degC",
    )
    missing = [name for name in required if not hasattr(measurement, name)]
    if missing:
        raise ValueError(f"Missing fields {missing} in {path}")
    arrays = {
        name: np.asarray(getattr(measurement, name), dtype=float).reshape(-1)
        for name in required
    }
    if len({array.size for array in arrays.values()}) != 1:
        raise ValueError(f"Inconsistent channel lengths in {path}")
    if not all(np.all(np.isfinite(arrays[name])) for name in ("Time", "Current", "Voltage", "Ah")):
        raise ValueError(f"Non-finite electrical data in {path}")
    if np.any(np.diff(arrays["Time"]) < 0.0):
        raise ValueError(f"Decreasing elapsed time in {path}")
    return MatProfile(
        path=path,
        time_s=arrays["Time"],
        current_a=arrays["Current"],
        voltage_v=arrays["Voltage"],
        ah=arrays["Ah"],
        battery_temp_c=arrays["Battery_Temp_degC"],
        chamber_temp_c=arrays["Chamber_Temp_degC"],
    )


def run_boundaries(state: np.ndarray) -> np.ndarray:
    return np.r_[0, np.flatnonzero(np.diff(state) != 0) + 1, state.size]


def maximum_rest_duration_h(profile: MatProfile) -> float:
    rest = np.abs(profile.current_a) <= ACTIVE_THRESHOLD_A
    edges = run_boundaries(rest.astype(int))
    durations = [
        float(profile.time_s[end - 1] - profile.time_s[start]) / 3600.0
        for start, end in zip(edges[:-1], edges[1:])
        if rest[start] and end - start >= 2
    ]
    return max(durations, default=0.0)


def candidate_branches(profile: MatProfile) -> list[CandidateBranch]:
    state = np.zeros(profile.current_a.size, dtype=int)
    state[profile.current_a > ACTIVE_THRESHOLD_A] = 1
    state[profile.current_a < -ACTIVE_THRESHOLD_A] = -1
    branches: list[CandidateBranch] = []
    max_rest_h = maximum_rest_duration_h(profile)
    for start, stop in zip(run_boundaries(state)[:-1], run_boundaries(state)[1:]):
        sign = int(state[start])
        if sign == 0 or stop - start < 2:
            continue
        indices = np.arange(start, stop)
        duration_h = float(profile.time_s[stop - 1] - profile.time_s[start]) / 3600.0
        magnitudes = np.abs(profile.current_a[indices])
        transferred_ah = abs(
            float(
                np.trapezoid(
                    profile.current_a[indices], profile.time_s[indices]
                )
                / 3600.0
            )
        )
        median_current = float(np.median(magnitudes))
        if not (
            duration_h >= CANDIDATE_MIN_DURATION_H
            and transferred_ah >= CANDIDATE_MIN_CAPACITY_AH
            and median_current <= CANDIDATE_MAX_MEDIAN_CURRENT_A
        ):
            continue
        current_mean = float(np.mean(magnitudes))
        p5, p95 = np.percentile(magnitudes, [5.0, 95.0])
        current_cv = float(np.std(magnitudes) / current_mean)
        ah_change = abs(float(profile.ah[stop - 1] - profile.ah[start]))
        ah_rel = abs(transferred_ah - ah_change) / transferred_ah
        tb_min, tb_max, tb_mean, _ = finite_stats(profile.battery_temp_c[indices])
        tc_min, tc_max, tc_mean, tc_fraction = finite_stats(profile.chamber_temp_c[indices])
        constant_current = bool(
            current_cv <= CONSTANT_CURRENT_MAX_CV
            and p95 / p5 <= CONSTANT_CURRENT_MAX_P95_P5_RATIO
        )
        spans_most_soc = bool(transferred_ah / Q_CELL_AH >= MOST_SOC_FRACTION)
        controlled_temperature = bool(
            tc_fraction >= 0.99
            and tc_max - tc_min <= 2.0
            and tb_max - tb_min <= 3.0
        )
        ah_counter_usable = bool(np.isfinite(ah_rel) and ah_rel <= 0.01)
        usable = bool(
            constant_current
            and spans_most_soc
            and controlled_temperature
            and ah_counter_usable
            and profile.voltage_v[indices].max() - profile.voltage_v[indices].min() >= 1.0
        )
        branches.append(
            CandidateBranch(
                path=profile.path,
                start=int(start),
                end=int(stop - 1),
                direction="charge" if sign > 0 else "discharge",
                nominal_temperature=nominal_temperature(profile.path),
                duration_h=duration_h,
                transferred_ah=transferred_ah,
                ah_counter_change=ah_change,
                ah_relative_difference=ah_rel,
                current_min_a=float(np.min(magnitudes)),
                current_median_a=median_current,
                current_mean_a=current_mean,
                current_max_a=float(np.max(magnitudes)),
                current_p5_a=float(p5),
                current_p95_a=float(p95),
                current_cv=current_cv,
                voltage_min_v=float(np.min(profile.voltage_v[indices])),
                voltage_max_v=float(np.max(profile.voltage_v[indices])),
                battery_temp_min_c=tb_min,
                battery_temp_max_c=tb_max,
                battery_temp_mean_c=tb_mean,
                chamber_temp_min_c=tc_min,
                chamber_temp_max_c=tc_max,
                chamber_temp_mean_c=tc_mean,
                chamber_finite_fraction=tc_fraction,
                max_rest_h=max_rest_h,
                spans_most_soc=spans_most_soc,
                has_long_rest=max_rest_h >= LONG_REST_MIN_H,
                ah_counter_usable=ah_counter_usable,
                constant_current=constant_current,
                controlled_temperature=controlled_temperature,
                usable_low_rate=usable,
            )
        )
    return branches


def fmt(value: float, digits: int = 3) -> str:
    return "missing" if not np.isfinite(value) else f"{value:.{digits}f}"


def yes_no(value: bool) -> str:
    return "yes" if value else "no"


def source_line(path: Path, exact_text: str) -> int:
    for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
        if exact_text in line:
            return number
    raise ValueError(f"Expected source expression not found in {path}: {exact_text}")


def rejection_reason(branch: CandidateBranch) -> str:
    reasons = []
    if not branch.constant_current:
        reasons.append("variable-current/CC–CV segment")
    if not branch.spans_most_soc:
        reasons.append("insufficient capacity coverage")
    if not branch.controlled_temperature:
        reasons.append("missing or changing temperature")
    if branch.voltage_max_v - branch.voltage_min_v < 1.0:
        reasons.append("limited voltage/SOC interval")
    return "; ".join(reasons) or "usable"


def main() -> None:
    if not DATASET_ROOT.is_dir() or not REFERENCE_PATH.is_file() or not README_PATH.is_file():
        raise FileNotFoundError("Panasonic dataset root, reference MAT, or README is missing.")

    mat_paths = sorted(DATASET_ROOT.rglob("*.mat"))
    csv_paths = sorted(DATASET_ROOT.rglob("*.csv"))
    non_eis_csv = [path for path in csv_paths if "eis" not in str(path.parent).lower()]
    if non_eis_csv:
        raise ValueError(f"Unexpected non-EIS CSV files: {non_eis_csv}")

    candidates: list[CandidateBranch] = []
    load_errors: list[str] = []
    for path in mat_paths:
        try:
            candidates.extend(candidate_branches(load_profile(path)))
        except Exception as exc:  # retained in report rather than silently skipped
            load_errors.append(f"{path}: {type(exc).__name__}: {exc}")

    if load_errors:
        raise RuntimeError("MAT scan was incomplete:\n" + "\n".join(load_errors))
    candidates.sort(key=lambda item: (str(item.path), item.start))
    usable = [item for item in candidates if item.usable_low_rate]
    usable_discharge = [item for item in usable if item.direction == "discharge"]
    usable_charge = [item for item in usable if item.direction == "charge"]
    if len(candidates) != 5 or len(usable) != 2:
        raise RuntimeError(
            f"Audit expectations changed: {len(candidates)} candidates, {len(usable)} usable."
        )
    if {item.path for item in usable} != {REFERENCE_PATH}:
        raise RuntimeError("A usable low-rate curve outside the reference file was found.")

    # A derivative cannot be inferred from charge/discharge direction mixing.
    distinct_discharge_temperatures = {item.nominal_temperature for item in usable_discharge}
    distinct_charge_temperatures = {item.nominal_temperature for item in usable_charge}
    identifiable = bool(
        len(distinct_discharge_temperatures) >= 2
        or len(distinct_charge_temperatures) >= 2
    )
    if identifiable:
        raise RuntimeError(
            "The current dataset unexpectedly passed the entropic-identifiability gate; "
            "coefficient implementation requires explicit review."
        )

    thermal_path = REPOSITORY_ROOT / "models" / "thermal.py"
    coupled_path = REPOSITORY_ROOT / "models" / "electrothermal_spm.py"
    grouped_path = REPOSITORY_ROOT / "models" / "grouped_spm.py"
    irr_line = source_line(
        thermal_path,
        "irreversible_w = current_a * (equilibrium_voltage_v - terminal_voltage_v)",
    )
    rev_line = source_line(
        thermal_path,
        "reversible_w = -current_a * temperature_k * coefficient",
    )
    total_line = source_line(thermal_path, "total_w = irreversible_w + reversible_w")
    wrapper_line = source_line(coupled_path, "terminal_voltage_v=terminal_voltage_v")
    ueq_line = source_line(coupled_path, "u_eq_v = electrical.U_p_v[index] - electrical.U_n_v[index]")
    surface_line = source_line(grouped_path, "u_n = np.asarray(U_n(x_n_surface), dtype=float)")

    core_rows = []
    temp_rows = []
    for identifier, branch in enumerate(candidates, 1):
        disposition = "usable low-rate reference" if branch.usable_low_rate else rejection_reason(branch)
        core_rows.append(
            f"| C{identifier} | `{branch.path}` | {branch.nominal_temperature} | "
            f"{branch.direction} | {branch.start}–{branch.end} | "
            f"{branch.current_min_a:.5f}/{branch.current_median_a:.5f}/{branch.current_max_a:.5f} | "
            f"{branch.duration_h:.6f} | {branch.transferred_ah:.9f} | "
            f"{branch.ah_counter_change:.9f} | {branch.voltage_min_v:.5f}–{branch.voltage_max_v:.5f} | "
            f"{yes_no(branch.spans_most_soc)} | {yes_no(branch.has_long_rest)} ({branch.max_rest_h:.3f} h max) | "
            f"{yes_no(branch.ah_counter_usable)} | {disposition} |"
        )
        temp_rows.append(
            f"| C{identifier} | {fmt(branch.chamber_temp_min_c)}–{fmt(branch.chamber_temp_max_c)} "
            f"(mean {fmt(branch.chamber_temp_mean_c)}) | "
            f"{fmt(branch.battery_temp_min_c)}–{fmt(branch.battery_temp_max_c)} "
            f"(mean {fmt(branch.battery_temp_mean_c)}) | "
            f"{100.0 * branch.chamber_finite_fraction:.1f}% | "
            f"{branch.current_p5_a:.5f}–{branch.current_p95_a:.5f}; CV={branch.current_cv:.3f} |"
        )

    usable_rows = []
    for branch in usable:
        coverage = 100.0 * branch.transferred_ah / Q_CELL_AH
        rests = (
            "initial 0.067 h, intermediate 1.000 h, final recorded 14.586 h"
            if branch.path == REFERENCE_PATH
            else f"maximum {branch.max_rest_h:.3f} h"
        )
        usable_rows.append(
            f"| `{branch.path}` | {branch.direction} | {branch.nominal_temperature} | "
            f"{branch.chamber_temp_min_c:.3f}–{branch.chamber_temp_max_c:.3f} | "
            f"{branch.battery_temp_min_c:.3f}–{branch.battery_temp_max_c:.3f} | "
            f"{branch.current_min_a:.5f}–{branch.current_max_a:.5f} A ({branch.direction}) | "
            f"{branch.duration_h:.6f} | {branch.transferred_ah:.9f} ({coverage:.3f}% Qcell) | "
            f"{branch.voltage_min_v:.5f}–{branch.voltage_max_v:.5f} | yes | {rests} | yes |"
        )

    report = f"""# Panasonic entropic-coefficient and thermal-data audit

## Decision

**PANASONIC ENTROPIC COEFFICIENT NOT IDENTIFIABLE FROM AVAILABLE DATA**

The local Panasonic tree contains one genuine C/20 experiment at one controlled temperature: the reference file has a 25 °C discharge and a subsequent 25 °C charge. Direction mixing is prohibited, so these two branches do not provide temperature replication. No second comparable low-rate discharge or charge curve exists at another controlled temperature. Consequently, no `dUeq/dT(SOC)` values were estimated, regressed, smoothed, plotted, or implemented.

No electrochemical model was run, no parameter was fitted, and no OCP or frozen electrical file was modified.

## Search coverage and reproducible candidate definition

- Dataset root: `{DATASET_ROOT}`
- Dataset README inspected: `{README_PATH}`
- MAT files inspected from their signals: **{len(mat_paths)} of {len(mat_paths)}**, with zero load failures.
- CSV files found: **{len(csv_paths)}**; all are inside `EIS` folders and are impedance spectra, not low-rate capacity curves.
- Reference C/20 file: `{REFERENCE_PATH}`
- Signals required from each MAT: `Time`, `Current`, `Voltage`, `Ah`, `Battery_Temp_degC`, and `Chamber_Temp_degC`.
- Active sign threshold: |I| > {ACTIVE_THRESHOLD_A:.3f} A. Contiguous positive and negative current runs were examined independently.
- Deliberately broad candidate gate: at least {CANDIDATE_MIN_DURATION_H:.1f} h, at least {CANDIDATE_MIN_CAPACITY_AH:.2f} Ah transferred, and median |I| no greater than {CANDIDATE_MAX_MEDIAN_CURRENT_A:.2f} A.
- A curve is marked usable only when it is near-constant current, covers at least {100.0 * MOST_SOC_FRACTION:.0f}% of frozen `Q_cell_Ah = {Q_CELL_AH:.12f} Ah`, has controlled temperature, a consistent Ah counter, and a broad voltage range.

The dataset README independently describes a C/20 discharge/charge only once, before HPPC/EIS/drive-cycle testing at 25, 10, 0, -10, and -20 °C. The conclusions below nevertheless come from the electrical and temperature signals rather than filenames alone.

## All signal-level low-current candidates

Indices are zero-based and inclusive. Current magnitude is shown as minimum/median/maximum over the active branch.

| ID | Exact path | Folder temperature | Direction | Indices | abs(I) min/median/max (A) | Duration (h) | Integrated Ah | Ah-counter Δ | Voltage range (V) | Most SOC? | Long rest in file? | Ah usable? | Disposition |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|---|
{chr(10).join(core_rows)}

### Candidate temperature and current-stability detail

| ID | Chamber °C range (mean) | Battery °C range (mean) | Chamber samples finite | abs(I) P5–P95 (A); coefficient of variation |
|---|---:|---:|---:|---:|
{chr(10).join(temp_rows)}

The three rejected candidates are not hidden low-rate OCV curves. Their current ranges extend into ordinary charge rates and then taper, their current coefficients of variation exceed 1.2, and they cover only restricted charge intervals. Two have no usable chamber-temperature channel; all three undergo substantial battery-temperature changes. Comparing them with the 25 °C C/20 branches would mix rate, direction/protocol, temperature history, and SOC coverage.

## Usable low-rate curves

| Exact path | Direction | Nominal temperature | Chamber °C | Battery °C | Current | Duration (h) | Transferred capacity | Voltage range (V) | Most SOC? | Rest context | Ah usable? |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|
{chr(10).join(usable_rows)}

The discharge is the previously used reference. Its chamber is exactly 25 °C and its current-integrated capacity is 2.994979139 Ah. The charge is also approximately C/20 at 25 °C, but it reaches only 2.613916708 Ah (87.277% of frozen cell capacity) before the 4.2 V cutoff. The branches are useful for the already documented direction/capacity comparison, not for a temperature derivative.

Long rest files at 0, -10, -20, and temperature-rise conditions were also encountered. They are isolated rest/conditioning records near limited states of charge, often with missing chamber values or large temperature transients. They do not supply full-SOC low-rate voltage curves and cannot replace the missing multi-temperature discharge series.

## Common-SOC and derivative gate

For any future discharge comparison, the required primary mapping remains

`SOC = 1 - Q_discharge_since_start / Q_cell_Ah`,

using actual transferred charge and `Q_cell_Ah = {Q_CELL_AH:.12f} Ah`. Charge would use its own direction-consistent coulomb count but must not be mixed into a discharge temperature regression. Interpolation would be allowed only on the intersection of measured SOC supports, with no endpoint stretching and no extrapolation.

Here there is only one usable discharge temperature and one usable charge temperature, both 25 °C. Therefore:

- number of temperatures available for a discharge-only regression: **1**;
- temperature span: **0 K**;
- common multi-temperature SOC interval: **none**;
- per-SOC slope, R², and temperature count: **not computed**.

The acceptance conditions fail before regression. A slope from these files would necessarily mix charge/discharge direction or unrelated dynamic/variable-temperature experiments and would be dominated by protocol and capacity-axis differences. No coefficient is invented and no literature coefficient is substituted.

## Heat-generation implementation audit

The existing implementation already declares discharge-positive model current (`I > 0` discharge, `I < 0` charge) and contains the requested signs and SI-unit structure:

- `{thermal_path}`, line {irr_line}: `Qdot_irr = I * (Ueq - V)`.
- `{thermal_path}`, line {rev_line}: `Qdot_rev = -I * T * dUeq_dT`.
- `{thermal_path}`, line {total_line}: `Qdot_gen = Qdot_irr + Qdot_rev`.

The generic heat function can evaluate either irreversible quantity depending on the voltage supplied, but the caller must label them separately:

```text
Qdot_irr_measured = I * (Ueq_average - V_measured)
Qdot_irr_model    = I * (Ueq_average - V_model)
Qdot_rev          = -I * T * dUeq_dT(SOC)
Qdot_gen          = Qdot_irr + Qdot_rev
```

The coupled wrapper currently passes its predicted terminal voltage at `{coupled_path}`, line {wrapper_line}; its exported `heat_irreversible_w` is therefore the **fully predictive model value**, not the measured-voltage calibration diagnostic. There is no separate measured-voltage heat output, so a later thermal-calibration script must calculate and retain `Qdot_irr_measured` under that explicit name rather than silently replacing the model value.

There is also an important state-definition finding. The current wrapper forms `Ueq` at `{coupled_path}`, line {ueq_line} from OCP arrays that are evaluated using **surface stoichiometries** (`{grouped_path}`, line {surface_line}). The thermal specification in this task instead requires equilibrium voltage from average electrode state:

```text
Ueq_average = U_p(x_p_avg) - U_n(x_n_avg)
```

That distinction must be explicit in the later Panasonic thermal calibration/validation interface. It was not changed here because the electrical architecture is frozen and this task authorizes an audit, not a model-equation change.

Reversible heat remains disabled when no coefficient is supplied. Because the identifiability gate failed, no Panasonic entropic module, coefficient CSV, coefficient plot, or entropic tests were created.

## Final scope statement

This audit stops at the required gate. `Cth`, `Gth`, electrical parameters, OCPs, stoichiometric windows, Arrhenius terms, and validation profiles remain untouched.
"""

    REPORT_PATH.parent.mkdir(parents=True, exist_ok=True)
    REPORT_PATH.write_text(report, encoding="utf-8")
    print(f"Scanned {len(mat_paths)} MAT files and {len(csv_paths)} EIS CSV files.")
    print(f"Signal-level candidates: {len(candidates)}; usable low-rate branches: {len(usable)}.")
    print("PANASONIC ENTROPIC COEFFICIENT NOT IDENTIFIABLE FROM AVAILABLE DATA")
    print(f"Wrote: {REPORT_PATH}")


if __name__ == "__main__":
    main()
