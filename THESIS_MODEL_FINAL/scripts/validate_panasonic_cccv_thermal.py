"""No-refit Panasonic 25 C CC-CV charge thermal validation.

Validation A uses measured voltage and measured battery temperature in the
heat calculation. Validation B uses the frozen electrical model voltage and
the recursively predicted thermal state. The frozen 1C-only one-node thermal
parameters are used throughout; no fitting or correction is performed.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from models.ocp_panasonic_qi_frozen import Q_CELL_AH, U_n, U_p  # noqa: E402
from models.panasonic_entropic import SUPPORTED_SOC_RANGE, dUeq_dT_panasonic  # noqa: E402
from models.panasonic_thermal_parameters import (  # noqa: E402
    C_TH_J_PER_K,
    G_TH_W_PER_K,
    TAU_TH_S,
)
from scripts.calibrate_panasonic_joint_1c_hppc_thermal import (  # noqa: E402
    exact_one_node,
    temperature_metrics,
)
from scripts.validate_panasonic_cccv_charge import (  # noqa: E402
    FARADAY_C_PER_MOL,
    MANIFEST_PATH,
    PARAMETERS,
    Q_N_C,
    Q_P_C,
    R_GAS_J_PER_MOL_K,
    RAW_SOURCE_PATH,
    REST_THRESHOLD_A,
    forward_charge,
    identify_phases,
    initial_stoichiometries,
    load_mat_channels,
    preceding_discharge_capacity_ah,
    remove_duplicate_times,
)


OUTPUT_DIRECTORY = REPOSITORY_ROOT / "thesis_data"
REPORT_PATH = OUTPUT_DIRECTORY / "panasonic_cccv_thermal_validation.md"
METRICS_PATH = OUTPUT_DIRECTORY / "panasonic_cccv_thermal_validation_metrics.csv"
TEMPERATURE_PLOT_PATH = OUTPUT_DIRECTORY / "panasonic_cccv_thermal_validation.png"
RESIDUAL_PLOT_PATH = OUTPUT_DIRECTORY / "panasonic_cccv_thermal_residual.png"
HEAT_PLOT_PATH = OUTPUT_DIRECTORY / "panasonic_cccv_heat_components.png"

LOWER_SOC, UPPER_SOC = SUPPORTED_SOC_RANGE
VALIDATION_A_CLASSIFICATION = "CHARGE-SIDE THERMAL TRANSFER LIMITED"
VALIDATION_B_CLASSIFICATION = "COUPLED CHARGE ELECTROTHERMAL TRANSFER FAILED"


@dataclass(frozen=True)
class RetainedData:
    source_columns: tuple[str, ...]
    source_raw_samples: int
    duplicate_rows_removed: int
    boundary_records: tuple[dict[str, object], ...]
    source_time_s: np.ndarray
    time_s: np.ndarray
    soc: np.ndarray
    model_current_a: np.ndarray
    measured_voltage_v: np.ndarray
    measured_temperature_k: np.ndarray
    ambient_temperature_k: np.ndarray
    exported_ah: np.ndarray
    x_n_avg: np.ndarray
    x_p_avg: np.ndarray
    x_n_surface: np.ndarray
    x_p_surface: np.ndarray
    measured_t_model_voltage_v: np.ndarray
    cc_mask: np.ndarray
    cv_mask: np.ndarray
    cv_start_retained_index: int
    source_cv_start_time_s: float
    source_cv_start_soc: float
    full_initial_soc: float
    full_final_soc: float
    full_charge_throughput_ah: float


@dataclass(frozen=True)
class CoupledResult:
    temperature_k: np.ndarray
    voltage_v: np.ndarray
    q_irr_w: np.ndarray
    q_rev_w: np.ndarray
    q_total_w: np.ndarray


def insert_soc_boundaries(
    arrays: dict[str, np.ndarray], initial_soc: float
) -> tuple[dict[str, np.ndarray], tuple[dict[str, object], ...]]:
    """Insert exact, linearly interpolated SOC boundary records."""

    output = {name: np.asarray(values, dtype=float).copy() for name, values in arrays.items()}
    records: list[dict[str, object]] = []
    for target in (LOWER_SOC, UPPER_SOC):
        soc = initial_soc + output["Ah"] / Q_CELL_AH
        exact = np.flatnonzero(np.isclose(soc, target, rtol=0.0, atol=2.0e-14))
        if exact.size:
            records.append(
                {
                    "target_soc": target,
                    "interpolated": False,
                    "left_index_before_insertion": int(exact[0]),
                    "right_index_before_insertion": int(exact[0]),
                    "fraction": 0.0,
                    "source_time_s": float(output["Time"][exact[0]]),
                }
            )
            continue
        right = int(np.searchsorted(soc, target, side="right"))
        left = right - 1
        if left < 0 or right >= soc.size or not soc[left] < target < soc[right]:
            raise ValueError(f"Cannot interpolate SOC boundary {target:g}.")
        fraction = float((target - soc[left]) / (soc[right] - soc[left]))
        boundary_values: dict[str, float] = {}
        for name, values in output.items():
            boundary_values[name] = float(values[left] + fraction * (values[right] - values[left]))
        boundary_values["Ah"] = float((target - initial_soc) * Q_CELL_AH)
        for name, values in output.items():
            output[name] = np.insert(values, right, boundary_values[name])
        records.append(
            {
                "target_soc": target,
                "interpolated": True,
                "left_index_before_insertion": left,
                "right_index_before_insertion": right,
                "fraction": fraction,
                "source_time_s": boundary_values["Time"],
            }
        )
    return output, tuple(records)


def prepare_data() -> RetainedData:
    manifest = MANIFEST_PATH.read_text(encoding="utf-8")
    if RAW_SOURCE_PATH.name not in manifest:
        raise RuntimeError("The selected CC-CV source is not recorded in the manifest.")
    if not np.isclose(Q_CELL_AH, 2.994979138519, rtol=0.0, atol=5.0e-13):
        raise RuntimeError("Frozen Panasonic cell capacity changed.")
    if not (
        np.isclose(C_TH_J_PER_K, 117.52317800099644, rtol=0.0, atol=5.0e-12)
        and np.isclose(G_TH_W_PER_K, 0.10911208974639978, rtol=0.0, atol=5.0e-13)
    ):
        raise RuntimeError("Frozen primary 1C-only thermal parameters changed.")

    source_columns, raw = load_mat_channels(RAW_SOURCE_PATH)
    source_raw_samples = int(raw["Time"].size)
    cleaned, duplicate_indices = remove_duplicate_times(raw)
    preceding_capacity_ah, _ = preceding_discharge_capacity_ah()
    initial_soc = 1.0 - preceding_capacity_ah / Q_CELL_AH
    augmented, boundary_records = insert_soc_boundaries(cleaned, initial_soc)
    source_soc = initial_soc + augmented["Ah"] / Q_CELL_AH
    if np.any(np.diff(source_soc) < -2.0e-13) or np.any((source_soc < 0.0) | (source_soc > 1.0)):
        raise ValueError("CC-CV SOC is nonmonotonic or leaves [0,1]; no clipping.")

    phases = identify_phases(augmented["Current"], augmented["Voltage"], augmented["Ah"])
    initial_x_n, initial_x_p = initial_stoichiometries(initial_soc)
    electrical = forward_charge(
        augmented["Time"],
        augmented["Current"],
        augmented["Ah"],
        augmented["Battery_Temp_degC"],
        initial_x_n,
        initial_x_p,
    )

    supported = (source_soc >= LOWER_SOC) & (source_soc <= UPPER_SOC)
    indices = np.flatnonzero(supported)
    if indices.size < 3 or not np.array_equal(indices, np.arange(indices[0], indices[-1] + 1)):
        raise ValueError("Entropy-supported charge interval is not contiguous.")
    start, end = int(indices[0]), int(indices[-1])
    if not np.isclose(source_soc[start], LOWER_SOC, rtol=0.0, atol=2.0e-13):
        raise RuntimeError("Lower supported SOC boundary is not exact.")
    if not np.isclose(source_soc[end], UPPER_SOC, rtol=0.0, atol=2.0e-13):
        raise RuntimeError("Upper supported SOC boundary is not exact.")

    sl = slice(start, end + 1)
    cc_mask = phases.cc[sl]
    cv_mask = phases.cv[sl]
    if not np.all(cc_mask | cv_mask):
        raise ValueError("Retained supported interval contains samples outside CC/CV charge.")
    cv_indices = np.flatnonzero(cv_mask)
    if cv_indices.size == 0:
        raise ValueError("CC-to-CV transition is not inside the supported interval.")
    model_current = -augmented["Current"][sl]
    if np.any(model_current[cc_mask] >= 0.0):
        raise RuntimeError("CC current sign conversion failed: I_model must be negative.")

    return RetainedData(
        source_columns=source_columns,
        source_raw_samples=source_raw_samples,
        duplicate_rows_removed=int(duplicate_indices.size),
        boundary_records=boundary_records,
        source_time_s=augmented["Time"][sl],
        time_s=augmented["Time"][sl] - augmented["Time"][start],
        soc=source_soc[sl],
        model_current_a=model_current,
        measured_voltage_v=augmented["Voltage"][sl],
        measured_temperature_k=augmented["Battery_Temp_degC"][sl] + 273.15,
        ambient_temperature_k=augmented["Chamber_Temp_degC"][sl] + 273.15,
        exported_ah=augmented["Ah"][sl],
        x_n_avg=electrical.x_n_avg[sl],
        x_p_avg=electrical.x_p_avg[sl],
        x_n_surface=electrical.x_n_surface[sl],
        x_p_surface=electrical.x_p_surface[sl],
        measured_t_model_voltage_v=electrical.voltage_v[sl],
        cc_mask=cc_mask,
        cv_mask=cv_mask,
        cv_start_retained_index=int(cv_indices[0]),
        source_cv_start_time_s=float(augmented["Time"][phases.cv_start_index]),
        source_cv_start_soc=float(source_soc[phases.cv_start_index]),
        full_initial_soc=float(initial_soc),
        full_final_soc=float(source_soc[-1]),
        full_charge_throughput_ah=float(
            augmented["Ah"][phases.last_charge_index] - augmented["Ah"][0]
        ),
    )


def frozen_voltage(
    model_current_a: np.ndarray | float,
    temperature_k: np.ndarray | float,
    x_n_surface: np.ndarray | float,
    x_p_surface: np.ndarray | float,
) -> np.ndarray:
    current = np.asarray(model_current_a, dtype=float)
    temperature = np.asarray(temperature_k, dtype=float)
    x_n = np.asarray(x_n_surface, dtype=float)
    x_p = np.asarray(x_p_surface, dtype=float)
    kinetic_factor = 2.0 * R_GAS_J_PER_MOL_K * temperature / FARADAY_C_PER_MOL
    eta_n = kinetic_factor * np.arcsinh(
        current
        / (6.0 * Q_N_C * PARAMETERS["d_n"] * np.sqrt(x_n * (1.0 - x_n)))
    )
    eta_p = kinetic_factor * np.arcsinh(
        current
        / (6.0 * Q_P_C * PARAMETERS["d_p"] * np.sqrt(x_p * (1.0 - x_p)))
    )
    return np.asarray(U_p(x_p) - U_n(x_n) + eta_p - eta_n - PARAMETERS["R0"] * current)


def one_node_step(
    temperature_k: float,
    ambient_k: float,
    heat_w: float,
    dt_s: float,
) -> float:
    x = G_TH_W_PER_K * dt_s / C_TH_J_PER_K
    approach = -np.expm1(-x)
    phi = approach / x if x != 0.0 else 1.0
    return float(
        (1.0 - approach) * temperature_k
        + approach * ambient_k
        + heat_w * dt_s * phi / C_TH_J_PER_K
    )


def validation_b(data: RetainedData, ueq_average_v: np.ndarray, dudt_v_per_k: np.ndarray) -> CoupledResult:
    count = data.time_s.size
    temperature = np.empty(count)
    voltage = np.empty(count)
    q_irr = np.empty(count)
    q_rev = np.empty(count)
    q_total = np.empty(count)
    temperature[0] = data.measured_temperature_k[0]
    for index in range(count):
        voltage[index] = float(
            frozen_voltage(
                data.model_current_a[index],
                temperature[index],
                data.x_n_surface[index],
                data.x_p_surface[index],
            )
        )
        q_irr[index] = data.model_current_a[index] * (
            ueq_average_v[index] - voltage[index]
        )
        q_rev[index] = (
            -data.model_current_a[index] * temperature[index] * dudt_v_per_k[index]
        )
        q_total[index] = q_irr[index] + q_rev[index]
        if index + 1 < count:
            temperature[index + 1] = one_node_step(
                temperature[index],
                data.ambient_temperature_k[index],
                q_total[index],
                data.time_s[index + 1] - data.time_s[index],
            )
    if not all(np.all(np.isfinite(values)) for values in (temperature, voltage, q_irr, q_rev, q_total)):
        raise ValueError("Non-finite fully predictive electrothermal trajectory.")
    return CoupledResult(temperature, voltage, q_irr, q_rev, q_total)


def subset_metrics(data: RetainedData, predicted_k: np.ndarray, mask: np.ndarray) -> dict[str, float]:
    indices = np.flatnonzero(mask)
    values = temperature_metrics(data.measured_temperature_k[mask], predicted_k[mask])
    values.update(
        {
            "sample_count": int(indices.size),
            "duration_s": float(data.time_s[indices[-1]] - data.time_s[indices[0]]),
            "soc_start": float(data.soc[indices[0]]),
            "soc_end": float(data.soc[indices[-1]]),
            "current_min_a": float(np.min(data.model_current_a[mask])),
            "current_max_a": float(np.max(data.model_current_a[mask])),
        }
    )
    return values


def zoh_energy(power_w: np.ndarray, time_s: np.ndarray) -> float:
    return float(np.sum(power_w[:-1] * np.diff(time_s)))


def heat_audit(data: RetainedData, q_irr: np.ndarray, q_rev: np.ndarray) -> dict[str, float | bool]:
    q_total = q_irr + q_rev
    duration = float(data.time_s[-1] - data.time_s[0])
    e_irr = zoh_energy(q_irr, data.time_s)
    e_rev = zoh_energy(q_rev, data.time_s)
    e_total = zoh_energy(q_total, data.time_s)
    dt = np.diff(data.time_s)
    interval_dudt = np.asarray(dUeq_dT_panasonic(data.soc[:-1]))
    sign_expected = np.sign(q_rev[:-1]) == np.sign(interval_dudt)
    nonzero = (q_rev[:-1] != 0.0) | (interval_dudt != 0.0)
    return {
        "e_irr_j": e_irr,
        "e_rev_j": e_rev,
        "e_total_j": e_total,
        "mean_q_irr_w": e_irr / duration,
        "mean_q_rev_w": e_rev / duration,
        "mean_q_total_w": e_total / duration,
        "signed_reversible_total_energy_fraction": e_rev / e_total,
        "absolute_reversible_heat_fraction": zoh_energy(np.abs(q_rev), data.time_s)
        / zoh_energy(np.abs(q_total), data.time_s),
        "time_fraction_q_rev_exothermic": float(np.sum(dt[q_rev[:-1] > 0.0]) / duration),
        "time_fraction_q_rev_endothermic": float(np.sum(dt[q_rev[:-1] < 0.0]) / duration),
        "charge_reversible_sign_rule_verified": bool(np.all(sign_expected[nonzero])),
    }


def electrical_heat_difference(
    data: RetainedData,
    q_irr_measured: np.ndarray,
    q_irr_model: np.ndarray,
    model_voltage_v: np.ndarray,
) -> dict[str, float]:
    difference = q_irr_measured - q_irr_model
    dt = np.diff(data.time_s)
    duration = float(data.time_s[-1] - data.time_s[0])
    return {
        "mean_difference_w": float(np.sum(difference[:-1] * dt) / duration),
        "rms_difference_w": float(np.sqrt(np.sum(difference[:-1] ** 2 * dt) / duration)),
        "integrated_difference_j": float(np.sum(difference[:-1] * dt)),
        "identity_max_abs_error_w": float(
            np.max(
                np.abs(
                    difference
                    - data.model_current_a
                    * (model_voltage_v - data.measured_voltage_v)
                )
            )
        ),
    }


def append_rows(rows: list[dict[str, object]], scope: str, subset: str, values: dict[str, object], n_samples: int) -> None:
    for metric, value in values.items():
        unit = "1"
        if metric.endswith("_k2"):
            unit = "K^2"
        elif metric.endswith("_k"):
            unit = "K"
        elif metric.endswith("_c"):
            unit = "degC"
        elif metric.endswith("_v"):
            unit = "V"
        elif metric.endswith("_w"):
            unit = "W"
        elif metric.endswith("_j"):
            unit = "J"
        elif metric.endswith("_s"):
            unit = "s"
        elif metric.endswith("_a"):
            unit = "A"
        elif metric.endswith("_ah"):
            unit = "Ah"
        rows.append(
            {
                "scope": scope,
                "subset": subset,
                "metric": metric,
                "value": value,
                "unit": unit,
                "n_samples": n_samples,
            }
        )


def write_metrics(
    data: RetainedData,
    audit: dict[str, object],
    metrics: dict[str, dict[str, float]],
    phases: dict[str, dict[str, dict[str, float]]],
    heats: dict[str, dict[str, float | bool]],
    difference: dict[str, float],
) -> None:
    rows: list[dict[str, object]] = []
    append_rows(rows, "data_audit", "retained_interval", audit, data.time_s.size)
    for name, values in metrics.items():
        append_rows(rows, "temperature_metrics", name, values, data.time_s.size)
    for path, path_phases in phases.items():
        for phase, values in path_phases.items():
            append_rows(rows, "phase_metrics", f"{path}_{phase}", values, int(values["sample_count"]))
    for name, values in heats.items():
        append_rows(rows, "heat_audit", name, values, data.time_s.size)
    append_rows(rows, "electrical_error_to_heat", "validation_b", difference, data.time_s.size)
    append_rows(
        rows,
        "classification",
        "validation_a",
        {"classification": VALIDATION_A_CLASSIFICATION},
        data.time_s.size,
    )
    append_rows(
        rows,
        "classification",
        "validation_b",
        {"classification": VALIDATION_B_CLASSIFICATION},
        data.time_s.size,
    )
    with METRICS_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def metric_table_row(name: str, values: dict[str, float]) -> str:
    return (
        f"| {name} | {values['rmse_k']:.9f} | {values['mae_k']:.9f} | "
        f"{values['mean_signed_error_k']:+.9f} | {values['max_abs_error_k']:.9f} | "
        f"{values['r2']:.9f} | {values['predicted_final_temperature_c']:.9f} |"
    )


def phase_table_row(name: str, values: dict[str, float]) -> str:
    return (
        f"| {name} | {int(values['sample_count'])} | {values['duration_s']:.6f} | "
        f"{values['soc_start']:.9f}-{values['soc_end']:.9f} | "
        f"{values['current_min_a']:.6f} to {values['current_max_a']:.6f} | "
        f"{values['rmse_k']:.9f} | {values['mae_k']:.9f} | "
        f"{values['mean_signed_error_k']:+.9f} | {values['max_abs_error_k']:.9f} |"
    )


def write_report(
    data: RetainedData,
    audit: dict[str, object],
    metrics: dict[str, dict[str, float]],
    phases: dict[str, dict[str, dict[str, float]]],
    heats: dict[str, dict[str, float | bool]],
    difference: dict[str, float],
) -> None:
    a = metrics["validation_a"]
    b = metrics["validation_b"]
    isolation = metrics["model_voltage_measured_t"]
    lines = [
        "# Panasonic 25 C CC-CV charge thermal validation",
        "",
        "## Strict no-refit scope",
        "",
        f"The frozen primary 1C-only one-node values `Cth={C_TH_J_PER_K:.12f} J/K`, `Gth={G_TH_W_PER_K:.12f} W/K`, and `tau={TAU_TH_S:.12f} s` were used unchanged. No joint 1C+HPPC parameter, optimizer, offset, correction, or architecture change was used.",
        "Validation B feeds its recursively predicted temperature into both the frozen Butler-Volmer voltage calculation and reversible heat after the shared measured initial condition; measured battery temperature is not used afterward, and no Arrhenius dependence is introduced.",
        "",
        "## Dataset and supported interval",
        "",
        f"- Exact source: `{RAW_SOURCE_PATH}`",
        f"- Raw columns: `{', '.join(data.source_columns)}`",
        f"- Raw samples: {data.source_raw_samples}; duplicate rows removed using the previously audited rule: {data.duplicate_rows_removed}.",
        f"- Complete charge SOC: {data.full_initial_soc:.12f}-{data.full_final_soc:.12f}. Retained entropy-supported SOC: **{data.soc[0]:.12f}-{data.soc[-1]:.12f}**, with both endpoints inserted by exact linear interpolation.",
        f"- Complete audited charge throughput: {data.full_charge_throughput_ah:.12f} Ah; retained supported-interval throughput: {audit['charge_throughput_ah']:.12f} Ah.",
        f"- Source start/end times: {data.source_time_s[0]:.9f}-{data.source_time_s[-1]:.9f} s; retained duration: {audit['duration_s']:.9f} s; samples: {data.time_s.size}.",
        f"- Model-current range: {audit['current_min_a']:.6f} to {audit['current_max_a']:.6f} A; CC values are numerically negative after the single loader-boundary sign conversion.",
        f"- Measured voltage: {audit['voltage_min_v']:.6f}-{audit['voltage_max_v']:.6f} V; battery temperature: {audit['battery_temperature_min_c']:.6f}-{audit['battery_temperature_max_c']:.6f} degC; chamber temperature: {audit['chamber_temperature_min_c']:.6f}-{audit['chamber_temperature_max_c']:.6f} degC.",
        f"- CC-to-CV transition begins at source time {data.source_cv_start_time_s:.9f} s and SOC {data.source_cv_start_soc:.12f}. It is inside the supported interval; SOC=0.90 occurs during CV, so no later CV or final-rest data were admitted.",
        "",
        "## Overall thermal validation",
        "",
        "| Path | RMSE [K] | MAE [K] | Bias [K] | Max abs [K] | R2 | Final predicted [degC] |",
        "|---|---:|---:|---:|---:|---:|---:|",
        metric_table_row("Validation A: measured voltage/measured T heat", a),
        metric_table_row("Validation B: fully predictive electrothermal", b),
        metric_table_row("Model-voltage / measured-T reversible diagnostic", isolation),
        "",
        f"Measured initial/final temperatures are {a['measured_final_temperature_c']-a['measured_temperature_rise_k']:.9f}/{a['measured_final_temperature_c']:.9f} degC. Validation A initial/final predictions are {a['predicted_final_temperature_c']-a['predicted_temperature_rise_k']:.9f}/{a['predicted_final_temperature_c']:.9f} degC; measured/predicted rises are {a['measured_temperature_rise_k']:+.9f}/{a['predicted_temperature_rise_k']:+.9f} K. Validation B initial/final predictions are {b['predicted_final_temperature_c']-b['predicted_temperature_rise_k']:.9f}/{b['predicted_final_temperature_c']:.9f} degC; its predicted rise is {b['predicted_temperature_rise_k']:+.9f} K.",
        "",
        "## CC/CV diagnostics",
        "",
        "| Path/phase | Samples | Duration [s] | SOC range | Model current [A] | RMSE [K] | MAE [K] | Bias [K] | Max abs [K] |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for path in ("validation_a", "validation_b", "model_voltage_measured_t"):
        for phase in ("cc", "cv"):
            lines.append(phase_table_row(f"{path} / {phase.upper()}", phases[path][phase]))
    measured_heat = heats["measured_voltage"]
    model_heat = heats["fully_predictive"]
    reversible_effect = heats["measured_voltage_reversible_temperature_effect"]
    lines.extend(
        [
            "",
            "## Charge-side heat audit",
            "",
            "| Heat path | Eirr [J] | Erev [J] | Etotal [J] | Mean Qirr [W] | Mean Qrev [W] | Mean total [W] |",
            "|---|---:|---:|---:|---:|---:|---:|",
            f"| Measured-voltage / measured-T | {measured_heat['e_irr_j']:.9f} | {measured_heat['e_rev_j']:.9f} | {measured_heat['e_total_j']:.9f} | {measured_heat['mean_q_irr_w']:.9f} | {measured_heat['mean_q_rev_w']:.9f} | {measured_heat['mean_q_total_w']:.9f} |",
            f"| Fully predictive | {model_heat['e_irr_j']:.9f} | {model_heat['e_rev_j']:.9f} | {model_heat['e_total_j']:.9f} | {model_heat['mean_q_irr_w']:.9f} | {model_heat['mean_q_rev_w']:.9f} | {model_heat['mean_q_total_w']:.9f} |",
            "",
            f"For the primary measured-voltage path, signed `Erev/Etotal={measured_heat['signed_reversible_total_energy_fraction']:+.9f}` and integral-absolute reversible fraction is `{measured_heat['absolute_reversible_heat_fraction']:.9f}`. Reversible heat is exothermic for {100.0*float(measured_heat['time_fraction_q_rev_exothermic']):.6f}% and endothermic for {100.0*float(measured_heat['time_fraction_q_rev_endothermic']):.6f}% of retained time.",
            f"The charge sign rule was verified numerically: **{measured_heat['charge_reversible_sign_rule_verified']}**. With I<0, positive dUeq/dT gives positive/exothermic reversible heat and negative dUeq/dT gives negative/endothermic heat.",
            f"Relative to an irreversible-only measured-voltage propagation, reversible heat changes the final predicted temperature by `{reversible_effect['final_temperature_delta_k']:+.9f} K`, with maximum absolute trajectory change `{reversible_effect['max_abs_temperature_delta_k']:.9f} K` and RMS change `{reversible_effect['rms_temperature_delta_k']:.9f} K`.",
            "",
            "## Electrical-error-to-heat diagnostic",
            "",
            f"`Qirr_measured-Qirr_model = I*(Vmodel-Vmeasured)` has mean `{difference['mean_difference_w']:+.9f} W`, RMS `{difference['rms_difference_w']:.9f} W`, and time integral `{difference['integrated_difference_j']:+.9f} J`. Numerical identity error is `{difference['identity_max_abs_error_w']:.3e} W`.",
            "The frozen electrical model predominantly underpredicts charge voltage. Because model current is negative, `I*(Vmodel-Vmeasured)` is therefore positive: the measured-voltage path supplies more irreversible heat than the model-voltage path. This is consistent with the fully predictive trajectory being substantially cooler.",
            "",
            "## Classifications and thesis interpretation",
            "",
            f"### Validation A: {VALIDATION_A_CLASSIFICATION}",
            "",
            f"### Validation B: {VALIDATION_B_CLASSIFICATION}",
            "",
            f"1. **Does the discharge-calibrated thermal pair transfer to charge with measured-voltage heat?** **Only partially.** Validation A has RMSE `{a['rmse_k']:.6f} K` and R2 `{a['r2']:.6f}`, but its near-zero overall bias is cancellation between CC underprediction and CV overprediction. It overpredicts the final temperature by `{a['predicted_final_temperature_c']-a['measured_final_temperature_c']:+.6f} K` and the temperature rise by `{a['predicted_temperature_rise_k']-a['measured_temperature_rise_k']:+.6f} K`.",
            f"2. **Is failure primarily thermal, electrical, or mixed?** **Mixed overall, with electrical-voltage heat dominant in Validation B.** Validation A already shows charge/CV thermal-transfer limitation. Replacing measured voltage with model voltage changes RMSE to `{isolation['rmse_k']:.6f} K`; adding predicted-temperature feedback changes it only to `{b['rmse_k']:.6f} K`. Thus nearly all additional coupled failure precedes thermal feedback and follows the `{difference['integrated_difference_j']:+.3f} J` electrical heat deficit.",
            f"3. **Does error differ materially between CC and CV?** **Yes.** Validation A changes from CC RMSE/bias `{phases['validation_a']['cc']['rmse_k']:.6f}/{phases['validation_a']['cc']['mean_signed_error_k']:+.6f} K` to CV `{phases['validation_a']['cv']['rmse_k']:.6f}/{phases['validation_a']['cv']['mean_signed_error_k']:+.6f} K`, including a bias-sign reversal. Validation B remains negatively biased in both phases.",
            f"4. **Does reversible heat materially affect prediction?** **Yes, but it is not the dominant total heat source.** It contributes `{100.0*float(measured_heat['signed_reversible_total_energy_fraction']):.3f}%` of signed total energy and shifts the final Validation-A prediction by `{reversible_effect['final_temperature_delta_k']:+.6f} K` (maximum trajectory effect `{reversible_effect['max_abs_temperature_delta_k']:.6f} K`).",
            "5. **Does this support transfer across current direction?** **It limits that claim.** The frozen pair retains useful CC shape information with measured-voltage heat, but the CV residual structure and the failed fully predictive path do not support general charge-direction transfer.",
            "",
            "No thermal/electrical parameter, OCP, entropy coefficient, SOC mapping, or frozen model equation was modified.",
        ]
    )
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def save_plots(
    data: RetainedData,
    prediction_a: np.ndarray,
    prediction_b: np.ndarray,
    q_irr_measured: np.ndarray,
    q_irr_model: np.ndarray,
    q_rev_measured: np.ndarray,
    q_total_measured: np.ndarray,
    q_total_model: np.ndarray,
) -> None:
    minutes = data.time_s / 60.0
    cv_minutes = minutes[data.cv_start_retained_index]

    fig, axis = plt.subplots(figsize=(11.2, 6.4), constrained_layout=True)
    axis.plot(minutes, data.measured_temperature_k - 273.15, color="black", marker="o", ms=2.5, lw=1.1, label="Measured battery")
    axis.plot(minutes, prediction_a - 273.15, color="#0072B2", lw=1.4, label="Validation A")
    axis.plot(minutes, prediction_b - 273.15, color="#D55E00", lw=1.3, label="Validation B")
    axis.plot(minutes, data.ambient_temperature_k - 273.15, color="#009E73", lw=1.0, label="Chamber")
    axis.axvline(cv_minutes, color="#CC79A7", ls="--", label="CC to CV")
    axis.set(xlabel="Elapsed time from SOC=0.10 [min]", ylabel="Temperature [degC]", title="Panasonic 25 C CC-CV thermal validation (SOC 0.10-0.90, no refit)")
    axis.grid(alpha=0.25)
    axis.legend()
    fig.savefig(TEMPERATURE_PLOT_PATH, dpi=180)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(11.2, 5.6), constrained_layout=True)
    axis.plot(minutes, prediction_a - data.measured_temperature_k, color="#0072B2", marker="o", ms=2.2, lw=1.0, label="Validation A residual")
    axis.plot(minutes, prediction_b - data.measured_temperature_k, color="#D55E00", marker="o", ms=2.2, lw=1.0, label="Validation B residual")
    axis.axhline(0.0, color="black", lw=0.8)
    axis.axvline(cv_minutes, color="#CC79A7", ls="--", label="CC to CV")
    axis.set(xlabel="Elapsed time from SOC=0.10 [min]", ylabel="Predicted - measured [K]", title="CC-CV thermal residuals")
    axis.grid(alpha=0.25)
    axis.legend()
    fig.savefig(RESIDUAL_PLOT_PATH, dpi=180)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(11.2, 6.2), constrained_layout=True)
    axis.plot(minutes, q_irr_measured, color="#0072B2", lw=1.2, label="Qirr measured voltage")
    axis.plot(minutes, q_irr_model, color="#56B4E9", lw=1.1, ls="--", label="Qirr model voltage")
    axis.plot(minutes, q_rev_measured, color="#009E73", lw=1.1, label="Qrev measured T")
    axis.plot(minutes, q_total_measured, color="black", lw=1.3, label="Total measured-voltage heat")
    axis.plot(minutes, q_total_model, color="#D55E00", lw=1.2, label="Total model-voltage heat")
    axis.axvline(cv_minutes, color="#CC79A7", ls="--", label="CC to CV")
    axis.axhline(0.0, color="gray", lw=0.7)
    axis.set(xlabel="Elapsed time from SOC=0.10 [min]", ylabel="Heat rate [W]", title="CC-CV charge heat components")
    axis.grid(alpha=0.25)
    axis.legend(ncol=2)
    fig.savefig(HEAT_PLOT_PATH, dpi=180)
    plt.close(fig)


def main() -> None:
    data = prepare_data()
    ueq_average = np.asarray(U_p(data.x_p_avg)) - np.asarray(U_n(data.x_n_avg))
    dudt = np.asarray(dUeq_dT_panasonic(data.soc))

    q_irr_measured = data.model_current_a * (ueq_average - data.measured_voltage_v)
    q_rev_measured = -data.model_current_a * data.measured_temperature_k * dudt
    q_total_measured = q_irr_measured + q_rev_measured
    prediction_a = exact_one_node(
        data.time_s,
        q_total_measured,
        data.ambient_temperature_k,
        float(data.measured_temperature_k[0]),
        C_TH_J_PER_K,
        G_TH_W_PER_K,
    )
    prediction_a_irreversible_only = exact_one_node(
        data.time_s,
        q_irr_measured,
        data.ambient_temperature_k,
        float(data.measured_temperature_k[0]),
        C_TH_J_PER_K,
        G_TH_W_PER_K,
    )

    coupled = validation_b(data, ueq_average, dudt)
    measured_t_model_voltage = frozen_voltage(
        data.model_current_a,
        data.measured_temperature_k,
        data.x_n_surface,
        data.x_p_surface,
    )
    if np.max(np.abs(measured_t_model_voltage - data.measured_t_model_voltage_v)) > 2.0e-12:
        raise RuntimeError("Frozen electrical voltage reproduction failed.")
    q_irr_model_measured_t = data.model_current_a * (ueq_average - measured_t_model_voltage)
    q_total_isolation = q_irr_model_measured_t + q_rev_measured
    prediction_isolation = exact_one_node(
        data.time_s,
        q_total_isolation,
        data.ambient_temperature_k,
        float(data.measured_temperature_k[0]),
        C_TH_J_PER_K,
        G_TH_W_PER_K,
    )

    metrics = {
        "validation_a": temperature_metrics(data.measured_temperature_k, prediction_a),
        "validation_b": temperature_metrics(data.measured_temperature_k, coupled.temperature_k),
        "model_voltage_measured_t": temperature_metrics(data.measured_temperature_k, prediction_isolation),
    }
    phases = {
        name: {
            "cc": subset_metrics(data, prediction, data.cc_mask),
            "cv": subset_metrics(data, prediction, data.cv_mask),
        }
        for name, prediction in (
            ("validation_a", prediction_a),
            ("validation_b", coupled.temperature_k),
            ("model_voltage_measured_t", prediction_isolation),
        )
    }
    heats = {
        "measured_voltage": heat_audit(data, q_irr_measured, q_rev_measured),
        "fully_predictive": heat_audit(data, coupled.q_irr_w, coupled.q_rev_w),
        "measured_voltage_reversible_temperature_effect": {
            "final_temperature_delta_k": float(
                prediction_a[-1] - prediction_a_irreversible_only[-1]
            ),
            "max_abs_temperature_delta_k": float(
                np.max(np.abs(prediction_a - prediction_a_irreversible_only))
            ),
            "mean_temperature_delta_k": float(
                np.mean(prediction_a - prediction_a_irreversible_only)
            ),
            "rms_temperature_delta_k": float(
                np.sqrt(np.mean((prediction_a - prediction_a_irreversible_only) ** 2))
            ),
        },
    }
    difference = electrical_heat_difference(
        data, q_irr_measured, coupled.q_irr_w, coupled.voltage_v
    )
    audit: dict[str, object] = {
        "source_path": str(RAW_SOURCE_PATH),
        "source_start_time_s": float(data.source_time_s[0]),
        "source_end_time_s": float(data.source_time_s[-1]),
        "duration_s": float(data.time_s[-1]),
        "sample_count": int(data.time_s.size),
        "soc_start": float(data.soc[0]),
        "soc_end": float(data.soc[-1]),
        "current_min_a": float(np.min(data.model_current_a)),
        "current_max_a": float(np.max(data.model_current_a)),
        "cc_current_max_a": float(np.max(data.model_current_a[data.cc_mask])),
        "voltage_min_v": float(np.min(data.measured_voltage_v)),
        "voltage_max_v": float(np.max(data.measured_voltage_v)),
        "battery_temperature_min_c": float(np.min(data.measured_temperature_k) - 273.15),
        "battery_temperature_max_c": float(np.max(data.measured_temperature_k) - 273.15),
        "chamber_temperature_min_c": float(np.min(data.ambient_temperature_k) - 273.15),
        "chamber_temperature_max_c": float(np.max(data.ambient_temperature_k) - 273.15),
        "charge_throughput_ah": float(data.exported_ah[-1] - data.exported_ah[0]),
        "full_charge_throughput_ah": data.full_charge_throughput_ah,
        "cv_start_source_time_s": data.source_cv_start_time_s,
        "cv_start_soc": data.source_cv_start_soc,
        "cc_to_cv_inside_supported_interval": True,
        "soc_0_90_occurs_during_cv": bool(data.cv_mask[-1]),
        "cc_model_current_strictly_negative": bool(np.all(data.model_current_a[data.cc_mask] < 0.0)),
        "entropy_extrapolation_count": 0,
        "soc_clipping_count": 0,
    }

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    write_metrics(data, audit, metrics, phases, heats, difference)
    write_report(data, audit, metrics, phases, heats, difference)
    save_plots(
        data,
        prediction_a,
        coupled.temperature_k,
        q_irr_measured,
        coupled.q_irr_w,
        q_rev_measured,
        q_total_measured,
        coupled.q_total_w,
    )
    print(f"Validation A: RMSE={metrics['validation_a']['rmse_k']:.9f} K, bias={metrics['validation_a']['mean_signed_error_k']:+.9f} K, R2={metrics['validation_a']['r2']:.9f}")
    print(f"Validation B: RMSE={metrics['validation_b']['rmse_k']:.9f} K, bias={metrics['validation_b']['mean_signed_error_k']:+.9f} K, R2={metrics['validation_b']['r2']:.9f}")
    print(f"Retained samples={data.time_s.size}; CC-to-CV SOC={data.source_cv_start_soc:.12f}")


if __name__ == "__main__":
    main()
