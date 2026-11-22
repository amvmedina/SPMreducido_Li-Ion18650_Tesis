"""Forward-only Panasonic NCR18650PF UDDS thermal validation.

The script evaluates the frozen 1C-derived lumped thermal parameters on the
independent 25 C UDDS experiment.  It performs no optimization, clipping,
entropy extrapolation, parameter release, or model modification.
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
from scipy.stats import spearmanr


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from models.grouped_spm import (  # noqa: E402
    FARADAY_C_PER_MOL,
    GAS_CONSTANT_J_PER_MOL_K,
)
from models.ocp_panasonic_qi import (  # noqa: E402
    MCMB_SUPPORTED_STOICHIOMETRY_RANGE,
    NCA_WORKING_STOICHIOMETRY_RANGE,
)
from models.ocp_panasonic_qi_frozen import (  # noqa: E402
    Q_CELL_AH,
    U_n,
    U_p,
    X_N_START,
    X_P_START,
)
from models.panasonic_entropic import (  # noqa: E402
    SUPPORTED_SOC_RANGE,
    dUeq_dT_panasonic,
)
from models.panasonic_thermal_parameters import (  # noqa: E402
    C_TH_J_PER_K,
    G_TH_W_PER_K,
    TAU_TH_S,
)
from scripts.audit_panasonic_hppc_identifiability import (  # noqa: E402
    Q_N_C as IMPLEMENTED_Q_N_C,
    Q_P_C as IMPLEMENTED_Q_P_C,
    simulate_local_group,
)
from scripts.validate_panasonic_udds_forward import (  # noqa: E402
    MANIFEST_PATH,
    NOMINAL_PARAMETERS,
    Q_N_C,
    Q_P_C,
    RAW_SOURCE_PATH,
    load_udds,
    remove_exact_duplicate_times,
)


OUTPUT_DIRECTORY = REPOSITORY_ROOT / "thesis_data"
REPORT_PATH = OUTPUT_DIRECTORY / "panasonic_udds_thermal_validation.md"
METRICS_PATH = OUTPUT_DIRECTORY / "panasonic_udds_thermal_validation_metrics.csv"
TEMPERATURE_PLOT_PATH = OUTPUT_DIRECTORY / "panasonic_udds_thermal_validation.png"
RESIDUAL_TIME_PATH = OUTPUT_DIRECTORY / "panasonic_udds_thermal_residual_time.png"
RESIDUAL_SOC_PATH = OUTPUT_DIRECTORY / "panasonic_udds_thermal_residual_soc.png"
HEAT_PLOT_PATH = OUTPUT_DIRECTORY / "panasonic_udds_heat_components.png"

NEAR_REST_CURRENT_A = 0.05
N_LARGEST_CURRENT_STEPS = 10


@dataclass(frozen=True)
class ValidationData:
    time_s: np.ndarray
    source_time_s: np.ndarray
    soc: np.ndarray
    current_a: np.ndarray
    measured_voltage_v: np.ndarray
    measured_temperature_k: np.ndarray
    ambient_temperature_k: np.ndarray
    x_n_avg: np.ndarray
    x_n_surface: np.ndarray
    x_p_avg: np.ndarray
    x_p_surface: np.ndarray
    ueq_average_v: np.ndarray
    dudt_v_per_k: np.ndarray
    model_voltage_measured_t_v: np.ndarray


def exact_thermal_trajectory(
    time_s: np.ndarray,
    heat_w: np.ndarray,
    ambient_temperature_k: np.ndarray,
    initial_temperature_k: float,
) -> np.ndarray:
    """Exact left-sample-ZOH update of the frozen lumped thermal equation."""

    if not (time_s.shape == heat_w.shape == ambient_temperature_k.shape):
        raise ValueError("Thermal inputs have different shapes.")
    if np.any(np.diff(time_s) <= 0.0):
        raise ValueError("Thermal timestamps must be strictly increasing.")
    temperature = np.empty_like(time_s)
    temperature[0] = initial_temperature_k
    for index, dt_s in enumerate(np.diff(time_s)):
        exponent = G_TH_W_PER_K * dt_s / C_TH_J_PER_K
        approach = -np.expm1(-exponent)
        phi = approach / exponent if exponent != 0.0 else 1.0
        temperature[index + 1] = (
            temperature[index]
            + approach * (ambient_temperature_k[index] - temperature[index])
            + heat_w[index] * dt_s * phi / C_TH_J_PER_K
        )
    return temperature


def _crossing_left_index(soc: np.ndarray, target: float, start: int = 0) -> int:
    candidates = np.flatnonzero(
        (np.arange(soc.size - 1) >= start)
        & (soc[:-1] > target)
        & (soc[1:] <= target)
    )
    if candidates.size != 1:
        raise ValueError(
            f"Expected one downward SOC={target:.2f} crossing after index {start}; "
            f"found {candidates.size}."
        )
    return int(candidates[0])


def _interpolate_at_soc(
    values: np.ndarray, soc: np.ndarray, left_index: int, target: float
) -> float:
    denominator = soc[left_index + 1] - soc[left_index]
    if denominator == 0.0:
        raise ValueError("Cannot interpolate an SOC boundary across equal SOC values.")
    fraction = (target - soc[left_index]) / denominator
    if not 0.0 <= fraction <= 1.0:
        raise ValueError("SOC boundary interpolation fraction is outside [0, 1].")
    return float(values[left_index] + fraction * (values[left_index + 1] - values[left_index]))


def restrict_to_continuous_entropy_interval(
    soc: np.ndarray, **signals: np.ndarray
) -> tuple[np.ndarray, dict[str, np.ndarray], tuple[int, int]]:
    """Return the unique chronological 0.90-to-0.10 interval.

    Boundary records are linearly interpolated. All original records strictly
    between the two crossings are retained, including regenerative excursions.
    """

    lower, upper = SUPPORTED_SOC_RANGE
    high_left = _crossing_left_index(soc, upper)
    low_left = _crossing_left_index(soc, lower, start=high_left + 1)
    interior = slice(high_left + 1, low_left + 1)
    selected_soc = np.concatenate(([upper], soc[interior], [lower]))
    if np.any((selected_soc < lower) | (selected_soc > upper)):
        raise ValueError("Continuous interval leaves entropy support; no clipping.")

    selected: dict[str, np.ndarray] = {}
    for name, raw in signals.items():
        values = np.asarray(raw, dtype=float)
        if values.shape != soc.shape:
            raise ValueError(f"Signal {name!r} does not match the SOC shape.")
        selected[name] = np.concatenate(
            (
                [_interpolate_at_soc(values, soc, high_left, upper)],
                values[interior],
                [_interpolate_at_soc(values, soc, low_left, lower)],
            )
        )
    if np.any(np.diff(selected["time_s"]) <= 0.0):
        raise ValueError("Restricted timestamps are not strictly increasing.")
    return selected_soc, selected, (high_left, low_left)


def terminal_voltage_at_temperature(
    current_a: np.ndarray,
    temperature_k: np.ndarray,
    x_n_surface: np.ndarray,
    x_p_surface: np.ndarray,
) -> np.ndarray:
    """Evaluate the frozen voltage equation at a prescribed temperature."""

    base_voltage, kinetic_per_k = terminal_voltage_components(
        current_a, x_n_surface, x_p_surface
    )
    voltage = base_voltage + temperature_k * kinetic_per_k
    if not np.all(np.isfinite(voltage)):
        raise ValueError("Non-finite frozen electrical voltage.")
    return voltage


def terminal_voltage_components(
    current_a: np.ndarray,
    x_n_surface: np.ndarray,
    x_p_surface: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return temperature-independent base voltage and kinetic V/K factor."""

    d_n = NOMINAL_PARAMETERS["d_n"]
    d_p = NOMINAL_PARAMETERS["d_p"]
    r0 = NOMINAL_PARAMETERS["R0"]
    eta_n_argument = current_a / (
        6.0 * Q_N_C * d_n * np.sqrt(x_n_surface * (1.0 - x_n_surface))
    )
    eta_p_argument = current_a / (
        6.0 * Q_P_C * d_p * np.sqrt(x_p_surface * (1.0 - x_p_surface))
    )
    kinetic_per_k = (
        2.0
        * GAS_CONSTANT_J_PER_MOL_K
        / FARADAY_C_PER_MOL
        * (np.arcsinh(eta_p_argument) - np.arcsinh(eta_n_argument))
    )
    base_voltage = (
        np.asarray(U_p(x_p_surface))
        - np.asarray(U_n(x_n_surface))
        - r0 * current_a
    )
    return base_voltage, kinetic_per_k


def prepare_data() -> tuple[ValidationData, dict[str, object]]:
    columns, raw = load_udds(RAW_SOURCE_PATH)
    raw_samples = raw["Time"].size
    cleaned, duplicate_indices = remove_exact_duplicate_times(raw)
    time_s = cleaned["Time"]
    source_soc = 1.0 + cleaned["Ah"] / Q_CELL_AH
    if np.any((source_soc < 0.0) | (source_soc > 1.0)):
        raise ValueError("Raw-Ah SOC leaves [0, 1]; no clipping.")

    current_a = -cleaned["Current"]  # sole raw-to-model sign conversion
    measured_temperature_k = cleaned["Battery_Temp_degC"] + 273.15
    electrical = simulate_local_group(
        time_s=time_s,
        discharge_positive_current_a=current_a,
        temperature_k=measured_temperature_k,
        initial_x_n=X_N_START,
        initial_x_p=X_P_START,
        dynamic_parameters=NOMINAL_PARAMETERS,
    )
    reproduced_voltage = terminal_voltage_at_temperature(
        current_a,
        measured_temperature_k,
        electrical.x_n_surface,
        electrical.x_p_surface,
    )
    voltage_reproduction_error = float(
        np.max(np.abs(reproduced_voltage - electrical.voltage_v))
    )
    if voltage_reproduction_error > 5.0e-13:
        raise RuntimeError("The locally coupled voltage equation changed.")

    soc, selected, crossings = restrict_to_continuous_entropy_interval(
        source_soc,
        time_s=time_s,
        current_a=current_a,
        measured_voltage_v=cleaned["Voltage"],
        measured_temperature_k=measured_temperature_k,
        ambient_temperature_k=cleaned["Chamber_Temp_degC"] + 273.15,
        x_n_avg=electrical.x_n_avg,
        x_n_surface=electrical.x_n_surface,
        x_p_avg=electrical.x_p_avg,
        x_p_surface=electrical.x_p_surface,
    )
    source_time_s = selected["time_s"].copy()
    selected["time_s"] -= selected["time_s"][0]

    xn_lo, xn_hi = MCMB_SUPPORTED_STOICHIOMETRY_RANGE
    xp_lo, xp_hi = NCA_WORKING_STOICHIOMETRY_RANGE
    domain_violation = (
        (selected["x_n_avg"] < xn_lo)
        | (selected["x_n_avg"] > xn_hi)
        | (selected["x_n_surface"] < xn_lo)
        | (selected["x_n_surface"] > xn_hi)
        | (selected["x_p_avg"] < xp_lo)
        | (selected["x_p_avg"] > xp_hi)
        | (selected["x_p_surface"] < xp_lo)
        | (selected["x_p_surface"] > xp_hi)
    )
    if np.any(domain_violation):
        raise ValueError("Frozen average/surface state left OCP support; no clipping.")

    ueq_average = np.asarray(U_p(selected["x_p_avg"])) - np.asarray(
        U_n(selected["x_n_avg"])
    )
    dudt = np.asarray(dUeq_dT_panasonic(soc))
    model_voltage_measured_t = terminal_voltage_at_temperature(
        selected["current_a"],
        selected["measured_temperature_k"],
        selected["x_n_surface"],
        selected["x_p_surface"],
    )
    data = ValidationData(
        time_s=selected["time_s"],
        source_time_s=source_time_s,
        soc=soc,
        current_a=selected["current_a"],
        measured_voltage_v=selected["measured_voltage_v"],
        measured_temperature_k=selected["measured_temperature_k"],
        ambient_temperature_k=selected["ambient_temperature_k"],
        x_n_avg=selected["x_n_avg"],
        x_n_surface=selected["x_n_surface"],
        x_p_avg=selected["x_p_avg"],
        x_p_surface=selected["x_p_surface"],
        ueq_average_v=ueq_average,
        dudt_v_per_k=dudt,
        model_voltage_measured_t_v=model_voltage_measured_t,
    )
    audit = {
        "columns": columns,
        "raw_samples": raw_samples,
        "duplicate_rows_removed": int(duplicate_indices.size),
        "high_crossing_left_index": crossings[0],
        "low_crossing_left_index": crossings[1],
        "ocp_domain_violations": int(np.count_nonzero(domain_violation)),
        "clipping_count": 0,
        "electrical_voltage_reproduction_max_abs_v": voltage_reproduction_error,
    }
    return data, audit


def fully_predictive_trajectory(
    data: ValidationData,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Recursively couple frozen voltage, heat, and thermal state."""

    count = data.time_s.size
    temperature = np.empty(count)
    voltage = np.empty(count)
    q_irr = np.empty(count)
    q_rev = np.empty(count)
    q_total = np.empty(count)
    temperature[0] = data.measured_temperature_k[0]
    base_voltage, kinetic_per_k = terminal_voltage_components(
        data.current_a, data.x_n_surface, data.x_p_surface
    )

    for index in range(count):
        voltage[index] = base_voltage[index] + temperature[index] * kinetic_per_k[index]
        q_irr[index] = data.current_a[index] * (
            data.ueq_average_v[index] - voltage[index]
        )
        q_rev[index] = (
            -data.current_a[index] * temperature[index] * data.dudt_v_per_k[index]
        )
        q_total[index] = q_irr[index] + q_rev[index]
        if index == count - 1:
            continue
        dt_s = data.time_s[index + 1] - data.time_s[index]
        exponent = G_TH_W_PER_K * dt_s / C_TH_J_PER_K
        approach = -np.expm1(-exponent)
        phi = approach / exponent if exponent != 0.0 else 1.0
        temperature[index + 1] = (
            temperature[index]
            + approach * (data.ambient_temperature_k[index] - temperature[index])
            + q_total[index] * dt_s * phi / C_TH_J_PER_K
        )
    return temperature, voltage, q_irr, q_rev, q_total


def temperature_metrics(measured: np.ndarray, predicted: np.ndarray) -> dict[str, float]:
    residual = predicted - measured
    denominator = float(np.sum((measured - np.mean(measured)) ** 2))
    return {
        "rmse_k": float(np.sqrt(np.mean(residual**2))),
        "mae_k": float(np.mean(np.abs(residual))),
        "mean_signed_error_k": float(np.mean(residual)),
        "max_abs_error_k": float(np.max(np.abs(residual))),
        "r2": float(1.0 - np.sum(residual**2) / denominator),
        "measured_final_temperature_c": float(measured[-1] - 273.15),
        "predicted_final_temperature_c": float(predicted[-1] - 273.15),
        "measured_temperature_rise_k": float(measured[-1] - measured[0]),
        "predicted_temperature_rise_k": float(predicted[-1] - predicted[0]),
    }


def soc_masks(soc: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "high": soc >= 2.0 / 3.0,
        "middle": (soc >= 1.0 / 3.0) & (soc < 2.0 / 3.0),
        "low": soc < 1.0 / 3.0,
    }


def left_zoh_energy(time_s: np.ndarray, power_w: np.ndarray) -> float:
    return float(np.sum(power_w[:-1] * np.diff(time_s)))


def time_weighted_mean(time_s: np.ndarray, values: np.ndarray) -> float:
    return left_zoh_energy(time_s, values) / float(time_s[-1] - time_s[0])


def time_fraction(time_s: np.ndarray, condition: np.ndarray) -> float:
    return float(np.sum(np.diff(time_s)[condition[:-1]]) / (time_s[-1] - time_s[0]))


def heat_audit(
    time_s: np.ndarray, q_irr: np.ndarray, q_rev: np.ndarray, q_total: np.ndarray
) -> dict[str, float]:
    e_irr = left_zoh_energy(time_s, q_irr)
    e_rev = left_zoh_energy(time_s, q_rev)
    e_total = left_zoh_energy(time_s, q_total)
    e_rev_abs = left_zoh_energy(time_s, np.abs(q_rev))
    e_total_abs_components = left_zoh_energy(time_s, np.abs(q_irr) + np.abs(q_rev))
    return {
        "irreversible_heat_energy_j": e_irr,
        "reversible_heat_energy_j": e_rev,
        "total_heat_energy_j": e_total,
        "mean_qdot_irr_w": time_weighted_mean(time_s, q_irr),
        "mean_qdot_rev_w": time_weighted_mean(time_s, q_rev),
        "mean_qdot_total_w": time_weighted_mean(time_s, q_total),
        "signed_reversible_fraction": e_rev / e_total,
        "absolute_reversible_component_fraction": e_rev_abs / e_total_abs_components,
        "reversible_exothermic_time_fraction": time_fraction(time_s, q_rev > 0.0),
        "reversible_endothermic_time_fraction": time_fraction(time_s, q_rev < 0.0),
    }


def sign_context(
    data: ValidationData, q_rev_a: np.ndarray, q_rev_b: np.ndarray
) -> tuple[dict[str, float | int], list[dict[str, object]]]:
    expected_sign = np.sign(-data.current_a * data.dudt_v_per_k)
    active = expected_sign != 0.0
    mismatch_a = int(np.count_nonzero(np.sign(q_rev_a[active]) != expected_sign[active]))
    mismatch_b = int(np.count_nonzero(np.sign(q_rev_b[active]) != expected_sign[active]))
    summary: dict[str, float | int] = {
        "validation_a_reversible_sign_mismatches": mismatch_a,
        "validation_b_reversible_sign_mismatches": mismatch_b,
    }
    contexts: list[dict[str, object]] = []
    cases = (
        ("discharge_dudt_positive", (data.current_a > 0.0) & (data.dudt_v_per_k > 0.0), "endothermic"),
        ("discharge_dudt_negative", (data.current_a > 0.0) & (data.dudt_v_per_k < 0.0), "exothermic"),
        ("regeneration_dudt_positive", (data.current_a < 0.0) & (data.dudt_v_per_k > 0.0), "exothermic"),
        ("regeneration_dudt_negative", (data.current_a < 0.0) & (data.dudt_v_per_k < 0.0), "endothermic"),
    )
    for name, mask, expected in cases:
        contexts.append(
            {
                "case": name,
                "expected": expected,
                "samples": int(np.count_nonzero(mask)),
                "time_fraction": time_fraction(data.time_s, mask),
                "mean_qrev_a_w": float(np.mean(q_rev_a[mask])) if np.any(mask) else float("nan"),
                "mean_qrev_b_w": float(np.mean(q_rev_b[mask])) if np.any(mask) else float("nan"),
            }
        )
    return summary, contexts


def dynamic_diagnostics(
    data: ValidationData,
    predicted_a: np.ndarray,
    predicted_b: np.ndarray,
) -> tuple[dict[str, float], list[dict[str, float | int]]]:
    residual_a = predicted_a - data.measured_temperature_k
    residual_b = predicted_b - data.measured_temperature_k
    current_magnitude = np.abs(data.current_a)
    loaded = current_magnitude > NEAR_REST_CURRENT_A
    regenerative = data.current_a < -NEAR_REST_CURRENT_A
    discharge = data.current_a > NEAR_REST_CURRENT_A
    diagnostics = {
        "validation_a_abs_error_vs_soc_depletion_spearman": float(
            spearmanr(np.abs(residual_a), 1.0 - data.soc).statistic
        ),
        "validation_b_abs_error_vs_soc_depletion_spearman": float(
            spearmanr(np.abs(residual_b), 1.0 - data.soc).statistic
        ),
        "validation_a_abs_error_vs_current_magnitude_loaded_spearman": float(
            spearmanr(np.abs(residual_a[loaded]), current_magnitude[loaded]).statistic
        ),
        "validation_b_abs_error_vs_current_magnitude_loaded_spearman": float(
            spearmanr(np.abs(residual_b[loaded]), current_magnitude[loaded]).statistic
        ),
        "validation_a_regenerative_rmse_k": float(np.sqrt(np.mean(residual_a[regenerative] ** 2))),
        "validation_b_regenerative_rmse_k": float(np.sqrt(np.mean(residual_b[regenerative] ** 2))),
        "validation_a_regenerative_bias_k": float(np.mean(residual_a[regenerative])),
        "validation_b_regenerative_bias_k": float(np.mean(residual_b[regenerative])),
        "validation_a_discharge_rmse_k": float(np.sqrt(np.mean(residual_a[discharge] ** 2))),
        "validation_b_discharge_rmse_k": float(np.sqrt(np.mean(residual_b[discharge] ** 2))),
    }

    # Exclude the two synthetic boundary-to-source transitions from event ranking.
    step_indices = np.arange(2, data.current_a.size - 1)
    changes = np.abs(data.current_a[step_indices] - data.current_a[step_indices - 1])
    ranked = step_indices[np.argsort(changes)[-N_LARGEST_CURRENT_STEPS:][::-1]]
    events: list[dict[str, float | int]] = []
    for rank, index in enumerate(ranked, start=1):
        events.append(
            {
                "rank": rank,
                "index": int(index),
                "source_time_s": float(data.source_time_s[index]),
                "elapsed_time_s": float(data.time_s[index]),
                "soc": float(data.soc[index]),
                "current_before_a": float(data.current_a[index - 1]),
                "current_after_a": float(data.current_a[index]),
                "delta_current_a": float(data.current_a[index] - data.current_a[index - 1]),
                "measured_temperature_c": float(data.measured_temperature_k[index] - 273.15),
                "validation_a_predicted_temperature_c": float(predicted_a[index] - 273.15),
                "validation_a_error_k": float(residual_a[index]),
                "validation_b_predicted_temperature_c": float(predicted_b[index] - 273.15),
                "validation_b_error_k": float(residual_b[index]),
            }
        )
    diagnostics["validation_a_large_step_rmse_k"] = float(
        np.sqrt(np.mean(np.asarray([event["validation_a_error_k"] for event in events]) ** 2))
    )
    diagnostics["validation_b_large_step_rmse_k"] = float(
        np.sqrt(np.mean(np.asarray([event["validation_b_error_k"] for event in events]) ** 2))
    )
    return diagnostics, events


def classify(metrics: dict[str, float], measured: np.ndarray, coupled: bool) -> str:
    """Evidence-scale classification; thresholds are relative to measured excursion."""

    excursion = float(np.max(measured) - np.min(measured))
    relative_rmse = metrics["rmse_k"] / excursion
    relative_bias = abs(metrics["mean_signed_error_k"]) / excursion
    if relative_rmse <= 0.15 and relative_bias <= 0.10 and metrics["r2"] >= 0.80:
        return "COUPLED ELECTROTHERMAL TRANSFER GOOD" if coupled else "THERMAL PARAMETERS TRANSFER WELL"
    if relative_rmse <= 0.50 and metrics["r2"] > 0.0:
        return "COUPLED ELECTROTHERMAL TRANSFER LIMITED" if coupled else "THERMAL TRANSFER LIMITED"
    return "COUPLED ELECTROTHERMAL TRANSFER FAILED" if coupled else "THERMAL TRANSFER FAILED"


def save_temperature_plot(
    data: ValidationData, predicted_a: np.ndarray, predicted_b: np.ndarray, isolation: np.ndarray
) -> None:
    shown = np.arange(0, data.time_s.size, max(1, data.time_s.size // 45000))
    hours = data.time_s / 3600.0
    fig, ax = plt.subplots(figsize=(12.0, 5.7), constrained_layout=True)
    ax.plot(hours[shown], data.measured_temperature_k[shown] - 273.15, color="black", lw=1.2, label="Measured battery")
    ax.plot(hours[shown], predicted_a[shown] - 273.15, color="#0072B2", lw=1.0, label="Validation A: measured-voltage heat")
    ax.plot(hours[shown], predicted_b[shown] - 273.15, color="#D55E00", lw=1.0, label="Validation B: fully predictive")
    ax.plot(hours[shown], isolation[shown] - 273.15, color="#CC79A7", lw=0.8, ls="--", label="Model-voltage / measured-T diagnostic")
    ax.plot(hours[shown], data.ambient_temperature_k[shown] - 273.15, color="#009E73", lw=0.9, label="Measured chamber")
    ax.set(xlabel="Elapsed time in entropy-supported interval [h]", ylabel="Temperature [degC]", title="Panasonic 25 C UDDS thermal validation (no refit)")
    ax.grid(alpha=0.22)
    ax.legend(loc="best", ncol=2)
    fig.savefig(TEMPERATURE_PLOT_PATH, dpi=180)
    plt.close(fig)


def save_residual_plots(
    data: ValidationData, residual_a: np.ndarray, residual_b: np.ndarray, residual_isolation: np.ndarray
) -> None:
    shown = np.arange(0, data.time_s.size, max(1, data.time_s.size // 45000))
    fig, ax = plt.subplots(figsize=(12.0, 5.0), constrained_layout=True)
    ax.plot(data.time_s[shown] / 3600.0, residual_a[shown], lw=0.9, color="#0072B2", label="Validation A")
    ax.plot(data.time_s[shown] / 3600.0, residual_b[shown], lw=0.9, color="#D55E00", label="Validation B")
    ax.plot(data.time_s[shown] / 3600.0, residual_isolation[shown], lw=0.7, ls="--", color="#CC79A7", label="Isolation diagnostic")
    ax.axhline(0.0, color="black", lw=0.7)
    ax.set(xlabel="Elapsed time [h]", ylabel="Predicted - measured [K]", title="UDDS thermal residual versus time")
    ax.grid(alpha=0.22)
    ax.legend()
    fig.savefig(RESIDUAL_TIME_PATH, dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.0, 5.2), constrained_layout=True)
    ax.scatter(100.0 * data.soc[shown], residual_a[shown], s=2, alpha=0.15, color="#0072B2", rasterized=True, label="Validation A")
    ax.scatter(100.0 * data.soc[shown], residual_b[shown], s=2, alpha=0.15, color="#D55E00", rasterized=True, label="Validation B")
    ax.axhline(0.0, color="black", lw=0.7)
    ax.axvline(100.0 / 3.0, color="0.45", ls="--", lw=0.8)
    ax.axvline(200.0 / 3.0, color="0.45", ls="--", lw=0.8)
    ax.invert_xaxis()
    ax.set(xlabel="Cell SOC [%]", ylabel="Predicted - measured [K]", title="UDDS thermal residual versus SOC")
    ax.grid(alpha=0.22)
    ax.legend()
    fig.savefig(RESIDUAL_SOC_PATH, dpi=180)
    plt.close(fig)


def save_heat_plot(
    data: ValidationData,
    heat_a: tuple[np.ndarray, np.ndarray, np.ndarray],
    heat_b: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> None:
    shown = np.arange(0, data.time_s.size, max(1, data.time_s.size // 45000))
    hours = data.time_s / 3600.0
    fig, axes = plt.subplots(2, 1, figsize=(12.0, 7.2), sharex=True, constrained_layout=True)
    for ax, title, heat in zip(axes, ("Validation A: measured-voltage / measured-T heat", "Validation B: model-voltage / predicted-T heat"), (heat_a, heat_b), strict=True):
        ax.plot(hours[shown], heat[0][shown], lw=0.65, color="#D55E00", label="Irreversible")
        ax.plot(hours[shown], heat[1][shown], lw=0.65, color="#0072B2", label="Reversible")
        ax.plot(hours[shown], heat[2][shown], lw=0.8, color="black", label="Total")
        ax.axhline(0.0, color="0.35", lw=0.6)
        ax.set(ylabel="Heat rate [W]", title=title)
        ax.grid(alpha=0.2)
        ax.legend(loc="upper right", ncol=3)
    axes[-1].set_xlabel("Elapsed time [h]")
    fig.savefig(HEAT_PLOT_PATH, dpi=180)
    plt.close(fig)


def add_rows(
    rows: list[dict[str, object]], validation: str, subset: str, values: dict[str, object], n: int
) -> None:
    units = {
        "rmse_k": "K", "mae_k": "K", "mean_signed_error_k": "K", "max_abs_error_k": "K",
        "r2": "1", "measured_final_temperature_c": "degC", "predicted_final_temperature_c": "degC",
        "measured_temperature_rise_k": "K", "predicted_temperature_rise_k": "K",
    }
    for metric, value in values.items():
        unit = units.get(metric, "1")
        if metric.endswith("_j"):
            unit = "J"
        elif metric.endswith("_w"):
            unit = "W"
        elif metric.endswith("_k"):
            unit = "K"
        elif metric.endswith("_s"):
            unit = "s"
        elif metric.endswith("_a"):
            unit = "A"
        elif metric.endswith("_v"):
            unit = "V"
        elif metric.endswith("_ah"):
            unit = "Ah"
        elif metric.endswith("_temperature_c"):
            unit = "degC"
        rows.append({"validation": validation, "subset": subset, "metric": metric, "value": value, "unit": unit, "n_samples": n})


def write_metrics(
    data: ValidationData,
    metrics_all: dict[str, dict[str, float]],
    regional: dict[str, dict[str, dict[str, float]]],
    heat: dict[str, dict[str, float]],
    diagnostics: dict[str, float],
    events: list[dict[str, float | int]],
    audit: dict[str, object],
) -> None:
    rows: list[dict[str, object]] = []
    for validation, values in metrics_all.items():
        add_rows(rows, validation, "overall", values, data.time_s.size)
    masks = soc_masks(data.soc)
    for validation, groups in regional.items():
        for region, values in groups.items():
            add_rows(rows, validation, f"soc_{region}", values, int(np.count_nonzero(masks[region])))
    for validation, values in heat.items():
        add_rows(rows, validation, "heat_components", values, data.time_s.size)
    add_rows(rows, "dynamic_diagnostic", "overall", diagnostics, data.time_s.size)
    scalar_audit = {key: value for key, value in audit.items() if key != "columns"}
    add_rows(rows, "dataset_audit", "supported_interval", scalar_audit, data.time_s.size)
    for event in events:
        add_rows(rows, "large_current_step", f"rank_{event['rank']}", event, 1)
    with METRICS_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _metric_table(values: dict[str, float]) -> str:
    return (
        f"| {values['rmse_k']:.9f} | {values['mae_k']:.9f} | "
        f"{values['mean_signed_error_k']:+.9f} | {values['max_abs_error_k']:.9f} | "
        f"{values['r2']:.9f} | {values['measured_final_temperature_c']:.9f} | "
        f"{values['predicted_final_temperature_c']:.9f} | "
        f"{values['measured_temperature_rise_k']:+.9f} | "
        f"{values['predicted_temperature_rise_k']:+.9f} |"
    )


def write_report(
    data: ValidationData,
    metrics_all: dict[str, dict[str, float]],
    regional: dict[str, dict[str, dict[str, float]]],
    heat: dict[str, dict[str, float]],
    diagnostics: dict[str, float],
    events: list[dict[str, float | int]],
    sign_summary: dict[str, float | int],
    sign_cases: list[dict[str, object]],
    audit: dict[str, object],
    classification_a: str,
    classification_b: str,
) -> None:
    a = metrics_all["validation_a"]
    b = metrics_all["validation_b"]
    iso = metrics_all["isolation_diagnostic"]
    lines = [
        "# Panasonic 25 C UDDS thermal validation",
        "",
        f"## {classification_a}",
        "",
        f"## {classification_b}",
        "",
        "This is an independent forward-only validation of the frozen 1C-derived effective lumped thermal parameters. No electrical, OCP, entropy, or thermal parameter was fitted or changed; there is no clipping, extrapolation, baseline correction, or temperature-offset fit.",
        "",
        "## Dataset and supported interval",
        "",
        f"- Exact source: `{RAW_SOURCE_PATH}`",
        f"- Manifest: `{MANIFEST_PATH}`",
        f"- Columns: `{', '.join(audit['columns'])}`",
        f"- Raw samples: {audit['raw_samples']:,}; one exact duplicate final record was removed before interval construction.",
        f"- Entropy-supported interval: exactly SOC {audit['start_soc']:.12f} to {audit['end_soc']:.12f}, created with linear interpolation at both boundaries. Source times are {audit['source_start_time_s']:.9f} to {audit['source_end_time_s']:.9f} s.",
        f"- Duration: {audit['duration_s']:.9f} s ({audit['duration_s']/3600.0:.6f} h); samples including the two exact boundary records: {audit['interval_samples']:,}.",
        f"- Actual positive sample interval: median {audit['dt_median_s']:.9f} s, mean {audit['dt_mean_s']:.9f} s, range {audit['dt_min_s']:.9f}-{audit['dt_max_s']:.9f} s.",
        f"- Internal current (discharge positive, regeneration negative): {audit['current_min_a']:.6f} to {audit['current_max_a']:.6f} A. Raw current was sign-converted once at the model boundary and was otherwise unchanged.",
        f"- Measured voltage: {audit['voltage_min_v']:.6f} to {audit['voltage_max_v']:.6f} V.",
        f"- Battery temperature: {audit['battery_temperature_min_c']:.6f} to {audit['battery_temperature_max_c']:.6f} degC; chamber temperature: {audit['chamber_temperature_min_c']:.6f} to {audit['chamber_temperature_max_c']:.6f} degC.",
        f"- Gross discharge throughput: {audit['gross_discharge_throughput_ah']:.9f} Ah; regenerative charge throughput: {audit['regenerative_charge_throughput_ah']:.9f} Ah; net discharge: {audit['net_discharge_throughput_ah']:.9f} Ah.",
        f"- OCP-domain violations: {audit['ocp_domain_violations']}; entropy-domain violations: {audit['entropy_domain_violations']}; clipping count: {audit['clipping_count']}.",
        "",
        "Actual nonuniform timestamps and all intervening regenerative excursions are preserved. Heat and ambient are held at the left sample for each exact analytical thermal-state update.",
        "",
        "## Frozen configuration",
        "",
        f"- `Cth = {C_TH_J_PER_K:.12f} J/K`, `Gth = {G_TH_W_PER_K:.12f} W/K`, `tau_th = {TAU_TH_S:.12f} s`.",
        f"- `Q_cell = {Q_CELL_AH:.12f} Ah`, `Q_n = {Q_N_C:.12f} C`, `Q_p = {Q_P_C:.12f} C`.",
        f"- Dynamic electrical values: `R0={NOMINAL_PARAMETERS['R0']:.13g} ohm`, `alpha_p={NOMINAL_PARAMETERS['alpha_p']:.12g} s`, `d_n={NOMINAL_PARAMETERS['d_n']:.12g} 1/s`, with fixed `alpha_n={NOMINAL_PARAMETERS['alpha_n']:.12g} s`, `d_p={NOMINAL_PARAMETERS['d_p']:.12g} 1/s`.",
        "- Frozen Guo-White MCMB negative OCP, frozen 13-knot positive OCP, nominal effective stoichiometric windows, and verified full-cell entropy table were used unchanged.",
        f"- Every thermal path starts at the exact interpolated measured battery temperature, {data.measured_temperature_k[0] - 273.15:.9f} degC; no initial offset is fitted.",
        f"- The temperature-aware local voltage calculation reproduces the existing frozen electrical implementation to a maximum absolute difference of {audit['electrical_voltage_reproduction_max_abs_v']:.3e} V before thermal coupling.",
        "",
        "| Frozen electrical state | Minimum | Maximum | Supported OCP interval |",
        "|---|---:|---:|---:|",
        f"| `x_n_avg` | {audit['x_n_avg_min']:.12f} | {audit['x_n_avg_max']:.12f} | [{MCMB_SUPPORTED_STOICHIOMETRY_RANGE[0]:.2f}, {MCMB_SUPPORTED_STOICHIOMETRY_RANGE[1]:.2f}] |",
        f"| `x_n_surf` | {audit['x_n_surface_min']:.12f} | {audit['x_n_surface_max']:.12f} | [{MCMB_SUPPORTED_STOICHIOMETRY_RANGE[0]:.2f}, {MCMB_SUPPORTED_STOICHIOMETRY_RANGE[1]:.2f}] |",
        f"| `x_p_avg` | {audit['x_p_avg_min']:.12f} | {audit['x_p_avg_max']:.12f} | [{NCA_WORKING_STOICHIOMETRY_RANGE[0]:.2f}, {NCA_WORKING_STOICHIOMETRY_RANGE[1]:.2f}] |",
        f"| `x_p_surf` | {audit['x_p_surface_min']:.12f} | {audit['x_p_surface_max']:.12f} | [{NCA_WORKING_STOICHIOMETRY_RANGE[0]:.2f}, {NCA_WORKING_STOICHIOMETRY_RANGE[1]:.2f}] |",
        "",
        "## Overall temperature metrics",
        "",
        "Residual is predicted minus measured battery temperature. Validation A uses measured terminal voltage and measured temperature in reversible heat. Validation B recursively uses frozen-model voltage and predicted temperature, including predicted temperature in the frozen Butler-Volmer terms. The isolation diagnostic uses model voltage but measured temperature in reversible heat.",
        "",
        "| Path | RMSE [K] | MAE [K] | Bias [K] | Max abs [K] | R2 | Measured final [degC] | Predicted final [degC] | Measured rise [K] | Predicted rise [K] |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        f"| Validation A | {_metric_table(a)[2:]}",
        f"| Validation B | {_metric_table(b)[2:]}",
        f"| MODEL-VOLTAGE / MEASURED-T REVERSIBLE DIAGNOSTIC | {_metric_table(iso)[2:]}",
        "",
        "## Metrics by SOC region",
        "",
        "| Path | Region | Samples | RMSE [K] | MAE [K] | Bias [K] | Max abs [K] | R2 |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    masks = soc_masks(data.soc)
    for validation, label in (("validation_a", "Validation A"), ("validation_b", "Validation B")):
        for region, interval in (("high", "SOC >= 2/3"), ("middle", "1/3 <= SOC < 2/3"), ("low", "0.10 <= SOC < 1/3")):
            values = regional[validation][region]
            lines.append(
                f"| {label} | {region} ({interval}) | {np.count_nonzero(masks[region]):,} | {values['rmse_k']:.9f} | {values['mae_k']:.9f} | {values['mean_signed_error_k']:+.9f} | {values['max_abs_error_k']:.9f} | {values['r2']:.9f} |"
            )
    lines.extend(["", "## Heat-component audit", ""])
    lines.extend([
        "The energy integral is consistent with the thermal propagator's left-sample zero-order hold. `Signed reversible fraction` is E_rev/E_total; the absolute-component fraction is integral(|Qrev|)/integral(|Qirr|+|Qrev|).",
        "",
        "| Path | Eirr [J] | Erev [J] | Etotal [J] | Mean Qirr [W] | Mean Qrev [W] | Mean total [W] | Signed reversible fraction | Absolute reversible fraction | Exothermic time fraction | Endothermic time fraction |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for key, label in (("validation_a", "Measured-heat path"), ("validation_b", "Model-heat path")):
        h = heat[key]
        lines.append(
            f"| {label} | {h['irreversible_heat_energy_j']:.9f} | {h['reversible_heat_energy_j']:+.9f} | {h['total_heat_energy_j']:.9f} | {h['mean_qdot_irr_w']:.9f} | {h['mean_qdot_rev_w']:+.9f} | {h['mean_qdot_total_w']:.9f} | {h['signed_reversible_fraction']:+.9f} | {h['absolute_reversible_component_fraction']:.9f} | {h['reversible_exothermic_time_fraction']:.9f} | {h['reversible_endothermic_time_fraction']:.9f} |"
        )
    lines.extend([
        "",
        f"The explicit sign check found {sign_summary['validation_a_reversible_sign_mismatches']} Validation-A and {sign_summary['validation_b_reversible_sign_mismatches']} Validation-B mismatches against `sign(-I*dUeq/dT)`. Thus reversible heat changes sign consistently with both current direction and entropy-coefficient sign.",
        "",
        "| Current / entropy-coefficient case | Expected Qrev sign | Samples | Time fraction | Mean Qrev A [W] | Mean Qrev B [W] |",
        "|---|---|---:|---:|---:|---:|",
    ])
    for case in sign_cases:
        lines.append(
            f"| {str(case['case']).replace('_', ' ')} | {case['expected']} | {case['samples']:,} | {case['time_fraction']:.9f} | {case['mean_qrev_a_w']:+.9f} | {case['mean_qrev_b_w']:+.9f} |"
        )
    lines.extend([
        "",
        "## Dynamic diagnostics",
        "",
        f"Absolute temperature error versus SOC depletion has Spearman rho {diagnostics['validation_a_abs_error_vs_soc_depletion_spearman']:+.6f} for A and {diagnostics['validation_b_abs_error_vs_soc_depletion_spearman']:+.6f} for B. On loaded samples, absolute error versus current magnitude has rho {diagnostics['validation_a_abs_error_vs_current_magnitude_loaded_spearman']:+.6f} for A and {diagnostics['validation_b_abs_error_vs_current_magnitude_loaded_spearman']:+.6f} for B.",
        "",
        f"Regenerative-current RMSE is {diagnostics['validation_a_regenerative_rmse_k']:.6f} K (A) and {diagnostics['validation_b_regenerative_rmse_k']:.6f} K (B), with biases {diagnostics['validation_a_regenerative_bias_k']:+.6f} and {diagnostics['validation_b_regenerative_bias_k']:+.6f} K. Discharge-current RMSE is {diagnostics['validation_a_discharge_rmse_k']:.6f} K (A) and {diagnostics['validation_b_discharge_rmse_k']:.6f} K (B).",
        "",
        f"At the ten largest real adjacent-sample current changes, temperature-error RMSE is {diagnostics['validation_a_large_step_rmse_k']:.6f} K (A) and {diagnostics['validation_b_large_step_rmse_k']:.6f} K (B). Thermal inertia prevents instantaneous temperature jumps; these rows diagnose whether accumulated errors are preferentially present at major electrical transitions.",
        "",
        "| Rank | Source time [s] | Interval elapsed [s] | SOC | I before [A] | I after [A] | delta I [A] | Measured T [degC] | Predicted A [degC] | Error A [K] | Predicted B [degC] | Error B [K] |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for event in events:
        lines.append(
            f"| {event['rank']} | {event['source_time_s']:.6f} | {event['elapsed_time_s']:.6f} | {event['soc']:.9f} | {event['current_before_a']:+.6f} | {event['current_after_a']:+.6f} | {event['delta_current_a']:+.6f} | {event['measured_temperature_c']:.6f} | {event['validation_a_predicted_temperature_c']:.6f} | {event['validation_a_error_k']:+.6f} | {event['validation_b_predicted_temperature_c']:.6f} | {event['validation_b_error_k']:+.6f} |"
        )
    excursion = float(np.ptp(data.measured_temperature_k))
    lines.extend([
        "",
        "## Interpretation",
        "",
        f"The observed battery-temperature excursion within the supported interval is {excursion:.6f} K. Validation A has RMSE {a['rmse_k']:.6f} K ({100.0*a['rmse_k']/excursion:.2f}% of that excursion), bias {a['mean_signed_error_k']:+.6f} K, and R2 {a['r2']:.6f}. **Classification: {classification_a}.**",
        "",
        f"Validation B has RMSE {b['rmse_k']:.6f} K ({100.0*b['rmse_k']/excursion:.2f}% of the excursion), bias {b['mean_signed_error_k']:+.6f} K, and R2 {b['r2']:.6f}. The isolation diagnostic RMSE is {iso['rmse_k']:.6f} K; its difference from B quantifies the incremental reversible-temperature/electrical-temperature feedback effect, while its difference from A reflects model-voltage heat. **Classification: {classification_b}.**",
        "",
        f"Validation A's error is a broad negative bias rather than a low-SOC-only failure: regional RMSE is {regional['validation_a']['high']['rmse_k']:.3f}/{regional['validation_a']['middle']['rmse_k']:.3f}/{regional['validation_a']['low']['rmse_k']:.3f} K (high/middle/low), and the SOC-depletion correlation is only {diagnostics['validation_a_abs_error_vs_soc_depletion_spearman']:+.3f}. Validation B develops clearer SOC-dependent drift, rising to {regional['validation_b']['low']['rmse_k']:.3f} K at low SOC with rho {diagnostics['validation_b_abs_error_vs_soc_depletion_spearman']:+.3f}.",
        "",
        f"Current-magnitude correlations are weak ({diagnostics['validation_a_abs_error_vs_current_magnitude_loaded_spearman']:+.3f} A, {diagnostics['validation_b_abs_error_vs_current_magnitude_loaded_spearman']:+.3f} B), regenerative and discharge RMSE are similar, and largest-step errors are not elevated above overall RMSE. The discrepancy is therefore dominated by persistent heat/temperature underprediction, not isolated current-step or regeneration events. Model-voltage irreversible energy is {100.0 * (heat['validation_b']['irreversible_heat_energy_j'] / heat['validation_a']['irreversible_heat_energy_j'] - 1.0):+.2f}% relative to the measured-voltage path; the nearly identical B and isolation trajectories show that model-voltage heat, not reversible-temperature feedback, causes almost all of the additional coupled-path error.",
        "",
        "The classifications use error scale relative to the measured thermal excursion together with bias, R2, regional residual structure, and dynamic diagnostics; they are not based on a fitted or task-specific acceptance target.",
        "",
        "No repair, refit, entropy/OCP change, Arrhenius term, CC-CV thermal validation, or own-TRC analysis was performed.",
    ])
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    if not MANIFEST_PATH.is_file() or RAW_SOURCE_PATH.name not in MANIFEST_PATH.read_text(encoding="utf-8"):
        raise RuntimeError("The selected UDDS raw file is not recorded in the manifest.")
    if not np.isclose(Q_CELL_AH, 2.994979138519, rtol=0.0, atol=5.0e-13):
        raise RuntimeError("Frozen Panasonic capacity changed.")
    if not (
        np.isclose(IMPLEMENTED_Q_N_C, Q_N_C, rtol=0.0, atol=5.0e-9)
        and np.isclose(IMPLEMENTED_Q_P_C, Q_P_C, rtol=0.0, atol=5.0e-9)
    ):
        raise RuntimeError("Implemented and frozen electrode capacities differ.")

    data, audit = prepare_data()
    qirr_a = data.current_a * (data.ueq_average_v - data.measured_voltage_v)
    qrev_a = -data.current_a * data.measured_temperature_k * data.dudt_v_per_k
    qtotal_a = qirr_a + qrev_a
    predicted_a = exact_thermal_trajectory(
        data.time_s,
        qtotal_a,
        data.ambient_temperature_k,
        data.measured_temperature_k[0],
    )

    predicted_b, model_voltage_b, qirr_b, qrev_b, qtotal_b = fully_predictive_trajectory(data)
    qirr_isolation = data.current_a * (
        data.ueq_average_v - data.model_voltage_measured_t_v
    )
    qtotal_isolation = qirr_isolation + qrev_a
    predicted_isolation = exact_thermal_trajectory(
        data.time_s,
        qtotal_isolation,
        data.ambient_temperature_k,
        data.measured_temperature_k[0],
    )

    arrays = (
        qirr_a, qrev_a, qtotal_a, predicted_a, model_voltage_b,
        qirr_b, qrev_b, qtotal_b, predicted_b, predicted_isolation,
    )
    if not all(np.all(np.isfinite(array)) for array in arrays):
        raise ValueError("A validation result is non-finite.")

    metrics_all = {
        "validation_a": temperature_metrics(data.measured_temperature_k, predicted_a),
        "validation_b": temperature_metrics(data.measured_temperature_k, predicted_b),
        "isolation_diagnostic": temperature_metrics(data.measured_temperature_k, predicted_isolation),
    }
    masks = soc_masks(data.soc)
    regional = {
        validation: {
            region: temperature_metrics(data.measured_temperature_k[mask], predicted[mask])
            for region, mask in masks.items()
        }
        for validation, predicted in (
            ("validation_a", predicted_a),
            ("validation_b", predicted_b),
        )
    }
    heat = {
        "validation_a": heat_audit(data.time_s, qirr_a, qrev_a, qtotal_a),
        "validation_b": heat_audit(data.time_s, qirr_b, qrev_b, qtotal_b),
    }
    sign_summary, sign_cases = sign_context(data, qrev_a, qrev_b)
    diagnostics, events = dynamic_diagnostics(data, predicted_a, predicted_b)

    dt_s = np.diff(data.time_s)
    discharge_interval_ah = np.maximum(data.current_a[:-1], 0.0) * dt_s / 3600.0
    regeneration_interval_ah = np.maximum(-data.current_a[:-1], 0.0) * dt_s / 3600.0
    audit.update(
        {
            "start_soc": float(data.soc[0]),
            "end_soc": float(data.soc[-1]),
            "source_start_time_s": float(data.source_time_s[0]),
            "source_end_time_s": float(data.source_time_s[-1]),
            "duration_s": float(data.time_s[-1]),
            "interval_samples": int(data.time_s.size),
            "dt_min_s": float(np.min(dt_s)),
            "dt_median_s": float(np.median(dt_s)),
            "dt_mean_s": float(np.mean(dt_s)),
            "dt_max_s": float(np.max(dt_s)),
            "current_min_a": float(np.min(data.current_a)),
            "current_max_a": float(np.max(data.current_a)),
            "voltage_min_v": float(np.min(data.measured_voltage_v)),
            "voltage_max_v": float(np.max(data.measured_voltage_v)),
            "battery_temperature_min_c": float(np.min(data.measured_temperature_k) - 273.15),
            "battery_temperature_max_c": float(np.max(data.measured_temperature_k) - 273.15),
            "chamber_temperature_min_c": float(np.min(data.ambient_temperature_k) - 273.15),
            "chamber_temperature_max_c": float(np.max(data.ambient_temperature_k) - 273.15),
            "gross_discharge_throughput_ah": float(np.sum(discharge_interval_ah)),
            "regenerative_charge_throughput_ah": float(np.sum(regeneration_interval_ah)),
            "net_discharge_throughput_ah": float(np.sum(discharge_interval_ah) - np.sum(regeneration_interval_ah)),
            "entropy_domain_violations": int(np.count_nonzero((data.soc < SUPPORTED_SOC_RANGE[0]) | (data.soc > SUPPORTED_SOC_RANGE[1]))),
            "x_n_avg_min": float(np.min(data.x_n_avg)),
            "x_n_avg_max": float(np.max(data.x_n_avg)),
            "x_n_surface_min": float(np.min(data.x_n_surface)),
            "x_n_surface_max": float(np.max(data.x_n_surface)),
            "x_p_avg_min": float(np.min(data.x_p_avg)),
            "x_p_avg_max": float(np.max(data.x_p_avg)),
            "x_p_surface_min": float(np.min(data.x_p_surface)),
            "x_p_surface_max": float(np.max(data.x_p_surface)),
        }
    )
    audit.update(sign_summary)

    classification_a = classify(metrics_all["validation_a"], data.measured_temperature_k, coupled=False)
    classification_b = classify(metrics_all["validation_b"], data.measured_temperature_k, coupled=True)

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    save_temperature_plot(data, predicted_a, predicted_b, predicted_isolation)
    save_residual_plots(
        data,
        predicted_a - data.measured_temperature_k,
        predicted_b - data.measured_temperature_k,
        predicted_isolation - data.measured_temperature_k,
    )
    save_heat_plot(data, (qirr_a, qrev_a, qtotal_a), (qirr_b, qrev_b, qtotal_b))
    write_metrics(data, metrics_all, regional, heat, diagnostics, events, audit)
    write_report(
        data, metrics_all, regional, heat, diagnostics, events,
        sign_summary, sign_cases, audit, classification_a, classification_b,
    )

    print("Panasonic UDDS thermal validation complete (no refit).")
    print(f"Supported interval: SOC {data.soc[0]:.12f} to {data.soc[-1]:.12f}; samples={data.time_s.size:,}; duration={data.time_s[-1]:.6f} s")
    print(f"Validation A: RMSE={metrics_all['validation_a']['rmse_k']:.9f} K; bias={metrics_all['validation_a']['mean_signed_error_k']:+.9f} K; R2={metrics_all['validation_a']['r2']:.9f}; {classification_a}")
    print(f"Validation B: RMSE={metrics_all['validation_b']['rmse_k']:.9f} K; bias={metrics_all['validation_b']['mean_signed_error_k']:+.9f} K; R2={metrics_all['validation_b']['r2']:.9f}; {classification_b}")
    print(f"Outputs: {REPORT_PATH.relative_to(REPOSITORY_ROOT)}, {METRICS_PATH.relative_to(REPOSITORY_ROOT)}, and four figures")


if __name__ == "__main__":
    main()
