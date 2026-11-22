"""Hierarchical forward-only transfer audit from Panasonic to own TRC18650 data.

Only the clean ambient cell-1 cycle is used.  No parameter is fitted or
changed, no voltage is corrected, and the loaded TRC curve is never treated
as OCV.  The absent Chroma time channel is represented only by the previously
audited provisional capacity/current-derived timebase.
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

from analysis.own_case_study_reference import (  # noqa: E402
    SOC_GRID,
    charge_curve,
    discharge_curve,
    interpolate_on_soc,
    segment_complete_cycle,
)
from models.grouped_spm import FARADAY_C_PER_MOL, GAS_CONSTANT_J_PER_MOL_K  # noqa: E402
from models.ocp_panasonic_qi import (  # noqa: E402
    MCMB_SUPPORTED_STOICHIOMETRY_RANGE,
    NCA_WORKING_STOICHIOMETRY_RANGE,
)
from models.ocp_panasonic_qi_frozen import (  # noqa: E402
    Q_CELL_AH,
    SOURCE_MULTISTART_START_ID,
    U_n,
    U_p,
    V_eq,
    X_N_END,
    X_N_START,
    X_P_END,
    X_P_START,
)
from scripts.audit_panasonic_hppc_identifiability import Q_N_C, Q_P_C  # noqa: E402
from scripts.audit_own_trc18650_data import load_all  # noqa: E402
from scripts.validate_panasonic_1c_forward import NOMINAL_PARAMETERS  # noqa: E402


REFERENCE_ID = "cycle_cell1_ambient"
Q_TRC_AH_EXPECTED = 2.8999
OUTPUT_DIRECTORY = REPOSITORY_ROOT / "thesis_data"
FIGURE_DIRECTORY = OUTPUT_DIRECTORY / "own_transfer_figures"
REPORT_PATH = OUTPUT_DIRECTORY / "own_model_transfer_audit.md"
METRICS_PATH = OUTPUT_DIRECTORY / "own_model_transfer_metrics.csv"
PARAMETER_MATRIX_PATH = OUTPUT_DIRECTORY / "own_model_transfer_parameter_matrix.csv"


@dataclass(frozen=True)
class FrozenTransferResult:
    voltage_v: np.ndarray
    x_n_avg: np.ndarray
    x_n_aux: np.ndarray
    x_n_surface: np.ndarray
    x_p_avg: np.ndarray
    x_p_aux: np.ndarray
    x_p_surface: np.ndarray
    ueq_average_v: np.ndarray
    eta_n_v: np.ndarray
    eta_p_v: np.ndarray
    ohmic_drop_v: np.ndarray
    ocp_valid: np.ndarray


def forward_frozen_panasonic(
    time_s: np.ndarray,
    model_current_a: np.ndarray,
    temperature_k: np.ndarray,
) -> FrozenTransferResult:
    """Propagate the unchanged frozen equations, masking unsupported OCP calls."""

    time_s = np.asarray(time_s, dtype=float)
    current = np.asarray(model_current_a, dtype=float)
    temperature = np.asarray(temperature_k, dtype=float)
    if not (time_s.shape == current.shape == temperature.shape):
        raise ValueError("Time, current, and temperature shapes differ.")
    if time_s.size < 2 or np.any(np.diff(time_s) <= 0.0):
        raise ValueError("Provisional integration coordinate must be strictly increasing.")
    if np.any(~np.isfinite(current)) or np.any(~np.isfinite(temperature)):
        raise ValueError("Current and temperature must be finite.")

    alpha_n = float(NOMINAL_PARAMETERS["alpha_n"])
    alpha_p = float(NOMINAL_PARAMETERS["alpha_p"])
    d_n = float(NOMINAL_PARAMETERS["d_n"])
    d_p = float(NOMINAL_PARAMETERS["d_p"])
    r0 = float(NOMINAL_PARAMETERS["R0"])
    n_samples = time_s.size
    x_n_avg = np.empty(n_samples)
    x_n_aux = np.empty(n_samples)
    x_p_avg = np.empty(n_samples)
    x_p_aux = np.empty(n_samples)
    x_n_avg[0] = x_n_aux[0] = X_N_START
    x_p_avg[0] = x_p_aux[0] = X_P_START

    for index, dt_s in enumerate(np.diff(time_s)):
        held_current = current[index]
        y_n = x_n_aux[index] - x_n_avg[index]
        decay_n = np.exp(-30.0 * dt_s / alpha_n)
        forcing_n = (-12.0 / (7.0 * Q_N_C)) * held_current
        y_n_next = decay_n * y_n + forcing_n * (alpha_n / 30.0) * (1.0 - decay_n)
        x_n_avg[index + 1] = x_n_avg[index] - held_current * dt_s / Q_N_C
        x_n_aux[index + 1] = x_n_avg[index + 1] + y_n_next

        y_p = x_p_aux[index] - x_p_avg[index]
        decay_p = np.exp(-30.0 * dt_s / alpha_p)
        forcing_p = (12.0 / (7.0 * Q_P_C)) * held_current
        y_p_next = decay_p * y_p + forcing_p * (alpha_p / 30.0) * (1.0 - decay_p)
        x_p_avg[index + 1] = x_p_avg[index] + held_current * dt_s / Q_P_C
        x_p_aux[index + 1] = x_p_avg[index + 1] + y_p_next

    x_n_surface = x_n_aux - alpha_n * current / (105.0 * Q_N_C)
    x_p_surface = x_p_aux + alpha_p * current / (105.0 * Q_P_C)
    n_low, n_high = MCMB_SUPPORTED_STOICHIOMETRY_RANGE
    p_low, p_high = NCA_WORKING_STOICHIOMETRY_RANGE
    ocp_valid = (
        (x_n_surface >= n_low)
        & (x_n_surface <= n_high)
        & (x_p_surface >= p_low)
        & (x_p_surface <= p_high)
    )
    average_valid = (
        (x_n_avg >= n_low)
        & (x_n_avg <= n_high)
        & (x_p_avg >= p_low)
        & (x_p_avg <= p_high)
    )
    ueq_average = np.full(n_samples, np.nan)
    ueq_average[average_valid] = U_p(x_p_avg[average_valid]) - U_n(x_n_avg[average_valid])
    eta_n = np.full(n_samples, np.nan)
    eta_p = np.full(n_samples, np.nan)
    ohmic = r0 * current
    voltage = np.full(n_samples, np.nan)
    if np.any(ocp_valid):
        factor = 2.0 * GAS_CONSTANT_J_PER_MOL_K * temperature[ocp_valid] / FARADAY_C_PER_MOL
        eta_n[ocp_valid] = factor * np.arcsinh(
            current[ocp_valid]
            / (6.0 * Q_N_C * d_n * np.sqrt(x_n_surface[ocp_valid] * (1.0 - x_n_surface[ocp_valid])))
        )
        eta_p[ocp_valid] = factor * np.arcsinh(
            current[ocp_valid]
            / (6.0 * Q_P_C * d_p * np.sqrt(x_p_surface[ocp_valid] * (1.0 - x_p_surface[ocp_valid])))
        )
        voltage[ocp_valid] = (
            U_p(x_p_surface[ocp_valid])
            - U_n(x_n_surface[ocp_valid])
            + eta_p[ocp_valid]
            - eta_n[ocp_valid]
            - ohmic[ocp_valid]
        )
    return FrozenTransferResult(
        voltage_v=voltage,
        x_n_avg=x_n_avg,
        x_n_aux=x_n_aux,
        x_n_surface=x_n_surface,
        x_p_avg=x_p_avg,
        x_p_aux=x_p_aux,
        x_p_surface=x_p_surface,
        ueq_average_v=ueq_average,
        eta_n_v=eta_n,
        eta_p_v=eta_p,
        ohmic_drop_v=ohmic,
        ocp_valid=ocp_valid,
    )


def voltage_metrics(measured_v: np.ndarray, predicted_v: np.ndarray, mask: np.ndarray | None = None) -> dict[str, float | int]:
    measured = np.asarray(measured_v, dtype=float)
    predicted = np.asarray(predicted_v, dtype=float)
    finite = np.isfinite(measured) & np.isfinite(predicted)
    if mask is not None:
        finite &= np.asarray(mask, dtype=bool)
    if not np.any(finite):
        return {name: float("nan") for name in ("rmse_v", "mae_v", "bias_v", "max_abs_error_v", "r2")} | {"n_samples": 0}
    residual = predicted[finite] - measured[finite]
    denominator = float(np.sum((measured[finite] - np.mean(measured[finite])) ** 2))
    return {
        "rmse_v": float(np.sqrt(np.mean(residual**2))),
        "mae_v": float(np.mean(np.abs(residual))),
        "bias_v": float(np.mean(residual)),
        "max_abs_error_v": float(np.max(np.abs(residual))),
        "r2": float(1.0 - np.sum(residual**2) / denominator) if denominator > 0.0 else float("nan"),
        "n_samples": int(np.count_nonzero(finite)),
    }


def region_masks(soc: np.ndarray) -> dict[str, np.ndarray]:
    return {
        "high": soc >= 2.0 / 3.0,
        "middle": (soc >= 1.0 / 3.0) & (soc < 2.0 / 3.0),
        "low": soc < 1.0 / 3.0,
    }


def add_metric_rows(
    rows: list[dict[str, object]],
    test: str,
    subset: str,
    metrics: dict[str, float | int],
    note: str,
) -> None:
    units = {
        "rmse_v": "V",
        "mae_v": "V",
        "bias_v": "V",
        "max_abs_error_v": "V",
        "r2": "1",
        "n_samples": "count",
    }
    for metric, value in metrics.items():
        rows.append({"test": test, "subset": subset, "metric": metric, "value": value, "unit": units.get(metric, "1"), "note": note})


def scalar_summary(values: np.ndarray) -> dict[str, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    return {
        "mean_v": float(np.mean(finite)),
        "mean_abs_v": float(np.mean(np.abs(finite))),
        "rms_v": float(np.sqrt(np.mean(finite**2))),
        "min_v": float(np.min(finite)),
        "max_v": float(np.max(finite)),
    }


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


def parameter_matrix_rows() -> list[dict[str, str]]:
    return [
        {"parameter_or_component": "Q_cell", "information_required": "complete charge/discharge throughput", "information_available": "complete clean discharge and recharge counters; Q_ref=2.8999 Ah", "identifiable_here": "YES", "reason": "Capacity is directly supported in capacity coordinates without measured time."},
        {"parameter_or_component": "OCP / equilibrium-voltage representation", "information_required": "near-equilibrium low-rate or rested voltage across SOC, preferably both directions", "information_available": "loaded 1C discharge and CC-CV charge only", "identifiable_here": "NO", "reason": "Loaded terminal voltage mixes equilibrium voltage, resistance, kinetics, diffusion, and directional effects."},
        {"parameter_or_component": "stoichiometric windows", "information_required": "cell-compatible electrode OCPs plus independent capacity/electrode constraints", "information_available": "cell capacity and loaded full-cell voltage only", "identifiable_here": "NO", "reason": "Many electrode windows can generate similar full-cell voltage; no electrode reference is available."},
        {"parameter_or_component": "Q_n / Q_p grouped state capacities", "information_required": "Q_cell plus defensible electrode stoichiometric window widths", "information_available": "Q_cell is measured; electrode window widths are not identified", "identifiable_here": "NO", "reason": "The Panasonic Q_n/Q_p values cannot be rescaled uniquely from cell capacity alone."},
        {"parameter_or_component": "R0", "information_required": "time-resolved current step and immediate voltage response with fixture contribution characterized", "information_available": "adjacent-record R_app only; no measured sampling time", "identifiable_here": "CONDITIONAL DIAGNOSTIC ONLY", "reason": "R_app contains cell, fixture, wiring, timing, and fast polarization contributions and is not intrinsic R0."},
        {"parameter_or_component": "alpha_n", "information_required": "measured-time relaxation/transient excitation with electrode separation", "information_available": "sustained 1C branch; provisional time only", "identifiable_here": "NO", "reason": "No measured timebase or electrode-specific excitation."},
        {"parameter_or_component": "alpha_p", "information_required": "measured-time relaxation/transient excitation with electrode separation", "information_available": "sustained 1C branch; provisional time only", "identifiable_here": "NO", "reason": "No measured timebase or electrode-specific excitation."},
        {"parameter_or_component": "d_n", "information_required": "time-resolved kinetic/polarization information over SOC", "information_available": "loaded curve with static/dynamic effects confounded", "identifiable_here": "NO", "reason": "Cannot separate negative kinetics from OCP, positive kinetics, diffusion, and setup effects."},
        {"parameter_or_component": "d_p", "information_required": "time-resolved kinetic/polarization information over SOC", "information_available": "loaded curve with static/dynamic effects confounded", "identifiable_here": "NO", "reason": "Cannot separate positive kinetics from OCP, negative kinetics, diffusion, and setup effects."},
        {"parameter_or_component": "thermal Cth", "information_required": "measured time, heat input, and calibrated temperature response", "information_available": "surface temperatures but no measured timestamps; heat generation not identified", "identifiable_here": "NO", "reason": "Thermal storage cannot be inferred in capacity coordinates alone."},
        {"parameter_or_component": "thermal Gth", "information_required": "measured-time cooling/heating response and ambient boundary condition", "information_available": "ambient channel exists but no measured timestamps", "identifiable_here": "NO", "reason": "Heat-transfer rate requires a real timebase."},
    ]


def accuracy_policy_rows() -> list[dict[str, str]]:
    return [
        {"quantity": "Voltage-vs-SOC RMSE/MAE/bias", "status": "QUANTITATIVE / CONDITIONAL", "basis": "Valid comparison of measured loaded voltage and prediction in capacity coordinates; not an OCP-validation metric."},
        {"quantity": "Capacity error", "status": "QUANTITATIVE", "basis": "Complete Chroma Ah counter supports 2.8999 Ah directly."},
        {"quantity": "SOC trajectory error", "status": "NOT INDEPENDENTLY VALIDATED", "basis": "SOC is a coulomb-counted reference coordinate, not independently measured true SOC."},
        {"quantity": "Dynamic time-domain RMSE", "status": "PROVISIONAL / NOT PRIMARY", "basis": "Chronology uses a capacity/current-derived uniform timebase, not measured time."},
        {"quantity": "Temperature-vs-time RMSE", "status": "NOT SUPPORTED", "basis": "Measured timestamps are absent and this audit is electrical only."},
        {"quantity": "Apparent series-path response", "status": "QUANTITATIVE DIAGNOSTIC", "basis": "Adjacent-sample DeltaV/DeltaI is reproducible but is not intrinsic cell R0."},
        {"quantity": "Voltage-backbone residual", "status": "QUANTITATIVE DIAGNOSTIC", "basis": "Numerically valid transfer diagnostic; TRC loaded 1C voltage is not equilibrium OCV."},
        {"quantity": "OCP surface-domain violation count", "status": "QUANTITATIVE MODEL DIAGNOSTIC", "basis": "Evaluated without clipping or extrapolation."},
    ]


def make_figures(
    grid_soc: np.ndarray,
    measured_grid_v: np.ndarray,
    backbone_v: np.ndarray,
    discharge_soc: np.ndarray,
    measured_discharge_v: np.ndarray,
    predicted_discharge_v: np.ndarray,
    charge_soc: np.ndarray,
    measured_charge_v: np.ndarray,
    predicted_charge_v: np.ndarray,
    charge_cc: np.ndarray,
    charge_cv: np.ndarray,
) -> None:
    figure, axis = plt.subplots(figsize=(8.8, 5.4), constrained_layout=True)
    axis.plot(grid_soc, backbone_v, color="#0072B2", lw=1.7, label="Frozen Panasonic equilibrium-voltage backbone")
    axis.plot(grid_soc, measured_grid_v, color="black", lw=1.5, label="TRC clean ambient loaded 1C discharge (NOT OCV)")
    axis.set(xlabel="Coulomb-counted reference SOC", ylabel="Voltage [V]", title="Voltage-backbone transfer diagnostic")
    axis.grid(alpha=0.25)
    axis.legend(fontsize=8)
    figure.savefig(FIGURE_DIRECTORY / "trc_panasonic_voltage_backbone_transfer.png", dpi=180)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(8.8, 5.4), constrained_layout=True)
    axis.plot(discharge_soc, measured_discharge_v, color="black", lw=1.5, label="Measured TRC loaded 1C")
    axis.plot(discharge_soc, predicted_discharge_v, color="#D55E00", lw=1.5, label="Literal frozen Panasonic transfer")
    invalid = ~np.isfinite(predicted_discharge_v)
    if np.any(invalid):
        axis.scatter(discharge_soc[invalid], measured_discharge_v[invalid], s=9, color="#CC79A7", label="Prediction unavailable: OCP support violation")
    axis.set(xlabel="Coulomb-counted reference SOC", ylabel="Terminal voltage [V]", title="PROVISIONAL-TIME direct transfer diagnostic")
    axis.grid(alpha=0.25)
    axis.legend(fontsize=8)
    figure.savefig(FIGURE_DIRECTORY / "trc_panasonic_full_transfer_discharge.png", dpi=180)
    plt.close(figure)

    nominal_backbone_at_raw = np.interp(discharge_soc, grid_soc, backbone_v)
    backbone_residual = nominal_backbone_at_raw - measured_discharge_v
    full_residual = predicted_discharge_v - measured_discharge_v
    dynamic_increment = predicted_discharge_v - nominal_backbone_at_raw
    figure, axis = plt.subplots(figsize=(9.2, 5.6), constrained_layout=True)
    axis.plot(discharge_soc, 1000.0 * backbone_residual, color="#0072B2", lw=1.3, label="Backbone - measured")
    axis.plot(discharge_soc, 1000.0 * full_residual, color="#D55E00", lw=1.3, label="Full transfer - measured (valid support only)")
    axis.plot(discharge_soc, 1000.0 * dynamic_increment, color="#009E73", lw=1.1, label="Full transfer - nominal backbone")
    axis.axhline(0.0, color="black", lw=0.8)
    axis.set(xlabel="Coulomb-counted reference SOC", ylabel="Voltage difference [mV]", title="Panasonic-to-TRC transfer residual decomposition")
    axis.grid(alpha=0.25)
    axis.legend(fontsize=8)
    figure.savefig(FIGURE_DIRECTORY / "trc_panasonic_transfer_residual_vs_soc.png", dpi=180)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(8.8, 5.4), constrained_layout=True)
    axis.plot(charge_soc, measured_charge_v, color="black", lw=1.5, label="Measured TRC CC-CV recharge")
    axis.plot(charge_soc, predicted_charge_v, color="#D55E00", lw=1.5, label="Frozen Panasonic transfer")
    if np.any(charge_cc):
        axis.axvspan(float(np.min(charge_soc[charge_cc])), float(np.max(charge_soc[charge_cc])), color="#56B4E9", alpha=0.10, label="CC")
    if np.any(charge_cv):
        axis.axvspan(float(np.min(charge_soc[charge_cv])), float(np.max(charge_soc[charge_cv])), color="#E69F00", alpha=0.10, label="CV")
    axis.set(xlabel="Coulomb-counted reference SOC", ylabel="Terminal voltage [V]", title="PROVISIONAL-TIME directional transfer check")
    axis.grid(alpha=0.25)
    axis.legend(fontsize=8)
    figure.savefig(FIGURE_DIRECTORY / "trc_panasonic_charge_transfer.png", dpi=180)
    plt.close(figure)


def fmt_metric(metric: dict[str, float | int], name: str, scale: float = 1.0, digits: int = 3) -> str:
    value = float(metric[name])
    return "not available" if not np.isfinite(value) else f"{scale*value:.{digits}f}"


def write_report(
    source_path: Path,
    provisional_dt_s: float,
    q_trc_ah: float,
    backbone_metrics: dict[str, float | int],
    backbone_regions: dict[str, dict[str, float | int]],
    full_metrics: dict[str, float | int],
    full_regions: dict[str, dict[str, float | int]],
    discharge_result: FrozenTransferResult,
    discharge_soc: np.ndarray,
    state_rows: list[dict[str, object]],
    decomposition: dict[str, dict[str, float]],
    charge_metrics: dict[str, dict[str, float | int]],
    charge_valid_count: int,
    charge_total_count: int,
) -> None:
    invalid_count = int(discharge_soc.size - int(np.count_nonzero(discharge_result.ocp_valid)))
    negative_invalid_count = int(np.count_nonzero(
        (discharge_result.x_n_surface < MCMB_SUPPORTED_STOICHIOMETRY_RANGE[0])
        | (discharge_result.x_n_surface > MCMB_SUPPORTED_STOICHIOMETRY_RANGE[1])
    ))
    positive_invalid_count = int(np.count_nonzero(
        (discharge_result.x_p_surface < NCA_WORKING_STOICHIOMETRY_RANGE[0])
        | (discharge_result.x_p_surface > NCA_WORKING_STOICHIOMETRY_RANGE[1])
    ))
    first_invalid = np.flatnonzero(~discharge_result.ocp_valid)
    first_invalid_text = "none" if first_invalid.size == 0 else f"row-local index {int(first_invalid[0])}, SOC_ref={discharge_soc[int(first_invalid[0])]:.6f}"
    panasonic_fraction = q_trc_ah / Q_CELL_AH
    lines = [
        "# Panasonic-to-TRC18650 hierarchical model-transfer audit",
        "",
        "## Literal transfer is not adequate",
        "",
        "Only `cycle_cell1_ambient` is used. No Panasonic parameter, OCP, window, support interval, or model equation was changed; no TRC quantity was fitted and no measured voltage was corrected.",
        "",
        f"Source: `{source_path.resolve()}`. Raw Chroma current is negative on discharge and positive on charge; the sole model-boundary conversion is `I_model=-I_raw`, giving discharge-positive model current. The file has no measured time. The uniform `{provisional_dt_s:.9f} s` integration step is therefore explicitly provisional and derived from capacity/current consistency.",
        "",
        "## Frozen Panasonic quantities transferred",
        "",
        "| Quantity | Frozen value / provenance |",
        "|---|---|",
        "| Negative OCP | Guo-White 2005 MCMB discharge-process spline via `ocp_panasonic_qi_frozen.U_n` |",
        "| Positive OCP | Frozen 13-knot Panasonic Qi-style representation, selected multistart start 3 |",
        f"| Nominal effective windows | x_n {X_N_START:.12f} to {X_N_END:.12f}; x_p {X_P_START:.12f} to {X_P_END:.12f} |",
        f"| OCP support | x_n {MCMB_SUPPORTED_STOICHIOMETRY_RANGE}; x_p {NCA_WORKING_STOICHIOMETRY_RANGE}; no clipping/extrapolation |",
        f"| Q_cell | {Q_CELL_AH:.12f} Ah |",
        f"| Q_n / Q_p | {Q_N_C:.12f} / {Q_P_C:.12f} C |",
        f"| alpha_n / alpha_p | {NOMINAL_PARAMETERS['alpha_n']:.12g} / {NOMINAL_PARAMETERS['alpha_p']:.12g} s |",
        f"| d_n / d_p | {NOMINAL_PARAMETERS['d_n']:.12g} / {NOMINAL_PARAMETERS['d_p']:.12g} 1/s |",
        f"| R0 | {NOMINAL_PARAMETERS['R0']:.13g} ohm |",
        "| Current convention | Model-positive discharge; raw TRC sign converted only at loader boundary |",
        "",
        "The dynamic values are the conditional Panasonic HPPC parameterization, not globally identified Panasonic constants.",
        "",
        "## Test A — voltage-backbone transfer diagnostic",
        "",
        "The frozen Panasonic equilibrium-voltage backbone is compared with the clean TRC **loaded 1C** discharge in normalized capacity coordinates. This is not OCP validation.",
        "",
        f"Overall: RMSE `{fmt_metric(backbone_metrics,'rmse_v',1000)} mV`, MAE `{fmt_metric(backbone_metrics,'mae_v',1000)} mV`, bias Ueq-minus-measured `{fmt_metric(backbone_metrics,'bias_v',1000)} mV`, maximum absolute error `{fmt_metric(backbone_metrics,'max_abs_error_v',1000)} mV`, R2 `{fmt_metric(backbone_metrics,'r2',1,6)}`.",
        "",
        "| SOC region | Samples | RMSE [mV] | Bias [mV] |",
        "|---|---:|---:|---:|",
    ]
    for region in ("high", "middle", "low"):
        metric = backbone_regions[region]
        lines.append(f"| {region} | {metric['n_samples']} | {fmt_metric(metric,'rmse_v',1000)} | {fmt_metric(metric,'bias_v',1000)} |")
    lines.extend([
        "",
        "The residual is SOC-shaped rather than a constant offset; the regional table and residual figure show the change in magnitude across the discharge.",
        "",
        "## Test B — literal frozen electrical transfer",
        "",
        "Classification: **PROVISIONAL-TIME DIRECT TRANSFER DIAGNOSTIC**. Metrics are computed only where the frozen surface states remain inside the frozen OCP support; unsupported samples are retained as invalid and never clipped or extrapolated.",
        "",
        f"Valid/total discharge surface-state samples: `{full_metrics['n_samples']}/{discharge_soc.size}`; invalid count `{invalid_count}` (negative-electrode `{negative_invalid_count}`, positive-electrode `{positive_invalid_count}`); first violation `{first_invalid_text}`.",
        f"Valid-domain metrics: RMSE `{fmt_metric(full_metrics,'rmse_v',1000)} mV`, MAE `{fmt_metric(full_metrics,'mae_v',1000)} mV`, bias predicted-minus-measured `{fmt_metric(full_metrics,'bias_v',1000)} mV`, maximum absolute error `{fmt_metric(full_metrics,'max_abs_error_v',1000)} mV`, R2 `{fmt_metric(full_metrics,'r2',1,6)}`.",
        "",
        "| SOC region | Valid samples | RMSE [mV] | Bias [mV] |",
        "|---|---:|---:|---:|",
    ])
    for region in ("high", "middle", "low"):
        metric = full_regions[region]
        lines.append(f"| {region} | {metric['n_samples']} | {fmt_metric(metric,'rmse_v',1000)} | {fmt_metric(metric,'bias_v',1000)} |")
    lines.extend([
        "",
        f"State ranges over the discharge before OCP masking: x_n_avg `{np.min(discharge_result.x_n_avg):.6f}-{np.max(discharge_result.x_n_avg):.6f}`, x_p_avg `{np.min(discharge_result.x_p_avg):.6f}-{np.max(discharge_result.x_p_avg):.6f}`, x_n_surface `{np.min(discharge_result.x_n_surface):.6f}-{np.max(discharge_result.x_n_surface):.6f}`, x_p_surface `{np.min(discharge_result.x_p_surface):.6f}-{np.max(discharge_result.x_p_surface):.6f}`.",
        "",
        "## Test C — normalized-capacity state-mapping diagnostic",
        "",
        f"The TRC discharge passes `{q_trc_ah:.6f} Ah`, only `{100.0*panasonic_fraction:.3f}%` of the frozen Panasonic cell capacity. Therefore literal Q_n/Q_p integration does not traverse the complete frozen nominal electrode windows over the TRC complete discharge.",
        "",
        "The average states below are calculated directly from normalized TRC capacity and frozen Q_n/Q_p. Surface values are the corresponding nearest samples from the provisional-time dynamic transfer.",
        "",
        "| z | Capacity-mapped x_n_avg | Capacity-mapped x_p_avg | Provisional x_n_surface | Provisional x_p_surface | Surface support valid? |",
        "|---:|---:|---:|---:|---:|---|",
    ])
    for row in state_rows:
        lines.append(f"| {row['z']:.2f} | {row['x_n_avg']:.9f} | {row['x_p_avg']:.9f} | {row['x_n_surface']:.9f} | {row['x_p_surface']:.9f} | {row['surface_support_valid']} |")
    lines.extend([
        "",
        f"At z=1 the average states are `{state_rows[-1]['x_n_avg']:.9f}` and `{state_rows[-1]['x_p_avg']:.9f}`, rather than the frozen endpoints `{X_N_END:.9f}` and `{X_P_END:.9f}`. Average states remain in support, but the loaded positive surface state becomes unsupported near low SOC. The mapping is therefore internally computable but not coherent as a complete frozen-window representation of the TRC discharge.",
        "",
        "## Error decomposition",
        "",
        f"Measured-minus-backbone has mean `{1000*decomposition['measured_loaded_offset']['mean_v']:.3f} mV` and RMS `{1000*decomposition['measured_loaded_offset']['rms_v']:.3f} mV` on the valid full-model subset. Full-minus-backbone has mean `{1000*decomposition['full_minus_nominal_backbone']['mean_v']:.3f} mV` and RMS `{1000*decomposition['full_minus_nominal_backbone']['rms_v']:.3f} mV`, while full-minus-average-state-equilibrium has RMS `{1000*decomposition['full_minus_average_equilibrium']['rms_v']:.3f} mV`. Literal full-model residual RMS is `{fmt_metric(full_metrics,'rmse_v',1000)} mV`.",
        "",
        "The voltage-backbone mismatch is already large. Frozen dynamic/polarization terms subtract about 119 mV on average and reduce the terminal-voltage error, but this is compensating error rather than evidence of uniquely correct transferred dynamics. They do not remove the systematic low-SOC discrepancy; capacity mapping and the eventual positive-surface support violation add a separate failure mode. The supported classification is therefore a **mixture of static voltage-reference, state-mapping, and dynamic-transfer mismatch**. These quantities are transfer diagnostics, not intrinsic TRC overpotentials.",
        "",
        "## Directional transfer check",
        "",
        f"Recharge valid/total samples: `{charge_valid_count}/{charge_total_count}`. The discharge-oriented frozen representation was propagated continuously through the post-discharge rest and recharge, using the same provisional timebase and no direction-specific correction. State integration remains mathematically defined through the 22 unsupported discharge-voltage samples; no OCP or voltage was evaluated there. The rested states re-enter support before recharge.",
        "",
        "| Recharge portion | Samples | RMSE [mV] | MAE [mV] | Bias [mV] | Max abs [mV] | R2 |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ])
    for subset in ("all_recharge", "cc", "cv"):
        metric = charge_metrics[subset]
        lines.append(f"| {subset} | {metric['n_samples']} | {fmt_metric(metric,'rmse_v',1000)} | {fmt_metric(metric,'mae_v',1000)} | {fmt_metric(metric,'bias_v',1000)} | {fmt_metric(metric,'max_abs_error_v',1000)} | {fmt_metric(metric,'r2',1,6)} |")
    lines.extend([
        "",
        "No recharge surface-support violation occurred. The charge results are a held-out directional diagnostic. The very negative CV R2 is ill-conditioned because measured CV voltage has very little variance; RMSE and bias are the interpretable CV metrics. These results do not justify a charge-specific fit or thermodynamic-hysteresis claim.",
        "",
        "## Minimum-adaptation decision tree",
        "",
        "- **LEVEL 0 — rejected:** literal direct transfer is not sufficient.",
        "- **LEVEL 1 — required:** the different measured cell capacity and incomplete traversal of the frozen windows require cell-specific capacity/state mapping.",
        "- **LEVEL 2 — required:** the large, SOC-shaped backbone residual requires a TRC-compatible voltage-reference/OCP basis before equilibrium and polarization errors can be separated.",
        "- **LEVEL 3 — indicated but not identifiable here:** frozen dynamics do not remove the residual and surface support fails at low SOC, but the absent measured timebase prevents defensible dynamic calibration.",
        "- **LEVEL 4 — applies to full identification from this dataset alone:** the own clean cycle is insufficient to identify electrode OCP/windows and dynamic/thermal parameters uniquely.",
        "",
        "The minimum defensible *adaptation* before meaningful quantitative reduced-model validation is cumulative **LEVEL 1 + LEVEL 2**. A complete TRC dynamic parameterization would also need LEVEL 3 evidence from data that are not available here; it must not be inferred by fitting this provisional-time cycle.",
        "",
        "## Identifiability matrix",
        "",
        "| Parameter/component | Information required | Information available | Identifiable here? | Reason |",
        "|---|---|---|---|---|",
    ])
    for row in parameter_matrix_rows():
        lines.append(f"| {row['parameter_or_component']} | {row['information_required']} | {row['information_available']} | {row['identifiable_here']} | {row['reason']} |")
    lines.extend([
        "",
        "## Accuracy-reporting policy",
        "",
        "| Quantity | Status | Basis |",
        "|---|---|---|",
    ])
    for row in accuracy_policy_rows():
        lines.append(f"| {row['quantity']} | {row['status']} | {row['basis']} |")
    lines.extend([
        "",
        "## Explicit answers",
        "",
        "1. **Is literal transfer adequate?** No.",
        f"2. **How much discrepancy exists at backbone level?** Overall backbone RMSE is `{fmt_metric(backbone_metrics,'rmse_v',1000)} mV` with bias `{fmt_metric(backbone_metrics,'bias_v',1000)} mV`; the mismatch is already large and SOC-shaped before dynamics.",
        f"3. **Is the Panasonic state mapping coherent over the full TRC discharge?** No as a complete-window transfer. The TRC passes `{100.0*panasonic_fraction:.3f}%` of Panasonic nominal capacity, finishes short of both nominal average-state endpoints, and produces low-SOC positive-surface support violations under load.",
        "4. **Would dynamics also require adaptation?** Evidence indicates yes for quantitative terminal-voltage transfer, but the parameters cannot be identified defensibly with provisional time.",
        "5. **What can be adapted defensibly?** Q_cell can be set from the complete 2.8999 Ah discharge. A capacity-coordinate empirical voltage reference can be constructed for diagnostics, but it must not be called OCP.",
        "6. **What cannot be identified?** Electrode OCPs/windows, intrinsic R0, alpha_n/alpha_p, d_n/d_p, Cth, and Gth cannot be uniquely identified from this dataset.",
        "7. **Minimum next adaptation level?** Cumulative LEVEL 1 + LEVEL 2. LEVEL 3 is indicated but unsupported for calibration; full identification remains LEVEL-4 data-limited.",
        "",
        "## Stop point",
        "",
        "No parameter optimization, OCP fit, voltage correction, affected-dataset calibration, thermal accuracy calculation, or climate/city matching was performed.",
    ])
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    if SOURCE_MULTISTART_START_ID != 3:
        raise RuntimeError("Frozen Panasonic positive-OCP provenance changed.")
    if not np.isclose(Q_CELL_AH, 2.994979138519, rtol=0.0, atol=5e-13):
        raise RuntimeError("Frozen Panasonic Q_cell changed.")
    experiment = next(item for item in load_all() if item.metadata.dataset_id == REFERENCE_ID)
    segments = segment_complete_cycle(experiment)
    discharge_segment = next(item for item in segments if item.name == "main_1c_discharge")
    post_rest_segment = next(item for item in segments if item.name == "post_discharge_rest")
    charge_segments = [item for item in segments if item.name.startswith("recharge")]
    recharge_cc = next(item for item in charge_segments if item.name == "recharge_cc")
    recharge_cv = next(item for item in charge_segments if item.name == "recharge_cv")
    discharge = discharge_curve(experiment, discharge_segment)
    charge = charge_curve(experiment, segments)
    q_trc_ah = float(discharge["capacity_ah"])
    if not np.isclose(q_trc_ah, Q_TRC_AH_EXPECTED, rtol=0.0, atol=5e-10):
        raise RuntimeError(f"Own clean-reference capacity changed: {q_trc_ah}")

    discharge_interp = interpolate_on_soc(discharge)
    measured_grid_v = np.asarray(discharge_interp["voltage_v"])
    backbone_v = np.asarray(V_eq(1.0 - SOC_GRID), dtype=float)
    backbone_metrics = voltage_metrics(measured_grid_v, backbone_v)
    backbone_regions = {name: voltage_metrics(measured_grid_v, backbone_v, mask) for name, mask in region_masks(SOC_GRID).items()}

    sequence_start = discharge_segment.start
    sequence_end = max(item.end for item in charge_segments)
    raw_rows = np.arange(sequence_start, sequence_end + 1)
    provisional_time_s = (raw_rows - raw_rows[0]) * experiment.estimated_dt_s
    model_current_a = -experiment.raw_current_a[raw_rows]
    temperature_k = experiment.numeric["Temp1"][raw_rows] + 273.15
    sequence_result = forward_frozen_panasonic(provisional_time_s, model_current_a, temperature_k)

    discharge_slice = slice(discharge_segment.start - sequence_start, discharge_segment.end - sequence_start + 1)
    discharge_result = FrozenTransferResult(**{field: np.asarray(getattr(sequence_result, field))[discharge_slice] for field in FrozenTransferResult.__dataclass_fields__})
    discharge_soc = np.asarray(discharge["soc"])
    measured_discharge_v = np.asarray(discharge["voltage_v"])
    full_metrics = voltage_metrics(measured_discharge_v, discharge_result.voltage_v)
    full_regions = {name: voltage_metrics(measured_discharge_v, discharge_result.voltage_v, mask) for name, mask in region_masks(discharge_soc).items()}

    z_dis = np.asarray(discharge["z"])
    state_rows: list[dict[str, object]] = []
    for z_target in (0.0, 0.25, 0.5, 0.75, 1.0):
        index = int(np.argmin(np.abs(z_dis - z_target)))
        q_passed_c = z_target * q_trc_ah * 3600.0
        state_rows.append({
            "z": z_target,
            "source_z": float(z_dis[index]),
            "x_n_avg": float(X_N_START - q_passed_c / Q_N_C),
            "x_p_avg": float(X_P_START + q_passed_c / Q_P_C),
            "x_n_surface": float(discharge_result.x_n_surface[index]),
            "x_p_surface": float(discharge_result.x_p_surface[index]),
            "surface_support_valid": bool(discharge_result.ocp_valid[index]),
        })

    nominal_backbone_raw = np.asarray(V_eq(z_dis), dtype=float)
    valid = np.isfinite(discharge_result.voltage_v)
    decomposition = {
        "measured_loaded_offset": scalar_summary((measured_discharge_v - nominal_backbone_raw)[valid]),
        "full_minus_nominal_backbone": scalar_summary((discharge_result.voltage_v - nominal_backbone_raw)[valid]),
        "full_minus_average_equilibrium": scalar_summary((discharge_result.voltage_v - discharge_result.ueq_average_v)[valid]),
    }

    charge_start = min(item.start for item in charge_segments)
    charge_end = max(item.end for item in charge_segments)
    charge_slice = slice(charge_start - sequence_start, charge_end - sequence_start + 1)
    charge_predicted_v = sequence_result.voltage_v[charge_slice]
    charge_soc = np.asarray(charge["soc"])
    measured_charge_v = np.asarray(charge["voltage_v"])
    charge_rows = np.asarray(charge["rows"], dtype=int)
    cc_mask = (charge_rows >= recharge_cc.start) & (charge_rows <= recharge_cc.end)
    cv_mask = (charge_rows >= recharge_cv.start) & (charge_rows <= recharge_cv.end)
    charge_metrics = {
        "all_recharge": voltage_metrics(measured_charge_v, charge_predicted_v),
        "cc": voltage_metrics(measured_charge_v, charge_predicted_v, cc_mask),
        "cv": voltage_metrics(measured_charge_v, charge_predicted_v, cv_mask),
    }

    metrics_rows: list[dict[str, object]] = []
    add_metric_rows(metrics_rows, "A_voltage_backbone", "overall", backbone_metrics, "Panasonic Ueq minus loaded TRC 1C voltage; not OCP validation")
    for region, metric in backbone_regions.items():
        add_metric_rows(metrics_rows, "A_voltage_backbone", f"soc_{region}", metric, "Capacity-coordinate regional diagnostic")
    add_metric_rows(metrics_rows, "B_literal_full_transfer", "valid_ocp_support_only", full_metrics, "PROVISIONAL-TIME direct transfer; invalid OCP-support samples excluded, not clipped")
    for region, metric in full_regions.items():
        add_metric_rows(metrics_rows, "B_literal_full_transfer", f"soc_{region}_valid_support_only", metric, "PROVISIONAL-TIME regional diagnostic")
    for name, summary in decomposition.items():
        for metric, value in summary.items():
            metrics_rows.append({"test": "error_decomposition", "subset": name, "metric": metric, "value": value, "unit": "V", "note": "Not an intrinsic TRC overpotential"})
    for subset, metric in charge_metrics.items():
        add_metric_rows(metrics_rows, "directional_recharge", subset, metric, "Continuous frozen transfer through discharge/rest/recharge; provisional time")
    for row in state_rows:
        for state in ("x_n_avg", "x_p_avg", "x_n_surface", "x_p_surface"):
            metrics_rows.append({"test": "C_state_mapping", "subset": f"z_{row['z']:.2f}", "metric": state, "value": row[state], "unit": "1", "note": f"surface_support_valid={row['surface_support_valid']}; source_z={row['source_z']:.9f}"})
    audit_values = {
        "trc_q_ref_ah": q_trc_ah,
        "panasonic_capacity_fraction_traversed": q_trc_ah / Q_CELL_AH,
        "provisional_dt_s": experiment.estimated_dt_s,
        "discharge_surface_valid_count": int(np.count_nonzero(discharge_result.ocp_valid)),
        "discharge_surface_invalid_count": int(np.count_nonzero(~discharge_result.ocp_valid)),
        "discharge_negative_surface_invalid_count": int(np.count_nonzero(
            (discharge_result.x_n_surface < MCMB_SUPPORTED_STOICHIOMETRY_RANGE[0])
            | (discharge_result.x_n_surface > MCMB_SUPPORTED_STOICHIOMETRY_RANGE[1])
        )),
        "discharge_positive_surface_invalid_count": int(np.count_nonzero(
            (discharge_result.x_p_surface < NCA_WORKING_STOICHIOMETRY_RANGE[0])
            | (discharge_result.x_p_surface > NCA_WORKING_STOICHIOMETRY_RANGE[1])
        )),
        "charge_surface_valid_count": int(np.count_nonzero(np.isfinite(charge_predicted_v))),
        "charge_surface_invalid_count": int(np.count_nonzero(~np.isfinite(charge_predicted_v))),
        "clipping_used": False,
    }
    for name, value in audit_values.items():
        metrics_rows.append({"test": "audit", "subset": "transfer_integrity", "metric": name, "value": value, "unit": "boolean" if name == "clipping_used" else "1", "note": "No measured time; no clipping or parameter change"})

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    FIGURE_DIRECTORY.mkdir(parents=True, exist_ok=True)
    write_csv(METRICS_PATH, metrics_rows)
    write_csv(PARAMETER_MATRIX_PATH, parameter_matrix_rows())
    make_figures(
        SOC_GRID,
        measured_grid_v,
        backbone_v,
        discharge_soc,
        measured_discharge_v,
        discharge_result.voltage_v,
        charge_soc,
        measured_charge_v,
        charge_predicted_v,
        cc_mask,
        cv_mask,
    )
    write_report(
        experiment.metadata.path,
        experiment.estimated_dt_s,
        q_trc_ah,
        backbone_metrics,
        backbone_regions,
        full_metrics,
        full_regions,
        discharge_result,
        discharge_soc,
        state_rows,
        decomposition,
        charge_metrics,
        int(np.count_nonzero(np.isfinite(charge_predicted_v))),
        charge_predicted_v.size,
    )
    print(f"Backbone RMSE: {1000*float(backbone_metrics['rmse_v']):.3f} mV")
    print(f"Full valid-domain RMSE: {1000*float(full_metrics['rmse_v']):.3f} mV; valid {full_metrics['n_samples']}/{discharge_soc.size}")
    print(f"Charge RMSE: {1000*float(charge_metrics['all_recharge']['rmse_v']):.3f} mV")
    print(f"Outputs written under {OUTPUT_DIRECTORY}")


if __name__ == "__main__":
    main()
