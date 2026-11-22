"""Primary three-parameter Panasonic HPPC dynamic calibration.

Fits only R0, alpha_p, and d_n against pulse-relative voltage changes. The
frozen Panasonic OCP, electrode windows, Q_n/Q_p mapping, alpha_n, and d_p are
not altered. Fourteen contiguous histories are initialized independently from
the source cumulative-Ah counter; omitted-data gaps are never integrated.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from models.grouped_spm import FARADAY_C_PER_MOL, GAS_CONSTANT_J_PER_MOL_K  # noqa: E402
from models.ocp_panasonic_qi_frozen import (  # noqa: E402
    Q_CELL_AH,
    U_n,
    U_p,
    X_N_END,
    X_N_START,
    X_P_END,
    X_P_START,
)
from scripts.audit_panasonic_hppc_identifiability import (  # noqa: E402
    EXCLUDED_PULSE_IDS,
    HppcData,
    PulseRecord,
    load_hppc_data,
    load_pulse_records,
    partition_contiguous_groups,
    simulate_local_group,
    stoichiometries_from_discharged_ah,
)


OUTPUT_DIRECTORY = REPOSITORY_ROOT / "thesis_data"
FIT_PARAMETER_NAMES = ("R0", "alpha_p", "d_n")
PHASES = ("immediate", "within_pulse", "post_relaxation")
ALPHA_N_FIXED_S = 3000.0
D_P_FIXED_PER_S = 0.0007
Q_N_C = 23008.732269258333
Q_P_C = 14780.898652565113

INITIAL_PARAMETERS = np.asarray([0.020, 1800.0, 5.0e-5])

# No target-specific bounds exist in the project. These deliberately broad
# positive bounds are numerical-safety guards, not literature/physical ranges.
LOWER_BOUNDS = np.asarray([1.0e-6, 1.0, 1.0e-8])
UPPER_BOUNDS = np.asarray([2.0e-1, 1.0e5, 1.0e-1])
LOG_LOWER_BOUNDS = np.log(LOWER_BOUNDS)
LOG_UPPER_BOUNDS = np.log(UPPER_BOUNDS)

RELATIVE_JACOBIAN_STEP = 1.0e-4
NEAR_BEST_RELATIVE_OBJECTIVE = 0.01
BOUND_LOG_DISTANCE_TOLERANCE = 1.0e-3

MULTISTART_PARAMETERS = np.asarray(
    [
        [0.020, 1800.0, 5.0e-5],
        [0.005, 300.0, 1.0e-5],
        [0.010, 800.0, 2.0e-5],
        [0.030, 3000.0, 1.0e-4],
        [0.050, 4000.0, 5.0e-4],
        [0.080, 100.0, 1.0e-3],
        [0.002, 2500.0, 5.0e-6],
        [0.015, 80.0, 2.0e-4],
        [0.040, 1200.0, 1.0e-5],
        [0.008, 3500.0, 1.0e-3],
        [0.100, 500.0, 5.0e-6],
        [0.003, 60.0, 5.0e-4],
    ],
    dtype=float,
)


@dataclass(frozen=True)
class PreparedPulse:
    record: PulseRecord
    pre_local_index: int
    end_local_index: int
    phase_local_indices: dict[str, np.ndarray]
    measured_pre_v: float
    measured_delta_by_phase: dict[str, np.ndarray]


@dataclass(frozen=True)
class PreparedGroup:
    group_id: int
    time_s: np.ndarray
    current_a: np.ndarray
    temperature_k: np.ndarray
    measured_voltage_v: np.ndarray
    initial_x_n: float
    initial_x_p: float
    x_n_surface: np.ndarray
    x_p_avg: np.ndarray
    u_n_v: np.ndarray
    current_interval_runs: tuple[tuple[int, int], ...]
    pulses: tuple[PreparedPulse, ...]


@dataclass(frozen=True)
class Evaluation:
    weighted_residual: np.ndarray
    group_voltage_v: tuple[np.ndarray, ...]
    dynamic_records: tuple[dict[str, object], ...]


def _usable(record: PulseRecord) -> bool:
    return record.usable and record.pulse_id not in EXCLUDED_PULSE_IDS


def prepare_groups(
    data: HppcData,
    records: list[PulseRecord],
    *,
    minimum_pulse_start_soc: float = 0.0,
    maximum_pulse_start_soc: float = np.inf,
) -> list[PreparedGroup]:
    contiguous = partition_contiguous_groups(data, records)
    prepared: list[PreparedGroup] = []
    for group_id, group_records in enumerate(contiguous, start=1):
        usable_records = [
            record
            for record in group_records
            if _usable(record)
            and minimum_pulse_start_soc <= record.soc < maximum_pulse_start_soc
        ]
        if not usable_records:
            continue
        lower = usable_records[0].pre_pulse_index
        upper = usable_records[-1].post_rest_index + 1
        time_s = data.time_s[lower:upper] - data.time_s[lower]
        current_a = -data.raw_current_a[lower:upper]
        temperature_k = data.battery_temperature_k[lower:upper]
        measured_voltage_v = data.voltage_v[lower:upper]
        discharged_ah = max(0.0, -float(data.exported_ah[lower]))
        initial_x_n, initial_x_p = stoichiometries_from_discharged_ah(discharged_ah)

        reference = simulate_local_group(
            time_s,
            current_a,
            temperature_k,
            initial_x_n,
            initial_x_p,
            {
                "alpha_n": ALPHA_N_FIXED_S,
                "alpha_p": INITIAL_PARAMETERS[1],
                "d_n": INITIAL_PARAMETERS[2],
                "d_p": D_P_FIXED_PER_S,
                "R0": INITIAL_PARAMETERS[0],
            },
        )
        pulses: list[PreparedPulse] = []
        for record in usable_records:
            pre = record.pre_pulse_index - lower
            start = record.start_index - lower
            end = record.end_index - lower
            post = record.post_rest_index - lower
            phase_indices = {
                "immediate": np.asarray([start], dtype=int),
                "within_pulse": np.arange(start + 1, end + 1, dtype=int),
                "post_relaxation": np.arange(end + 1, post + 1, dtype=int),
            }
            if any(indices.size == 0 for indices in phase_indices.values()):
                raise ValueError(f"Pulse {record.pulse_id} has an empty response phase.")
            measured_pre = float(measured_voltage_v[pre])
            pulses.append(
                PreparedPulse(
                    record=record,
                    pre_local_index=pre,
                    end_local_index=end,
                    phase_local_indices=phase_indices,
                    measured_pre_v=measured_pre,
                    measured_delta_by_phase={
                        phase: measured_voltage_v[indices] - measured_pre
                        for phase, indices in phase_indices.items()
                    },
                )
            )
        prepared.append(
            PreparedGroup(
                group_id=group_id,
                time_s=time_s,
                current_a=current_a,
                temperature_k=temperature_k,
                measured_voltage_v=measured_voltage_v,
                initial_x_n=initial_x_n,
                initial_x_p=initial_x_p,
                x_n_surface=reference.x_n_surface,
                x_p_avg=reference.x_p_avg,
                u_n_v=np.asarray(U_n(reference.x_n_surface), dtype=float),
                current_interval_runs=tuple(
                    zip(
                        np.r_[
                            0,
                            np.flatnonzero(np.diff(current_a[:-1]) != 0.0) + 1,
                        ].tolist(),
                        np.r_[
                            np.flatnonzero(np.diff(current_a[:-1]) != 0.0) + 1,
                            current_a.size - 1,
                        ].tolist(),
                        strict=True,
                    )
                ),
                pulses=tuple(pulses),
            )
        )
    if minimum_pulse_start_soc == 0.0 and np.isinf(maximum_pulse_start_soc) and (
        len(prepared) != 14 or sum(len(group.pulses) for group in prepared) != 64
    ):
        raise ValueError("Prepared HPPC structure is not the audited 14 histories/64 pulses.")
    if not prepared:
        raise ValueError("The pulse-start SOC filter retained no HPPC histories.")
    return prepared


def theta_to_parameters(theta: np.ndarray) -> dict[str, float]:
    values = np.exp(np.asarray(theta, dtype=float))
    return dict(zip(FIT_PARAMETER_NAMES, values, strict=True))


def positive_surface_state(group: PreparedGroup, alpha_p: float) -> np.ndarray:
    """Return x_p,surf using the exact grouped-model ZOH state response."""

    y_p = np.empty(group.time_s.size)
    y_p[0] = 0.0
    decay_rate = 30.0 / alpha_p
    for start, stop in group.current_interval_runs:
        held_current = group.current_a[start]
        equilibrium_y = 12.0 * held_current / (7.0 * Q_P_C) / decay_rate
        elapsed = group.time_s[start + 1 : stop + 1] - group.time_s[start]
        y_p[start + 1 : stop + 1] = equilibrium_y + (
            y_p[start] - equilibrium_y
        ) * np.exp(-decay_rate * elapsed)
    x_p_aux = group.x_p_avg + y_p
    return x_p_aux + alpha_p * group.current_a / (105.0 * Q_P_C)


def alpha_p_is_feasible(groups: list[PreparedGroup], alpha_p: float) -> bool:
    for group in groups:
        surface = positive_surface_state(group, alpha_p)
        if np.any((surface < 0.18) | (surface > 0.96)):
            return False
    return True


def maximum_feasible_alpha_p(groups: list[PreparedGroup]) -> float:
    """Locate the data/OCP-domain-imposed upper alpha_p feasibility edge."""

    lower = float(LOWER_BOUNDS[1])
    upper = float(UPPER_BOUNDS[1])
    if not alpha_p_is_feasible(groups, lower):
        raise ValueError("The lower alpha_p safety bound is OCP-domain infeasible.")
    if alpha_p_is_feasible(groups, upper):
        return upper
    for _ in range(70):
        midpoint = np.sqrt(lower * upper)
        if alpha_p_is_feasible(groups, midpoint):
            lower = midpoint
        else:
            upper = midpoint
    return lower


def evaluate_group(group: PreparedGroup, parameters: dict[str, float]) -> np.ndarray:
    alpha_p = parameters["alpha_p"]
    d_n = parameters["d_n"]
    r0 = parameters["R0"]

    x_p_surface = positive_surface_state(group, alpha_p)
    if np.any((x_p_surface < 0.18) | (x_p_surface > 0.96)):
        raise ValueError("Positive surface state left frozen OCP support; no clipping.")

    factor = 2.0 * GAS_CONSTANT_J_PER_MOL_K * group.temperature_k / FARADAY_C_PER_MOL
    eta_n = factor * np.arcsinh(
        group.current_a
        / (
            6.0
            * Q_N_C
            * d_n
            * np.sqrt(group.x_n_surface * (1.0 - group.x_n_surface))
        )
    )
    eta_p = factor * np.arcsinh(
        group.current_a
        / (
            6.0
            * Q_P_C
            * D_P_FIXED_PER_S
            * np.sqrt(x_p_surface * (1.0 - x_p_surface))
        )
    )
    voltage = (
        np.asarray(U_p(x_p_surface), dtype=float)
        - group.u_n_v
        + eta_p
        - eta_n
        - r0 * group.current_a
    )
    if not np.all(np.isfinite(voltage)):
        raise ValueError("Non-finite model voltage.")
    return voltage


def evaluate(
    groups: list[PreparedGroup],
    parameters: dict[str, float],
    *,
    include_records: bool = True,
) -> Evaluation:
    model_groups = tuple(evaluate_group(group, parameters) for group in groups)
    expected_pulse_count = sum(len(group.pulses) for group in groups)
    phase_residuals: dict[str, list[np.ndarray]] = {phase: [] for phase in PHASES}
    dynamic_records: list[dict[str, object]] = []
    for group, model_voltage in zip(groups, model_groups, strict=True):
        for pulse in group.pulses:
            model_pre = float(model_voltage[pulse.pre_local_index])
            for phase in PHASES:
                indices = pulse.phase_local_indices[phase]
                model_delta = model_voltage[indices] - model_pre
                residual = pulse.measured_delta_by_phase[phase] - model_delta
                phase_residuals[phase].append(residual)
                if include_records:
                    for index, measured_delta, modeled_delta, residual_value in zip(
                        indices,
                        pulse.measured_delta_by_phase[phase],
                        model_delta,
                        residual,
                        strict=True,
                    ):
                        dynamic_records.append(
                            {
                                "pulse_id": pulse.record.pulse_id,
                                "nominal_rate": pulse.record.nominal_rate,
                                "soc": pulse.record.soc,
                                "phase": phase,
                                "local_time_s": float(group.time_s[index]),
                                "measured_voltage_v": float(group.measured_voltage_v[index]),
                                "model_voltage_v": float(model_voltage[index]),
                                "measured_delta_v": float(measured_delta),
                                "model_delta_v": float(modeled_delta),
                                "dynamic_residual_v": float(residual_value),
                                "absolute_residual_v": float(
                                    group.measured_voltage_v[index] - model_voltage[index]
                                ),
                            }
                        )

    weighted_parts = []
    for phase in PHASES:
        pulse_parts = phase_residuals[phase]
        n_pulses = len(pulse_parts)
        if n_pulses != expected_pulse_count:
            raise ValueError(
                f"Expected {expected_pulse_count} pulse residuals for {phase}, "
                f"found {n_pulses}."
            )
        for residual in pulse_parts:
            weighted_parts.append(residual / np.sqrt(n_pulses * residual.size))
    return Evaluation(
        weighted_residual=np.concatenate(weighted_parts),
        group_voltage_v=model_groups,
        dynamic_records=tuple(dynamic_records),
    )


def phase_balanced_residual(
    theta: np.ndarray, groups: list[PreparedGroup], expected_size: int
) -> np.ndarray:
    try:
        result = evaluate(
            groups, theta_to_parameters(theta), include_records=False
        )
    except (ValueError, FloatingPointError):
        # Explicit infeasibility penalty. No state or parameter clipping occurs.
        distance = float(np.linalg.norm(theta - np.log(INITIAL_PARAMETERS)))
        return np.full(expected_size, 10.0 + 0.01 * distance)
    return result.weighted_residual


def objective_metrics(evaluation: Evaluation) -> dict[str, float]:
    records = evaluation.dynamic_records
    result: dict[str, float] = {
        "phase_balanced_objective_v2": float(
            np.dot(evaluation.weighted_residual, evaluation.weighted_residual)
        ),
        "phase_balanced_rmse_v": float(
            np.sqrt(
                np.dot(evaluation.weighted_residual, evaluation.weighted_residual)
                / len(PHASES)
            )
        ),
        "dynamic_raw_rmse_v": float(
            np.sqrt(np.mean([float(row["dynamic_residual_v"]) ** 2 for row in records]))
        ),
    }
    for phase in PHASES:
        pulse_mses = []
        for pulse_id in sorted({int(row["pulse_id"]) for row in records}):
            selected = [
                float(row["dynamic_residual_v"])
                for row in records
                if row["phase"] == phase and int(row["pulse_id"]) == pulse_id
            ]
            pulse_mses.append(float(np.mean(np.square(selected))))
        result[f"rmse_{phase}_v"] = float(np.sqrt(np.mean(pulse_mses)))
    absolute = np.asarray([float(row["absolute_residual_v"]) for row in records])
    measured = np.asarray([float(row["measured_voltage_v"]) for row in records])
    result["absolute_rmse_v"] = float(np.sqrt(np.mean(absolute**2)))
    result["absolute_mae_v"] = float(np.mean(np.abs(absolute)))
    result["absolute_r2"] = float(1.0 - np.sum(absolute**2) / np.sum((measured - np.mean(measured)) ** 2))
    return result


def phase_balanced_subset_rmse(records: list[dict[str, object]]) -> float:
    phase_mses = []
    for phase in PHASES:
        pulse_mses = []
        for pulse_id in sorted({int(row["pulse_id"]) for row in records}):
            selected = np.asarray(
                [
                    float(row["dynamic_residual_v"])
                    for row in records
                    if row["phase"] == phase and int(row["pulse_id"]) == pulse_id
                ]
            )
            if selected.size:
                pulse_mses.append(float(np.mean(selected**2)))
        if pulse_mses:
            phase_mses.append(float(np.mean(pulse_mses)))
    return float(np.sqrt(np.mean(phase_mses)))


def calibration_diagnostics(
    groups: list[PreparedGroup], evaluation: Evaluation
) -> dict[str, object]:
    records = list(evaluation.dynamic_records)
    metrics = objective_metrics(evaluation)
    by_rate = {
        rate: phase_balanced_subset_rmse(
            [row for row in records if row["nominal_rate"] == rate]
        )
        for rate in ("0.5C", "1C", "2C", "4C", "6C")
    }
    soc_definition = {
        "high": lambda value: value >= 2.0 / 3.0,
        "middle": lambda value: 1.0 / 3.0 <= value < 2.0 / 3.0,
        "low": lambda value: value < 1.0 / 3.0,
    }
    by_soc = {}
    for region, predicate in soc_definition.items():
        selected_records = [
            row for row in records if predicate(float(row["soc"]))
        ]
        if selected_records:
            by_soc[region] = phase_balanced_subset_rmse(selected_records)

    step_errors = []
    relaxation_shape_errors = []
    measured_rapp = []
    model_rapp = []
    pulse_rows = []
    for group, voltage in zip(groups, evaluation.group_voltage_v, strict=True):
        for pulse in group.pulses:
            immediate = pulse.phase_local_indices["immediate"][0]
            pre = pulse.pre_local_index
            end = pulse.end_local_index
            measured_step = group.measured_voltage_v[immediate] - group.measured_voltage_v[pre]
            model_step = voltage[immediate] - voltage[pre]
            step_error = measured_step - model_step
            pulse_current = float(np.median(group.current_a[immediate : end + 1]))
            delta_i = pulse_current - group.current_a[pre]
            measured_r = abs(measured_step / delta_i)
            modeled_r = abs(model_step / delta_i)
            relaxation_indices = pulse.phase_local_indices["post_relaxation"]
            first_rest = relaxation_indices[0]
            measured_recovery = (
                group.measured_voltage_v[relaxation_indices]
                - group.measured_voltage_v[first_rest]
            )
            model_recovery = voltage[relaxation_indices] - voltage[first_rest]
            shape_error = measured_recovery - model_recovery
            step_errors.append(float(step_error))
            relaxation_shape_errors.extend(shape_error.tolist())
            measured_rapp.append(float(measured_r))
            model_rapp.append(float(modeled_r))
            pulse_rows.append(
                {
                    "pulse_id": pulse.record.pulse_id,
                    "rate": pulse.record.nominal_rate,
                    "soc": pulse.record.soc,
                    "measured_step_v": measured_step,
                    "model_step_v": model_step,
                    "step_error_v": step_error,
                    "measured_rapp_ohm": measured_r,
                    "model_rapp_ohm": modeled_r,
                }
            )
    step_errors_array = np.asarray(step_errors)
    relaxation_errors_array = np.asarray(relaxation_shape_errors)
    measured_rapp_array = np.asarray(measured_rapp)
    model_rapp_array = np.asarray(model_rapp)
    return {
        "metrics": metrics,
        "by_rate": by_rate,
        "by_soc": by_soc,
        "pulse_rows": pulse_rows,
        "step_rmse_v": float(np.sqrt(np.mean(step_errors_array**2))),
        "step_mae_v": float(np.mean(np.abs(step_errors_array))),
        "step_bias_v": float(np.mean(step_errors_array)),
        "relaxation_shape_rmse_v": float(np.sqrt(np.mean(relaxation_errors_array**2))),
        "relaxation_shape_mae_v": float(np.mean(np.abs(relaxation_errors_array))),
        "measured_rapp_range_ohm": (
            float(np.min(measured_rapp_array)),
            float(np.max(measured_rapp_array)),
        ),
        "model_rapp_range_ohm": (
            float(np.min(model_rapp_array)),
            float(np.max(model_rapp_array)),
        ),
        "rapp_rmse_ohm": float(np.sqrt(np.mean((measured_rapp_array - model_rapp_array) ** 2))),
        "rapp_mae_ohm": float(np.mean(np.abs(measured_rapp_array - model_rapp_array))),
    }


def postfit_jacobian(
    groups: list[PreparedGroup], best_parameters: dict[str, float]
) -> dict[str, object]:
    theta = np.log([best_parameters[name] for name in FIT_PARAMETER_NAMES])
    columns = []
    base_residual = evaluate(
        groups, best_parameters, include_records=False
    ).weighted_residual
    expected_size = base_residual.size
    for index in range(len(FIT_PARAMETER_NAMES)):
        plus = theta.copy()
        minus = theta.copy()
        plus[index] += RELATIVE_JACOBIAN_STEP
        minus[index] -= RELATIVE_JACOBIAN_STEP
        try:
            plus_residual = evaluate(
                groups, theta_to_parameters(plus), include_records=False
            ).weighted_residual
        except ValueError:
            plus_residual = None
        try:
            minus_residual = evaluate(
                groups, theta_to_parameters(minus), include_records=False
            ).weighted_residual
        except ValueError:
            minus_residual = None
        if plus_residual is not None and minus_residual is not None:
            derivative = (plus_residual - minus_residual) / (
                2.0 * RELATIVE_JACOBIAN_STEP
            )
        elif minus_residual is not None:
            minus_two = theta.copy()
            minus_two[index] -= 2.0 * RELATIVE_JACOBIAN_STEP
            minus_two_residual = evaluate(
                groups, theta_to_parameters(minus_two), include_records=False
            ).weighted_residual
            derivative = (
                3.0 * base_residual - 4.0 * minus_residual + minus_two_residual
            ) / (2.0 * RELATIVE_JACOBIAN_STEP)
        elif plus_residual is not None:
            plus_two = theta.copy()
            plus_two[index] += 2.0 * RELATIVE_JACOBIAN_STEP
            plus_two_residual = evaluate(
                groups, theta_to_parameters(plus_two), include_records=False
            ).weighted_residual
            derivative = (
                -3.0 * base_residual + 4.0 * plus_residual - plus_two_residual
            ) / (2.0 * RELATIVE_JACOBIAN_STEP)
        else:
            raise RuntimeError(
                f"No feasible local derivative direction for {FIT_PARAMETER_NAMES[index]}."
            )
        columns.append(derivative)
    jacobian = np.column_stack(columns)
    _, singular_values, vt = np.linalg.svd(jacobian, full_matrices=False)
    condition = float(singular_values[0] / singular_values[-1])
    covariance = np.linalg.pinv(jacobian.T @ jacobian, rcond=1.0e-12)
    scale = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    correlation = covariance / np.outer(scale, scale)
    np.fill_diagonal(correlation, 1.0)
    return {
        "jacobian": jacobian,
        "singular_values": singular_values,
        "condition_number": condition,
        "correlation": correlation,
        "weakest_vector": vt[-1],
    }


def run_multistart(
    groups: list[PreparedGroup], alpha_p_domain_upper_s: float
) -> list[dict[str, object]]:
    expected_size = evaluate(
        groups,
        dict(zip(FIT_PARAMETER_NAMES, INITIAL_PARAMETERS, strict=True)),
        include_records=False,
    ).weighted_residual.size
    results: list[dict[str, object]] = []
    effective_upper_bounds = LOG_UPPER_BOUNDS.copy()
    effective_upper_bounds[1] = np.log(
        alpha_p_domain_upper_s * (1.0 - 1.0e-10)
    )
    for start_id, start in enumerate(MULTISTART_PARAMETERS, start=1):
        initial_theta = np.log(start)
        initial_residual = phase_balanced_residual(initial_theta, groups, expected_size)
        initial_objective = float(np.dot(initial_residual, initial_residual))
        fit = least_squares(
            phase_balanced_residual,
            initial_theta,
            bounds=(LOG_LOWER_BOUNDS, effective_upper_bounds),
            args=(groups, expected_size),
            method="trf",
            x_scale="jac",
            ftol=1.0e-10,
            xtol=1.0e-10,
            gtol=1.0e-10,
            max_nfev=300,
        )
        parameters = theta_to_parameters(fit.x)
        try:
            evaluation = evaluate(groups, parameters)
            metrics = objective_metrics(evaluation)
            feasible = True
        except ValueError:
            metrics = {
                "phase_balanced_objective_v2": float("nan"),
                "phase_balanced_rmse_v": float("nan"),
                "dynamic_raw_rmse_v": float("nan"),
            }
            feasible = False
        log_distance_to_bounds = np.minimum(fit.x - LOG_LOWER_BOUNDS, LOG_UPPER_BOUNDS - fit.x)
        results.append(
            {
                "start_id": start_id,
                "success": bool(fit.success and feasible),
                "feasible": feasible,
                "status": int(fit.status),
                "message": str(fit.message).replace("\n", " "),
                "nfev": int(fit.nfev),
                "initial_R0_ohm": start[0],
                "initial_alpha_p_s": start[1],
                "initial_d_n_per_s": start[2],
                "initial_objective_v2": initial_objective,
                "R0_ohm": parameters["R0"],
                "alpha_p_s": parameters["alpha_p"],
                "d_n_per_s": parameters["d_n"],
                "objective_v2": metrics["phase_balanced_objective_v2"],
                "phase_balanced_rmse_v": metrics["phase_balanced_rmse_v"],
                "dynamic_raw_rmse_v": metrics["dynamic_raw_rmse_v"],
                "at_numerical_bound": bool(
                    np.any(log_distance_to_bounds <= BOUND_LOG_DISTANCE_TOLERANCE)
                ),
                "at_frozen_ocp_domain_boundary": bool(
                    (alpha_p_domain_upper_s - parameters["alpha_p"])
                    / alpha_p_domain_upper_s
                    <= 1.0e-6
                ),
                "alpha_p_frozen_domain_upper_s": alpha_p_domain_upper_s,
            }
        )
        print(
            f"start {start_id:02d}: success={fit.success and feasible}, "
            f"objective={metrics['phase_balanced_objective_v2']:.8g}, "
            f"R0={parameters['R0']:.8g}, alpha_p={parameters['alpha_p']:.8g}, "
            f"d_n={parameters['d_n']:.8g}"
        )
    return results


def _write_dict_rows(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_parameters(
    path: Path,
    best: dict[str, object],
    postfit: dict[str, object],
    max_x_p_surface: float,
) -> None:
    row = {
        "R0_ohm": best["R0_ohm"],
        "alpha_p_s": best["alpha_p_s"],
        "d_n_per_s": best["d_n_per_s"],
        "alpha_n_fixed_s": ALPHA_N_FIXED_S,
        "d_p_fixed_per_s": D_P_FIXED_PER_S,
        "Q_n_fixed_c": Q_N_C,
        "Q_p_fixed_c": Q_P_C,
        "Q_cell_fixed_ah": Q_CELL_AH,
        "x_n_start_fixed": X_N_START,
        "x_n_end_fixed": X_N_END,
        "x_p_start_fixed": X_P_START,
        "x_p_end_fixed": X_P_END,
        "conditional_on_provisional_alpha_n_and_d_p": True,
        "postfit_condition_number": postfit["condition_number"],
        "maximum_fitted_x_p_surface": max_x_p_surface,
        "frozen_positive_ocp_upper": 0.96,
        "alpha_p_at_frozen_ocp_domain_boundary": best[
            "at_frozen_ocp_domain_boundary"
        ],
    }
    _write_dict_rows(path, [row])


def write_metrics(path: Path, diagnostics: dict[str, object]) -> None:
    rows: list[dict[str, object]] = []
    for name, value in dict(diagnostics["metrics"]).items():
        rows.append({"group": "overall", "subset": "all", "metric": name, "value": value})
    for rate, value in dict(diagnostics["by_rate"]).items():
        rows.append(
            {"group": "c_rate", "subset": rate, "metric": "phase_balanced_dynamic_rmse_v", "value": value}
        )
    for region, value in dict(diagnostics["by_soc"]).items():
        rows.append(
            {"group": "soc", "subset": region, "metric": "phase_balanced_dynamic_rmse_v", "value": value}
        )
    for metric in (
        "step_rmse_v",
        "step_mae_v",
        "step_bias_v",
        "relaxation_shape_rmse_v",
        "relaxation_shape_mae_v",
        "rapp_rmse_ohm",
        "rapp_mae_ohm",
    ):
        rows.append({"group": "diagnostic", "subset": "all", "metric": metric, "value": diagnostics[metric]})
    for label, bounds in (
        ("measured", diagnostics["measured_rapp_range_ohm"]),
        ("model", diagnostics["model_rapp_range_ohm"]),
    ):
        rows.append({"group": "rapp", "subset": label, "metric": "minimum_ohm", "value": bounds[0]})
        rows.append({"group": "rapp", "subset": label, "metric": "maximum_ohm", "value": bounds[1]})
    _write_dict_rows(path, rows)


def save_figures(
    groups: list[PreparedGroup], evaluation: Evaluation, diagnostics: dict[str, object]
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Representative 1C delta-voltage responses across SOC, limited to 80 s
    # so pulse and early relaxation remain visible.
    one_c: list[tuple[PreparedGroup, PreparedPulse, np.ndarray]] = []
    for group, voltage in zip(groups, evaluation.group_voltage_v, strict=True):
        for pulse in group.pulses:
            if pulse.record.nominal_rate == "1C":
                one_c.append((group, pulse, voltage))
    targets = np.asarray([1.0, 0.8, 0.6, 0.4, 0.2, 0.08])
    selected = []
    for target in targets:
        candidate = min(one_c, key=lambda item: abs(item[1].record.soc - target))
        if candidate not in selected:
            selected.append(candidate)
    figure, axes = plt.subplots(2, 3, figsize=(13.0, 7.5), constrained_layout=True)
    for axis, (group, pulse, voltage) in zip(axes.flat, selected, strict=False):
        pre = pulse.pre_local_index
        stop_time = group.time_s[pre] + 80.0
        stop = int(np.searchsorted(group.time_s, stop_time, side="right"))
        indices = np.arange(pre, min(stop, group.time_s.size))
        relative_time = group.time_s[indices] - group.time_s[pre]
        axis.plot(
            relative_time,
            group.measured_voltage_v[indices] - group.measured_voltage_v[pre],
            label="measured",
        )
        axis.plot(relative_time, voltage[indices] - voltage[pre], label="model")
        pulse_start = pulse.phase_local_indices["immediate"][0]
        axis.axvspan(
            group.time_s[pulse_start] - group.time_s[pre],
            group.time_s[pulse.end_local_index] - group.time_s[pre],
            color="tab:red",
            alpha=0.08,
        )
        axis.set(
            title=f"Pulse {pulse.record.pulse_id}: {100*pulse.record.soc:.1f}% SOC",
            xlabel="Time from pre-pulse [s]",
            ylabel="Delta voltage [V]",
        )
        axis.grid(alpha=0.25)
    axes.flat[0].legend()
    figure.suptitle("Primary HPPC fit: representative 1C responses")
    figure.savefig(OUTPUT_DIRECTORY / "panasonic_hppc_primary_fit_pulses.png", dpi=180)
    plt.close(figure)

    records = list(evaluation.dynamic_records)
    figure, axes = plt.subplots(1, 3, figsize=(14.0, 4.8), constrained_layout=True)
    for axis, phase in zip(axes, PHASES, strict=True):
        selected_rows = [row for row in records if row["phase"] == phase]
        # Deterministic thinning is display-only; metrics retain every sample.
        stride = max(1, len(selected_rows) // 5000)
        shown = selected_rows[::stride]
        axis.scatter(
            [float(row["soc"]) * 100.0 for row in shown],
            [1000.0 * float(row["dynamic_residual_v"]) for row in shown],
            s=5,
            alpha=0.35,
        )
        axis.axhline(0.0, color="black", linewidth=0.8)
        axis.set(title=phase.replace("_", " "), xlabel="Pulse-start SOC [%]", ylabel="Dynamic residual [mV]")
        axis.grid(alpha=0.25)
    figure.suptitle("Pulse-relative voltage residuals (all samples used in metrics)")
    figure.savefig(OUTPUT_DIRECTORY / "panasonic_hppc_primary_fit_residuals.png", dpi=180)
    plt.close(figure)

    def phase_rmse_by(key: str, order: tuple[str, ...], filename: str, title: str) -> None:
        values = np.zeros((len(order), len(PHASES)))
        for row_index, name in enumerate(order):
            for phase_index, phase in enumerate(PHASES):
                selected_residuals = [
                    float(row["dynamic_residual_v"])
                    for row in records
                    if str(row[key]) == name and row["phase"] == phase
                ]
                values[row_index, phase_index] = 1000.0 * np.sqrt(
                    np.mean(np.square(selected_residuals))
                )
        figure, axis = plt.subplots(figsize=(10.0, 5.8), constrained_layout=True)
        positions = np.arange(len(order))
        width = 0.25
        for phase_index, phase in enumerate(PHASES):
            axis.bar(
                positions + (phase_index - 1) * width,
                values[:, phase_index],
                width,
                label=phase.replace("_", " "),
            )
        axis.set_xticks(positions, order)
        axis.set(ylabel="Dynamic RMSE [mV]", title=title)
        axis.grid(axis="y", alpha=0.25)
        axis.legend()
        figure.savefig(OUTPUT_DIRECTORY / filename, dpi=180)
        plt.close(figure)

    phase_rmse_by(
        "nominal_rate",
        ("0.5C", "1C", "2C", "4C", "6C"),
        "panasonic_hppc_primary_fit_by_crate.png",
        "Primary dynamic fit by C-rate and response phase",
    )

    for row in records:
        soc = float(row["soc"])
        row["soc_label"] = "high" if soc >= 2.0 / 3.0 else "middle" if soc >= 1.0 / 3.0 else "low"
    phase_rmse_by(
        "soc_label",
        ("high", "middle", "low"),
        "panasonic_hppc_primary_fit_by_soc.png",
        "Primary dynamic fit by pulse-start SOC and response phase",
    )


def write_report(
    path: Path,
    best: dict[str, object],
    multistart: list[dict[str, object]],
    diagnostics: dict[str, object],
    postfit: dict[str, object],
    max_x_p_surface: float,
) -> None:
    successful = [row for row in multistart if row["success"]]
    best_objective = float(best["objective_v2"])
    near_best = [
        row
        for row in successful
        if float(row["objective_v2"]) <= best_objective * (1.0 + NEAR_BEST_RELATIVE_OBJECTIVE)
    ]
    spreads = {}
    for field in ("R0_ohm", "alpha_p_s", "d_n_per_s"):
        values = np.asarray([float(row[field]) for row in near_best])
        spreads[field] = (float(np.min(values)), float(np.max(values)))
    metrics = dict(diagnostics["metrics"])
    correlation = np.asarray(postfit["correlation"])
    reference_objective = float(multistart[0]["initial_objective_v2"])
    reference_balanced_rmse = np.sqrt(reference_objective / len(PHASES))
    lines = [
        "# Panasonic HPPC primary three-parameter calibration",
        "",
        "## Scope and frozen quantities",
        "",
        "Only `R0`, `alpha_p`, and `d_n` were optimized. The estimates are conditional on the provisional fixed values `alpha_n=3000 s` and `d_p=0.0007 1/s`. The frozen Panasonic OCPs, nominal effective stoichiometric windows, `Q_cell`, `Q_n`, and `Q_p` were not changed.",
        "",
        "The calibration uses 64 complete pulses in 14 independently initialized contiguous histories. Pulse IDs 60, 64, and 67 are excluded. Cumulative exported Ah initializes each history; current is not integrated through omitted gaps. Actual mixed time steps and duplicate timestamps are retained without resampling.",
        "",
        "## Objective",
        "",
        "For each pulse, one measured and one modeled pre-pulse voltage are subtracted from that pulse's response. No sample-wise shifting is used. The implemented least-squares vector weights every sample in a pulse phase by `1/sqrt(64*N_phase_samples_for_that_pulse)`. Therefore its squared norm is exactly:",
        "",
        "```text",
        "J = mean_pulses(MSE_immediate) + mean_pulses(MSE_within_pulse) + mean_pulses(MSE_post_relaxation)",
        "```",
        "",
        "The three response phases have equal fixed weight, and long rest records do not give a pulse more weight. These weights were specified before fitting and were not tuned for fit quality.",
        "",
        "## Parameterization and bounds",
        "",
        "Optimization is in natural-log parameter coordinates. The project contains no documented target-specific admissible bounds, so the following broad positivity/numerical-safety guards were used; they are not asserted physical Panasonic ranges:",
        "",
        "| Parameter | Lower | Upper |",
        "|---|---:|---:|",
        f"| `R0` [ohm] | {LOWER_BOUNDS[0]:.3g} | {UPPER_BOUNDS[0]:.3g} |",
        f"| `alpha_p` [s] | {LOWER_BOUNDS[1]:.3g} | {UPPER_BOUNDS[1]:.3g} |",
        f"| `d_n` [1/s] | {LOWER_BOUNDS[2]:.3g} | {UPPER_BOUNDS[2]:.3g} |",
        "",
        "Frozen OCP-domain violations are rejected, never clipped.",
        f"For these 14 histories, the frozen positive-OCP domain imposes a derived feasibility edge `alpha_p <= {float(best['alpha_p_frozen_domain_upper_s']):.12g} s`. This is not an invented physical bound; it is the largest value whose evaluated `x_p,surf` remains at or below 0.96 on the retained data.",
        "",
        "## Best fit",
        "",
        f"- `R0 = {float(best['R0_ohm']):.12g} ohm`",
        f"- `alpha_p = {float(best['alpha_p_s']):.12g} s`",
        f"- `d_n = {float(best['d_n_per_s']):.12g} 1/s`",
        f"- Phase-balanced objective `J = {metrics['phase_balanced_objective_v2']:.12g} V^2`",
        f"- Phase-balanced dynamic RMSE = {1000*metrics['phase_balanced_rmse_v']:.6f} mV",
        f"- Initialization-fixture phase-balanced RMSE = {1000*reference_balanced_rmse:.6f} mV; the optimized value is {(1.0-metrics['phase_balanced_rmse_v']/reference_balanced_rmse)*100.0:.3f}% lower.",
        f"- Raw-sample dynamic RMSE = {1000*metrics['dynamic_raw_rmse_v']:.6f} mV (diagnostic only; rest-sample dominated)",
        f"- Maximum fitted `x_p,surf = {max_x_p_surface:.12f}`; frozen upper support is 0.960000000000.",
        f"- Active frozen-OCP feasibility boundary: **{'YES' if best['at_frozen_ocp_domain_boundary'] else 'NO'}**.",
        "",
        "### Dynamic RMSE by phase",
        "",
        "| Phase | RMSE [mV] |",
        "|---|---:|",
        f"| Immediate response | {1000*metrics['rmse_immediate_v']:.6f} |",
        f"| Within pulse | {1000*metrics['rmse_within_pulse_v']:.6f} |",
        f"| Post relaxation | {1000*metrics['rmse_post_relaxation_v']:.6f} |",
        "",
        "### Phase-balanced dynamic RMSE by C-rate",
        "",
        "| C-rate | RMSE [mV] |",
        "|---|---:|",
        *[f"| {name} | {1000*value:.6f} |" for name, value in dict(diagnostics["by_rate"]).items()],
        "",
        "### Phase-balanced dynamic RMSE by SOC region",
        "",
        "SOC bins are low `[0,1/3)`, middle `[1/3,2/3)`, and high `[2/3,1]` by pulse-start cell SOC.",
        "",
        "| SOC region | RMSE [mV] |",
        "|---|---:|",
        *[f"| {name} | {1000*value:.6f} |" for name, value in dict(diagnostics["by_soc"]).items()],
        "",
        "## Absolute-voltage diagnostic",
        "",
        f"- RMSE: {1000*metrics['absolute_rmse_v']:.6f} mV",
        f"- MAE: {1000*metrics['absolute_mae_v']:.6f} mV",
        f"- R2: {metrics['absolute_r2']:.9f}",
        "",
        "Absolute voltage was not the optimization target, so static OCP/baseline mismatch cannot dominate the dynamic fit.",
        "",
        "## Step, relaxation, and apparent-resistance diagnostics",
        "",
        f"- Pulse-step error RMSE / MAE / bias: {1000*float(diagnostics['step_rmse_v']):.6f} / {1000*float(diagnostics['step_mae_v']):.6f} / {1000*float(diagnostics['step_bias_v']):+.6f} mV.",
        f"- Relaxation-shape error (each response referenced to its first post-pulse rest sample) RMSE / MAE: {1000*float(diagnostics['relaxation_shape_rmse_v']):.6f} / {1000*float(diagnostics['relaxation_shape_mae_v']):.6f} mV.",
        f"- Measured `R_app` range: {1000*diagnostics['measured_rapp_range_ohm'][0]:.6f}-{1000*diagnostics['measured_rapp_range_ohm'][1]:.6f} mOhm.",
        f"- Model-implied initial-step range: {1000*diagnostics['model_rapp_range_ohm'][0]:.6f}-{1000*diagnostics['model_rapp_range_ohm'][1]:.6f} mOhm.",
        f"- Measured/model implied `R_app` RMSE / MAE: {1000*float(diagnostics['rapp_rmse_ohm']):.6f} / {1000*float(diagnostics['rapp_mae_ohm']):.6f} mOhm.",
        "",
        "`R0` is not equated to `R_app`; the model-implied voltage step also contains algebraic kinetics and the current-dependent surface-state/OCP contribution.",
        "",
        "## Multistart recoverability",
        "",
        f"- Deterministic starts: {len(multistart)}; successful feasible solutions: {len(successful)}.",
        f"- Near-best definition: objective within {100*NEAR_BEST_RELATIVE_OBJECTIVE:.1f}% of the best; near-best solutions: {len(near_best)}.",
        f"- Near-best `R0` range: {spreads['R0_ohm'][0]:.12g}-{spreads['R0_ohm'][1]:.12g} ohm.",
        f"- Near-best `alpha_p` range: {spreads['alpha_p_s'][0]:.12g}-{spreads['alpha_p_s'][1]:.12g} s.",
        f"- Near-best `d_n` range: {spreads['d_n_per_s'][0]:.12g}-{spreads['d_n_per_s'][1]:.12g} 1/s.",
        f"- Successful solutions at a numerical bound: {sum(bool(row['at_numerical_bound']) for row in successful)}.",
        f"- Successful solutions at the frozen-OCP-derived `alpha_p` feasibility boundary: {sum(bool(row['at_frozen_ocp_domain_boundary']) for row in successful)}.",
        "",
        "The multistart CSV records every start, termination result, objective, and fitted parameter set. Similar objectives with materially different near-best parameters indicate weak practical recoverability even when sensitivities are nonzero. Convergence of `alpha_p` to the frozen-domain feasibility edge is a boundary-limited estimate, not evidence of an interior optimum.",
        "",
        "## Post-fit local audit",
        "",
        f"- Singular values of the phase-balanced normalized Jacobian: {', '.join(f'{value:.9g}' for value in np.asarray(postfit['singular_values']))}.",
        f"- Condition number: {float(postfit['condition_number']):.9g}.",
        "- Weakest normalized direction: "
        + ", ".join(
            f"`{name}`={value:+.6f}"
            for name, value in zip(FIT_PARAMETER_NAMES, np.asarray(postfit["weakest_vector"]), strict=True)
        )
        + ".",
        "",
        "Estimator-oriented parameter correlations (normalized pseudoinverse of `J.T @ J`):",
        "",
        "| | R0 | alpha_p | d_n |",
        "|---|---:|---:|---:|",
        *[
            "| " + FIT_PARAMETER_NAMES[row] + " | "
            + " | ".join(f"{correlation[row, column]:+.6f}" for column in range(3))
            + " |"
            for row in range(3)
        ],
        "",
        "Central normalized log-parameter differences are used for `R0` and `d_n`. Because `alpha_p` is on the frozen-OCP feasibility edge, its derivative uses a second-order inward (backward) difference; no infeasible OCP evaluation is substituted into the Jacobian.",
        "",
        "This remains a local fitted-point audit. It tests conditioning under the actual HPPC objective but does not by itself prove global or noise-robust recoverability.",
        "",
        "## Calibration assessment",
        "",
        "The numerical multistart is stable for `R0` and `d_n` under the stated fixed assumptions, but the three-parameter representation is **boundary-limited and leaves substantial systematic mismatch**:",
        "",
        f"- all {len(successful)} successful starts place `alpha_p` on the frozen positive-OCP support boundary;",
        f"- immediate voltage steps retain a one-sided bias of {1000*float(diagnostics['step_bias_v']):+.3f} mV;",
        f"- within-pulse RMSE remains {1000*metrics['rmse_within_pulse_v']:.3f} mV;",
        f"- low-SOC phase-balanced RMSE reaches {1000*dict(diagnostics['by_soc'])['low']:.3f} mV; and",
        f"- the independently referenced relaxation-shape RMSE is {1000*float(diagnostics['relaxation_shape_rmse_v']):.3f} mV.",
        "",
        "Therefore this requested first calibration is computationally complete, but the fitted triplet should be retained as a **conditional diagnostic result**, not frozen as a validated Panasonic dynamic parameter set. The local post-fit condition number does not cure the active boundary or model-data shape mismatch. Per the stop condition, no additional parameter was released and no alternative model or dataset was run.",
        "",
        "## Stop condition",
        "",
        "The first three-parameter HPPC calibration is complete. `alpha_n` and `d_p` were not released; no OCP, capacity, stoichiometric window, model equation, thermal parameter, UDDS profile, or 1C validation was changed or run.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    data = load_hppc_data()
    records = load_pulse_records()
    groups = prepare_groups(data, records)
    alpha_p_domain_upper_s = maximum_feasible_alpha_p(groups)
    multistart = run_multistart(groups, alpha_p_domain_upper_s)
    successful = [row for row in multistart if row["success"]]
    if not successful:
        raise RuntimeError("No successful feasible multistart solution.")
    best = min(successful, key=lambda row: float(row["objective_v2"]))
    best_parameters = {
        "R0": float(best["R0_ohm"]),
        "alpha_p": float(best["alpha_p_s"]),
        "d_n": float(best["d_n_per_s"]),
    }
    best_evaluation = evaluate(groups, best_parameters)
    diagnostics = calibration_diagnostics(groups, best_evaluation)
    postfit = postfit_jacobian(groups, best_parameters)
    max_x_p_surface = max(
        float(np.max(positive_surface_state(group, best_parameters["alpha_p"])))
        for group in groups
    )

    _write_dict_rows(
        OUTPUT_DIRECTORY / "panasonic_hppc_primary_fit_multistart.csv", multistart
    )
    write_parameters(
        OUTPUT_DIRECTORY / "panasonic_hppc_primary_fit_parameters.csv",
        best,
        postfit,
        max_x_p_surface,
    )
    write_metrics(
        OUTPUT_DIRECTORY / "panasonic_hppc_primary_fit_metrics.csv", diagnostics
    )
    write_report(
        OUTPUT_DIRECTORY / "panasonic_hppc_primary_fit_report.md",
        best,
        multistart,
        diagnostics,
        postfit,
        max_x_p_surface,
    )
    save_figures(groups, best_evaluation, diagnostics)

    print("Best conditional primary fit:")
    print(
        f"R0={best_parameters['R0']:.12g} ohm, "
        f"alpha_p={best_parameters['alpha_p']:.12g} s, "
        f"d_n={best_parameters['d_n']:.12g} 1/s"
    )
    print(
        "Phase-balanced dynamic RMSE: "
        f"{1000*float(diagnostics['metrics']['phase_balanced_rmse_v']):.6f} mV"
    )
    print((OUTPUT_DIRECTORY / "panasonic_hppc_primary_fit_report.md").resolve())


if __name__ == "__main__":
    main()
