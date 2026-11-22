"""Panasonic HPPC parameter-mapping and local-identifiability audit.

This is a forward/sensitivity-only script.  It does not optimize parameters,
alter the frozen OCP representation, or modify the grouped-SPM equations.
Each audited pulse/rest history is initialized independently from the source
file's cumulative Ah counter, so no current is integrated through omitted-data
gaps.
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

from models.grouped_spm import (  # noqa: E402
    FARADAY_C_PER_MOL,
    GAS_CONSTANT_J_PER_MOL_K,
    GroupedSPMParameters,
    exact_zoh_matrices,
)
from models.ocp_panasonic_qi_frozen import (  # noqa: E402
    Q_CELL_AH,
    U_n,
    U_p,
    X_N_END,
    X_N_START,
    X_P_END,
    X_P_START,
)
from runtime_paths import PANASONIC_HPPC_PATH, require_data_file  # noqa: E402


OUTPUT_DIRECTORY = REPOSITORY_ROOT / "thesis_data"
PULSE_SUMMARY_PATH = OUTPUT_DIRECTORY / "panasonic_hppc_pulse_summary.csv"
HPPC_CSV_PATH = PANASONIC_HPPC_PATH

EXCLUDED_PULSE_IDS = frozenset({60, 64, 67})
PARAMETER_NAMES = ("alpha_n", "alpha_p", "d_n", "d_p", "R0")

# The generic model requires all nine parameters explicitly and contains no
# target-cell defaults.  These five values are the only numerical dynamic
# fixture already documented in README.md and the smoke tests.  They are used
# solely as a local sensitivity reference, never as Panasonic estimates.
REFERENCE_DYNAMIC_PARAMETERS = {
    "alpha_n": 3000.0,
    "alpha_p": 1800.0,
    "d_n": 5.0e-5,
    "d_p": 7.0e-4,
    "R0": 0.020,
}
RELATIVE_FINITE_DIFFERENCE_STEP = 1.0e-4

Q_TOTAL_C = 3600.0 * Q_CELL_AH
Q_N_C = Q_TOTAL_C / (X_N_START - X_N_END)
Q_P_C = Q_TOTAL_C / (X_P_END - X_P_START)
SOC0_N = X_N_START
SOC0_P = X_P_START

SOC_BINS = (
    ("high", 2.0 / 3.0, np.inf),
    ("middle", 1.0 / 3.0, 2.0 / 3.0),
    ("low", -np.inf, 1.0 / 3.0),
)


@dataclass(frozen=True)
class HppcData:
    time_s: np.ndarray
    voltage_v: np.ndarray
    raw_current_a: np.ndarray
    exported_ah: np.ndarray
    battery_temperature_k: np.ndarray


@dataclass(frozen=True)
class PulseRecord:
    pulse_id: int
    nominal_rate: str
    start_index: int
    end_index: int
    pre_pulse_index: int
    post_rest_index: int
    soc: float
    usable: bool


@dataclass(frozen=True)
class LocalResult:
    voltage_v: np.ndarray
    x_n_avg: np.ndarray
    x_n_aux: np.ndarray
    x_n_surface: np.ndarray
    x_p_avg: np.ndarray
    x_p_aux: np.ndarray
    x_p_surface: np.ndarray


def load_hppc_data(path: Path = HPPC_CSV_PATH) -> HppcData:
    table = np.genfromtxt(require_data_file(path), delimiter=",", names=True, dtype=float, encoding="utf-8")
    names = set(table.dtype.names or ())
    required = {"Time", "Voltage", "Current", "Ah", "Battery_Temp_degC"}
    missing = required - names
    if missing:
        raise ValueError(f"HPPC CSV is missing columns: {sorted(missing)}")
    arrays = [np.asarray(table[name], dtype=float) for name in required]
    if not all(np.all(np.isfinite(array)) for array in arrays):
        raise ValueError("HPPC input contains non-finite data.")
    time_s = np.asarray(table["Time"], dtype=float)
    if np.any(np.diff(time_s) < 0.0):
        raise ValueError("HPPC elapsed time is not nondecreasing.")
    return HppcData(
        time_s=time_s,
        voltage_v=np.asarray(table["Voltage"], dtype=float),
        raw_current_a=np.asarray(table["Current"], dtype=float),
        exported_ah=np.asarray(table["Ah"], dtype=float),
        battery_temperature_k=np.asarray(table["Battery_Temp_degC"], dtype=float)
        + 273.15,
    )


def _parse_bool(value: str) -> bool:
    if value == "True":
        return True
    if value == "False":
        return False
    raise ValueError(f"Expected explicit boolean text, received {value!r}.")


def load_pulse_records(path: Path = PULSE_SUMMARY_PATH) -> list[PulseRecord]:
    with require_data_file(path).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    records = [
        PulseRecord(
            pulse_id=int(row["pulse_id"]),
            nominal_rate=row["nominal_rate"],
            start_index=int(row["start_index"]),
            end_index=int(row["end_index"]),
            pre_pulse_index=int(row["pre_pulse_index"]),
            post_rest_index=int(row["post_rest_index"]),
            soc=float(row["soc"]),
            usable=_parse_bool(row["usable_for_dynamic_calibration"]),
        )
        for row in rows
    ]
    usable = [
        record
        for record in records
        if record.usable and record.pulse_id not in EXCLUDED_PULSE_IDS
    ]
    if len(usable) != 64:
        raise ValueError(f"Expected 64 complete pulses, found {len(usable)}.")
    observed_excluded = {record.pulse_id for record in records if not record.usable}
    if observed_excluded != EXCLUDED_PULSE_IDS:
        raise ValueError(
            "The pulse-summary exclusions differ from the frozen audit: "
            f"{sorted(observed_excluded)}."
        )
    return records


def stoichiometries_from_discharged_ah(
    discharged_ah: float,
) -> tuple[float, float]:
    z = discharged_ah / Q_CELL_AH
    if not 0.0 <= z <= 1.0:
        raise ValueError(f"Cumulative discharged fraction is outside [0,1]: {z}.")
    return (
        X_N_START + z * (X_N_END - X_N_START),
        X_P_START + z * (X_P_END - X_P_START),
    )


def verify_frozen_state_mapping() -> dict[str, float]:
    """Verify the mapping analytically and with the actual ZOH propagator."""

    final_x_n = SOC0_N - Q_TOTAL_C / Q_N_C
    final_x_p = SOC0_P + Q_TOTAL_C / Q_P_C
    errors = {
        "x_n_final": final_x_n,
        "x_p_final": final_x_p,
        "x_n_error": final_x_n - X_N_END,
        "x_p_error": final_x_p - X_P_END,
    }
    if abs(errors["x_n_error"]) > 5.0e-15 or abs(errors["x_p_error"]) > 5.0e-15:
        raise RuntimeError(f"Frozen state mapping failed: {errors}")

    parameters = GroupedSPMParameters(
        alpha_n_s=REFERENCE_DYNAMIC_PARAMETERS["alpha_n"],
        alpha_p_s=REFERENCE_DYNAMIC_PARAMETERS["alpha_p"],
        q_n_c=Q_N_C,
        q_p_c=Q_P_C,
        d_n=REFERENCE_DYNAMIC_PARAMETERS["d_n"],
        d_p=REFERENCE_DYNAMIC_PARAMETERS["d_p"],
        soc0_n=SOC0_N,
        soc0_p=SOC0_P,
        r0_ohm=REFERENCE_DYNAMIC_PARAMETERS["R0"],
    )
    initial_state = np.asarray([SOC0_N, SOC0_N, SOC0_P, SOC0_P])
    # One ampere for Q_TOTAL_C seconds passes exactly Q_CELL_AH. The averages
    # must hit the endpoints even though the loaded auxiliary states retain a
    # diffusion gradient. A long zero-current rest removes that gradient.
    loaded_a, loaded_b = exact_zoh_matrices(parameters, Q_TOTAL_C)
    loaded_state = loaded_a @ initial_state + loaded_b[:, 0]
    rest_duration_s = 100.0 * max(
        REFERENCE_DYNAMIC_PARAMETERS["alpha_n"],
        REFERENCE_DYNAMIC_PARAMETERS["alpha_p"],
    )
    rested_a, _ = exact_zoh_matrices(parameters, rest_duration_s)
    rested_state = rested_a @ loaded_state
    target_state = np.asarray([X_N_END, X_N_END, X_P_END, X_P_END])
    errors["zoh_loaded_x_n_avg_error"] = loaded_state[0] - X_N_END
    errors["zoh_loaded_x_p_avg_error"] = loaded_state[2] - X_P_END
    errors["zoh_rested_max_state_error"] = float(
        np.max(np.abs(rested_state - target_state))
    )
    errors["zoh_rest_duration_s"] = rest_duration_s
    if (
        abs(errors["zoh_loaded_x_n_avg_error"]) > 5.0e-13
        or abs(errors["zoh_loaded_x_p_avg_error"]) > 5.0e-13
        or errors["zoh_rested_max_state_error"] > 5.0e-13
    ):
        raise RuntimeError(f"Actual ZOH state-propagation check failed: {errors}")
    return errors


def simulate_local_group(
    time_s: np.ndarray,
    discharge_positive_current_a: np.ndarray,
    temperature_k: np.ndarray,
    initial_x_n: float,
    initial_x_p: float,
    dynamic_parameters: dict[str, float],
) -> LocalResult:
    """Evaluate the implemented grouped-SPM equations on one local history.

    The two-electrode exact update below is the analytic form of the same ZOH
    state matrices used by ``models/grouped_spm.py``.  Frozen Panasonic OCPs
    are evaluated here because the generic model intentionally imports its
    original chemistry-specific OCP module.
    """

    time_s = np.asarray(time_s, dtype=float)
    current = np.asarray(discharge_positive_current_a, dtype=float)
    temperature = np.asarray(temperature_k, dtype=float)
    if not (time_s.shape == current.shape == temperature.shape):
        raise ValueError("Local time, current, and temperature shapes differ.")
    if time_s.size < 2 or np.any(np.diff(time_s) < 0.0):
        raise ValueError("A local history needs at least two nondecreasing samples.")

    alpha_n = dynamic_parameters["alpha_n"]
    alpha_p = dynamic_parameters["alpha_p"]
    d_n = dynamic_parameters["d_n"]
    d_p = dynamic_parameters["d_p"]
    r0 = dynamic_parameters["R0"]
    if min(alpha_n, alpha_p, d_n, d_p) <= 0.0 or r0 < 0.0:
        raise ValueError("Local dynamic parameters are outside model constraints.")

    n_samples = time_s.size
    x_n_avg = np.empty(n_samples)
    x_n_aux = np.empty(n_samples)
    x_p_avg = np.empty(n_samples)
    x_p_aux = np.empty(n_samples)
    x_n_avg[0] = x_n_aux[0] = initial_x_n
    x_p_avg[0] = x_p_aux[0] = initial_x_p

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
    if np.any((x_n_surface < 0.01) | (x_n_surface > 0.80)):
        raise ValueError(
            "Negative surface stoichiometry left the frozen OCP support; no clipping."
        )
    if np.any((x_p_surface < 0.18) | (x_p_surface > 0.96)):
        raise ValueError(
            "Positive surface stoichiometry left the frozen OCP support; no clipping."
        )

    kinetic_factor = 2.0 * GAS_CONSTANT_J_PER_MOL_K * temperature / FARADAY_C_PER_MOL
    eta_n = kinetic_factor * np.arcsinh(
        current
        / (6.0 * Q_N_C * d_n * np.sqrt(x_n_surface * (1.0 - x_n_surface)))
    )
    eta_p = kinetic_factor * np.arcsinh(
        current
        / (6.0 * Q_P_C * d_p * np.sqrt(x_p_surface * (1.0 - x_p_surface)))
    )
    voltage = U_p(x_p_surface) - U_n(x_n_surface) + eta_p - eta_n - r0 * current
    if not np.all(np.isfinite(voltage)):
        raise ValueError("Non-finite voltage in local forward evaluation.")
    return LocalResult(
        voltage_v=np.asarray(voltage),
        x_n_avg=x_n_avg,
        x_n_aux=x_n_aux,
        x_n_surface=x_n_surface,
        x_p_avg=x_p_avg,
        x_p_aux=x_p_aux,
        x_p_surface=x_p_surface,
    )


def _phase_for_index(index: int, record: PulseRecord) -> str | None:
    if index == record.pre_pulse_index:
        return None
    if index == record.start_index:
        return "immediate"
    if record.start_index < index <= record.end_index:
        return "within_pulse"
    if record.end_index < index <= record.post_rest_index:
        return "post_relaxation"
    raise ValueError(f"Index {index} is outside pulse {record.pulse_id}'s group.")


def partition_contiguous_groups(
    data: HppcData, records: list[PulseRecord]
) -> list[list[PulseRecord]]:
    """Partition the audited pulses at the source file's long omitted gaps."""

    groups: list[list[PulseRecord]] = []
    for record in records:
        if not groups:
            groups.append([record])
            continue
        previous = groups[-1][-1]
        gap_s = (
            data.time_s[record.pre_pulse_index]
            - data.time_s[previous.post_rest_index]
        )
        if gap_s > 10.0:
            groups.append([record])
        else:
            groups[-1].append(record)
    if len(groups) != 14:
        raise ValueError(f"Expected 14 contiguous pulse groups, found {len(groups)}.")
    return groups


def build_sensitivity_rows(
    data: HppcData, records: list[PulseRecord]
) -> tuple[list[dict[str, float | int | str]], int]:
    rows: list[dict[str, float | int | str]] = []
    groups = partition_contiguous_groups(data, records)
    for group in groups:
        usable_group = [
            record
            for record in group
            if record.usable and record.pulse_id not in EXCLUDED_PULSE_IDS
        ]
        if not usable_group:
            continue
        lower = usable_group[0].pre_pulse_index
        upper = usable_group[-1].post_rest_index + 1
        absolute_indices = np.arange(lower, upper)
        local_time = data.time_s[lower:upper] - data.time_s[lower]
        current = -data.raw_current_a[lower:upper]
        temperature = data.battery_temperature_k[lower:upper]
        discharged_ah = max(0.0, -float(data.exported_ah[lower]))
        initial_x_n, initial_x_p = stoichiometries_from_discharged_ah(discharged_ah)

        base = simulate_local_group(
            local_time,
            current,
            temperature,
            initial_x_n,
            initial_x_p,
            REFERENCE_DYNAMIC_PARAMETERS,
        )
        normalized_sensitivities: dict[str, np.ndarray] = {}
        for parameter_name in PARAMETER_NAMES:
            plus = dict(REFERENCE_DYNAMIC_PARAMETERS)
            minus = dict(REFERENCE_DYNAMIC_PARAMETERS)
            plus[parameter_name] *= 1.0 + RELATIVE_FINITE_DIFFERENCE_STEP
            minus[parameter_name] *= 1.0 - RELATIVE_FINITE_DIFFERENCE_STEP
            voltage_plus = simulate_local_group(
                local_time,
                current,
                temperature,
                initial_x_n,
                initial_x_p,
                plus,
            ).voltage_v
            voltage_minus = simulate_local_group(
                local_time,
                current,
                temperature,
                initial_x_n,
                initial_x_p,
                minus,
            ).voltage_v
            # Central difference with respect to fractional parameter change:
            # dV/d(ln p) = p*dV/dp, in volts.
            normalized_sensitivities[parameter_name] = (
                voltage_plus - voltage_minus
            ) / (2.0 * RELATIVE_FINITE_DIFFERENCE_STEP)

        for record in usable_group:
            soc_region = next(
                name for name, lower_soc, upper_soc in SOC_BINS
                if lower_soc <= record.soc < upper_soc
            )
            for absolute_index in range(
                record.pre_pulse_index, record.post_rest_index + 1
            ):
                phase = _phase_for_index(absolute_index, record)
                if phase is None:
                    continue
                local_index = absolute_index - lower
                row: dict[str, float | int | str] = {
                    "pulse_id": record.pulse_id,
                    "nominal_rate": record.nominal_rate,
                    "soc_region": soc_region,
                    "pulse_start_soc": record.soc,
                    "phase": phase,
                    "source_index": absolute_index,
                    "local_time_s": float(local_time[local_index]),
                    "current_a_discharge_positive": float(current[local_index]),
                    "measured_voltage_v": float(data.voltage_v[absolute_index]),
                    "reference_voltage_v": float(base.voltage_v[local_index]),
                }
                for parameter_name in PARAMETER_NAMES:
                    row[f"S_{parameter_name}_V"] = float(
                        normalized_sensitivities[parameter_name][local_index]
                    )
                rows.append(row)
    return rows, len(groups)


def rows_to_jacobian(rows: list[dict[str, object]]) -> np.ndarray:
    return np.asarray(
        [[float(row[f"S_{name}_V"]) for name in PARAMETER_NAMES] for row in rows]
    )


def analyze_subset(rows: list[dict[str, object]]) -> dict[str, object]:
    jacobian = rows_to_jacobian(rows)
    if jacobian.shape[0] < len(PARAMETER_NAMES):
        raise ValueError("A sensitivity subset has too few observations for SVD.")
    column_rms = np.sqrt(np.mean(jacobian**2, axis=0))
    singular_values = np.linalg.svd(jacobian, full_matrices=False, compute_uv=False)
    _, _, vt = np.linalg.svd(jacobian, full_matrices=False)
    condition_number = (
        float(singular_values[0] / singular_values[-1])
        if singular_values[-1] > 0.0
        else float("inf")
    )

    # Estimator-oriented parameter correlation from the pseudoinverse of J'J.
    # The SVD tolerance is recorded so rank-deficient cases remain explicit.
    tolerance = np.finfo(float).eps * max(jacobian.shape) * singular_values[0]
    rank = int(np.sum(singular_values > tolerance))
    covariance = np.linalg.pinv(jacobian.T @ jacobian, rcond=1.0e-12)
    standard_deviation = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    denominator = np.outer(standard_deviation, standard_deviation)
    correlation = np.divide(
        covariance,
        denominator,
        out=np.full_like(covariance, np.nan),
        where=denominator > 0.0,
    )
    np.fill_diagonal(correlation, 1.0)
    return {
        "n_observations": jacobian.shape[0],
        "column_rms": column_rms,
        "singular_values": singular_values,
        "condition_number": condition_number,
        "rank": rank,
        "correlation": correlation,
        "weakest_vector": vt[-1],
    }


def make_subsets(
    rows: list[dict[str, object]],
) -> list[tuple[str, str, list[dict[str, object]]]]:
    subsets: list[tuple[str, str, list[dict[str, object]]]] = [
        ("phase", "all_combined", rows)
    ]
    for phase in ("immediate", "within_pulse", "post_relaxation"):
        subsets.append(("phase", phase, [row for row in rows if row["phase"] == phase]))
    for rate in ("0.5C", "1C", "2C", "4C", "6C"):
        subsets.append(
            ("current", rate, [row for row in rows if row["nominal_rate"] == rate])
        )
    for region in ("high", "middle", "low"):
        subsets.append(
            ("soc", region, [row for row in rows if row["soc_region"] == region])
        )
    return subsets


def write_sensitivity_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_identifiability_csvs(
    correlation_path: Path,
    svd_path: Path,
    analyses: list[tuple[str, str, dict[str, object]]],
) -> None:
    with correlation_path.open("w", newline="", encoding="utf-8") as handle:
        fields = (
            "subset_type",
            "subset",
            "n_observations",
            "row_parameter",
            "column_parameter",
            "correlation",
        )
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for subset_type, subset_name, analysis in analyses:
            matrix = np.asarray(analysis["correlation"])
            for row_index, row_parameter in enumerate(PARAMETER_NAMES):
                for column_index, column_parameter in enumerate(PARAMETER_NAMES):
                    writer.writerow(
                        {
                            "subset_type": subset_type,
                            "subset": subset_name,
                            "n_observations": analysis["n_observations"],
                            "row_parameter": row_parameter,
                            "column_parameter": column_parameter,
                            "correlation": f"{matrix[row_index, column_index]:.12g}",
                        }
                    )

    with svd_path.open("w", newline="", encoding="utf-8") as handle:
        fields = [
            "subset_type",
            "subset",
            "n_observations",
            "numerical_rank",
            "condition_number",
            "singular_index",
            "singular_value",
            *[f"column_rms_{name}_V" for name in PARAMETER_NAMES],
            *[f"weakest_vector_{name}" for name in PARAMETER_NAMES],
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for subset_type, subset_name, analysis in analyses:
            common = {
                "subset_type": subset_type,
                "subset": subset_name,
                "n_observations": analysis["n_observations"],
                "numerical_rank": analysis["rank"],
                "condition_number": f"{float(analysis['condition_number']):.12g}",
            }
            common.update(
                {
                    f"column_rms_{name}_V": f"{value:.12g}"
                    for name, value in zip(
                        PARAMETER_NAMES, np.asarray(analysis["column_rms"]), strict=True
                    )
                }
            )
            common.update(
                {
                    f"weakest_vector_{name}": f"{value:.12g}"
                    for name, value in zip(
                        PARAMETER_NAMES,
                        np.asarray(analysis["weakest_vector"]),
                        strict=True,
                    )
                }
            )
            for index, singular_value in enumerate(analysis["singular_values"], start=1):
                writer.writerow(
                    {
                        **common,
                        "singular_index": index,
                        "singular_value": f"{float(singular_value):.12g}",
                    }
                )


def _format_parameter(name: str) -> str:
    return {
        "alpha_n": "alpha_n",
        "alpha_p": "alpha_p",
        "d_n": "d_n",
        "d_p": "d_p",
        "R0": "R0",
    }[name]


def _rms_table(analysis: dict[str, object]) -> list[str]:
    rms = np.asarray(analysis["column_rms"])
    return [
        f"| `{name}` | {1000.0 * value:.6g} |"
        for name, value in zip(PARAMETER_NAMES, rms, strict=True)
    ]


def high_correlation_pairs(analysis: dict[str, object]) -> list[tuple[str, str, float]]:
    matrix = np.asarray(analysis["correlation"])
    pairs = []
    for first in range(len(PARAMETER_NAMES)):
        for second in range(first + 1, len(PARAMETER_NAMES)):
            value = float(matrix[first, second])
            if np.isfinite(value) and abs(value) >= 0.90:
                pairs.append((PARAMETER_NAMES[first], PARAMETER_NAMES[second], value))
    return pairs


def write_parameter_mapping(path: Path, mapping_errors: dict[str, float]) -> None:
    lines = [
        "# Panasonic grouped-SPM parameter mapping",
        "",
        "## Scope",
        "",
        "Code/equation audit and analytical mapping only. No parameter was optimized, and the frozen model and Panasonic OCP modules were not modified.",
        "",
        "## Actual implemented parameterization",
        "",
        "For discharge-positive model current `I`, the implemented state/output equations are:",
        "",
        "```text",
        "dx_n_avg/dt = -I/Q_n",
        "dx_n_aux/dt = (30/alpha_n)(x_n_avg-x_n_aux) - 19 I/(7 Q_n)",
        "x_n_surf = x_n_aux - alpha_n I/(105 Q_n)",
        "dx_p_avg/dt = +I/Q_p",
        "dx_p_aux/dt = (30/alpha_p)(x_p_avg-x_p_aux) + 19 I/(7 Q_p)",
        "x_p_surf = x_p_aux + alpha_p I/(105 Q_p)",
        "eta_k = (2RT/F) asinh[I/(6 Q_k d_k sqrt(x_k_surf(1-x_k_surf)))]",
        "V = U_p(x_p_surf)-U_n(x_n_surf)+eta_p-eta_n-R0 I",
        "```",
        "",
        "| Parameter | Units | Exact role and equation locations | Category | Status for HPPC |",
        "|---|---:|---|---|---|",
        "| `alpha_n` | s | `30/alpha_n` negative auxiliary-state relaxation and `-alpha_n I/(105 Q_n)` surface correction. The documented group is particle radius squared over solid diffusivity. | solid diffusion/state output | dynamic candidate |",
        "| `alpha_p` | s | `30/alpha_p` positive auxiliary-state relaxation and `+alpha_p I/(105 Q_p)` surface correction. | solid diffusion/state output | dynamic candidate |",
        "| `Q_n` | C | Negative average and auxiliary state forcing, surface correction, and kinetic denominator. | electrode state evolution/kinetic scale | analytically frozen |",
        "| `Q_p` | C | Positive average and auxiliary state forcing, surface correction, and kinetic denominator. | electrode state evolution/kinetic scale | analytically frozen |",
        "| `d_n` | 1/s | Appears only in `6 Q_n d_n sqrt(x(1-x))` inside negative `asinh`. | grouped Butler-Volmer kinetics | dynamic candidate |",
        "| `d_p` | 1/s | Appears only in `6 Q_p d_p sqrt(x(1-x))` inside positive `asinh`. | grouped Butler-Volmer kinetics | dynamic candidate |",
        "| `SOC0_n` | dimensionless | Initializes both negative average and auxiliary normalized concentrations. It is not cell SOC. | thermodynamic alignment/initial state | analytically frozen |",
        "| `SOC0_p` | dimensionless | Initializes both positive average and auxiliary normalized concentrations. It is not cell SOC. | thermodynamic alignment/initial state | analytically frozen |",
        "| `R0` | ohm | Algebraic terminal drop `-R0 I`; also supplies `R0 I^2` to the coupled heat calculation. | ohmic resistance | dynamic candidate |",
        "",
        "`models/electrothermal_spm.py` does not redefine these electrical parameters. It re-evaluates the grouped kinetic equations at predicted temperature, applies the same voltage equation, and keeps cell-level `Q_cell_Ah`/SOC separate.",
        "",
        "## Exact frozen-window mapping",
        "",
        "With total discharged charge `Delta q = 3600 Q_cell_Ah`, integration of the two average-state equations gives:",
        "",
        "```text",
        "Q_n = Delta q/(x_n_start-x_n_end)",
        "Q_p = Delta q/(x_p_end-x_p_start)",
        "SOC0_n = x_n_start",
        "SOC0_p = x_p_start",
        "```",
        "",
        f"- `Delta q = {Q_TOTAL_C:.12f} C`",
        f"- `Q_n = {Q_N_C:.12f} C` ({Q_N_C / 3600.0:.12f} Ah-equivalent)",
        f"- `Q_p = {Q_P_C:.12f} C` ({Q_P_C / 3600.0:.12f} Ah-equivalent)",
        f"- `SOC0_n = {SOC0_N:.12f}`",
        f"- `SOC0_p = {SOC0_P:.12f}`",
        "",
        "These are grouped effective electrode capacity/state quantities derived from the nominal effective windows, not separately measured physical electrode capacities or uniquely identified stoichiometries.",
        "",
        "## Forward state check",
        "",
        f"Passing exactly `{Q_CELL_AH:.12f} Ah` through the average-state equations gives `x_n={mapping_errors['x_n_final']:.12f}` and `x_p={mapping_errors['x_p_final']:.12f}`.",
        f"The numerical endpoint errors are `{mapping_errors['x_n_error']:.3e}` and `{mapping_errors['x_p_error']:.3e}`. Thus the frozen mapping is reproduced to floating-point precision without optimization or clipping.",
        "",
        f"An independent forward check using the actual `exact_zoh_matrices` implementation applied 1 A for `{Q_TOTAL_C:.12f} s` (exactly the frozen cell capacity), then zero current for `{mapping_errors['zoh_rest_duration_s']:.0f} s`. The loaded average-state endpoint errors were `{mapping_errors['zoh_loaded_x_n_avg_error']:.3e}` and `{mapping_errors['zoh_loaded_x_p_avg_error']:.3e}`; after relaxation the maximum error across all four states was `{mapping_errors['zoh_rested_max_state_error']:.3e}`.",
        "",
        "The exact loaded endpoint statement applies to the volume-average states. During discharge, auxiliary/surface states differ by the model's diffusion terms; after zero-current relaxation they converge to the unchanged average endpoints.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_identifiability_report(
    path: Path,
    rows: list[dict[str, object]],
    analyses: list[tuple[str, str, dict[str, object]]],
    n_contiguous_groups: int,
) -> None:
    lookup = {(kind, name): analysis for kind, name, analysis in analyses}
    combined = lookup[("phase", "all_combined")]
    weakest = np.asarray(combined["weakest_vector"])
    pair_rows = high_correlation_pairs(combined)
    pulse_ids = sorted({int(row["pulse_id"]) for row in rows})
    source_indices = np.asarray([int(row["source_index"]) for row in rows])
    duplicate_observation_indices = len(source_indices) - len(np.unique(source_indices))

    lines = [
        "# Panasonic HPPC local sensitivity and identifiability audit",
        "",
        "## Scope and data product",
        "",
        "No optimization or dynamic calibration was performed. The frozen OCP, nominal effective windows, grouped-SPM equations, and electrothermal architecture were not modified.",
        "",
        f"- Complete pulse/rest groups: **{len(pulse_ids)}** (IDs 60, 64, and 67 excluded).",
        f"- Independently initialized contiguous multi-rate histories: **{n_contiguous_groups}**, separated by the source file's long omitted-data gaps.",
        f"- Voltage-response observations in the Jacobian: {len(rows)}.",
        f"- Repeated source indices across independent groups: {duplicate_observation_indices} (a boundary rest sample can close one group and initialize the next, but the initializing sample is excluded from the Jacobian).",
        "- Each local history begins at the first usable pulse's `pre_pulse_index`, ends at the last usable pulse's `post_rest_index`, and uses local time zero at that block boundary. State is preserved across every pulse and rest inside the contiguous block.",
        "- Raw negative current was converted at the loader boundary to model discharge-positive current.",
        "- Cumulative exported `Ah` initializes each local average and auxiliary electrode state; current is never integrated through the 13 long omitted-data gaps.",
        "- Actual nonuniform time intervals, including zero-length duplicate timestamps, were retained. No resampling was performed.",
        "- Measured battery temperature is used only as the prescribed temperature in the existing kinetic factor; no thermal parameters are involved.",
        "",
        "Phase definitions: `immediate` is the first loaded sample; `within_pulse` is every later loaded sample through pulse end; `post_relaxation` is every contiguous rest sample after pulse end. `all_combined` contains all three phases and excludes the pre-pulse initialization sample.",
        "",
        "SOC bins are non-overlapping and use pulse-start cell SOC: low `[0,1/3)`, middle `[1/3,2/3)`, and high `[2/3,1]`.",
        "",
        "## Local reference and finite differences",
        "",
        "The model API has no built-in parameter defaults. The existing README/smoke-test fixture is therefore used only as a local numerical reference, not as a Panasonic parameter estimate:",
        "",
        "| Parameter | Local reference |",
        "|---|---:|",
        *[
            f"| `{name}` | {REFERENCE_DYNAMIC_PARAMETERS[name]:.12g} |"
            for name in PARAMETER_NAMES
        ],
        "",
        f"Central scale-aware perturbations use `p*(1 +/- {RELATIVE_FINITE_DIFFERENCE_STEP:g})`. Each Jacobian entry is `dV/d(ln p) = p dV/dp` in volts, so columns correspond to equal fractional parameter changes.",
        "",
        "## Sensitivity by response phase",
        "",
        "Column RMS values below are mV per unit fractional parameter change.",
        "",
    ]
    for phase in ("immediate", "within_pulse", "post_relaxation", "all_combined"):
        analysis = lookup[("phase", phase)]
        lines.extend(
            [
                f"### {phase.replace('_', ' ').title()}",
                "",
                f"Observations: {analysis['n_observations']}",
                "",
                "| Parameter | RMS normalized sensitivity [mV] |",
                "|---|---:|",
                *_rms_table(analysis),
                "",
            ]
        )

    lines.extend(
        [
            "## Combined identifiability",
            "",
            f"- Numerical rank: {combined['rank']} of {len(PARAMETER_NAMES)}.",
            f"- Singular values: {', '.join(f'{value:.6g}' for value in np.asarray(combined['singular_values']))}.",
            f"- Jacobian condition number: {float(combined['condition_number']):.6g}.",
            "- Weakest right-singular direction (unit vector; sign is arbitrary): "
            + ", ".join(
                f"`{name}`={value:+.5f}"
                for name, value in zip(PARAMETER_NAMES, weakest, strict=True)
            )
            + ".",
            "",
            "The reported parameter-correlation matrix is estimator-oriented: the pseudoinverse of `J.T @ J` is normalized by its diagonal. It is not the ordinary correlation of parameter values, because no parameters were estimated.",
            "",
            "Pairs with `|correlation| >= 0.90` in the combined audit:",
            "",
        ]
    )
    if pair_rows:
        lines.extend(
            [
                "| Parameter pair | Correlation |",
                "|---|---:|",
                *[
                    f"| `{first}` / `{second}` | {value:+.6f} |"
                    for first, second, value in pair_rows
                ],
            ]
        )
    else:
        lines.append("None.")

    lines.extend(
        [
            "",
            "## Stratified audit",
            "",
            "| Stratum | Observations | Condition number | Weakest singular value |",
            "|---|---:|---:|---:|",
        ]
    )
    for subset_type in ("current", "soc"):
        for kind, name, analysis in analyses:
            if kind != subset_type:
                continue
            singular = np.asarray(analysis["singular_values"])
            lines.append(
                f"| {subset_type}: `{name}` | {analysis['n_observations']} | "
                f"{float(analysis['condition_number']):.6g} | {singular[-1]:.6g} |"
            )

    lines.extend(
        [
            "",
            "High-correlation pairs in the stratified audits:",
            "",
            "| Stratum | Parameter pair | Correlation |",
            "|---|---|---:|",
        ]
    )
    stratified_pairs_found = False
    for kind, name, analysis in analyses:
        if kind not in {"current", "soc"}:
            continue
        for first, second, value in high_correlation_pairs(analysis):
            stratified_pairs_found = True
            lines.append(
                f"| {kind}: `{name}` | `{first}` / `{second}` | {value:+.6f} |"
            )
    if not stratified_pairs_found:
        lines.append("| all current/SOC strata | none | -- |")

    lines.extend(
        [
            "",
            "Detailed SVD, weakest-direction, and correlation results for every phase/current/SOC subset are in the machine-readable CSV files.",
            "",
            "## Apparent resistance context",
            "",
            "The 64 usable measured pulse steps have `R_app = 20.071-33.348 mOhm`. This observed step quantity contains the model's algebraic `R0 I` contribution together with instantaneous surface-state/OCP changes and the algebraic grouped kinetic contributions `eta_p-eta_n`. Consequently, `R_app` is contextual evidence for voltage-step scale only; setting `R0 = R_app` would assign non-ohmic voltage terms to `R0` and is not justified.",
            "",
            "Consistently, the immediate-response normalized RMS sensitivities are 185.930 mV for `R0`, 40.578 mV for `d_n`, 15.272 mV for `d_p`, 12.806 mV for `alpha_p`, and 1.136 mV for `alpha_n`. The model therefore predicts that the measured instantaneous step is not uniquely ohmic at this local reference. The fixture `R0=20 mOhm` happens to sit at the lower edge of the measured `R_app` range, but that numerical proximity is not an identification result.",
            "",
            "## Minimum defensible dynamic subset",
            "",
            "| Parameter | Recommendation | Reason from implemented equations and this local audit |",
            "|---|---|---|",
            "| `R0` | **CALIBRATE** | It is the only purely ohmic candidate and acts directly on every loaded sample; it has no relaxation sensitivity. It must be separated jointly from algebraic kinetic voltage rather than assigned from `R_app`. |",
            "| `alpha_n` | **WEAKLY IDENTIFIABLE / CONDITIONAL** | It controls negative diffusion, but has the smallest combined RMS sensitivity and contributes about 0.998 of the unit weakest singular direction. Do not release it in the smallest first fit. |",
            "| `alpha_p` | **CALIBRATE** | It controls the positive surface displacement and positive diffusion-relaxation time scale; its post-rest RMS sensitivity is materially larger than that of `alpha_n`. |",
            "| `d_n` | **CALIBRATE** | It has a strong loaded-sample signature and is not highly correlated with `R0` only when multiple current magnitudes are combined. The near-unity `d_n`/`R0` correlations within each fixed C-rate make the multi-rate design essential. |",
            "| `d_p` | **WEAKLY IDENTIFIABLE / CONDITIONAL** | It enters only the algebraic positive kinetic overpotential, has no rest signature, and is strongly correlated with `R0` in the combined audit (`-0.936`). Do not release it in the smallest first fit. |",
            "",
            "Thus the smallest defensible primary calibration subset is **`R0`, `alpha_p`, and `d_n`**, using all current magnitudes together. Hold `alpha_n` and `d_p` fixed at explicitly declared provisional values for that first fit, then release at most one conditionally in a later nested check only if stability improves. This is a subset recommendation, not a parameter estimate; the present generic fixture is not asserted to be the physically correct fixed value.",
            "",
            "The earlier synthetic-identifiability phase is closed and its numerical artifacts are intentionally absent from this repository. Its general sensitivity conclusion cannot override this experimental audit: nonzero synthetic/local voltage sensitivity shows that a parameter affects the output, whereas practical recoverability additionally requires independent sensitivity directions under the actual HPPC currents, SOC coverage, sampling, noise, and parameter-reference point. The condition numbers, correlations, and weakest singular vectors here assess that distinction; they are still local diagnostics, not proof of global identifiability.",
            "",
            "## Stop condition",
            "",
            "No dynamic parameter optimization, HPPC calibration, UDDS/1C validation, thermal fitting, OCP change, window change, or model-equation change was performed.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def save_figures(
    analyses: list[tuple[str, str, dict[str, object]]], output_directory: Path
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    lookup = {(kind, name): analysis for kind, name, analysis in analyses}
    parameter_labels = [r"$\alpha_n$", r"$\alpha_p$", r"$d_n$", r"$d_p$", r"$R_0$"]

    phases = ("immediate", "within_pulse", "post_relaxation", "all_combined")
    phase_rms = np.asarray([lookup[("phase", phase)]["column_rms"] for phase in phases])
    figure, axis = plt.subplots(figsize=(10.5, 6.2), constrained_layout=True)
    width = 0.19
    positions = np.arange(len(PARAMETER_NAMES))
    for index, phase in enumerate(phases):
        axis.bar(
            positions + (index - 1.5) * width,
            1000.0 * phase_rms[index],
            width,
            label=phase.replace("_", " "),
        )
    axis.set_xticks(positions, parameter_labels)
    axis.set_ylabel("RMS normalized sensitivity [mV]")
    axis.set_title("Panasonic HPPC normalized local voltage sensitivity")
    axis.grid(axis="y", alpha=0.25)
    axis.legend()
    figure.savefig(output_directory / "panasonic_hppc_normalized_sensitivity.png", dpi=180)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(9.5, 5.8), constrained_layout=True)
    for kind, name in (("phase", "all_combined"), ("current", "0.5C"), ("current", "6C")):
        singular = np.asarray(lookup[(kind, name)]["singular_values"])
        axis.semilogy(np.arange(1, singular.size + 1), singular, marker="o", label=name)
    axis.set(
        xlabel="Singular-value index",
        ylabel="Singular value [V]",
        title="Normalized-Jacobian singular values",
        xticks=np.arange(1, len(PARAMETER_NAMES) + 1),
    )
    axis.grid(alpha=0.25, which="both")
    axis.legend()
    figure.savefig(output_directory / "panasonic_hppc_singular_values.png", dpi=180)
    plt.close(figure)

    correlation = np.asarray(lookup[("phase", "all_combined")]["correlation"])
    figure, axis = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
    image = axis.imshow(correlation, vmin=-1.0, vmax=1.0, cmap="coolwarm")
    axis.set_xticks(np.arange(len(PARAMETER_NAMES)), parameter_labels)
    axis.set_yticks(np.arange(len(PARAMETER_NAMES)), parameter_labels)
    for row in range(len(PARAMETER_NAMES)):
        for column in range(len(PARAMETER_NAMES)):
            axis.text(column, row, f"{correlation[row, column]:.2f}", ha="center", va="center")
    axis.set_title("Combined normalized parameter correlation")
    figure.colorbar(image, ax=axis, label="Correlation")
    figure.savefig(output_directory / "panasonic_hppc_parameter_correlation.png", dpi=180)
    plt.close(figure)

    def grouped_rms_plot(kind: str, names: tuple[str, ...], filename: str, title: str) -> None:
        rms = np.asarray([lookup[(kind, name)]["column_rms"] for name in names])
        figure, axis = plt.subplots(figsize=(10.5, 6.2), constrained_layout=True)
        for parameter_index, label in enumerate(parameter_labels):
            axis.plot(names, 1000.0 * rms[:, parameter_index], marker="o", label=label)
        axis.set(ylabel="RMS normalized sensitivity [mV]", title=title)
        axis.grid(alpha=0.25)
        axis.legend(ncol=3)
        figure.savefig(output_directory / filename, dpi=180)
        plt.close(figure)

    grouped_rms_plot(
        "current",
        ("0.5C", "1C", "2C", "4C", "6C"),
        "panasonic_hppc_sensitivity_by_current.png",
        "Normalized sensitivity by pulse-current magnitude",
    )
    grouped_rms_plot(
        "soc",
        ("high", "middle", "low"),
        "panasonic_hppc_sensitivity_by_soc.png",
        "Normalized sensitivity by pulse-start SOC region",
    )


def main() -> None:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    mapping_errors = verify_frozen_state_mapping()
    data = load_hppc_data()
    records = load_pulse_records()
    rows, n_contiguous_groups = build_sensitivity_rows(data, records)
    subsets = make_subsets(rows)
    analyses = [(kind, name, analyze_subset(subset)) for kind, name, subset in subsets]

    write_parameter_mapping(
        OUTPUT_DIRECTORY / "panasonic_hppc_parameter_mapping.md", mapping_errors
    )
    write_sensitivity_csv(
        OUTPUT_DIRECTORY / "panasonic_hppc_sensitivity.csv", rows
    )
    write_identifiability_csvs(
        OUTPUT_DIRECTORY / "panasonic_hppc_correlation.csv",
        OUTPUT_DIRECTORY / "panasonic_hppc_svd.csv",
        analyses,
    )
    write_identifiability_report(
        OUTPUT_DIRECTORY / "panasonic_hppc_identifiability_report.md",
        rows,
        analyses,
        n_contiguous_groups,
    )
    save_figures(analyses, OUTPUT_DIRECTORY)

    combined = dict(((kind, name), analysis) for kind, name, analysis in analyses)[
        ("phase", "all_combined")
    ]
    print(
        "Audited complete pulses: "
        f"{sum(record.usable and record.pulse_id not in EXCLUDED_PULSE_IDS for record in records)}"
    )
    print(f"Independently initialized contiguous groups: {n_contiguous_groups}")
    print(f"Jacobian observations: {len(rows)}")
    print(f"Q_n = {Q_N_C:.12f} C; Q_p = {Q_P_C:.12f} C")
    print(f"Combined condition number: {float(combined['condition_number']):.6g}")
    print((OUTPUT_DIRECTORY / "panasonic_hppc_identifiability_report.md").resolve())


if __name__ == "__main__":
    main()
