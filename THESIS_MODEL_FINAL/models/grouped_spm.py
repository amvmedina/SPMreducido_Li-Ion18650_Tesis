"""Reusable four-state grouped single-particle-model implementation.

The public current convention is positive for discharge. The implementation
uses exact zero-order-hold discretization on every supplied time interval.
No state, OCP, or voltage clipping is applied: an invalid state is reported and
outputs from the first domain violation onward are NaN.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Sequence

import numpy as np
from scipy.linalg import expm

from .ocp import U_n, U_p


FARADAY_C_PER_MOL = 96487.0
GAS_CONSTANT_J_PER_MOL_K = 8.314472
DEFAULT_TEMPERATURE_K = 298.15
PARAMETER_ORDER = (
    "alpha_n_s",
    "alpha_p_s",
    "q_n_c",
    "q_p_c",
    "d_n",
    "d_p",
    "soc0_n",
    "soc0_p",
    "r0_ohm",
)
CurrentConvention = Literal["discharge_positive", "guo_charge_positive"]


@dataclass(frozen=True)
class GroupedSPMParameters:
    """Nine grouped electrical parameters; values must be supplied explicitly."""

    alpha_n_s: float
    alpha_p_s: float
    q_n_c: float
    q_p_c: float
    d_n: float
    d_p: float
    soc0_n: float
    soc0_p: float
    r0_ohm: float

    def __post_init__(self) -> None:
        positive = (
            self.alpha_n_s,
            self.alpha_p_s,
            self.q_n_c,
            self.q_p_c,
            self.d_n,
            self.d_p,
        )
        if any(value <= 0.0 or not np.isfinite(value) for value in positive):
            raise ValueError("Diffusion, capacity, and kinetic groups must be positive.")
        if not 0.0 < self.soc0_n < 1.0 or not 0.0 < self.soc0_p < 1.0:
            raise ValueError("Initial normalized concentrations must be in (0, 1).")
        if self.r0_ohm < 0.0 or not np.isfinite(self.r0_ohm):
            raise ValueError("R0 must be finite and nonnegative.")

    def as_vector(self) -> np.ndarray:
        """Return parameters in ``PARAMETER_ORDER``."""

        return np.asarray([getattr(self, name) for name in PARAMETER_ORDER], dtype=float)

    @classmethod
    def from_vector(cls, values: Sequence[float]) -> "GroupedSPMParameters":
        """Build parameters from a nine-value vector in ``PARAMETER_ORDER``."""

        vector = np.asarray(values, dtype=float)
        if vector.shape != (len(PARAMETER_ORDER),):
            raise ValueError(f"Expected {len(PARAMETER_ORDER)} parameter values.")
        return cls(**dict(zip(PARAMETER_ORDER, vector, strict=True)))


@dataclass(frozen=True)
class DomainViolation:
    time_s: float
    current_a: float
    sample_index: int
    state_name: str
    electrode: str
    value: float


@dataclass(frozen=True)
class GroupedSPMResult:
    time_s: np.ndarray
    current_a: np.ndarray
    voltage_v: np.ndarray
    x_n_avg: np.ndarray
    x_n_aux: np.ndarray
    x_n_surface: np.ndarray
    x_p_avg: np.ndarray
    x_p_aux: np.ndarray
    x_p_surface: np.ndarray
    U_n_v: np.ndarray
    U_p_v: np.ndarray
    eta_n_v: np.ndarray
    eta_p_v: np.ndarray
    ohmic_drop_v: np.ndarray
    temperature_k: np.ndarray
    domain_valid: bool
    first_domain_violation: DomainViolation | None


def continuous_state_matrices(
    parameters: GroupedSPMParameters,
    *,
    current_convention: CurrentConvention = "discharge_positive",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return A, B, C, D for [x_n_avg, x_n_aux, x_p_avg, x_p_aux]."""

    if current_convention not in ("discharge_positive", "guo_charge_positive"):
        raise ValueError(f"Unsupported current convention: {current_convention}")
    p = parameters
    A = np.array(
        [
            [0.0, 0.0, 0.0, 0.0],
            [30.0 / p.alpha_n_s, -30.0 / p.alpha_n_s, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 30.0 / p.alpha_p_s, -30.0 / p.alpha_p_s],
        ],
        dtype=float,
    )
    B_guo = np.array(
        [
            1.0 / p.q_n_c,
            19.0 / (7.0 * p.q_n_c),
            -1.0 / p.q_p_c,
            -19.0 / (7.0 * p.q_p_c),
        ],
        dtype=float,
    )
    D_guo = np.array(
        [
            0.0,
            p.alpha_n_s / (105.0 * p.q_n_c),
            0.0,
            -p.alpha_p_s / (105.0 * p.q_p_c),
        ],
        dtype=float,
    )
    sign = -1.0 if current_convention == "discharge_positive" else 1.0
    return A, sign * B_guo[:, None], np.eye(4), sign * D_guo[:, None]


def exact_zoh_matrices(
    parameters: GroupedSPMParameters,
    dt_s: float,
    *,
    current_convention: CurrentConvention = "discharge_positive",
) -> tuple[np.ndarray, np.ndarray]:
    """Discretize one interval exactly for a zero-order-held current."""

    if dt_s < 0.0 or not np.isfinite(dt_s):
        raise ValueError("dt must be finite and nonnegative.")
    A, B, _, _ = continuous_state_matrices(
        parameters, current_convention=current_convention
    )
    if dt_s == 0.0:
        return np.eye(4), np.zeros((4, 1))
    augmented = np.zeros((5, 5), dtype=float)
    augmented[:4, :4] = A
    augmented[:4, 4:] = B
    discrete = expm(augmented * dt_s)
    return discrete[:4, :4], discrete[:4, 4:]


def surface_outputs(
    state: np.ndarray,
    current_a: float,
    parameters: GroupedSPMParameters,
    *,
    current_convention: CurrentConvention = "discharge_positive",
) -> tuple[float, float]:
    """Return negative and positive normalized surface concentrations."""

    _, _, C, D = continuous_state_matrices(
        parameters, current_convention=current_convention
    )
    output = C @ np.asarray(state, dtype=float) + D[:, 0] * float(current_a)
    return float(output[1]), float(output[3])


def grouped_overpotentials(
    current_a: float,
    x_n_surface: float,
    x_p_surface: float,
    parameters: GroupedSPMParameters,
    temperature_k: float = DEFAULT_TEMPERATURE_K,
) -> tuple[float, float]:
    """Return grouped negative- and positive-electrode overpotentials."""

    values = (x_n_surface, x_p_surface, temperature_k)
    if not all(np.isfinite(values)) or not 0.0 < x_n_surface < 1.0:
        raise ValueError("Negative surface concentration is outside (0, 1).")
    if not 0.0 < x_p_surface < 1.0:
        raise ValueError("Positive surface concentration is outside (0, 1).")
    if temperature_k <= 0.0:
        raise ValueError("Absolute temperature must be positive.")
    factor = 2.0 * GAS_CONSTANT_J_PER_MOL_K * temperature_k / FARADAY_C_PER_MOL
    denominator_n = (
        6.0
        * parameters.q_n_c
        * parameters.d_n
        * np.sqrt(x_n_surface * (1.0 - x_n_surface))
    )
    denominator_p = (
        6.0
        * parameters.q_p_c
        * parameters.d_p
        * np.sqrt(x_p_surface * (1.0 - x_p_surface))
    )
    eta_n = factor * np.arcsinh(float(current_a) / denominator_n)
    eta_p = factor * np.arcsinh(float(current_a) / denominator_p)
    return float(eta_n), float(eta_p)


def simulate_grouped_spm(
    time: np.ndarray,
    current: np.ndarray,
    parameters: GroupedSPMParameters,
    temperature: float | np.ndarray = DEFAULT_TEMPERATURE_K,
    current_convention: CurrentConvention = "discharge_positive",
) -> GroupedSPMResult:
    """Simulate the grouped SPM with exact interval-by-interval ZOH."""

    time_s = np.asarray(time, dtype=float)
    current_input = np.asarray(current, dtype=float)
    if time_s.ndim != 1 or time_s.size == 0:
        raise ValueError("time must be a nonempty one-dimensional array.")
    if current_input.shape != time_s.shape:
        raise ValueError("current must have the same shape as time.")
    if np.any(~np.isfinite(time_s)) or np.any(~np.isfinite(current_input)):
        raise ValueError("time and current must be finite.")
    if np.any(np.diff(time_s) < 0.0):
        raise ValueError("time must be nondecreasing.")
    if current_convention == "discharge_positive":
        current_model = current_input
    elif current_convention == "guo_charge_positive":
        current_model = -current_input
    else:
        raise ValueError(f"Unsupported current convention: {current_convention}")

    if np.isscalar(temperature):
        temperature_k = np.full(time_s.size, float(temperature))
    else:
        temperature_k = np.asarray(temperature, dtype=float)
        if temperature_k.shape != time_s.shape:
            raise ValueError("temperature array must have the same shape as time.")
    if np.any(~np.isfinite(temperature_k)) or np.any(temperature_k <= 0.0):
        raise ValueError("temperature must be finite and positive in kelvin.")

    names = (
        "voltage",
        "x_n_avg",
        "x_n_aux",
        "x_n_surface",
        "x_p_avg",
        "x_p_aux",
        "x_p_surface",
        "u_n",
        "u_p",
        "eta_n",
        "eta_p",
        "ohmic",
    )
    fields = {name: np.full(time_s.size, np.nan, dtype=float) for name in names}
    states = np.full((time_s.size, 4), np.nan, dtype=float)
    states[0] = [
        parameters.soc0_n,
        parameters.soc0_n,
        parameters.soc0_p,
        parameters.soc0_p,
    ]
    zoh_cache: dict[float, tuple[np.ndarray, np.ndarray]] = {}
    for index, dt_s_value in enumerate(np.diff(time_s)):
        dt_s = float(dt_s_value)
        matrices = zoh_cache.get(dt_s)
        if matrices is None:
            matrices = exact_zoh_matrices(parameters, dt_s)
            zoh_cache[dt_s] = matrices
        states[index + 1] = (
            matrices[0] @ states[index]
            + matrices[1][:, 0] * float(current_model[index])
        )

    _, _, _, surface_direct = continuous_state_matrices(parameters)
    surface_states = states + current_model[:, None] * surface_direct[:, 0]
    x_n_surface_all = surface_states[:, 1]
    x_p_surface_all = surface_states[:, 3]
    domain_values = np.column_stack(
        (
            states[:, 0],
            states[:, 1],
            x_n_surface_all,
            states[:, 2],
            states[:, 3],
            x_p_surface_all,
        )
    )
    domain_names = (
        ("x_n_avg", "negative"),
        ("x_n_aux", "negative"),
        ("x_n_surface", "negative"),
        ("x_p_avg", "positive"),
        ("x_p_aux", "positive"),
        ("x_p_surface", "positive"),
    )
    invalid_locations = np.argwhere(
        (~np.isfinite(domain_values))
        | (domain_values <= 0.0)
        | (domain_values >= 1.0)
    )
    violation = None
    valid_count = time_s.size
    if invalid_locations.size:
        sample_index, value_index = invalid_locations[0]
        state_name, electrode = domain_names[int(value_index)]
        violation = DomainViolation(
            time_s=float(time_s[sample_index]),
            current_a=float(current_model[sample_index]),
            sample_index=int(sample_index),
            state_name=state_name,
            electrode=electrode,
            value=float(domain_values[sample_index, value_index]),
        )
        valid_count = int(sample_index)

    if valid_count:
        valid_slice = slice(0, valid_count)
        model_current = current_model[valid_slice]
        x_n_surface = x_n_surface_all[valid_slice]
        x_p_surface = x_p_surface_all[valid_slice]
        u_n = np.asarray(U_n(x_n_surface), dtype=float)
        u_p = np.asarray(U_p(x_p_surface), dtype=float)
        kinetic_factor = (
            2.0
            * GAS_CONSTANT_J_PER_MOL_K
            * temperature_k[valid_slice]
            / FARADAY_C_PER_MOL
        )
        denominator_n = (
            6.0
            * parameters.q_n_c
            * parameters.d_n
            * np.sqrt(x_n_surface * (1.0 - x_n_surface))
        )
        denominator_p = (
            6.0
            * parameters.q_p_c
            * parameters.d_p
            * np.sqrt(x_p_surface * (1.0 - x_p_surface))
        )
        eta_n = kinetic_factor * np.arcsinh(model_current / denominator_n)
        eta_p = kinetic_factor * np.arcsinh(model_current / denominator_p)
        ohmic = parameters.r0_ohm * model_current
        voltage = u_p - u_n + eta_p - eta_n - ohmic
        evaluated = np.column_stack((u_n, u_p, eta_n, eta_p, ohmic, voltage))
        nonfinite = np.argwhere(~np.isfinite(evaluated))
        if nonfinite.size:
            local_index = int(nonfinite[0, 0])
            violation = DomainViolation(
                time_s=float(time_s[local_index]),
                current_a=float(current_model[local_index]),
                sample_index=local_index,
                state_name="voltage_evaluation",
                electrode="cell",
                value=float(voltage[local_index]),
            )
            valid_count = local_index
            valid_slice = slice(0, valid_count)
            u_n = u_n[valid_slice]
            u_p = u_p[valid_slice]
            eta_n = eta_n[valid_slice]
            eta_p = eta_p[valid_slice]
            ohmic = ohmic[valid_slice]
            voltage = voltage[valid_slice]

        fields["x_n_avg"][valid_slice] = states[valid_slice, 0]
        fields["x_n_aux"][valid_slice] = states[valid_slice, 1]
        fields["x_p_avg"][valid_slice] = states[valid_slice, 2]
        fields["x_p_aux"][valid_slice] = states[valid_slice, 3]
        fields["x_n_surface"][valid_slice] = x_n_surface_all[valid_slice]
        fields["x_p_surface"][valid_slice] = x_p_surface_all[valid_slice]
        fields["u_n"][valid_slice] = u_n
        fields["u_p"][valid_slice] = u_p
        fields["eta_n"][valid_slice] = eta_n
        fields["eta_p"][valid_slice] = eta_p
        fields["ohmic"][valid_slice] = ohmic
        fields["voltage"][valid_slice] = voltage

    return GroupedSPMResult(
        time_s=time_s.copy(),
        current_a=current_model.copy(),
        voltage_v=fields["voltage"],
        x_n_avg=fields["x_n_avg"],
        x_n_aux=fields["x_n_aux"],
        x_n_surface=fields["x_n_surface"],
        x_p_avg=fields["x_p_avg"],
        x_p_aux=fields["x_p_aux"],
        x_p_surface=fields["x_p_surface"],
        U_n_v=fields["u_n"],
        U_p_v=fields["u_p"],
        eta_n_v=fields["eta_n"],
        eta_p_v=fields["eta_p"],
        ohmic_drop_v=fields["ohmic"],
        temperature_k=temperature_k.copy(),
        domain_valid=violation is None,
        first_domain_violation=violation,
    )


def charge_conservation_error_c(
    result: GroupedSPMResult,
    parameters: GroupedSPMParameters,
) -> np.ndarray:
    """Return Q_n Delta x_n + Q_p Delta x_p in coulombs."""

    return (
        parameters.q_n_c * (result.x_n_avg - result.x_n_avg[0])
        + parameters.q_p_c * (result.x_p_avg - result.x_p_avg[0])
    )

