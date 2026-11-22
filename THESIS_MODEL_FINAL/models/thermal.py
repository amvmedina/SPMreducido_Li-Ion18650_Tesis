"""Reduced heat generation and first-order lumped thermal dynamics.

All heat rates are positive when they add heat to the cell.  The electrical
expressions use the repository-wide convention ``I > 0`` for discharge.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class LumpedThermalParameters:
    """Cell-level thermal parameters with explicit SI units."""

    c_th_j_per_k: float
    g_th_w_per_k: float
    t_initial_k: float
    t_ambient_k: float

    def __post_init__(self) -> None:
        values = (
            self.c_th_j_per_k,
            self.g_th_w_per_k,
            self.t_initial_k,
            self.t_ambient_k,
        )
        if not all(np.isfinite(values)):
            raise ValueError("Thermal parameters must be finite.")
        if self.c_th_j_per_k <= 0.0:
            raise ValueError("Thermal capacitance must be positive.")
        if self.g_th_w_per_k < 0.0:
            raise ValueError("Thermal conductance must be nonnegative.")
        if self.t_initial_k <= 0.0 or self.t_ambient_k <= 0.0:
            raise ValueError("Absolute temperatures must be positive.")


@dataclass(frozen=True)
class HeatGeneration:
    """Instantaneous cell heat-rate contributions in watts."""

    irreversible_w: float
    ohmic_w: float
    kinetic_w: float
    reversible_w: float
    total_w: float
    irreversible_nonnegative: bool


@dataclass(frozen=True)
class LumpedThermalResult:
    """Temperature trajectory produced by a prescribed heat-rate trajectory."""

    time_s: np.ndarray
    heat_generation_w: np.ndarray
    temperature_k: np.ndarray
    ambient_temperature_k: np.ndarray


def calculate_heat_generation(
    *,
    current_a: float,
    equilibrium_voltage_v: float,
    terminal_voltage_v: float,
    ohmic_drop_v: float,
    eta_n_v: float,
    eta_p_v: float,
    temperature_k: float,
    d_u_eq_d_t_v_per_k: float | None = None,
) -> HeatGeneration:
    """Calculate loss and optional entropic heat for discharge-positive current.

    The irreversible term is obtained directly from electrical power loss,
    ``I * (U_eq - V)``.  The decomposition follows the implemented voltage
    equation, so the kinetic term is ``I * (eta_n - eta_p)`` and the ohmic
    term is ``I * (R0 I)``.  If an entropic coefficient is supplied, the
    reversible contribution is ``-I*T*dU_eq/dT``; otherwise it is zero.
    """

    coefficient = 0.0 if d_u_eq_d_t_v_per_k is None else d_u_eq_d_t_v_per_k
    values = (
        current_a,
        equilibrium_voltage_v,
        terminal_voltage_v,
        ohmic_drop_v,
        eta_n_v,
        eta_p_v,
        temperature_k,
        coefficient,
    )
    if not all(np.isfinite(values)):
        raise ValueError("Heat-generation inputs must be finite.")
    if temperature_k <= 0.0:
        raise ValueError("Absolute temperature must be positive.")

    irreversible_w = current_a * (equilibrium_voltage_v - terminal_voltage_v)
    ohmic_w = current_a * ohmic_drop_v
    kinetic_w = current_a * (eta_n_v - eta_p_v)
    reversible_w = -current_a * temperature_k * coefficient
    total_w = irreversible_w + reversible_w
    tolerance_w = 1.0e-12 * max(1.0, abs(irreversible_w))
    return HeatGeneration(
        irreversible_w=float(irreversible_w),
        ohmic_w=float(ohmic_w),
        kinetic_w=float(kinetic_w),
        reversible_w=float(reversible_w),
        total_w=float(total_w),
        irreversible_nonnegative=bool(irreversible_w >= -tolerance_w),
    )


def lumped_thermal_step(
    temperature_k: float,
    heat_generation_w: float,
    dt_s: float,
    parameters: LumpedThermalParameters,
) -> float:
    """Advance one ZOH interval using the analytical first-order solution."""

    values = (temperature_k, heat_generation_w, dt_s)
    if not all(np.isfinite(values)):
        raise ValueError("Thermal step inputs must be finite.")
    if temperature_k <= 0.0:
        raise ValueError("Absolute temperature must be positive.")
    if dt_s < 0.0:
        raise ValueError("dt must be nonnegative.")
    if dt_s == 0.0:
        return float(temperature_k)

    if parameters.g_th_w_per_k == 0.0:
        return float(
            temperature_k
            + heat_generation_w * dt_s / parameters.c_th_j_per_k
        )

    decay_argument = (
        -parameters.g_th_w_per_k * dt_s / parameters.c_th_j_per_k
    )
    approach_fraction = -np.expm1(decay_argument)
    equilibrium_temperature_k = (
        parameters.t_ambient_k
        + heat_generation_w / parameters.g_th_w_per_k
    )
    return float(
        temperature_k
        + approach_fraction * (equilibrium_temperature_k - temperature_k)
    )


def simulate_lumped_thermal(
    time: np.ndarray,
    heat_generation: np.ndarray,
    parameters: LumpedThermalParameters,
) -> LumpedThermalResult:
    """Simulate temperature with left-sample ZOH heat generation."""

    time_s = np.asarray(time, dtype=float)
    heat_generation_w = np.asarray(heat_generation, dtype=float)
    if time_s.ndim != 1 or time_s.size == 0:
        raise ValueError("time must be a nonempty one-dimensional array.")
    if heat_generation_w.shape != time_s.shape:
        raise ValueError("heat_generation must have the same shape as time.")
    if np.any(~np.isfinite(time_s)) or np.any(~np.isfinite(heat_generation_w)):
        raise ValueError("time and heat_generation must be finite.")
    if np.any(np.diff(time_s) < 0.0):
        raise ValueError("time must be nondecreasing.")

    temperature_k = np.empty(time_s.size, dtype=float)
    temperature_k[0] = parameters.t_initial_k
    for index, dt_s in enumerate(np.diff(time_s)):
        temperature_k[index + 1] = lumped_thermal_step(
            temperature_k[index],
            heat_generation_w[index],
            float(dt_s),
            parameters,
        )

    return LumpedThermalResult(
        time_s=time_s.copy(),
        heat_generation_w=heat_generation_w.copy(),
        temperature_k=temperature_k,
        ambient_temperature_k=np.full(time_s.size, parameters.t_ambient_k),
    )
