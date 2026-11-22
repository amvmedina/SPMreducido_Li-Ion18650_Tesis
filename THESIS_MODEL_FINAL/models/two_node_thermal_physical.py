"""Clean physical two-node core/surface thermal model.

This module preserves the legacy physical parameterization while removing
all SOC/regime multipliers, heat clipping, sign inference, additive heat,
and temperature-bias corrections. Inputs are used exactly as supplied.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class TwoNodeThermalParameters:
    """Legacy-equivalent physical parameters, all strictly positive."""

    c_c_j_per_k: float
    c_s_j_per_k: float
    g_cs_w_per_k: float
    hA_w_per_k: float


@dataclass(frozen=True)
class TwoNodeThermalResult:
    time_s: np.ndarray
    q_gen_w: np.ndarray
    t_ambient_k: np.ndarray
    t_core_k: np.ndarray
    t_surface_k: np.ndarray


def validate_parameters(parameters: TwoNodeThermalParameters) -> None:
    values = np.asarray(
        [
            parameters.c_c_j_per_k,
            parameters.c_s_j_per_k,
            parameters.g_cs_w_per_k,
            parameters.hA_w_per_k,
        ],
        dtype=float,
    )
    if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("All two-node physical parameters must be finite and positive.")


def thermal_eigenvalues(parameters: TwoNodeThermalParameters) -> np.ndarray:
    """Return the two continuous-time state-matrix eigenvalues [1/s]."""

    validate_parameters(parameters)
    cc = parameters.c_c_j_per_k
    cs = parameters.c_s_j_per_k
    g = parameters.g_cs_w_per_k
    h = parameters.hA_w_per_k
    matrix = np.asarray([[-g / cc, g / cc], [g / cs, -(g + h) / cs]])
    return np.sort(np.linalg.eigvals(matrix).real)


def thermal_time_constants_s(parameters: TwoNodeThermalParameters) -> np.ndarray:
    """Return the fast then slow positive modal time constants [s]."""

    eigenvalues = thermal_eigenvalues(parameters)
    return np.sort(-1.0 / eigenvalues)


def simulate_two_node_physical(
    time_s: np.ndarray,
    q_gen_w: np.ndarray,
    t_ambient_k: np.ndarray | float,
    parameters: TwoNodeThermalParameters,
    initial_core_temperature_k: float,
    initial_surface_temperature_k: float,
) -> TwoNodeThermalResult:
    """Propagate the physical two-node model with exact ZOH interval updates.

    The equations are

    ``Cc*dTc/dt = Qgen - Gcs*(Tc-Ts)`` and
    ``Cs*dTs/dt = Gcs*(Tc-Ts) - hA*(Ts-Tamb)``.

    Heat and ambient are held at the left sample over each actual interval.
    There are deliberately no empirical correction inputs.
    """

    validate_parameters(parameters)
    time = np.asarray(time_s, dtype=float)
    heat = np.asarray(q_gen_w, dtype=float)
    if time.ndim != 1 or heat.ndim != 1 or time.shape != heat.shape:
        raise ValueError("time_s and q_gen_w must be matching one-dimensional arrays.")
    if time.size < 2 or not np.all(np.isfinite(time)) or np.any(np.diff(time) <= 0.0):
        raise ValueError("At least two finite, strictly increasing timestamps are required.")
    if not np.all(np.isfinite(heat)):
        raise ValueError("Heat generation must be finite.")
    if np.isscalar(t_ambient_k):
        ambient = np.full(time.shape, float(t_ambient_k))
    else:
        ambient = np.asarray(t_ambient_k, dtype=float)
    if ambient.shape != time.shape or not np.all(np.isfinite(ambient)):
        raise ValueError("Ambient temperature must be finite and scalar or time-shaped.")
    if not np.isfinite(initial_core_temperature_k) or not np.isfinite(
        initial_surface_temperature_k
    ):
        raise ValueError("Initial temperatures must be finite.")

    cc = parameters.c_c_j_per_k
    cs = parameters.c_s_j_per_k
    g = parameters.g_cs_w_per_k
    h = parameters.hA_w_per_k
    a00 = -g / cc
    a01 = g / cc
    a10 = g / cs
    a11 = -(g + h) / cs
    half_trace = 0.5 * (a00 + a11)
    half_difference = 0.5 * (a00 - a11)
    spectral_radius = float(np.sqrt(half_difference**2 + a01 * a10))
    lambda_slow = half_trace + spectral_radius
    lambda_fast = half_trace - spectral_radius
    eigenvalue_gap = lambda_slow - lambda_fast

    core = np.empty_like(time)
    surface = np.empty_like(time)
    core[0] = float(initial_core_temperature_k)
    surface[0] = float(initial_surface_temperature_k)

    for index, dt in enumerate(np.diff(time)):
        # With positive h and g, the constant-input equilibrium is explicit.
        surface_eq = ambient[index] + heat[index] / h
        core_eq = surface_eq + heat[index] / g
        slow_decay = np.exp(lambda_slow * dt)
        fast_decay = np.exp(lambda_fast * dt)
        m00 = (
            slow_decay * (a00 - lambda_fast)
            - fast_decay * (a00 - lambda_slow)
        ) / eigenvalue_gap
        m01 = a01 * (slow_decay - fast_decay) / eigenvalue_gap
        m10 = a10 * (slow_decay - fast_decay) / eigenvalue_gap
        m11 = (
            slow_decay * (a11 - lambda_fast)
            - fast_decay * (a11 - lambda_slow)
        ) / eigenvalue_gap
        delta_core = core[index] - core_eq
        delta_surface = surface[index] - surface_eq
        core[index + 1] = core_eq + m00 * delta_core + m01 * delta_surface
        surface[index + 1] = surface_eq + m10 * delta_core + m11 * delta_surface

    if not np.all(np.isfinite(core)) or not np.all(np.isfinite(surface)):
        raise ValueError("Two-node propagation produced non-finite temperatures.")
    return TwoNodeThermalResult(
        time_s=time,
        q_gen_w=heat,
        t_ambient_k=ambient,
        t_core_k=core,
        t_surface_k=surface,
    )


__all__ = [
    "TwoNodeThermalParameters",
    "TwoNodeThermalResult",
    "simulate_two_node_physical",
    "thermal_eigenvalues",
    "thermal_time_constants_s",
    "validate_parameters",
]
