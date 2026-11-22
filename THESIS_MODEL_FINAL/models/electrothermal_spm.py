"""Cell-level SOC and causal electrochemical-thermal SPM coupling."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from .grouped_spm import (
    CurrentConvention,
    DomainViolation,
    GroupedSPMParameters,
    grouped_overpotentials,
    simulate_grouped_spm,
)
from .thermal import (
    LumpedThermalParameters,
    calculate_heat_generation,
    lumped_thermal_step,
)


EntropicCoefficient = Callable[[float, float, float], float]


@dataclass(frozen=True)
class CellSOCParameters:
    """Independent cell-level coulomb-counting configuration."""

    q_cell_ah: float
    soc_cell_0: float

    def __post_init__(self) -> None:
        if not np.isfinite(self.q_cell_ah) or self.q_cell_ah <= 0.0:
            raise ValueError("Cell capacity must be finite and positive.")
        if not np.isfinite(self.soc_cell_0):
            raise ValueError("Initial cell SOC must be finite.")


@dataclass(frozen=True)
class SOCViolation:
    """First un-clipped cell SOC value outside the physical interval."""

    time_s: float
    current_a: float
    sample_index: int
    value: float


@dataclass(frozen=True)
class ElectrothermalSPMResult:
    """Explicit electrical, SOC, heat, thermal, and electrode outputs."""

    time_s: np.ndarray
    current_a: np.ndarray
    terminal_voltage_v: np.ndarray
    cell_soc: np.ndarray
    cell_soc_valid: np.ndarray
    cell_soc_domain_valid: bool
    first_soc_violation: SOCViolation | None
    temperature_k: np.ndarray
    heat_generation_total_w: np.ndarray
    heat_irreversible_w: np.ndarray
    heat_ohmic_w: np.ndarray
    heat_kinetic_w: np.ndarray
    heat_reversible_w: np.ndarray
    heat_irreversible_nonnegative: np.ndarray
    u_eq_v: np.ndarray
    U_n_v: np.ndarray
    U_p_v: np.ndarray
    eta_n_v: np.ndarray
    eta_p_v: np.ndarray
    ohmic_drop_v: np.ndarray
    x_n_avg: np.ndarray
    x_n_aux: np.ndarray
    x_n_surface: np.ndarray
    x_p_avg: np.ndarray
    x_p_aux: np.ndarray
    x_p_surface: np.ndarray
    electrochemical_domain_valid: bool
    first_electrochemical_domain_violation: DomainViolation | None
    reversible_heat_enabled: bool


def _cell_soc_trajectory(
    time_s: np.ndarray,
    discharge_positive_current_a: np.ndarray,
    parameters: CellSOCParameters,
) -> tuple[np.ndarray, np.ndarray, SOCViolation | None]:
    cell_soc = np.empty(time_s.size, dtype=float)
    cell_soc[0] = parameters.soc_cell_0
    scale_c = 3600.0 * parameters.q_cell_ah
    for index, dt_s in enumerate(np.diff(time_s)):
        cell_soc[index + 1] = (
            cell_soc[index]
            - discharge_positive_current_a[index] * float(dt_s) / scale_c
        )
    cell_soc_valid = np.isfinite(cell_soc) & (cell_soc >= 0.0) & (cell_soc <= 1.0)
    invalid = np.flatnonzero(~cell_soc_valid)
    violation = None
    if invalid.size:
        index = int(invalid[0])
        violation = SOCViolation(
            time_s=float(time_s[index]),
            current_a=float(discharge_positive_current_a[index]),
            sample_index=index,
            value=float(cell_soc[index]),
        )
    return cell_soc, cell_soc_valid, violation


def simulate_electrothermal_spm(
    time: np.ndarray,
    current: np.ndarray,
    electrochemical_parameters: GroupedSPMParameters,
    soc_parameters: CellSOCParameters,
    thermal_parameters: LumpedThermalParameters,
    *,
    current_convention: CurrentConvention = "discharge_positive",
    entropic_coefficient: float | EntropicCoefficient | None = None,
) -> ElectrothermalSPMResult:
    """Run the grouped SPM with predicted temperature fed back at each sample.

    A callable entropic coefficient receives ``(x_n_surface, x_p_surface,
    temperature_k)`` and must return ``dU_eq/dT`` in V/K.  No reversible heat
    is included when ``entropic_coefficient`` is ``None``.
    """

    # The solid-state dynamics are temperature independent in grouped_spm.py.
    # Run them once through the existing implementation, then causally evaluate
    # its kinetic equation at each predicted temperature.
    electrical = simulate_grouped_spm(
        time,
        current,
        electrochemical_parameters,
        temperature=thermal_parameters.t_initial_k,
        current_convention=current_convention,
    )
    time_s = electrical.time_s
    current_a = electrical.current_a
    cell_soc, cell_soc_valid, soc_violation = _cell_soc_trajectory(
        time_s, current_a, soc_parameters
    )

    temperature_k = np.full(time_s.size, np.nan, dtype=float)
    temperature_k[0] = thermal_parameters.t_initial_k
    float_fields = {
        name: np.full(time_s.size, np.nan, dtype=float)
        for name in (
            "voltage",
            "u_eq",
            "eta_n",
            "eta_p",
            "ohmic",
            "heat_total",
            "heat_irreversible",
            "heat_ohmic",
            "heat_kinetic",
            "heat_reversible",
        )
    }
    heat_nonnegative = np.zeros(time_s.size, dtype=bool)

    for index in range(time_s.size):
        state_values = (
            electrical.x_n_surface[index],
            electrical.x_p_surface[index],
            electrical.U_n_v[index],
            electrical.U_p_v[index],
            temperature_k[index],
        )
        if not all(np.isfinite(state_values)):
            break

        eta_n_v, eta_p_v = grouped_overpotentials(
            current_a[index],
            electrical.x_n_surface[index],
            electrical.x_p_surface[index],
            electrochemical_parameters,
            temperature_k[index],
        )
        u_eq_v = electrical.U_p_v[index] - electrical.U_n_v[index]
        ohmic_drop_v = electrochemical_parameters.r0_ohm * current_a[index]
        terminal_voltage_v = (
            u_eq_v + eta_p_v - eta_n_v - ohmic_drop_v
        )

        if entropic_coefficient is None:
            coefficient_v_per_k = None
        elif callable(entropic_coefficient):
            coefficient_v_per_k = float(
                entropic_coefficient(
                    electrical.x_n_surface[index],
                    electrical.x_p_surface[index],
                    temperature_k[index],
                )
            )
        else:
            coefficient_v_per_k = float(entropic_coefficient)

        heat = calculate_heat_generation(
            current_a=current_a[index],
            equilibrium_voltage_v=u_eq_v,
            terminal_voltage_v=terminal_voltage_v,
            ohmic_drop_v=ohmic_drop_v,
            eta_n_v=eta_n_v,
            eta_p_v=eta_p_v,
            temperature_k=temperature_k[index],
            d_u_eq_d_t_v_per_k=coefficient_v_per_k,
        )
        float_fields["voltage"][index] = terminal_voltage_v
        float_fields["u_eq"][index] = u_eq_v
        float_fields["eta_n"][index] = eta_n_v
        float_fields["eta_p"][index] = eta_p_v
        float_fields["ohmic"][index] = ohmic_drop_v
        float_fields["heat_total"][index] = heat.total_w
        float_fields["heat_irreversible"][index] = heat.irreversible_w
        float_fields["heat_ohmic"][index] = heat.ohmic_w
        float_fields["heat_kinetic"][index] = heat.kinetic_w
        float_fields["heat_reversible"][index] = heat.reversible_w
        heat_nonnegative[index] = heat.irreversible_nonnegative

        if index + 1 < time_s.size:
            temperature_k[index + 1] = lumped_thermal_step(
                temperature_k[index],
                heat.total_w,
                float(time_s[index + 1] - time_s[index]),
                thermal_parameters,
            )

    return ElectrothermalSPMResult(
        time_s=time_s.copy(),
        current_a=current_a.copy(),
        terminal_voltage_v=float_fields["voltage"],
        cell_soc=cell_soc,
        cell_soc_valid=cell_soc_valid,
        cell_soc_domain_valid=bool(np.all(cell_soc_valid)),
        first_soc_violation=soc_violation,
        temperature_k=temperature_k,
        heat_generation_total_w=float_fields["heat_total"],
        heat_irreversible_w=float_fields["heat_irreversible"],
        heat_ohmic_w=float_fields["heat_ohmic"],
        heat_kinetic_w=float_fields["heat_kinetic"],
        heat_reversible_w=float_fields["heat_reversible"],
        heat_irreversible_nonnegative=heat_nonnegative,
        u_eq_v=float_fields["u_eq"],
        U_n_v=electrical.U_n_v.copy(),
        U_p_v=electrical.U_p_v.copy(),
        eta_n_v=float_fields["eta_n"],
        eta_p_v=float_fields["eta_p"],
        ohmic_drop_v=float_fields["ohmic"],
        x_n_avg=electrical.x_n_avg.copy(),
        x_n_aux=electrical.x_n_aux.copy(),
        x_n_surface=electrical.x_n_surface.copy(),
        x_p_avg=electrical.x_p_avg.copy(),
        x_p_aux=electrical.x_p_aux.copy(),
        x_p_surface=electrical.x_p_surface.copy(),
        electrochemical_domain_valid=electrical.domain_valid,
        first_electrochemical_domain_violation=electrical.first_domain_violation,
        reversible_heat_enabled=entropic_coefficient is not None,
    )
