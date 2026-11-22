"""Public API for the reduced electrochemical-thermal cell model."""

from .electrothermal_spm import (
    CellSOCParameters,
    ElectrothermalSPMResult,
    EntropicCoefficient,
    SOCViolation,
    simulate_electrothermal_spm,
)

from .grouped_spm import (
    DEFAULT_TEMPERATURE_K,
    PARAMETER_ORDER,
    DomainViolation,
    GroupedSPMParameters,
    GroupedSPMResult,
    charge_conservation_error_c,
    continuous_state_matrices,
    exact_zoh_matrices,
    grouped_overpotentials,
    simulate_grouped_spm,
    surface_outputs,
)
from .ocp import U_n, U_p
from .thermal import (
    HeatGeneration,
    LumpedThermalParameters,
    LumpedThermalResult,
    calculate_heat_generation,
    lumped_thermal_step,
    simulate_lumped_thermal,
)

__all__ = [
    "DEFAULT_TEMPERATURE_K",
    "PARAMETER_ORDER",
    "CellSOCParameters",
    "DomainViolation",
    "ElectrothermalSPMResult",
    "EntropicCoefficient",
    "GroupedSPMParameters",
    "GroupedSPMResult",
    "HeatGeneration",
    "LumpedThermalParameters",
    "LumpedThermalResult",
    "SOCViolation",
    "U_n",
    "U_p",
    "charge_conservation_error_c",
    "calculate_heat_generation",
    "continuous_state_matrices",
    "exact_zoh_matrices",
    "grouped_overpotentials",
    "lumped_thermal_step",
    "simulate_electrothermal_spm",
    "simulate_grouped_spm",
    "simulate_lumped_thermal",
    "surface_outputs",
]
