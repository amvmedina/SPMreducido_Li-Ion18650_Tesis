"""Frozen effective Panasonic lumped thermal parameterization.

These values were fitted only to the NCR18650PF 25 °C 1C discharge over
0.10 <= cell SOC <= 0.90 using measured-voltage irreversible heat and the
verified full-cell entropy coefficient. They are effective lumped parameters,
not independently measured material properties or globally identified values.

The generic thermal model and its defaults are not modified by this module.
"""

from __future__ import annotations

from types import MappingProxyType


C_TH_J_PER_K = 117.52317800099644
G_TH_W_PER_K = 0.10911208974639978
TAU_TH_S = C_TH_J_PER_K / G_TH_W_PER_K

CALIBRATION_SOC_RANGE = (0.10, 0.90)
CALIBRATION_HEAT_INPUT = "measured_voltage_irreversible_plus_reversible"
CALIBRATION_SOURCE = "Panasonic NCR18650PF 25 C 1C discharge"

EFFECTIVE_LUMPED_THERMAL_PARAMETERS = MappingProxyType(
    {
        "Cth_J_per_K": C_TH_J_PER_K,
        "Gth_W_per_K": G_TH_W_PER_K,
        "tau_th_s": TAU_TH_S,
    }
)

assert 5.0 < C_TH_J_PER_K < 200.0
assert 1.0e-4 < G_TH_W_PER_K < 1.0


__all__ = [
    "CALIBRATION_HEAT_INPUT",
    "CALIBRATION_SOC_RANGE",
    "CALIBRATION_SOURCE",
    "C_TH_J_PER_K",
    "EFFECTIVE_LUMPED_THERMAL_PARAMETERS",
    "G_TH_W_PER_K",
    "TAU_TH_S",
]
