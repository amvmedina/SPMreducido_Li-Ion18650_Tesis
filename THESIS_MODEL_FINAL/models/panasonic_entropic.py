"""Verified full-cell entropy coefficient for Panasonic NCR18650PF.

Source: He et al. (2024) Supplementary Information, Table S2, which cites
Lu et al. (2019). The tabulation is a full-cell coefficient versus average
cell SOC, not separate positive- and negative-electrode coefficients.

All public coefficient values and function outputs use V/K. The literature
values were supplied in mV/K and are converted exactly by a factor of 1e-3.
Linear interpolation is supported only for 0.10 <= SOC <= 0.90. There is no
clipping or extrapolation.
"""

from __future__ import annotations

import numpy as np


SUPPORTED_SOC_RANGE = (0.10, 0.90)

# Ascending SOC is convenient for numpy.interp. These arrays contain exactly
# the values supplied from verified He et al. Table S2.
SOC_AVERAGE = np.asarray(
    [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90],
    dtype=float,
)
DU_EQ_DT_MV_PER_K = np.asarray(
    [0.069, -0.082, 0.058, 0.141, 0.169, 0.143, 0.065, 0.083, 0.083],
    dtype=float,
)
DU_EQ_DT_V_PER_K = DU_EQ_DT_MV_PER_K * 1.0e-3

SOC_AVERAGE.setflags(write=False)
DU_EQ_DT_MV_PER_K.setflags(write=False)
DU_EQ_DT_V_PER_K.setflags(write=False)


def dUeq_dT_panasonic(soc_average: float | np.ndarray) -> float | np.ndarray:
    """Return the full-cell Panasonic entropy coefficient in V/K.

    Parameters
    ----------
    soc_average:
        Average cell SOC, restricted to the inclusive literature domain
        ``0.10 <= SOC <= 0.90``.

    Raises
    ------
    ValueError
        If the input is empty, non-finite, or outside the supported domain.
        Values are never clipped and the function never extrapolates.
    """

    soc = np.asarray(soc_average, dtype=float)
    scalar_output = soc.ndim == 0
    if soc.size == 0 or not np.all(np.isfinite(soc)):
        raise ValueError("Average cell SOC must be non-empty and finite.")
    lower, upper = SUPPORTED_SOC_RANGE
    if np.any((soc < lower) | (soc > upper)):
        raise ValueError(
            "Average cell SOC is outside the literature-supported domain "
            "0.10 <= SOC <= 0.90; no clipping or extrapolation is allowed."
        )
    coefficient = np.interp(soc, SOC_AVERAGE, DU_EQ_DT_V_PER_K)
    if scalar_output:
        return float(coefficient)
    return coefficient


__all__ = [
    "DU_EQ_DT_MV_PER_K",
    "DU_EQ_DT_V_PER_K",
    "SOC_AVERAGE",
    "SUPPORTED_SOC_RANGE",
    "dUeq_dT_panasonic",
]
