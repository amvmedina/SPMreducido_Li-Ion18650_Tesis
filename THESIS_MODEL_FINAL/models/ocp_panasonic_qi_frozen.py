"""Frozen nominal Panasonic OCP parameterization selected from multistart 3.

The electrode endpoints below are *nominal effective stoichiometric windows*.
They are not uniquely identified physical electrode stoichiometries.

Negative OCP
------------
Fixed Guo--White (2005) MCMB discharge-process spline, provided by
``models.ocp_panasonic_qi.U_n``.

Positive OCP
------------
Qi-style monotonic linear interpolation using the 13 fixed coordinates and the
13 start-3 voltages recorded in ``panasonic_c20_qi_up_multistart.csv``.

This module performs forward evaluation only.  It contains no optimization,
dynamic, kinetic, resistance, diffusion, or thermal parameters.
"""

from __future__ import annotations

import numpy as np
from types import MappingProxyType

from models.ocp_panasonic_qi import (
    MCMB_SUPPORTED_STOICHIOMETRY_RANGE,
    NCA_WORKING_STOICHIOMETRY_RANGE,
    UP_KNOT_STOICHIOMETRIES as _UP_KNOT_STOICHIOMETRIES,
    U_n as _guo_white_mcmb_discharge_ocp,
    U_p as _interpolated_positive_ocp,
)


SOURCE_MULTISTART_START_ID = 3
Q_CELL_AH = 2.994979138519

# Nominal effective stoichiometric windows.  These are retained at the full
# precision recorded for successful multistart solution 3, not claimed as
# uniquely identified physical electrode stoichiometries.
X_N_START = 0.720451136448
X_N_END = 0.251849704066
X_P_START = 0.224027933455
X_P_END = 0.953477823621

UP_KNOT_STOICHIOMETRIES = np.array(_UP_KNOT_STOICHIOMETRIES, copy=True)
UP_KNOT_VOLTAGES_V = np.asarray(
    [
        4.3271226714,
        4.23395701159,
        4.16106937748,
        4.06723948864,
        3.99908742406,
        3.92819762226,
        3.83141441646,
        3.77104313682,
        3.72457004303,
        3.67243722837,
        3.5853445102,
        3.46621450891,
        3.34685441976,
    ],
    dtype=float,
)
UP_KNOT_STOICHIOMETRIES.setflags(write=False)
UP_KNOT_VOLTAGES_V.setflags(write=False)

NOMINAL_EFFECTIVE_STOICHIOMETRIC_WINDOWS = MappingProxyType({
    "x_n_start": X_N_START,
    "x_n_end": X_N_END,
    "x_p_start": X_P_START,
    "x_p_end": X_P_END,
})


def U_n(x_n: float | np.ndarray) -> float | np.ndarray:
    """Evaluate the fixed Guo--White MCMB discharge OCP [V vs Li/Li+]."""

    return _guo_white_mcmb_discharge_ocp(x_n)


def U_p(x_p: float | np.ndarray) -> float | np.ndarray:
    """Evaluate the frozen start-3 Qi-style positive OCP [V vs Li/Li+]."""

    return _interpolated_positive_ocp(x_p, UP_KNOT_VOLTAGES_V)


def electrode_stoichiometries(
    normalized_discharged_capacity: float | np.ndarray,
) -> tuple[float | np.ndarray, float | np.ndarray]:
    """Return the nominal effective electrode trajectories for 0 <= z <= 1."""

    z = np.asarray(normalized_discharged_capacity, dtype=float)
    scalar_output = z.ndim == 0
    if z.size == 0 or not np.all(np.isfinite(z)):
        raise ValueError("Normalized discharged capacity must be non-empty and finite.")
    if np.any((z < 0.0) | (z > 1.0)):
        raise ValueError("Normalized discharged capacity must satisfy 0 <= z <= 1; no clipping.")
    x_n = X_N_START + z * (X_N_END - X_N_START)
    x_p = X_P_START + z * (X_P_END - X_P_START)
    if scalar_output:
        return float(x_n), float(x_p)
    return x_n, x_p


def V_eq(normalized_discharged_capacity: float | np.ndarray) -> float | np.ndarray:
    """Evaluate the frozen quasi-equilibrium full-cell voltage [V]."""

    x_n, x_p = electrode_stoichiometries(normalized_discharged_capacity)
    return U_p(x_p) - U_n(x_n)


assert MCMB_SUPPORTED_STOICHIOMETRY_RANGE[0] <= X_N_END < X_N_START <= MCMB_SUPPORTED_STOICHIOMETRY_RANGE[1]
assert NCA_WORKING_STOICHIOMETRY_RANGE[0] <= X_P_START < X_P_END <= NCA_WORKING_STOICHIOMETRY_RANGE[1]
assert np.all(np.diff(UP_KNOT_VOLTAGES_V) <= 0.0)


__all__ = [
    "NOMINAL_EFFECTIVE_STOICHIOMETRIC_WINDOWS",
    "Q_CELL_AH",
    "SOURCE_MULTISTART_START_ID",
    "UP_KNOT_STOICHIOMETRIES",
    "UP_KNOT_VOLTAGES_V",
    "U_n",
    "U_p",
    "V_eq",
    "X_N_END",
    "X_N_START",
    "X_P_END",
    "X_P_START",
    "electrode_stoichiometries",
]
