"""OCP functions for the Panasonic Qi-style identification branch.

The negative electrode is the fixed Guo--White (2005) MCMB discharge spline.
The positive electrode is a monotonic piecewise-linear function on 13 fixed
stoichiometric knots; its voltages are supplied by the Panasonic C/20 fit.

No function in this module clips or extrapolates stoichiometry.
"""

from __future__ import annotations

import numpy as np


MCMB_SUPPORTED_STOICHIOMETRY_RANGE = (0.01, 0.80)
NCA_WORKING_STOICHIOMETRY_RANGE = (0.18, 0.96)
N_UP_KNOTS = 13
UP_KNOT_STOICHIOMETRIES = np.linspace(*NCA_WORKING_STOICHIOMETRY_RANGE, N_UP_KNOTS)
UP_VOLTAGE_BOUNDS_V = (2.5, 4.4)

GUO_MCMB_DISCHARGE_KNOTS = np.asarray(
    [0.47423, 0.16833, 0.078218, 0.015397, 0.0059697], dtype=float
)
GUO_MCMB_DISCHARGE_BASE_COEFFICIENTS = np.asarray(
    [-6.3949e-2, 1.2623, -2.3261, 1.2853], dtype=float
)
GUO_MCMB_DISCHARGE_SPLINE_COEFFICIENTS = np.asarray(
    [-3.7276, 1.4102e1, -1.1094e3, -1.0475e5, -9.4697e6], dtype=float
)


def _validated_stoichiometry(
    stoichiometry: float | np.ndarray,
    supported_range: tuple[float, float],
    electrode_name: str,
) -> tuple[np.ndarray, bool]:
    x = np.asarray(stoichiometry, dtype=float)
    scalar_output = x.ndim == 0
    if x.size == 0:
        raise ValueError("Stoichiometry must contain at least one value.")
    if not np.all(np.isfinite(x)):
        raise ValueError("Stoichiometry must be finite.")
    lower, upper = supported_range
    if np.any((x < lower) | (x > upper)):
        raise ValueError(
            f"{electrode_name} stoichiometry must remain inside the supported "
            f"interval {lower:.2f} <= x <= {upper:.2f}; observed range "
            f"[{float(np.min(x)):.9g}, {float(np.max(x)):.9g}]. No clipping or "
            "extrapolation is permitted."
        )
    return x, scalar_output


def U_n(x_n: float | np.ndarray) -> float | np.ndarray:
    """Guo--White MCMB discharge-process OCP [V vs Li/Li+]."""

    x, scalar_output = _validated_stoichiometry(
        x_n, MCMB_SUPPORTED_STOICHIOMETRY_RANGE, "MCMB negative-electrode"
    )
    a, b, c, d = GUO_MCMB_DISCHARGE_BASE_COEFFICIENTS
    result = a + b * x + c * x**2 + d * x**3
    for knot, coefficient in zip(
        GUO_MCMB_DISCHARGE_KNOTS,
        GUO_MCMB_DISCHARGE_SPLINE_COEFFICIENTS,
        strict=True,
    ):
        active = x <= knot
        result = result + active * coefficient * (x - knot) ** 3
    return float(result) if scalar_output else result


def validate_positive_knot_voltages(
    knot_voltages_v: np.ndarray | list[float] | tuple[float, ...],
) -> np.ndarray:
    """Validate the 13 bounded, non-increasing positive-OCP knot voltages."""

    voltages = np.asarray(knot_voltages_v, dtype=float)
    if voltages.shape != (N_UP_KNOTS,):
        raise ValueError(f"Exactly {N_UP_KNOTS} positive-OCP knot voltages are required.")
    if not np.all(np.isfinite(voltages)):
        raise ValueError("Positive-OCP knot voltages must be finite.")
    lower, upper = UP_VOLTAGE_BOUNDS_V
    if np.any((voltages < lower) | (voltages > upper)):
        raise ValueError(
            f"Positive-OCP knot voltages must remain inside {lower:.1f} to "
            f"{upper:.1f} V; no clipping is permitted."
        )
    if np.any(np.diff(voltages) > 1.0e-12):
        raise ValueError(
            "Positive-OCP knot voltages must be non-increasing as x_p increases."
        )
    return voltages


def U_p(
    x_p: float | np.ndarray,
    knot_voltages_v: np.ndarray | list[float] | tuple[float, ...],
) -> float | np.ndarray:
    """Evaluate the bounded, monotonic piecewise-linear positive OCP [V]."""

    x, scalar_output = _validated_stoichiometry(
        x_p, NCA_WORKING_STOICHIOMETRY_RANGE, "Positive-electrode"
    )
    voltages = validate_positive_knot_voltages(knot_voltages_v)
    result = np.interp(x, UP_KNOT_STOICHIOMETRIES, voltages)
    return float(result) if scalar_output else result


__all__ = [
    "GUO_MCMB_DISCHARGE_BASE_COEFFICIENTS",
    "GUO_MCMB_DISCHARGE_KNOTS",
    "GUO_MCMB_DISCHARGE_SPLINE_COEFFICIENTS",
    "MCMB_SUPPORTED_STOICHIOMETRY_RANGE",
    "NCA_WORKING_STOICHIOMETRY_RANGE",
    "N_UP_KNOTS",
    "UP_KNOT_STOICHIOMETRIES",
    "UP_VOLTAGE_BOUNDS_V",
    "U_n",
    "U_p",
    "validate_positive_knot_voltages",
]
