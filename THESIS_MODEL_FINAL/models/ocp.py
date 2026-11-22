"""Electrode open-circuit-potential functions used by the grouped SPM.

The equations reproduce the public FrankSuperG/CPG-SPMT implementation at
Git commit 2b396b93713a9ab908fca27f4e4ff3b9a9bf9343. They are chemistry-specific
and are not generic lithium-ion OCP functions. Their valid model domain is
0 < x < 1; the simulator checks this domain without clipping.
"""

from __future__ import annotations

import numpy as np


SOURCE_URL = "https://github.com/FrankSuperG/CPG-SPMT/blob/main/cpg_spmt.py"
SOURCE_BLOB_SHA = "2b396b93713a9ab908fca27f4e4ff3b9a9bf9343"


def U_p(x_p: float | np.ndarray) -> float | np.ndarray:
    """Return positive-electrode OCP in volts."""

    x_p_array = np.asarray(x_p, dtype=float)
    value = (
        1.1127 * np.exp(-1.9474 - 2000.6182 * x_p_array)
        + 0.0388 * np.exp(-7.5995 * x_p_array)
        + 0.1220 * np.exp(-52.9452 * x_p_array)
        - 0.6990 * np.tanh(15.7948 * (x_p_array - 1.0))
        - 0.6353 * x_p_array**87.1117
        + 2.7313
    )
    return float(value) if value.ndim == 0 else value


def U_n(x_n: float | np.ndarray) -> float | np.ndarray:
    """Return negative-electrode OCP in volts."""

    x_n_array = np.asarray(x_n, dtype=float)
    value = (
        55.26
        + 9.451 * np.exp(-146.4 * x_n_array)
        - 0.05217 * np.tanh(9.453 * x_n_array - 1.783)
        + 0.1051 * np.tanh(4.452 * x_n_array - 5.533)
        - 0.07638 * np.tanh(19.25 * x_n_array - 18.04)
        - 0.02072 * np.tanh(25.87 * x_n_array - 14.56)
        - 54.9 * np.tanh(2.106 * x_n_array + 26.37)
        + 0.01532 * np.tanh(63.36 * x_n_array - 56.77)
        - 0.111 * np.tanh(52.65 * x_n_array - 2.862)
        - 0.0496 * np.tanh(77.08 * x_n_array - 2.383)
    )
    return float(value) if value.ndim == 0 else value


ocp_pos = U_p
ocp_neg = U_n

