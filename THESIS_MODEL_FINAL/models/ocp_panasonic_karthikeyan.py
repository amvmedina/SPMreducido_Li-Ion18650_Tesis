"""Karthikeyan et al. (2008) electrode equilibrium potentials.

This standalone module implements the Redlich--Kister thermodynamic model for
the two literature electrodes selected for the Panasonic analysis:

* positive: LiNi0.8Co0.15Al0.05O2 (NCA);
* negative: mesocarbon microbead (MCMB) graphite.

Source
------
Deepak K. Karthikeyan, Godfrey Sikha, and Ralph E. White,
"Thermodynamic model development for lithium intercalation electrodes,"
Journal of Power Sources 185 (2008) 1398--1407, equations (23)--(24).

The stable polynomial basis below corresponds to the Redlich--Kister term
whose direct form contains ``(2*x - 1)**(k - 1)``.  A transcription with the
exponent ``1 - k`` is not a removable singularity: for k >= 2 it diverges at
x = 0.5.  Using ``t = 2*x - 1`` and ``x*(1-x) = (1-t**2)/4``, the finite basis
is

    (1 + k/2) * t**(k+1) - (k/2) * t**(k-1),  k >= 1,

while the k=0 basis is simply t.  This is algebraically equivalent to the
finite Redlich--Kister expression away from t=0 and is continuous at t=0.

The low-level evaluator is mathematically available for ``0 < x < 1``.  The
Panasonic-facing ``U_n`` and ``U_p`` wrappers additionally enforce the
stoichiometric spans supported by the experimental data in the paper.  No
stoichiometry clipping or literature-domain extrapolation is performed.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np


GAS_CONSTANT_J_PER_MOL_K = 8.314
FARADAY_C_PER_MOL = 96485.0
DEFAULT_TEMPERATURE_K = 298.15

# Literature-supported fit domains.  The MCMB interval is the experimentally
# supported Karthikeyan range selected for this work.  For NCA, the
# experimental markers and fitted curve in Karthikeyan et al. Figs. 2 and 6
# span approximately x=0.18 through x=0.96; those plotted endpoints are used
# as hard bounds rather than extending the fit to the mathematical endpoints.
MCMB_GRAPHITE_LITERATURE_STOICHIOMETRY_RANGE = (0.01, 0.80)
NCA_LITERATURE_STOICHIOMETRY_RANGE = (0.18, 0.96)

NCA_U0_V = 3.9662
NCA_COEFFICIENTS_J_PER_MOL = (
    -5.8992e4,
    -2.0881e4,
    -1.3273e4,
    -6.9538e3,
    -2.6023e4,
    1.0715e4,
)

MCMB_GRAPHITE_U0_V = -1.7203
MCMB_GRAPHITE_COEFFICIENTS_J_PER_MOL = (
    -0.35799e6,
    -0.35008e6,
    -0.35247e6,
    -0.35692e6,
    -0.38633e6,
    -0.35908e6,
    -0.28794e6,
    -0.14979e6,
    -0.39912e6,
    -0.96172e6,
    -0.63262e6,
)


def _validate_literature_stoichiometry(
    stoichiometry: float | np.ndarray,
    supported_range: tuple[float, float],
    electrode_name: str,
) -> None:
    """Reject values outside an electrode's inclusive literature fit range."""

    x = np.asarray(stoichiometry, dtype=float)
    if x.size == 0:
        raise ValueError("Stoichiometry must contain at least one value.")
    if not np.all(np.isfinite(x)):
        raise ValueError("Stoichiometry must be finite.")

    lower, upper = supported_range
    if np.any((x < lower) | (x > upper)):
        raise ValueError(
            f"{electrode_name} stoichiometry is outside the literature-supported "
            f"interval {lower:.2f} <= x <= {upper:.2f}; observed range "
            f"[{float(np.min(x)):.9g}, {float(np.max(x)):.9g}]. No extrapolation "
            "or clipping is permitted."
        )


def _validated_inputs(
    stoichiometry: float | np.ndarray,
    temperature_k: float | np.ndarray,
) -> tuple[np.ndarray, np.ndarray, bool]:
    x = np.asarray(stoichiometry, dtype=float)
    temperature = np.asarray(temperature_k, dtype=float)
    scalar_output = x.ndim == 0 and temperature.ndim == 0
    try:
        x, temperature = np.broadcast_arrays(x, temperature)
    except ValueError as error:
        raise ValueError("Stoichiometry and temperature are not broadcastable.") from error

    if not np.all(np.isfinite(x)):
        raise ValueError("Stoichiometry must be finite.")
    if np.any((x <= 0.0) | (x >= 1.0)):
        raise ValueError(
            "Electrode stoichiometry must remain strictly inside 0 < x < 1; "
            f"observed range [{float(np.min(x)):.9g}, {float(np.max(x)):.9g}]."
        )
    if not np.all(np.isfinite(temperature)) or np.any(temperature <= 0.0):
        raise ValueError("Absolute temperature must be finite and positive.")
    return x, temperature, scalar_output


def _stable_redlich_kister_sum(
    x: np.ndarray,
    coefficients_j_per_mol: tuple[float, ...],
) -> np.ndarray:
    """Return the finite Redlich--Kister excess chemical potential [J/mol]."""

    t = 2.0 * x - 1.0
    total = coefficients_j_per_mol[0] * t
    for k, coefficient in enumerate(coefficients_j_per_mol[1:], start=1):
        basis = (1.0 + 0.5 * k) * t ** (k + 1) - 0.5 * k * t ** (k - 1)
        total = total + coefficient * basis
    return total


def redlich_kister_equilibrium_potential(
    stoichiometry: float | np.ndarray,
    temperature_k: float | np.ndarray,
    *,
    u0_v: float,
    coefficients_j_per_mol: tuple[float, ...],
) -> float | np.ndarray:
    """Evaluate a Karthikeyan potential on the mathematical domain 0 < x < 1.

    This low-level function deliberately does not impose an electrode-specific
    literature fit range.  Panasonic-facing calculations must call ``U_n`` or
    ``U_p``, which reject literature-domain extrapolation.
    """

    if not np.isfinite(u0_v):
        raise ValueError("u0_v must be finite.")
    if not coefficients_j_per_mol or not np.all(
        np.isfinite(np.asarray(coefficients_j_per_mol, dtype=float))
    ):
        raise ValueError("At least one finite Redlich--Kister coefficient is required.")

    x, temperature, scalar_output = _validated_inputs(stoichiometry, temperature_k)
    ideal_v = (
        GAS_CONSTANT_J_PER_MOL_K
        * temperature
        / FARADAY_C_PER_MOL
        * np.log((1.0 - x) / x)
    )
    excess_v = _stable_redlich_kister_sum(x, coefficients_j_per_mol) / FARADAY_C_PER_MOL
    result = u0_v + ideal_v + excess_v
    return float(result) if scalar_output else result


def U_p(
    x_p: float | np.ndarray,
    temperature_k: float | np.ndarray = DEFAULT_TEMPERATURE_K,
) -> float | np.ndarray:
    """NCA equilibrium potential within Karthikeyan's fitted data range."""

    _validate_literature_stoichiometry(
        x_p,
        NCA_LITERATURE_STOICHIOMETRY_RANGE,
        "NCA positive-electrode",
    )

    return redlich_kister_equilibrium_potential(
        x_p,
        temperature_k,
        u0_v=NCA_U0_V,
        coefficients_j_per_mol=NCA_COEFFICIENTS_J_PER_MOL,
    )


def U_n(
    x_n: float | np.ndarray,
    temperature_k: float | np.ndarray = DEFAULT_TEMPERATURE_K,
) -> float | np.ndarray:
    """MCMB equilibrium potential within Karthikeyan's fitted data range."""

    _validate_literature_stoichiometry(
        x_n,
        MCMB_GRAPHITE_LITERATURE_STOICHIOMETRY_RANGE,
        "MCMB graphite negative-electrode",
    )

    return redlich_kister_equilibrium_potential(
        x_n,
        temperature_k,
        u0_v=MCMB_GRAPHITE_U0_V,
        coefficients_j_per_mol=MCMB_GRAPHITE_COEFFICIENTS_J_PER_MOL,
    )


def save_diagnostic_plot(output_path: str | Path) -> Path:
    """Save the requested two-panel OCP diagnostic plot and return its path."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    x_n = np.linspace(*MCMB_GRAPHITE_LITERATURE_STOICHIOMETRY_RANGE, 2000)
    x_p = np.linspace(*NCA_LITERATURE_STOICHIOMETRY_RANGE, 2000)

    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.5), constrained_layout=True)
    axes[0].plot(x_n, U_n(x_n), color="tab:blue", linewidth=1.6)
    axes[0].set(
        xlabel="MCMB graphite stoichiometry, $x_n$",
        ylabel="$U_n$ [V vs Li/Li$^+$]",
        title="MCMB: literature range $0.01 \\leq x_n \\leq 0.80$",
    )
    axes[1].plot(x_p, U_p(x_p), color="tab:red", linewidth=1.6)
    axes[1].set(
        xlabel="NCA stoichiometry, $x_p$",
        ylabel="$U_p$ [V vs Li/Li$^+$]",
        title="NCA: literature range $0.18 \\leq x_p \\leq 0.96$",
    )
    for axis in axes:
        axis.grid(alpha=0.25)
    figure.suptitle(
        "Karthikeyan et al. (2008) OCPs: literature-supported domains only, 298.15 K"
    )
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path.resolve()


if __name__ == "__main__":
    default_output = (
        Path(__file__).resolve().parents[1]
        / "thesis_data"
        / "panasonic_karthikeyan_ocp_diagnostic.png"
    )
    print(save_diagnostic_plot(default_output))


__all__ = [
    "DEFAULT_TEMPERATURE_K",
    "FARADAY_C_PER_MOL",
    "GAS_CONSTANT_J_PER_MOL_K",
    "MCMB_GRAPHITE_COEFFICIENTS_J_PER_MOL",
    "MCMB_GRAPHITE_LITERATURE_STOICHIOMETRY_RANGE",
    "MCMB_GRAPHITE_U0_V",
    "NCA_COEFFICIENTS_J_PER_MOL",
    "NCA_LITERATURE_STOICHIOMETRY_RANGE",
    "NCA_U0_V",
    "U_n",
    "U_p",
    "redlich_kister_equilibrium_potential",
    "save_diagnostic_plot",
]
