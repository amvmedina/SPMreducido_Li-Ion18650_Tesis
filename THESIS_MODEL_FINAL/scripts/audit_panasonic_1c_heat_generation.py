"""Forward-only Panasonic 1C heat-generation component audit.

The frozen electrical solution is reproduced without fitting. Heat generation
is evaluated only for 0.10 <= average cell SOC <= 0.90 using the verified
full-cell literature entropy coefficient. Irreversible heat based on measured
and model voltage is retained as two explicitly different quantities.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from models.ocp_panasonic_qi import (  # noqa: E402
    MCMB_SUPPORTED_STOICHIOMETRY_RANGE,
    NCA_WORKING_STOICHIOMETRY_RANGE,
)
from models.ocp_panasonic_qi_frozen import (  # noqa: E402
    Q_CELL_AH,
    U_n,
    U_p,
    X_N_END,
    X_N_START,
    X_P_END,
    X_P_START,
)
from models.panasonic_entropic import (  # noqa: E402
    DU_EQ_DT_MV_PER_K,
    DU_EQ_DT_V_PER_K,
    SOC_AVERAGE,
    SUPPORTED_SOC_RANGE,
    dUeq_dT_panasonic,
)
from scripts.audit_panasonic_hppc_identifiability import (  # noqa: E402
    Q_N_C as IMPLEMENTED_Q_N_C,
    Q_P_C as IMPLEMENTED_Q_P_C,
    simulate_local_group,
)
from scripts.validate_panasonic_1c_forward import (  # noqa: E402
    NOMINAL_PARAMETERS,
    Q_N_C,
    Q_P_C,
    RAW_SOURCE_PATH,
    load_raw_discharge,
    longest_contiguous_true,
)


OUTPUT_DIRECTORY = REPOSITORY_ROOT / "thesis_data"
COEFFICIENT_CSV_PATH = OUTPUT_DIRECTORY / "panasonic_entropic_coefficient_literature.csv"
METRICS_CSV_PATH = OUTPUT_DIRECTORY / "panasonic_1c_heat_generation_metrics.csv"
REPORT_PATH = OUTPUT_DIRECTORY / "panasonic_1c_heat_generation_audit.md"
TIME_PLOT_PATH = OUTPUT_DIRECTORY / "panasonic_1c_heat_components_time.png"
SOC_PLOT_PATH = OUTPUT_DIRECTORY / "panasonic_1c_heat_components_soc.png"

DISCHARGE_CURRENT_THRESHOLD_A = -0.5
EXPECTED_FULL_FORWARD_RMSE_V = 0.11436218800084591


@dataclass(frozen=True)
class RestrictedProfile:
    time_s: np.ndarray
    soc: np.ndarray
    current_a: np.ndarray
    temperature_k: np.ndarray
    measured_voltage_v: np.ndarray
    model_voltage_v: np.ndarray
    x_n_avg: np.ndarray
    x_p_avg: np.ndarray


def interpolate_descending(
    target: float,
    coordinate_descending: np.ndarray,
    values: np.ndarray,
) -> float:
    if target < coordinate_descending[-1] or target > coordinate_descending[0]:
        raise ValueError("Requested interpolation point is outside measured SOC support.")
    return float(np.interp(target, coordinate_descending[::-1], values[::-1]))


def restrict_to_exact_soc_domain(
    soc: np.ndarray,
    **fields: np.ndarray,
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Preserve actual interior samples and interpolate exact supported endpoints."""

    lower, upper = SUPPORTED_SOC_RANGE
    if np.any(np.diff(soc) >= 0.0):
        raise ValueError("The discharge SOC coordinate must be strictly decreasing.")
    if not (soc[0] >= upper and soc[-1] <= lower):
        raise ValueError("The 1C discharge does not span the full entropy domain.")
    interior = (soc < upper) & (soc > lower)
    restricted_soc = np.r_[upper, soc[interior], lower]
    restricted: dict[str, np.ndarray] = {}
    for name, values in fields.items():
        values = np.asarray(values, dtype=float)
        if values.shape != soc.shape:
            raise ValueError(f"Field {name} does not match the SOC shape.")
        restricted[name] = np.r_[
            interpolate_descending(upper, soc, values),
            values[interior],
            interpolate_descending(lower, soc, values),
        ]
    if np.any(np.diff(restricted["time_s"]) <= 0.0):
        raise ValueError("Restricted timestamps are not strictly increasing.")
    return restricted_soc, restricted


def time_weighted_metrics(values_w: np.ndarray, time_s: np.ndarray) -> dict[str, float]:
    duration_s = float(time_s[-1] - time_s[0])
    if duration_s <= 0.0:
        raise ValueError("Heat metrics require positive duration.")
    energy_j = float(np.trapezoid(values_w, time_s))
    return {
        "min_w": float(np.min(values_w)),
        "max_w": float(np.max(values_w)),
        "mean_w": energy_j / duration_s,
        "rms_w": float(np.sqrt(np.trapezoid(values_w**2, time_s) / duration_s)),
        "integrated_energy_j": energy_j,
    }


def zero_crossings() -> np.ndarray:
    roots: list[float] = []
    for left_soc, right_soc, left_y, right_y in zip(
        SOC_AVERAGE[:-1],
        SOC_AVERAGE[1:],
        DU_EQ_DT_V_PER_K[:-1],
        DU_EQ_DT_V_PER_K[1:],
    ):
        if left_y == 0.0:
            roots.append(float(left_soc))
        elif left_y * right_y < 0.0:
            roots.append(
                float(left_soc - left_y * (right_soc - left_soc) / (right_y - left_y))
            )
    return np.asarray(roots)


def write_coefficient_csv(path: Path) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=(
                "soc_average",
                "dUeq_dT_mV_per_K",
                "dUeq_dT_V_per_K",
                "source",
            ),
        )
        writer.writeheader()
        for soc, mv_per_k, v_per_k in zip(
            SOC_AVERAGE[::-1],
            DU_EQ_DT_MV_PER_K[::-1],
            DU_EQ_DT_V_PER_K[::-1],
        ):
            writer.writerow(
                {
                    "soc_average": f"{soc:.2f}",
                    "dUeq_dT_mV_per_K": f"{mv_per_k:.6f}",
                    "dUeq_dT_V_per_K": f"{v_per_k:.9f}",
                    "source": "He et al. 2024 Supplementary Table S2; cites Lu et al. 2019",
                }
            )


def write_metrics_csv(
    path: Path,
    metrics: dict[str, dict[str, float]],
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        fields = (
            "component",
            "minimum_w",
            "maximum_w",
            "time_weighted_mean_w",
            "time_weighted_rms_w",
            "time_integrated_energy_j",
        )
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for name, values in metrics.items():
            writer.writerow(
                {
                    "component": name,
                    "minimum_w": f"{values['min_w']:.12g}",
                    "maximum_w": f"{values['max_w']:.12g}",
                    "time_weighted_mean_w": f"{values['mean_w']:.12g}",
                    "time_weighted_rms_w": f"{values['rms_w']:.12g}",
                    "time_integrated_energy_j": f"{values['integrated_energy_j']:.12g}",
                }
            )


def make_component_plot(
    path: Path,
    x: np.ndarray,
    x_label: str,
    components: dict[str, np.ndarray],
    *,
    invert_x: bool,
    roots: np.ndarray,
) -> None:
    colors = {
        "Qdot_irr_measured": "#0072B2",
        "Qdot_irr_model": "#D55E00",
        "Qdot_rev": "#009E73",
        "Qdot_total_measured": "#56B4E9",
        "Qdot_total_model": "#CC79A7",
    }
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 7.5), sharex=True, constrained_layout=True)
    for name in ("Qdot_irr_measured", "Qdot_irr_model", "Qdot_rev"):
        axes[0].plot(x, components[name], label=name, color=colors[name], linewidth=1.7)
    for name in ("Qdot_total_measured", "Qdot_total_model"):
        axes[1].plot(x, components[name], label=name, color=colors[name], linewidth=1.7)
    axes[0].axhline(0.0, color="black", linewidth=0.8)
    axes[1].axhline(0.0, color="black", linewidth=0.8)
    axes[0].set_ylabel("Heat rate (W)")
    axes[1].set_ylabel("Total heat rate (W)")
    axes[1].set_xlabel(x_label)
    axes[0].grid(alpha=0.25)
    axes[1].grid(alpha=0.25)
    axes[0].legend(loc="best", ncol=3)
    axes[1].legend(loc="best", ncol=2)
    if invert_x:
        for root in roots:
            axes[0].axvline(root, color="0.45", linestyle="--", linewidth=0.9)
            axes[1].axvline(root, color="0.45", linestyle="--", linewidth=0.9)
        axes[0].invert_xaxis()
    fig.suptitle("Panasonic 25 °C 1C heat-generation audit (SOC 0.10–0.90)")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def metric_rows(metrics: dict[str, dict[str, float]]) -> str:
    return "\n".join(
        f"| `{name}` | {values['min_w']:.9f} | {values['max_w']:.9f} | "
        f"{values['mean_w']:.9f} | {values['rms_w']:.9f} | "
        f"{values['integrated_energy_j']:.6f} |"
        for name, values in metrics.items()
    )


def main() -> None:
    if not np.isclose(Q_CELL_AH, 2.994979138519, rtol=0.0, atol=5.0e-13):
        raise RuntimeError("Frozen Panasonic capacity changed.")
    if not (
        np.isclose(IMPLEMENTED_Q_N_C, Q_N_C, rtol=0.0, atol=5.0e-9)
        and np.isclose(IMPLEMENTED_Q_P_C, Q_P_C, rtol=0.0, atol=5.0e-9)
    ):
        raise RuntimeError("Frozen electrode capacities changed.")

    columns, raw = load_raw_discharge(RAW_SOURCE_PATH)
    indices = longest_contiguous_true(raw["Current"] <= DISCHARGE_CURRENT_THRESHOLD_A)
    time_full = raw["Time"][indices]
    measured_voltage_full = raw["Voltage"][indices]
    raw_current_full = raw["Current"][indices]
    current_full = -raw_current_full
    temperature_full_k = raw["Battery_Temp_degC"][indices] + 273.15
    exported_ah = raw["Ah"][indices]
    source_soc_full = 1.0 - (exported_ah[0] - exported_ah) / Q_CELL_AH
    if np.any(np.diff(time_full) <= 0.0) or np.any(np.diff(source_soc_full) >= 0.0):
        raise ValueError("The active 1C time/SOC coordinates are not strictly monotonic.")
    if np.any(current_full <= 0.0):
        raise ValueError("Model-boundary current conversion did not yield discharge-positive current.")

    initial_x_n = X_N_START
    initial_x_p = X_P_START
    result = simulate_local_group(
        time_s=time_full,
        discharge_positive_current_a=current_full,
        temperature_k=temperature_full_k,
        initial_x_n=initial_x_n,
        initial_x_p=initial_x_p,
        dynamic_parameters=NOMINAL_PARAMETERS,
    )
    model_voltage_full = result.voltage_v
    full_rmse = float(np.sqrt(np.mean((model_voltage_full - measured_voltage_full) ** 2)))
    if not np.isclose(full_rmse, EXPECTED_FULL_FORWARD_RMSE_V, rtol=0.0, atol=5.0e-13):
        raise RuntimeError(
            "Frozen 1C forward prediction no longer reproduces the recorded validation RMSE."
        )

    soc, selected = restrict_to_exact_soc_domain(
        source_soc_full,
        time_s=time_full,
        current_a=current_full,
        temperature_k=temperature_full_k,
        measured_voltage_v=measured_voltage_full,
        model_voltage_v=model_voltage_full,
        x_n_avg=result.x_n_avg,
        x_p_avg=result.x_p_avg,
    )
    profile = RestrictedProfile(soc=soc, **selected)
    if not (
        np.all((profile.x_n_avg >= MCMB_SUPPORTED_STOICHIOMETRY_RANGE[0]) & (profile.x_n_avg <= MCMB_SUPPORTED_STOICHIOMETRY_RANGE[1]))
        and np.all((profile.x_p_avg >= NCA_WORKING_STOICHIOMETRY_RANGE[0]) & (profile.x_p_avg <= NCA_WORKING_STOICHIOMETRY_RANGE[1]))
    ):
        raise ValueError("Average electrode state left frozen OCP support; no clipping.")

    ueq_average_v = np.asarray(U_p(profile.x_p_avg)) - np.asarray(U_n(profile.x_n_avg))
    coefficient_v_per_k = np.asarray(dUeq_dT_panasonic(profile.soc))
    qdot_irr_measured = profile.current_a * (ueq_average_v - profile.measured_voltage_v)
    qdot_irr_model = profile.current_a * (ueq_average_v - profile.model_voltage_v)
    qdot_rev = -profile.current_a * profile.temperature_k * coefficient_v_per_k
    qdot_total_measured = qdot_irr_measured + qdot_rev
    qdot_total_model = qdot_irr_model + qdot_rev
    qdot_irr_difference = qdot_irr_measured - qdot_irr_model

    components = {
        "Qdot_irr_measured": qdot_irr_measured,
        "Qdot_irr_model": qdot_irr_model,
        "Qdot_rev": qdot_rev,
        "Qdot_total_measured": qdot_total_measured,
        "Qdot_total_model": qdot_total_model,
    }
    if not all(np.all(np.isfinite(values)) for values in components.values()):
        raise ValueError("A heat component is non-finite.")
    metrics = {
        **{name: time_weighted_metrics(values, profile.time_s) for name, values in components.items()},
        "Qdot_irr_measured_minus_model": time_weighted_metrics(
            qdot_irr_difference, profile.time_s
        ),
    }

    roots = zero_crossings()
    expected_roots = np.asarray([0.1456953642384106, 0.25857142857142856])
    np.testing.assert_allclose(roots, expected_roots, rtol=0.0, atol=5.0e-15)
    epsilon = 1.0e-7
    for root in roots:
        coefficient_above = dUeq_dT_panasonic(root + epsilon)
        coefficient_below = dUeq_dT_panasonic(root - epsilon)
        if coefficient_above * coefficient_below >= 0.0:
            raise RuntimeError("Literature coefficient did not change sign at an expected root.")
    if not (
        dUeq_dT_panasonic(0.20) < 0.0
        and (-2.9 * 298.15 * dUeq_dT_panasonic(0.20)) > 0.0
        and dUeq_dT_panasonic(0.50) > 0.0
        and (-2.9 * 298.15 * dUeq_dT_panasonic(0.50)) < 0.0
    ):
        raise RuntimeError("Discharge reversible-heat sign convention check failed.")

    duration_s = float(profile.time_s[-1] - profile.time_s[0])
    e_rev = metrics["Qdot_rev"]["integrated_energy_j"]
    e_total_measured = metrics["Qdot_total_measured"]["integrated_energy_j"]
    e_total_model = metrics["Qdot_total_model"]["integrated_energy_j"]
    signed_fraction_measured = e_rev / e_total_measured
    signed_fraction_model = e_rev / e_total_model
    absolute_rev_energy = float(np.trapezoid(np.abs(qdot_rev), profile.time_s))
    absolute_total_measured = float(np.trapezoid(np.abs(qdot_total_measured), profile.time_s))
    absolute_total_model = float(np.trapezoid(np.abs(qdot_total_model), profile.time_s))

    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    write_coefficient_csv(COEFFICIENT_CSV_PATH)
    write_metrics_csv(METRICS_CSV_PATH, metrics)
    make_component_plot(
        TIME_PLOT_PATH,
        (profile.time_s - profile.time_s[0]) / 60.0,
        "Elapsed time within supported interval (min)",
        components,
        invert_x=False,
        roots=roots,
    )
    make_component_plot(
        SOC_PLOT_PATH,
        profile.soc,
        "Average cell SOC",
        components,
        invert_x=True,
        roots=roots,
    )

    report = f"""# Panasonic 25 °C 1C heat-generation audit

## Scope

This is a forward-only heat-component audit. No value was optimized; `Cth` and `Gth` were not fitted. The frozen Panasonic electrical parameters, OCP functions, nominal effective stoichiometric windows, capacities, and voltage equations were used unchanged.

- 1C source: `{RAW_SOURCE_PATH}`
- Source columns: `{', '.join(columns)}`
- Raw-current sign: negative on discharge; converted once at the model-loader boundary to the internal convention `I > 0` for discharge.
- Frozen cell capacity: `{Q_CELL_AH:.12f} Ah`.
- Full 1C forward voltage RMSE reproduced before the heat calculation: `{1000.0 * full_rmse:.9f} mV`.

## Verified literature coefficient

The full-cell entropy coefficient is from He et al. (2024) Supplementary Information, Table S2, which explicitly cites Lu et al. (2019). The specified cell is Panasonic NCR18650PF, Li(NiCoAl)O2/carbon, 2.9 Ah, 2.5–4.2 V. No literature search or coefficient fitting was performed here.

`models/panasonic_entropic.py` stores the supplied table in mV/K and converts it to V/K. `dUeq_dT_panasonic(SOC)` uses linear interpolation only on the inclusive literature domain `0.10 <= SOC <= 0.90`; it raises outside that domain and never clips or extrapolates.

## Data and state preparation

- The complete active discharge was propagated once with actual nonuniform timestamps, measured discharge current, measured battery surface temperature, and frozen dynamics.
- Cell SOC was aligned exactly as in the closed 1C validation using relative raw Ah and frozen `Q_cell_Ah`.
- The audit interval was restricted to exact boundary-interpolated SOC endpoints `{profile.soc[0]:.12f}` and `{profile.soc[-1]:.12f}`, preserving all {profile.soc.size - 2} original interior timestamps.
- Evaluation points: `{profile.soc.size}`; duration: `{duration_s:.9f} s` (`{duration_s / 60.0:.6f} min`).
- Measured battery temperature: `{profile.temperature_k.min() - 273.15:.6f}` to `{profile.temperature_k.max() - 273.15:.6f} °C` (`{profile.temperature_k.min():.6f}` to `{profile.temperature_k.max():.6f} K`).
- Discharge-positive current: `{profile.current_a.min():.6f}` to `{profile.current_a.max():.6f} A`.
- Average negative stoichiometry: `{profile.x_n_avg.min():.12f}` to `{profile.x_n_avg.max():.12f}`.
- Average positive stoichiometry: `{profile.x_p_avg.min():.12f}` to `{profile.x_p_avg.max():.12f}`.
- Average-state OCP domain violations: **0**; clipping: **none**.

For all heat calculations in this audit,

```text
Ueq_average = U_p(x_p_avg) - U_n(x_n_avg)
```

was evaluated from the frozen OCPs. Surface stoichiometries were not used for `Ueq_average`. Its range over the supported interval is `{ueq_average_v.min():.9f}` to `{ueq_average_v.max():.9f} V`.

## Heat definitions

```text
Qdot_irr_measured = I * (Ueq_average - V_measured)
Qdot_irr_model    = I * (Ueq_average - V_model)
Qdot_rev          = -I * T_K * dUeq_dT(SOC)
Qdot_total_measured = Qdot_irr_measured + Qdot_rev
Qdot_total_model    = Qdot_irr_model + Qdot_rev
```

The measured-voltage and model-voltage definitions remain separate throughout the CSV, plots, metrics, and energy accounting.

## Heat-component metrics

Means and RMS values are time-weighted using the actual timestamps. Energies are trapezoidal time integrals over exact SOC 0.90–0.10 support.

| Component | Minimum (W) | Maximum (W) | Mean (W) | RMS (W) | Integrated energy (J) |
|---|---:|---:|---:|---:|---:|
{metric_rows(metrics)}

### Integrated heat and reversible fractions

- Measured-voltage irreversible heat: `{metrics['Qdot_irr_measured']['integrated_energy_j']:.6f} J`.
- Model-voltage irreversible heat: `{metrics['Qdot_irr_model']['integrated_energy_j']:.6f} J`.
- Net reversible heat: `{e_rev:.6f} J`.
- Measured-voltage total heat: `{e_total_measured:.6f} J`.
- Model-voltage total heat: `{e_total_model:.6f} J`.
- Signed reversible/total energy fraction, measured-voltage total: `{signed_fraction_measured:.9f}` (`{100.0 * signed_fraction_measured:.6f}%`).
- Signed reversible/total energy fraction, model-voltage total: `{signed_fraction_model:.9f}` (`{100.0 * signed_fraction_model:.6f}%`).
- Absolute reversible/total integral ratio, measured-voltage total: `{absolute_rev_energy / absolute_total_measured:.9f}` (`{100.0 * absolute_rev_energy / absolute_total_measured:.6f}%`).
- Absolute reversible/total integral ratio, model-voltage total: `{absolute_rev_energy / absolute_total_model:.9f}` (`{100.0 * absolute_rev_energy / absolute_total_model:.6f}%`).

The signed fractions retain the endothermic sign of net reversible energy; the absolute ratios describe reversible heat magnitude without cancellation.

## Physical sign and sign-change audit

For discharge, all evaluated current and absolute temperature values are positive. Numerically:

- at SOC 0.50, `dUeq/dT = {1.0e3 * dUeq_dT_panasonic(0.50):.6f} mV/K` and `Qdot_rev < 0`: endothermic;
- at SOC 0.20, `dUeq/dT = {1.0e3 * dUeq_dT_panasonic(0.20):.6f} mV/K` and `Qdot_rev > 0`: exothermic.

Under the specified piecewise-linear interpolation, the exact zero crossings are:

1. `SOC = {roots[1]:.12f}` between the 0.30 and 0.20 table points;
2. `SOC = {roots[0]:.12f}` between the 0.20 and 0.10 table points.

In discharge order from SOC 0.90 toward 0.10:

- endothermic reversible contribution: `0.90 >= SOC > {roots[1]:.12f}`;
- exothermic reversible contribution: `{roots[1]:.12f} > SOC > {roots[0]:.12f}`;
- endothermic reversible contribution: `{roots[0]:.12f} > SOC >= 0.10`.

`Qdot_rev` is zero at the two crossings. Direct evaluations immediately on both sides of each root verify the sign reversals; no sign or monotonicity constraint was imposed.

## Measured-versus-model irreversible heat

The explicitly reported difference is

```text
Qdot_irr_measured - Qdot_irr_model = I * (V_model - V_measured).
```

It ranges from `{metrics['Qdot_irr_measured_minus_model']['min_w']:.9f}` to `{metrics['Qdot_irr_measured_minus_model']['max_w']:.9f} W`, has time-weighted mean `{metrics['Qdot_irr_measured_minus_model']['mean_w']:.9f} W`, RMS `{metrics['Qdot_irr_measured_minus_model']['rms_w']:.9f} W`, and integrates to `{metrics['Qdot_irr_measured_minus_model']['integrated_energy_j']:.6f} J`. This difference is electrical-model voltage error expressed as heat-rate discrepancy; it was not pushed into any thermal parameter.

## Stop condition

The verified entropy coefficient and forward-only 1C heat audit are complete. No `Cth`, `Gth`, electrical, OCP, Arrhenius, UDDS, or CC–CV thermal calibration/validation was performed.
"""
    REPORT_PATH.write_text(report, encoding="utf-8")

    print(f"Supported interval: SOC {profile.soc[-1]:.2f} to {profile.soc[0]:.2f}")
    print(f"Evaluation points: {profile.soc.size}; duration: {duration_s:.6f} s")
    print(f"Sign-change SOC values: {roots.tolist()}")
    for name, values in metrics.items():
        print(
            f"{name}: min={values['min_w']:.9f}, max={values['max_w']:.9f}, "
            f"mean={values['mean_w']:.9f}, rms={values['rms_w']:.9f}, "
            f"energy={values['integrated_energy_j']:.6f} J"
        )
    print(f"Wrote: {COEFFICIENT_CSV_PATH}")
    print(f"Wrote: {METRICS_CSV_PATH}")
    print(f"Wrote: {REPORT_PATH}")
    print(f"Wrote: {TIME_PLOT_PATH}")
    print(f"Wrote: {SOC_PLOT_PATH}")


if __name__ == "__main__":
    main()
