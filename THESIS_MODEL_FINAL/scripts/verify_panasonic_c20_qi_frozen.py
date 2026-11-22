"""Forward-only closure verification of the frozen Panasonic OCP parameters."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from models.ocp_panasonic_qi_frozen import (  # noqa: E402
    Q_CELL_AH,
    SOURCE_MULTISTART_START_ID,
    UP_KNOT_STOICHIOMETRIES,
    UP_KNOT_VOLTAGES_V,
    V_eq,
    X_N_END,
    X_N_START,
    X_P_END,
    X_P_START,
)


DATA_PATH = Path(
    r"C:\Users\amvel\Documents\Andes\Tesis\datasets\converted_csv"
    r"\c20_ocv_25C_converted.csv"
)
MULTISTART_PATH = REPOSITORY_ROOT / "thesis_data" / "panasonic_c20_qi_up_multistart.csv"
REPORT_PATH = REPOSITORY_ROOT / "thesis_data" / "panasonic_ocp_identification_closure.md"
DISCHARGE_START_INDEX = 6
DISCHARGE_END_INDEX_INCLUSIVE = 1246
FIT_Z_MIN = 0.05
FIT_Z_MAX = 0.95


def _source_record() -> dict[str, str]:
    with MULTISTART_PATH.open(newline="", encoding="utf-8") as handle:
        matches = [
            row
            for row in csv.DictReader(handle)
            if int(row["start_id"]) == SOURCE_MULTISTART_START_ID
        ]
    if len(matches) != 1:
        raise ValueError("Expected exactly one recorded multistart solution 3.")
    record = matches[0]
    if record["successful_solution"].lower() != "true" or record["feasible"].lower() != "true":
        raise ValueError("Recorded multistart solution 3 is not successful and feasible.")
    return record


def _verify_frozen_values(record: dict[str, str]) -> None:
    frozen_endpoints = np.asarray([X_N_START, X_N_END, X_P_START, X_P_END])
    recorded_endpoints = np.asarray(
        [float(record[name]) for name in ("x_n_start", "x_n_end", "x_p_start", "x_p_end")]
    )
    recorded_knots = np.asarray(
        [float(record[f"u_p_knot_{index:02d}_v"]) for index in range(13)]
    )
    np.testing.assert_allclose(frozen_endpoints, recorded_endpoints, rtol=0.0, atol=5.0e-13)
    np.testing.assert_allclose(UP_KNOT_VOLTAGES_V, recorded_knots, rtol=0.0, atol=5.0e-13)


def _load_c20() -> tuple[np.ndarray, np.ndarray]:
    table = np.genfromtxt(DATA_PATH, delimiter=",", names=True, dtype=float, encoding="utf-8")
    selection = slice(DISCHARGE_START_INDEX, DISCHARGE_END_INDEX_INCLUSIVE + 1)
    time_s = np.asarray(table["Time"][selection], dtype=float)
    current_a = np.asarray(table["Current"][selection], dtype=float)
    voltage_v = np.asarray(table["Voltage"][selection], dtype=float)
    increments_ah = (
        0.5
        * (np.abs(current_a[:-1]) + np.abs(current_a[1:]))
        * np.diff(time_s)
        / 3600.0
    )
    capacity_ah = np.concatenate(([0.0], np.cumsum(increments_ah)))
    if not np.isclose(capacity_ah[-1], Q_CELL_AH, rtol=0.0, atol=1.0e-9):
        raise ValueError("Panasonic capacity differs from the frozen capacity.")
    return capacity_ah / capacity_ah[-1], voltage_v


def _metrics(measured_v: np.ndarray, fitted_v: np.ndarray) -> tuple[float, float, float]:
    residual_v = measured_v - fitted_v
    rmse_v = float(np.sqrt(np.mean(residual_v**2)))
    mae_v = float(np.mean(np.abs(residual_v)))
    ss_res = float(np.sum(residual_v**2))
    ss_total = float(np.sum((measured_v - np.mean(measured_v)) ** 2))
    return rmse_v, mae_v, 1.0 - ss_res / ss_total


def main() -> None:
    record = _source_record()
    _verify_frozen_values(record)
    z, measured_v = _load_c20()
    fit_mask = (z >= FIT_Z_MIN) & (z <= FIT_Z_MAX)
    fitted_v = np.asarray(V_eq(z), dtype=float)
    rmse_v, mae_v, r_squared = _metrics(measured_v[fit_mask], fitted_v[fit_mask])
    recorded_rmse_v = float(record["fit_rmse_v"])
    rmse_difference_v = abs(rmse_v - recorded_rmse_v)
    if rmse_difference_v > 5.0e-12:
        raise AssertionError(
            f"Frozen RMSE {rmse_v:.12g} does not reproduce recorded RMSE "
            f"{recorded_rmse_v:.12g}."
        )

    lines = [
        "# Panasonic OCP-identification closure",
        "",
        "## Status: CLOSED",
        "",
        "- **Voltage/OCP representation accepted.**",
        "- **Electrode stoichiometric windows weakly identifiable.**",
        "- **Nominal solution selected for subsequent dynamic calibration.**",
        "- No dynamic calibration or new optimization was performed during closure.",
        "",
        "## Frozen source",
        "",
        f"- Multistart source: `{MULTISTART_PATH}`",
        f"- Selected successful solution: start `{SOURCE_MULTISTART_START_ID}`",
        "- Negative OCP: fixed Guo--White (2005) MCMB discharge spline",
        "- Positive OCP: 13-knot monotonic Qi-style representation from start 3",
        f"- Panasonic capacity: `{Q_CELL_AH:.12f} Ah`",
        "",
        "## Nominal effective stoichiometric windows",
        "",
        "These are nominal effective stoichiometric windows, not uniquely "
        "identified physical electrode stoichiometries.",
        "",
        "| Endpoint | Frozen value |",
        "| --- | ---: |",
        f"| `x_n_start` | {X_N_START:.12f} |",
        f"| `x_n_end` | {X_N_END:.12f} |",
        f"| `x_p_start` | {X_P_START:.12f} |",
        f"| `x_p_end` | {X_P_END:.12f} |",
        "",
        "## Frozen positive-OCP knots",
        "",
        "| Knot | `x_p` | `U_p` [V] |",
        "| ---: | ---: | ---: |",
    ]
    for index, (x_p, voltage_v) in enumerate(
        zip(UP_KNOT_STOICHIOMETRIES, UP_KNOT_VOLTAGES_V, strict=True)
    ):
        lines.append(f"| {index} | {x_p:.6f} | {voltage_v:.12f} |")
    lines.extend(
        [
            "",
            "## Forward-only C/20 verification",
            "",
            f"- Objective region: `0.05 <= z <= 0.95` ({int(np.sum(fit_mask))} samples)",
            f"- Reproduced RMSE: `{1000.0 * rmse_v:.9f} mV`",
            f"- Recorded start-3 RMSE: `{1000.0 * recorded_rmse_v:.9f} mV`",
            f"- Absolute RMSE difference: `{1000.0 * rmse_difference_v:.12g} mV`",
            f"- Reproduced MAE: `{1000.0 * mae_v:.9f} mV`",
            f"- Reproduced R²: `{r_squared:.12f}`",
            "",
            "The frozen forward evaluation reproduces the recorded start-3 fit "
            "within numerical precision. The Panasonic OCP-identification stage "
            "is closed; HPPC and dynamic calibration were not run.",
        ]
    )
    REPORT_PATH.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Frozen start-3 RMSE verification: {1000.0 * rmse_v:.9f} mV")
    print(REPORT_PATH.resolve())


if __name__ == "__main__":
    main()
