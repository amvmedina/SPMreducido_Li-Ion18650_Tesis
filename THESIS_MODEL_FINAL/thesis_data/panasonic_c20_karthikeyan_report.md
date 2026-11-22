# Panasonic C/20 Karthikeyan OCP compatibility gate

## Decision: **FAIL**

This gate adjusts only `x_n_start`, `x_n_end`, `x_p_start`, and `x_p_end`. It contains no resistance, kinetic, diffusion, overpotential, or thermal terms, and performs no clipping.

## Data preparation

- Source: `C:\Users\amvel\Documents\Andes\Tesis\datasets\converted_csv\c20_ocv_25C_converted.csv`
- Discharge source rows: 6 through 1246, inclusive
- Samples: 1241
- Raw current sign: negative throughout (-0.14536 to -0.14454 A)
- Voltage span: 4.17030 to 2.49948 V
- Integrated Panasonic C/20 capacity: 2.994979138519 Ah
- Fixed Panasonic capacity reference: 2.994979138519 Ah
- Normalization: `z = Q_discharge / Q_discharge_final`, with endpoints exactly 0 and 1

## Best feasible solution

| Metric | Value |
| --- | ---: |
| RMSE | 53.882 mV |
| MAE | 40.303 mV |
| R² | 0.96159673 |
| `x_n_start` | 0.238047229 |
| `x_n_end` | 0.028850400 |
| `x_p_start` | 0.401280996 |
| `x_p_end` | 0.883551177 |

Supported domains are `0.01 <= x_n <= 0.80` and `0.18 <= x_p <= 0.96`; the fitted linear trajectories remain inside these domains by construction and validation.

## Residual-shape evidence

- Beginning residual (measured - fitted): 128.640 mV
- Residual near `z=0.50`: -63.660 mV
- Residual near `z=0.95`: 92.622 mV
- Terminal residual at `z=1`: -496.771 mV
- Largest absolute point residual: 496.771 mV at `z=1.000000`

The alternating broad residual lobes and the terminal mismatch show a systematic voltage-shape incompatibility, not unstructured noise.

## Gate checks

| Check | Criterion | Result |
| --- | --- | --- |
| Voltage RMSE | <= 50 mV | FAIL |
| Voltage MAE | <= 40 mV | FAIL |
| R² | >= 0.990 | FAIL |
| Largest absolute 5%-bin mean residual | <= 75 mV; observed 93.117 mV | FAIL |
| Span of 5%-bin mean residuals | <= 125 mV; observed 155.453 mV | FAIL |
| Literature domains and physical direction | Required | PASS |
| Non-pathological electrode-window widths | >= 0.05 each | PASS |

## Multistart and identifiability

- Deterministic starts attempted: 40 (seed `20260912` for generated starts)
- Successful feasible solutions: 40
- Good solutions within 1.0 mV RMSE of the best: 2
- Maximum good-solution endpoint span: 0.561952771
- Materially different windows (endpoint span >= 0.05) with nearly identical error: **YES**

Good-solution endpoint spans:

| Endpoint | Span across good solutions |
| --- | ---: |
| `x_n_start` | 0.561952771 |
| `x_n_end` | 0.014473787 |
| `x_p_start` | 0.079990570 |
| `x_p_end` | 0.003976919 |

All good multistart solutions:

| Start | RMSE [mV] | MAE [mV] | R² | x_n_start | x_n_end | x_p_start | x_p_end |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 28 | 53.882 | 40.303 | 0.96159673 | 0.238047229 | 0.028850400 | 0.401280996 | 0.883551177 |
| 4 | 54.725 | 46.493 | 0.96038529 | 0.800000000 | 0.014376613 | 0.481271566 | 0.887528096 |

The complete record of every start and every successful solution is stored in `panasonic_c20_karthikeyan_multistart.csv`.

## Interpretation

The Karthikeyan NCA/MCMB pair fails this quasi-equilibrium compatibility gate. Failed checks: rmse, mae, r_squared, shape_absolute, shape_span.
