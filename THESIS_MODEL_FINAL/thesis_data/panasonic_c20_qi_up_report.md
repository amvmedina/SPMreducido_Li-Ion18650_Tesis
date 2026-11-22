# Panasonic NCR18650PF Qi-style OCP identification

## Decision: **ACCEPTABLE OCP IDENTIFICATION**

The negative OCP is the fixed Guo--White MCMB discharge spline. The positive OCP is represented by 13 fitted voltages at fixed stoichiometric knots with monotonic linear interpolation. Only the four stoichiometric endpoints and those 13 voltages were fitted.

## Data and objective

- Panasonic capacity: `2.994979138519 Ah`
- Discharge rows: 6--1246, inclusive
- Full samples: 1241
- Primary objective samples: 1116, from `z=0.050011` to `z=0.949451`
- No clipping, resistance, kinetics, diffusion, overpotential, or thermal terms

## Fit-region metrics: 5--95%

| Metric | Value |
| --- | ---: |
| RMSE | 0.877 mV |
| MAE | 0.702 mV |
| R² | 0.999984571 |
| Mean-residual bin range | -0.998 to 0.728 mV |
| Mean-residual bin span | 1.725 mV |

## Full-discharge evaluation: 0--100%

| Metric | Value |
| --- | ---: |
| RMSE | 23.334 mV |
| MAE | 4.035 mV |
| R² | 0.992798036 |
| Residual at `z=0` | 47.583 mV |
| Residual at `z=1` | -408.057 mV |

## Fitted stoichiometric windows

| Endpoint | Value |
| --- | ---: |
| `x_n_start` | 0.800000000 |
| `x_n_end` | 0.011762206 |
| `x_p_start` | 0.223425408 |
| `x_p_end` | 0.952707790 |

## Identified positive-OCP knots

| Knot | `x_p` | `U_p` [V] | Drop to next [mV] |
| ---: | ---: | ---: | ---: |
| 0 | 0.180000 | 4.257696238 | 29.516 |
| 1 | 0.245000 | 4.228180602 | 72.941 |
| 2 | 0.310000 | 4.155239554 | 90.727 |
| 3 | 0.375000 | 4.064512619 | 63.824 |
| 4 | 0.440000 | 4.000688694 | 68.216 |
| 5 | 0.505000 | 3.932472973 | 97.003 |
| 6 | 0.570000 | 3.835470174 | 60.419 |
| 7 | 0.635000 | 3.775051077 | 41.735 |
| 8 | 0.700000 | 3.733315698 | 39.056 |
| 9 | 0.765000 | 3.694259398 | 59.790 |
| 10 | 0.830000 | 3.634469207 | 86.532 |
| 11 | 0.895000 | 3.547936730 | 50.250 |
| 12 | 0.960000 | 3.497686884 | -- |

## Constraint and pathology audit

| Check | Finding |
| --- | --- |
| Physical directions and domains | PASS |
| Monotonic positive OCP | PASS |
| Stoichiometric endpoint near a domain bound | YES |
| Positive-OCP knot near 2.5 or 4.4 V | NO |
| Long flat knot run | NO (longest adjacent run: 0) |
| Broad systematic fit-region residual lobes | NO |

## Multistart identifiability

- Deterministic starts: 40
- Successful feasible solutions: 25
- Good solutions within 0.5 mV RMSE of best: 4
- Largest endpoint span across good solutions: 0.243592156
- Largest fit-region voltage difference among good solutions: 6.363 mV
- Materially different windows with nearly identical voltage error: **YES**

| Start | RMSE [mV] | MAE [mV] | R² | `x_n_start` | `x_n_end` | `x_p_start` | `x_p_end` |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0.877 | 0.702 | 0.999984571 | 0.800000000 | 0.011762206 | 0.223425408 | 0.952707790 |
| 2 | 0.877 | 0.702 | 0.999984571 | 0.800000000 | 0.011762619 | 0.223425027 | 0.952707775 |
| 3 | 1.086 | 0.828 | 0.999976353 | 0.720451136 | 0.251849704 | 0.224027933 | 0.953477824 |
| 6 | 1.348 | 1.100 | 0.999963568 | 0.556407844 | 0.198863391 | 0.180000000 | 0.869999314 |

Every successful result, its metrics, four endpoints, and 13 knot voltages is recorded in `panasonic_c20_qi_up_multistart.csv`.

## Decision checks

| Requirement | Result |
| --- | --- |
| `fit_rmse` | PASS |
| `fit_mae` | PASS |
| `fit_r_squared` | PASS |
| `residual_shape` | PASS |
| `domains` | PASS |
| `monotonic` | PASS |
| `no_long_flat_run` | PASS |
| `voltage_bounds_not_active` | PASS |

## Interpretation

The Guo-negative/Qi-style-positive construction reproduces the Panasonic 5--95% C/20 voltage shape without broad systematic residual mismatch and passes the OCP compatibility checks.
**Voltage/OCP representation acceptable, stoichiometric parameters weakly identifiable.** Materially different near-optimal windows exist, so the stoichiometric solution is not unique.
