# Panasonic OCP-identification closure

## Status: CLOSED

- **Voltage/OCP representation accepted.**
- **Electrode stoichiometric windows weakly identifiable.**
- **Nominal solution selected for subsequent dynamic calibration.**
- No dynamic calibration or new optimization was performed during closure.

## Frozen source

- Multistart source: `C:\Users\amvel\Documents\Andes\THESIS_MODEL_FINAL\thesis_data\panasonic_c20_qi_up_multistart.csv`
- Selected successful solution: start `3`
- Negative OCP: fixed Guo--White (2005) MCMB discharge spline
- Positive OCP: 13-knot monotonic Qi-style representation from start 3
- Panasonic capacity: `2.994979138519 Ah`

## Nominal effective stoichiometric windows

These are nominal effective stoichiometric windows, not uniquely identified physical electrode stoichiometries.

| Endpoint | Frozen value |
| --- | ---: |
| `x_n_start` | 0.720451136448 |
| `x_n_end` | 0.251849704066 |
| `x_p_start` | 0.224027933455 |
| `x_p_end` | 0.953477823621 |

## Frozen positive-OCP knots

| Knot | `x_p` | `U_p` [V] |
| ---: | ---: | ---: |
| 0 | 0.180000 | 4.327122671400 |
| 1 | 0.245000 | 4.233957011590 |
| 2 | 0.310000 | 4.161069377480 |
| 3 | 0.375000 | 4.067239488640 |
| 4 | 0.440000 | 3.999087424060 |
| 5 | 0.505000 | 3.928197622260 |
| 6 | 0.570000 | 3.831414416460 |
| 7 | 0.635000 | 3.771043136820 |
| 8 | 0.700000 | 3.724570043030 |
| 9 | 0.765000 | 3.672437228370 |
| 10 | 0.830000 | 3.585344510200 |
| 11 | 0.895000 | 3.466214508910 |
| 12 | 0.960000 | 3.346854419760 |

## Forward-only C/20 verification

- Objective region: `0.05 <= z <= 0.95` (1116 samples)
- Reproduced RMSE: `1.085749161 mV`
- Recorded start-3 RMSE: `1.085749161 mV`
- Absolute RMSE difference: `4.67464608689e-12 mV`
- Reproduced MAE: `0.828250698 mV`
- Reproduced R²: `0.999976353112`

The frozen forward evaluation reproduces the recorded start-3 fit within numerical precision. The Panasonic OCP-identification stage is closed; HPPC and dynamic calibration were not run.
