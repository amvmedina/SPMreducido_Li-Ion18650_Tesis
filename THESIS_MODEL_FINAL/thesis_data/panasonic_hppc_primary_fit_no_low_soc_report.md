# Panasonic HPPC diagnostic calibration without low-SOC pulses

## BOUNDARY WAS PRIMARILY LOW-SOC DRIVEN

This is the single requested diagnostic rerun. The sole calibration-data change is `pulse-start cell SOC >= 1/3`; all lower-SOC pulses are excluded. Frozen OCPs, capacities, nominal effective stoichiometric windows, model equations, phase weighting, broad numerical bounds, multistart definitions, and block initialization logic are unchanged.

- Retained pulses: **40** in 8 independently initialized contiguous histories.
- Fixed `alpha_n = 3000 s` and `d_p = 0.0007 1/s`.
- Actual mixed sampling is retained; no omitted gap is integrated and no resampling is used.

## Fitted conditional parameters

| Parameter | SOC >= 1/3 diagnostic | Original all-SOC fit |
|---|---:|---:|
| `R0` [ohm] | 0.0259594900586 | 0.0183337802186 |
| `alpha_p` [s] | 3204.22765511 | 11088.8154739 |
| `d_n` [1/s] | 0.000183029692927 | 7.58662286968e-05 |

The retained-data frozen-OCP feasibility edge is `alpha_p <= 23435.6278698 s`. `alpha_p` at that edge: **NO**.

## Dynamic fit metrics

| Metric | SOC >= 1/3 diagnostic [mV] | Original all-SOC [mV] |
|---|---:|---:|
| Phase-balanced RMSE | 36.835198 | 55.485121 |
| Immediate RMSE | 43.114901 | 50.117074 |
| Within-pulse RMSE | 44.261897 | 80.508042 |
| Post-relaxation RMSE | 15.889782 | 15.573373 |

Phase-balanced RMSE by retained SOC region:

| Region | SOC >= 1/3 diagnostic [mV] | Original all-SOC fit [mV] |
|---|---:|---:|
| High SOC | 37.760683 | 38.771898 |
| Middle SOC | 35.238749 | 37.786283 |

## Multistart recoverability

- Successful feasible starts: 12 of 12.
- Near-best solutions (within 1.0%): 12.
- Near-best `R0` spread: 0.0237820739005-0.0260982505475 ohm.
- Near-best `alpha_p` spread: 3204.22757834-6350.11220171 s.
- Near-best `d_n` spread: 0.000173561943181-0.000192513821843 1/s.
- Solutions at the retained-data frozen-OCP feasibility edge: 0.

All near-best solutions are interior with respect to the OCP support, but they occupy distinct `alpha_p` bands near 3204, 4847, and 6350 s while differing in objective by less than 1%. Thus the boundary result is removed, but `alpha_p` remains practically weak/non-unique under this diagnostic.

## Post-fit local audit

- Normalized-Jacobian singular values: 0.363731562, 0.0080305727, 0.00643585696.
- Condition number: **56.5164149**.
- Estimator-oriented correlations:

| | R0 | alpha_p | d_n |
|---|---:|---:|---:|
| R0 | +1.000000 | -0.552742 | +0.811047 |
| alpha_p | -0.552742 | +1.000000 | +0.018572 |
| d_n | +0.811047 | +0.018572 | +1.000000 |

The fitted-point Jacobian uses the same normalized finite differences as the original calibration. If the optimum is on the OCP feasibility edge, `alpha_p` uses the same second-order inward difference.

## Comparison and decision

Removing low-SOC pulses changes the phase-balanced RMSE from 55.485 to 36.835 mV. The original all-SOC `alpha_p` was boundary-limited at 11088.815474 s; the retained-data value is 3204.227655 s against its recalculated support edge of 23435.627870 s.

**Final classification: BOUNDARY WAS PRIMARILY LOW-SOC DRIVEN.**

No other calibration, validation, parameter release, OCP/support change, capacity/window change, or model-equation change was performed.
