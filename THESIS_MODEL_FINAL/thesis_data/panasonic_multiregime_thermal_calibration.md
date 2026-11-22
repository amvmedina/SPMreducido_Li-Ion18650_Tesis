# Panasonic common-parameter multiregime thermal calibration

## Scope and result

This was the single requested experiment. Only `Cth` and `Gth` were fitted; no model equation, OCP, entropy coefficient, electrical parameter, ambient offset, heat sign, or frozen parameter file was changed. UDDS remained fully withheld from calibration.

**Classification: MULTI-REGIME THERMAL CALIBRATION DOES NOT IMPROVE TRANSFER**

## Objective and preprocessing

The minimized objective was exactly

`J = (MSE_1C + MSE_HPPC + MSE_CC-CV) / 3 = 0.175331903462 K^2`.

- 1C: accepted SOC 0.90 to 0.10 interval; 299 samples.
- HPPC: 11 retained blocks and 81637 retained samples; no propagation across 13 omitted gaps; every block initialized at its measured starting temperature. The HPPC experiment contributes one retained-sample MSE, matching the preceding joint calibration convention.
- CC-CV: accepted SOC 0.10 to 0.90 interval; duplicate removal and existing CC/CV phase masks retained; 53 samples.
- Heat for all three calibration experiments: measured-voltage irreversible heat plus measured-temperature reversible heat.
- Bounds were numerical guards: `5 <= Cth <= 1000 J/K`, `1e-4 <= Gth <= 1 W/K`. Optimization was in log space with 64 deterministic starts.

## Best common parameters

| Quantity | Result |
|---|---:|
| Cth [J/K] | 76.5012809971 |
| Gth [W/K] | 0.123084630683 |
| tau = Cth/Gth [s] | 621.533984969 |
| Objective J [K²] | 0.175331903462 |

Relative to the accepted 1C-only pair, Cth changed by -34.905% and Gth by +12.806%.

## Calibration performance

| Experiment/phase | RMSE [K] | MAE [K] | Bias (predicted-measured) [K] |
|---|---:|---:|---:|
| 1C | 0.234333114 | 0.184241646 | -0.002883663 |
| HPPC | 0.639331724 | 0.582739544 | -0.580499194 |
| CC-CV total | 0.249677090 | 0.207237345 | -0.009403054 |
| CC phase | 0.226406250 | 0.185765877 | -0.058530743 |
| CV phase | 0.386164533 | 0.375430513 | 0.375430513 |

## Bounds, multistart stability, and practical conditioning

- Log-range distance to nearest bound: Cth 0.485144; Gth 0.227449. Neither parameter is within 1% of its log-bound range: True.
- Successful starts: 64 / 64.
- All successful Cth range: [76.4890785039, 76.5142070559] J/K; Gth range: [0.123081280944, 0.12308779925] W/K.
- Near-optimal (within 1.0e-06 relative objective) starts: 64; Cth range [76.4890785039, 76.5142070559], Gth range [0.123081280944, 0.12308779925].
- Numerical rank: 2.
- Singular values of the equal-experiment residual Jacobian in log parameters: 2.62033424645, 0.386765839955.
- Condition number: 6.77498883.
- Sensitivity-column correlation: 0.696269741.
- Covariance-derived log-parameter correlation: -0.697076336.

These are practical local-conditioning diagnostics, not a structural-identifiability claim. The accepted 1C-only analysis had condition number 6.120213126, sensitivity correlation 0.304113503, and covariance-derived parameter correlation -0.866785425.

## Strict withheld UDDS validation

Validation A used the same audited SOC 0.90 to 0.10 interval, chamber temperature, initialization, measured-voltage irreversible heat, and measured-temperature reversible heat, with no refit.

| Metric | Accepted 1C-only | New common pair |
|---|---:|---:|
| RMSE [K] | 0.526697652 | 0.572197773 |
| MAE [K] | 0.512960900 | 0.563603486 |
| Bias [K] | -0.512714512 | -0.563445905 |
| Maximum absolute error [K] | 0.875239690 | 0.848014389 |
| R² | -4.119410369 | -5.042122214 |
| Final measured temperature [°C] | 26.469807973 | 26.469807973 |
| Final predicted temperature [°C] | 26.112601243 | 26.018712689 |
| Measured temperature rise [K] | 0.636968211 | 0.636968211 |
| Predicted temperature rise [K] | 0.279761481 | 0.185872927 |

Relative to the accepted pair, the withheld UDDS RMSE worsened by 8.639%, and the absolute bias worsened by 9.895%.

## Secondary fully predictive UDDS diagnostic

| Metric | Result |
|---|---:|
| RMSE [K] | 0.715024880 |
| Bias [K] | -0.696713341 |
| Maximum absolute error [K] | 1.381721630 |

This path used frozen model voltage and recursively predicted temperature in reversible heat. It did not influence fitting. Its RMSE is 1.250 times the primary Validation-A RMSE.

## Compact comparison

The old pairs were not refitted. Their 1C, HPPC, and UDDS entries reproduce the accepted results. The previously unreported old-pair CC-CV entries are forward-only evaluations using the already-audited CC-CV Validation-A preprocessing.

| Case | Cth [J/K] | Gth [W/K] | tau [s] | 1C RMSE [K] | HPPC RMSE [K] | CC-CV A RMSE [K] | UDDS A RMSE [K] | UDDS bias [K] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Accepted 1C-only pair | 117.523178001 | 0.109112090 | 1077.086675 | 0.143622065 | 0.617273179 | 0.427780949 | 0.526697652 | -0.512714512 |
| Previous 1C+HPPC joint optimum | 200.000000000 | 0.078451555 | 2549.343991 | 0.254140709 | 0.512740397 | 0.650974767 | 0.378519120 | -0.353251679 |
| New 1C+HPPC+CC-CV common pair | 76.501280997 | 0.123084631 | 621.533985 | 0.234333114 | 0.639331724 | 0.249677090 | 0.572197773 | -0.563445905 |

## Required decisions

1. **Does multiregime calibration materially change Cth and/or Gth?** Yes: changes are Cth -34.905% and Gth +12.806%.
2. **Can one common pair represent 1C, HPPC, and CC-CV reasonably?** No: 1C and CC-CV remain below 0.30 K RMSE, but HPPC has 0.639331724 K RMSE and -0.580499194 K bias, both worse than required for accepting a common representation. The individual results remain visible and are not hidden by the mean objective.
3. **Is the new pair numerically well-defined?** Yes: rank, condition, bounds, and multistart convergence are reported above. This is a practical statement only.
4. **Does withheld UDDS Validation A improve relative to 1C-only?** No; RMSE changes from 0.526697652 to 0.572197773 K.
5. **Does systematic negative UDDS bias decrease?** No; bias changes from -0.512714512 to -0.563445905 K.
6. **Is the fully predictive path still primarily limited by electrical heat-generation error?** No clear evidence that it is the primary limitation: its RMSE is 1.250 times Validation A, so electrical heat adds error, but the Validation-A thermal-transfer discrepancy remains dominant. No repair was attempted.
7. **Is there sufficient evidence to replace the frozen 1C-only pair?** No.

The currently frozen thermal parameter file was not overwritten.
