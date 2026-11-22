# Panasonic HPPC primary three-parameter calibration

## Scope and frozen quantities

Only `R0`, `alpha_p`, and `d_n` were optimized. The estimates are conditional on the provisional fixed values `alpha_n=3000 s` and `d_p=0.0007 1/s`. The frozen Panasonic OCPs, nominal effective stoichiometric windows, `Q_cell`, `Q_n`, and `Q_p` were not changed.

The calibration uses 64 complete pulses in 14 independently initialized contiguous histories. Pulse IDs 60, 64, and 67 are excluded. Cumulative exported Ah initializes each history; current is not integrated through omitted gaps. Actual mixed time steps and duplicate timestamps are retained without resampling.

## Objective

For each pulse, one measured and one modeled pre-pulse voltage are subtracted from that pulse's response. No sample-wise shifting is used. The implemented least-squares vector weights every sample in a pulse phase by `1/sqrt(64*N_phase_samples_for_that_pulse)`. Therefore its squared norm is exactly:

```text
J = mean_pulses(MSE_immediate) + mean_pulses(MSE_within_pulse) + mean_pulses(MSE_post_relaxation)
```

The three response phases have equal fixed weight, and long rest records do not give a pulse more weight. These weights were specified before fitting and were not tuned for fit quality.

## Parameterization and bounds

Optimization is in natural-log parameter coordinates. The project contains no documented target-specific admissible bounds, so the following broad positivity/numerical-safety guards were used; they are not asserted physical Panasonic ranges:

| Parameter | Lower | Upper |
|---|---:|---:|
| `R0` [ohm] | 1e-06 | 0.2 |
| `alpha_p` [s] | 1 | 1e+05 |
| `d_n` [1/s] | 1e-08 | 0.1 |

Frozen OCP-domain violations are rejected, never clipped.
For these 14 histories, the frozen positive-OCP domain imposes a derived feasibility edge `alpha_p <= 11088.815475 s`. This is not an invented physical bound; it is the largest value whose evaluated `x_p,surf` remains at or below 0.96 on the retained data.

## Best fit

- `R0 = 0.0183337802186 ohm`
- `alpha_p = 11088.8154739 s`
- `d_n = 7.58662286968e-05 1/s`
- Phase-balanced objective `J = 0.00923579593227 V^2`
- Phase-balanced dynamic RMSE = 55.485121 mV
- Initialization-fixture phase-balanced RMSE = 67.277603 mV; the optimized value is 17.528% lower.
- Raw-sample dynamic RMSE = 23.429229 mV (diagnostic only; rest-sample dominated)
- Maximum fitted `x_p,surf = 0.959999999987`; frozen upper support is 0.960000000000.
- Active frozen-OCP feasibility boundary: **YES**.

### Dynamic RMSE by phase

| Phase | RMSE [mV] |
|---|---:|
| Immediate response | 50.117074 |
| Within pulse | 80.508042 |
| Post relaxation | 15.573373 |

### Phase-balanced dynamic RMSE by C-rate

| C-rate | RMSE [mV] |
|---|---:|
| 0.5C | 25.121243 |
| 1C | 54.446893 |
| 2C | 65.507099 |
| 4C | 57.342596 |
| 6C | 68.390138 |

### Phase-balanced dynamic RMSE by SOC region

SOC bins are low `[0,1/3)`, middle `[1/3,2/3)`, and high `[2/3,1]` by pulse-start cell SOC.

| SOC region | RMSE [mV] |
|---|---:|
| high | 38.771898 |
| middle | 37.786283 |
| low | 75.837480 |

## Absolute-voltage diagnostic

- RMSE: 36.997186 mV
- MAE: 22.907830 mV
- R2: 0.983905477

Absolute voltage was not the optimization target, so static OCP/baseline mismatch cannot dominate the dynamic fit.

## Step, relaxation, and apparent-resistance diagnostics

- Pulse-step error RMSE / MAE / bias: 50.117074 / 43.668889 / +43.668889 mV.
- Relaxation-shape error (each response referenced to its first post-pulse rest sample) RMSE / MAE: 145.431000 / 108.833505 mV.
- Measured `R_app` range: 20.071433-33.347788 mOhm.
- Model-implied initial-step range: 28.086614-39.022182 mOhm.
- Measured/model implied `R_app` RMSE / MAE: 9.350505 / 8.465919 mOhm.

`R0` is not equated to `R_app`; the model-implied voltage step also contains algebraic kinetics and the current-dependent surface-state/OCP contribution.

## Multistart recoverability

- Deterministic starts: 12; successful feasible solutions: 12.
- Near-best definition: objective within 1.0% of the best; near-best solutions: 12.
- Near-best `R0` range: 0.0183337787203-0.018333792713 ohm.
- Near-best `alpha_p` range: 11088.8154739-11088.8154739 s.
- Near-best `d_n` range: 7.58661971579e-05-7.58664865096e-05 1/s.
- Successful solutions at a numerical bound: 0.
- Successful solutions at the frozen-OCP-derived `alpha_p` feasibility boundary: 12.

The multistart CSV records every start, termination result, objective, and fitted parameter set. Similar objectives with materially different near-best parameters indicate weak practical recoverability even when sensitivities are nonzero. Convergence of `alpha_p` to the frozen-domain feasibility edge is a boundary-limited estimate, not evidence of an interior optimum.

## Post-fit local audit

- Singular values of the phase-balanced normalized Jacobian: 0.272556193, 0.0297469266, 0.018097652.
- Condition number: 15.0603069.
- Weakest normalized direction: `R0`=+0.187868, `alpha_p`=+0.020576, `d_n`=+0.981979.

Estimator-oriented parameter correlations (normalized pseudoinverse of `J.T @ J`):

| | R0 | alpha_p | d_n |
|---|---:|---:|---:|
| R0 | +1.000000 | -0.766720 | +0.539627 |
| alpha_p | -0.766720 | +1.000000 | +0.076011 |
| d_n | +0.539627 | +0.076011 | +1.000000 |

Central normalized log-parameter differences are used for `R0` and `d_n`. Because `alpha_p` is on the frozen-OCP feasibility edge, its derivative uses a second-order inward (backward) difference; no infeasible OCP evaluation is substituted into the Jacobian.

This remains a local fitted-point audit. It tests conditioning under the actual HPPC objective but does not by itself prove global or noise-robust recoverability.

## Calibration assessment

The numerical multistart is stable for `R0` and `d_n` under the stated fixed assumptions, but the three-parameter representation is **boundary-limited and leaves substantial systematic mismatch**:

- all 12 successful starts place `alpha_p` on the frozen positive-OCP support boundary;
- immediate voltage steps retain a one-sided bias of +43.669 mV;
- within-pulse RMSE remains 80.508 mV;
- low-SOC phase-balanced RMSE reaches 75.837 mV; and
- the independently referenced relaxation-shape RMSE is 145.431 mV.

Therefore this requested first calibration is computationally complete, but the fitted triplet should be retained as a **conditional diagnostic result**, not frozen as a validated Panasonic dynamic parameter set. The local post-fit condition number does not cure the active boundary or model-data shape mismatch. Per the stop condition, no additional parameter was released and no alternative model or dataset was run.

## Stop condition

The first three-parameter HPPC calibration is complete. `alpha_n` and `d_p` were not released; no OCP, capacity, stoichiometric window, model equation, thermal parameter, UDDS profile, or 1C validation was changed or run.
