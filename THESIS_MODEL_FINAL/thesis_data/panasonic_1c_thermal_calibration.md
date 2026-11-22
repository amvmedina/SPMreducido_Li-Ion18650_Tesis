# Panasonic 25 °C 1C primary thermal calibration

## Result: THERMAL PARAMETERS WELL CONDITIONED

The calibration fitted only `Cth` and `Gth`. They are retained as **effective lumped thermal parameters**, not uniquely measured material properties: the target is cell-surface temperature, the cell has one uniform model temperature, and all external heat-transfer detail is collapsed into `Gth`.

No electrical parameter, OCP, nominal effective stoichiometric window, capacity, entropy coefficient, or voltage equation was changed. No initial-temperature offset was fitted.

## Data and prescribed calibration heat

- Source: `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\C20 OCV and 1C discharge tests_start_of_tests\03-09-17_17.59 3349_Dis1C_1.mat`
- Exact entropy support used: `0.10 <= SOC <= 0.90`, with the same exact endpoint interpolation as the preceding heat audit.
- Points: 299; duration: 2974.928411972 s.
- Raw negative discharge current was converted once to model discharge-positive current: 2.899000–2.899820 A.
- Target surface temperature: measured `Battery_Temp_degC`, 26.863600–30.986661 °C.
- Measured time-varying ambient: `Chamber_Temp_degC`, 25.000000–26.000000 °C. It was not replaced by a fixed 25 °C value.
- Initial condition: measured and modeled temperature both 26.863600000 °C.
- Frozen electrical forward result reproduced: 114.362188001 mV full-discharge RMSE.

The primary prescribed heat was

```text
Ueq_average = U_p(x_p_avg) - U_n(x_n_avg)
Qdot_irr_measured = I * (Ueq_average - V_measured)
Qdot_rev = -I * T_measured_K * dUeq_dT(SOC)
Qdot_cal = Qdot_irr_measured + Qdot_rev
```

Measured temperature—not fitted model temperature—was used in `Qdot_rev`. `Qdot_irr_model` was not used in the fit.

## Optimization

The model was propagated exactly on every nonuniform interval under left-sample piecewise-constant `Qdot_cal` and measured `Tamb`. The implementation uses `expm1` and a finite small-argument ratio, retaining the adiabatic limiting form without Euler integration.

- Guard bounds: `5 <= Cth <= 200 J/K`; `1e-4 <= Gth <= 1 W/K`.
- Deterministic starts: 25 on a 5×5 grid distributed geometrically in parameter space (uniformly in log space).
- Successful starts: 25.
- Primary objective: mean squared unweighted temperature residual.
- Best objective: 0.020627297616 K².
- **Cth = 117.523178001 J/K**.
- **Gth = 0.109112089746 W/K**.
- **Thermal time constant = 1077.08667549 s (17.951444592 min)**.
- Frozen effective-parameter module: `models/panasonic_thermal_parameters.py`. It does not alter generic thermal defaults.
- Cth near a guard bound: **NO**; Gth near a guard bound: **NO**. “Near” means the closest 1% of the log-bound interval.

### Primary measured-voltage-heat fit

| Metric | Value |
|---|---:|
| RMSE | 0.143622065 K |
| MAE | 0.119569322 K |
| Mean signed error | -0.008530653 K |
| Maximum absolute error | 0.383292660 K |
| R² | 0.964702003 |

## Multistart stability

- Successful-solution Cth range: 117.523100528–117.523257798 J/K; spread 0.00015727 J/K.
- Successful-solution Gth range: 0.109112067422–0.109112118133 W/K; spread 5.07117e-08 W/K.
- Near-optimal definition: MSE no more than 0.0100% above best.
- Near-optimal starts: 25.
- Near-optimal Cth range: 117.523100528–117.523257798 J/K.
- Near-optimal Gth range: 0.109112067422–0.109112118133 W/K.

## Local conditioning audit

Finite central differences used step `1.0e-05` in `log(Cth)` and `log(Gth)`.

- Numerical rank: **2 of 2**.
- Singular values: `42.2173717424`, `6.89802313643` K per unit log parameter.
- Condition number: `6.12021312591`.
- Sensitivity-column correlation: `0.304113502931`.
- Local estimator parameter correlation from `(JᵀJ)^−1`: `-0.866785424705`.
- Sensitivity RMS: `log(Cth)` 0.839257377997 K; `log(Gth)` 2.32716089319 K.

This is a practical local-conditioning audit for one heating experiment, not a claim of structural identifiability.

The covariance-derived correlation is appreciable, but it does not produce practical solution ambiguity here: the Jacobian condition number is low, both singular values are resolved, and all 25 broadly separated starts collapse to the same narrow parameter ranges. That evidence supports the stated well-conditioned classification for this experiment while the parameters remain explicitly effective and experiment-conditional.

## Temperature diagnostics

- Initial measured/model temperature: 26.863600000/26.863600000 °C.
- Final measured/model temperature: 30.986661355/30.611425097 °C.
- Measured/predicted temperature rise: 4.123061355/3.747825097 K.
- Measured/predicted maximum temperature: 30.986661355/30.611425097 °C.
- Maximum absolute primary error: 0.383292660 K at elapsed time 48.968925627 min and SOC 0.109890294155.
- Pearson correlation of primary residual with SOC: -0.074680059.

### Residuals by discharge region

| Region | SOC interval | Samples | RMSE (K) | Bias (K) | Maximum absolute error (K) |
|---|---|---:|---:|---:|---:|
| Early/high | [0.633333, 0.90] | 100 | 0.181418378 | -0.035017144 | 0.353231791 |
| Middle | [0.366667, 0.633333) | 99 | 0.100961014 | 0.022006458 | 0.249041075 |
| Late/low | [0.10, 0.366667) | 100 | 0.136644784 | -0.012275900 | 0.383292660 |

These regional values and the residual plots document any early-, middle-, late-, or SOC-dependent structure without adding weights or correction terms.

The residual is not strongly monotonic with SOC: regional biases alternate from slightly negative (early) to slightly positive (middle) and back to slightly negative (late), and the residual/SOC correlation is small. The largest error is late in the discharge near SOC 0.110, but the regional RMSE is highest in the early/high-SOC third rather than progressively worsening toward low SOC.

## Secondary no-refit model-voltage heat diagnostic

After freezing the best values above, the same thermal parameters and measured-temperature reversible heat were propagated once with

```text
Qdot_total_model = I * (Ueq_average - V_model) + Qdot_rev.
```

No parameter was refitted.

| Metric | Value |
|---|---:|
| RMSE | 1.187950709 K |
| MAE | 1.012001174 K |
| Mean signed error | -1.011259502 K |
| Maximum absolute error | 2.758052735 K |
| R² | -1.414930121 |
| Final model temperature | 28.228608620 °C |
| Predicted temperature rise | 1.365008620 K |
| Maximum-error time | 49.582140200 min |
| Maximum-error SOC | 0.100000000000 |

Relative to the primary fit, switching to model-voltage irreversible heat changes RMSE by +1.044328644 K and final predicted temperature by -2.382816477 K. This loss is attributed to the already documented electrical-voltage transfer error expressed through heat input; it was not repaired thermally.

## Interpretation and stop condition

**A. Thermal-submodel calibration:** the fitted effective parameters use measured-voltage irreversible heat.

**B. Coupled electrothermal forward diagnostic:** the no-refit result uses model-voltage irreversible heat with exactly the same effective parameters.

The calibration stops here. No UDDS or CC–CV thermal fitting/validation, electrical refit, Arrhenius dependence, entropy modification, OCP change, or own-TRC data use was performed.
