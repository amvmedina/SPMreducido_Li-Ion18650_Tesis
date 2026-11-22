# Panasonic 25 C UDDS thermal validation

## THERMAL TRANSFER FAILED

## COUPLED ELECTROTHERMAL TRANSFER FAILED

This is an independent forward-only validation of the frozen 1C-derived effective lumped thermal parameters. No electrical, OCP, entropy, or thermal parameter was fitted or changed; there is no clipping, extrapolation, baseline correction, or temperature-offset fit.

## Dataset and supported interval

- Exact source: `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\Drive cycles\03-21-17_00.29 25degC_UDDS_Pan18650PF.mat`
- Manifest: `C:\Users\amvel\Documents\Andes\THESIS_MODEL_FINAL\experimental_data_manifest.md`
- Columns: `TimeStamp, Voltage, Current, Ah, Wh, Power, Battery_Temp_degC, Time, Chamber_Temp_degC`
- Raw samples: 224,187; one exact duplicate final record was removed before interval construction.
- Entropy-supported interval: exactly SOC 0.900000000000 to 0.100000000000, created with linear interpolation at both boundaries. Source times are 2643.105209092 to 22144.751588421 s.
- Duration: 19501.646379329 s (5.417124 h); samples including the two exact boundary records: 194,744.
- Actual positive sample interval: median 0.098007917 s, mean 0.100140423 s, range 0.003789892-2.323000878 s.
- Internal current (discharge positive, regeneration negative): -4.124970 to 6.283890 A. Raw current was sign-converted once at the model boundary and was otherwise unchanged.
- Measured voltage: 2.833602 to 4.144090 V.
- Battery temperature: 25.619490 to 27.110180 degC; chamber temperature: 25.000000 to 25.000000 degC.
- Gross discharge throughput: 3.526425306 Ah; regenerative charge throughput: 1.133866608 Ah; net discharge: 2.392558697 Ah.
- OCP-domain violations: 0; entropy-domain violations: 0; clipping count: 0.

Actual nonuniform timestamps and all intervening regenerative excursions are preserved. Heat and ambient are held at the left sample for each exact analytical thermal-state update.

## Frozen configuration

- `Cth = 117.523178000996 J/K`, `Gth = 0.109112089746 W/K`, `tau_th = 1077.086675492568 s`.
- `Q_cell = 2.994979138519 Ah`, `Q_n = 23008.732269258333 C`, `Q_p = 14780.898652565113 C`.
- Dynamic electrical values: `R0=0.0259594900586 ohm`, `alpha_p=3204.22765511 s`, `d_n=0.000183029692927 1/s`, with fixed `alpha_n=3000 s`, `d_p=0.0007 1/s`.
- Frozen Guo-White MCMB negative OCP, frozen 13-knot positive OCP, nominal effective stoichiometric windows, and verified full-cell entropy table were used unchanged.
- Every thermal path starts at the exact interpolated measured battery temperature, 25.832839762 degC; no initial offset is fitted.
- The temperature-aware local voltage calculation reproduces the existing frozen electrical implementation to a maximum absolute difference of 1.776e-15 V before thermal coupling.

| Frozen electrical state | Minimum | Maximum | Supported OCP interval |
|---|---:|---:|---:|
| `x_n_avg` | 0.299333920774 | 0.673679244000 | [0.01, 0.80] |
| `x_n_surf` | 0.289534080018 | 0.673909921179 | [0.01, 0.80] |
| `x_p_avg` | 0.296835546602 | 0.879561368680 | [0.18, 0.96] |
| `x_p_surf` | 0.296295406542 | 0.895765359060 | [0.18, 0.96] |

## Overall temperature metrics

Residual is predicted minus measured battery temperature. Validation A uses measured terminal voltage and measured temperature in reversible heat. Validation B recursively uses frozen-model voltage and predicted temperature, including predicted temperature in the frozen Butler-Volmer terms. The isolation diagnostic uses model voltage but measured temperature in reversible heat.

| Path | RMSE [K] | MAE [K] | Bias [K] | Max abs [K] | R2 | Measured final [degC] | Predicted final [degC] | Measured rise [K] | Predicted rise [K] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Validation A | 0.526697652 | 0.512960900 | -0.512714512 | 0.875239690 | -4.119410369 | 26.469807973 | 26.112601243 | +0.636968211 | +0.279761481 |
| Validation B | 0.674696405 | 0.651358329 | -0.651117626 | 1.394657444 | -7.400670890 | 26.469807973 | 25.499201217 | +0.636968211 | -0.333638545 |
| MODEL-VOLTAGE / MEASURED-T REVERSIBLE DIAGNOSTIC | 0.674816837 | 0.651492375 | -0.651251673 | 1.394446351 | -7.403670185 | 26.469807973 | 25.499231580 | +0.636968211 | -0.333608182 |

## Metrics by SOC region

| Path | Region | Samples | RMSE [K] | MAE [K] | Bias [K] | Max abs [K] | R2 |
|---|---|---:|---:|---:|---:|---:|---:|
| Validation A | high (SOC >= 2/3) | 58,753 | 0.474352541 | 0.455115352 | -0.454298668 | 0.799236827 | -18.982308216 |
| Validation A | middle (1/3 <= SOC < 2/3) | 81,634 | 0.552954126 | 0.547511595 | -0.547511595 | 0.850228281 | -23.054293141 |
| Validation A | low (0.10 <= SOC < 1/3) | 54,357 | 0.539884337 | 0.523595926 | -0.523595926 | 0.875239690 | -2.820720020 |
| Validation B | high (SOC >= 2/3) | 58,753 | 0.541744373 | 0.518782934 | -0.517985092 | 0.904487099 | -25.063456353 |
| Validation B | middle (1/3 <= SOC < 2/3) | 81,634 | 0.654606516 | 0.649082845 | -0.649082845 | 0.985943415 | -32.711264387 |
| Validation B | low (0.10 <= SOC < 1/3) | 54,357 | 0.818614461 | 0.798072812 | -0.798072812 | 1.394657444 | -7.784207440 |

## Heat-component audit

The energy integral is consistent with the thermal propagator's left-sample zero-order hold. `Signed reversible fraction` is E_rev/E_total; the absolute-component fraction is integral(|Qrev|)/integral(|Qirr|+|Qrev|).

| Path | Eirr [J] | Erev [J] | Etotal [J] | Mean Qirr [W] | Mean Qrev [W] | Mean total [W] | Signed reversible fraction | Absolute reversible fraction | Exothermic time fraction | Endothermic time fraction |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Measured-heat path | 1340.589785411 | -210.217865876 | 1130.371919535 | 0.068742390 | -0.010779493 | 0.057962897 | -0.185972300 | 0.256090571 | 0.292435367 | 0.706556355 |
| Model-heat path | 973.624975893 | -209.787063508 | 763.837912385 | 0.049925271 | -0.010757403 | 0.039167868 | -0.274648666 | 0.323993515 | 0.292435367 | 0.706556355 |

The explicit sign check found 0 Validation-A and 0 Validation-B mismatches against `sign(-I*dUeq/dT)`. Thus reversible heat changes sign consistently with both current direction and entropy-coefficient sign.

| Current / entropy-coefficient case | Expected Qrev sign | Samples | Time fraction | Mean Qrev A [W] | Mean Qrev B [W] |
|---|---|---:|---:|---:|---:|
| discharge dudt positive | endothermic | 132,146 | 0.678793007 | -0.025127421 | -0.025074729 |
| discharge dudt negative | exothermic | 18,975 | 0.097491870 | +0.011571747 | +0.011541433 |
| regeneration dudt positive | exothermic | 38,018 | 0.194943496 | +0.028124642 | +0.028064890 |
| regeneration dudt negative | endothermic | 5,415 | 0.027763348 | -0.012972255 | -0.012936914 |

## Dynamic diagnostics

Absolute temperature error versus SOC depletion has Spearman rho +0.175154 for A and +0.613152 for B. On loaded samples, absolute error versus current magnitude has rho +0.084430 for A and +0.105054 for B.

Regenerative-current RMSE is 0.530120 K (A) and 0.680711 K (B), with biases -0.516084 and -0.656561 K. Discharge-current RMSE is 0.526190 K (A) and 0.673583 K (B).

At the ten largest real adjacent-sample current changes, temperature-error RMSE is 0.494727 K (A) and 0.663944 K (B). Thermal inertia prevents instantaneous temperature jumps; these rows diagnose whether accumulated errors are preferentially present at major electrical transitions.

| Rank | Source time [s] | Interval elapsed [s] | SOC | I before [A] | I after [A] | delta I [A] | Measured T [degC] | Predicted A [degC] | Error A [K] | Predicted B [degC] | Error B [K] |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 22066.957998 | 19423.852789 | 0.101626464 | -0.100450 | -4.124970 | -4.024520 | 26.471310 | 26.059211 | -0.412099 | 25.481697 | -0.989613 |
| 2 | 20695.123999 | 18052.018790 | 0.166338100 | -0.097180 | -3.989400 | -3.892220 | 26.056610 | 25.755039 | -0.301571 | 25.486516 | -0.570094 |
| 3 | 17951.105992 | 15308.000783 | 0.287106887 | -0.093100 | -3.795850 | -3.702750 | 25.843660 | 25.415880 | -0.427780 | 25.293299 | -0.550361 |
| 4 | 16579.449994 | 13936.344785 | 0.345671569 | -0.092280 | -3.781150 | -3.688870 | 25.843660 | 25.335782 | -0.507878 | 25.235433 | -0.608227 |
| 5 | 15207.329997 | 12564.224788 | 0.403591839 | -0.092280 | -3.737050 | -3.644770 | 25.843660 | 25.287520 | -0.556140 | 25.203646 | -0.640014 |
| 6 | 19323.007995 | 16679.902786 | 0.227507140 | -0.095550 | -3.731340 | -3.635790 | 26.056610 | 25.556535 | -0.500075 | 25.395702 | -0.660908 |
| 7 | 13835.588001 | 11192.482792 | 0.460851003 | -0.090650 | -3.716640 | -3.625990 | 25.854860 | 25.275659 | -0.579201 | 25.198028 | -0.656832 |
| 8 | 12463.590997 | 9820.485788 | 0.517502482 | -0.089830 | -3.666820 | -3.576990 | 25.843660 | 25.291058 | -0.552602 | 25.214721 | -0.628939 |
| 9 | 11091.585996 | 8448.480787 | 0.573426077 | -0.088200 | -3.608840 | -3.520640 | 25.854860 | 25.326234 | -0.528626 | 25.231080 | -0.623780 |
| 10 | 9719.763996 | 7076.658787 | 0.628468197 | -0.087380 | -3.572900 | -3.485520 | 25.854860 | 25.338928 | -0.515932 | 25.247792 | -0.607068 |

## Interpretation

The observed battery-temperature excursion within the supported interval is 1.490690 K. Validation A has RMSE 0.526698 K (35.33% of that excursion), bias -0.512715 K, and R2 -4.119410. **Classification: THERMAL TRANSFER FAILED.**

Validation B has RMSE 0.674696 K (45.26% of the excursion), bias -0.651118 K, and R2 -7.400671. The isolation diagnostic RMSE is 0.674817 K; its difference from B quantifies the incremental reversible-temperature/electrical-temperature feedback effect, while its difference from A reflects model-voltage heat. **Classification: COUPLED ELECTROTHERMAL TRANSFER FAILED.**

Validation A's error is a broad negative bias rather than a low-SOC-only failure: regional RMSE is 0.474/0.553/0.540 K (high/middle/low), and the SOC-depletion correlation is only +0.175. Validation B develops clearer SOC-dependent drift, rising to 0.819 K at low SOC with rho +0.613.

Current-magnitude correlations are weak (+0.084 A, +0.105 B), regenerative and discharge RMSE are similar, and largest-step errors are not elevated above overall RMSE. The discrepancy is therefore dominated by persistent heat/temperature underprediction, not isolated current-step or regeneration events. Model-voltage irreversible energy is -27.37% relative to the measured-voltage path; the nearly identical B and isolation trajectories show that model-voltage heat, not reversible-temperature feedback, causes almost all of the additional coupled-path error.

The classifications use error scale relative to the measured thermal excursion together with bias, R2, regional residual structure, and dynamic diagnostics; they are not based on a fitted or task-specific acceptance target.

No repair, refit, entropy/OCP change, Arrhenius term, CC-CV thermal validation, or own-TRC analysis was performed.
