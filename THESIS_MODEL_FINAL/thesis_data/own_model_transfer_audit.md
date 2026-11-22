# Panasonic-to-TRC18650 hierarchical model-transfer audit

## Literal transfer is not adequate

Only `cycle_cell1_ambient` is used. No Panasonic parameter, OCP, window, support interval, or model equation was changed; no TRC quantity was fitted and no measured voltage was corrected.

Source: `C:\Users\amvel\Documents\Andes\Tesis\TRC18650experimental data\amvTRC18650LiIonCycle1C Tambient cell1.csv`. Raw Chroma current is negative on discharge and positive on charge; the sole model-boundary conversion is `I_model=-I_raw`, giving discharge-positive model current. The file has no measured time. The uniform `2.994702465 s` integration step is therefore explicitly provisional and derived from capacity/current consistency.

## Frozen Panasonic quantities transferred

| Quantity | Frozen value / provenance |
|---|---|
| Negative OCP | Guo-White 2005 MCMB discharge-process spline via `ocp_panasonic_qi_frozen.U_n` |
| Positive OCP | Frozen 13-knot Panasonic Qi-style representation, selected multistart start 3 |
| Nominal effective windows | x_n 0.720451136448 to 0.251849704066; x_p 0.224027933455 to 0.953477823621 |
| OCP support | x_n (0.01, 0.8); x_p (0.18, 0.96); no clipping/extrapolation |
| Q_cell | 2.994979138519 Ah |
| Q_n / Q_p | 23008.732269258333 / 14780.898652565113 C |
| alpha_n / alpha_p | 3000 / 3204.22765511 s |
| d_n / d_p | 0.000183029692927 / 0.0007 1/s |
| R0 | 0.0259594900586 ohm |
| Current convention | Model-positive discharge; raw TRC sign converted only at loader boundary |

The dynamic values are the conditional Panasonic HPPC parameterization, not globally identified Panasonic constants.

## Test A — voltage-backbone transfer diagnostic

The frozen Panasonic equilibrium-voltage backbone is compared with the clean TRC **loaded 1C** discharge in normalized capacity coordinates. This is not OCP validation.

Overall: RMSE `188.882 mV`, MAE `179.277 mV`, bias Ueq-minus-measured `179.277 mV`, maximum absolute error `356.685 mV`, R2 `0.610108`.

| SOC region | Samples | RMSE [mV] | Bias [mV] |
|---|---:|---:|---:|
| high | 334 | 146.479 | 144.493 |
| middle | 333 | 143.922 | 141.490 |
| low | 334 | 254.588 | 251.734 |

The residual is SOC-shaped rather than a constant offset; the regional table and residual figure show the change in magnitude across the discharge.

## Test B — literal frozen electrical transfer

Classification: **PROVISIONAL-TIME DIRECT TRANSFER DIAGNOSTIC**. Metrics are computed only where the frozen surface states remain inside the frozen OCP support; unsupported samples are retained as invalid and never clipped or extrapolated.

Valid/total discharge surface-state samples: `1182/1204`; invalid count `22` (negative-electrode `0`, positive-electrode `22`); first violation `row-local index 1182, SOC_ref=0.017001`.
Valid-domain metrics: RMSE `82.134 mV`, MAE `59.903 mV`, bias predicted-minus-measured `56.982 mV`, maximum absolute error `207.020 mV`, R2 `0.920997`.

| SOC region | Valid samples | RMSE [mV] | Bias [mV] |
|---|---:|---:|---:|
| high | 401 | 30.306 | 17.680 |
| middle | 401 | 44.215 | 26.140 |
| low | 380 | 133.983 | 131.003 |

State ranges over the discharge before OCP masking: x_n_avg `0.266403-0.720451`, x_p_avg `0.224028-0.930823`, x_n_surface `0.241198-0.716850`, x_p_surface `0.230015-0.972730`.

## Test C — normalized-capacity state-mapping diagnostic

The TRC discharge passes `2.899900 Ah`, only `96.825%` of the frozen Panasonic cell capacity. Therefore literal Q_n/Q_p integration does not traverse the complete frozen nominal electrode windows over the TRC complete discharge.

The average states below are calculated directly from normalized TRC capacity and frozen Q_n/Q_p. Surface values are the corresponding nearest samples from the provisional-time dynamic transfer.

| z | Capacity-mapped x_n_avg | Capacity-mapped x_p_avg | Provisional x_n_surface | Provisional x_p_surface | Surface support valid? |
|---:|---:|---:|---:|---:|---|
| 0.00 | 0.720451136 | 0.224027933 | 0.716850392 | 0.230014617 | True |
| 0.25 | 0.607019854 | 0.400601095 | 0.581641676 | 0.442773586 | True |
| 0.50 | 0.493588572 | 0.577174256 | 0.468408448 | 0.619042254 | True |
| 0.75 | 0.380157290 | 0.753747417 | 0.354803170 | 0.795886051 | True |
| 1.00 | 0.266726008 | 0.930320578 | 0.241197882 | 0.972730033 | False |

At z=1 the average states are `0.266726008` and `0.930320578`, rather than the frozen endpoints `0.251849704` and `0.953477824`. Average states remain in support, but the loaded positive surface state becomes unsupported near low SOC. The mapping is therefore internally computable but not coherent as a complete frozen-window representation of the TRC discharge.

## Error decomposition

Measured-minus-backbone has mean `-176.379 mV` and RMS `185.029 mV` on the valid full-model subset. Full-minus-backbone has mean `-119.397 mV` and RMS `119.882 mV`, while full-minus-average-state-equilibrium has RMS `135.195 mV`. Literal full-model residual RMS is `82.134 mV`.

The voltage-backbone mismatch is already large. Frozen dynamic/polarization terms subtract about 119 mV on average and reduce the terminal-voltage error, but this is compensating error rather than evidence of uniquely correct transferred dynamics. They do not remove the systematic low-SOC discrepancy; capacity mapping and the eventual positive-surface support violation add a separate failure mode. The supported classification is therefore a **mixture of static voltage-reference, state-mapping, and dynamic-transfer mismatch**. These quantities are transfer diagnostics, not intrinsic TRC overpotentials.

## Directional transfer check

Recharge valid/total samples: `1477/1477`. The discharge-oriented frozen representation was propagated continuously through the post-discharge rest and recharge, using the same provisional timebase and no direction-specific correction. State integration remains mathematically defined through the 22 unsupported discharge-voltage samples; no OCP or voltage was evaluated there. The rested states re-enter support before recharge.

| Recharge portion | Samples | RMSE [mV] | MAE [mV] | Bias [mV] | Max abs [mV] | R2 |
|---|---:|---:|---:|---:|---:|---:|
| all_recharge | 1477 | 40.092 | 36.744 | -34.933 | 147.100 | 0.976533 |
| cc | 1041 | 39.883 | 35.492 | -32.924 | 147.100 | 0.969330 |
| cv | 436 | 40.587 | 39.732 | -39.732 | 52.410 | -4259.955300 |

No recharge surface-support violation occurred. The charge results are a held-out directional diagnostic. The very negative CV R2 is ill-conditioned because measured CV voltage has very little variance; RMSE and bias are the interpretable CV metrics. These results do not justify a charge-specific fit or thermodynamic-hysteresis claim.

## Minimum-adaptation decision tree

- **LEVEL 0 — rejected:** literal direct transfer is not sufficient.
- **LEVEL 1 — required:** the different measured cell capacity and incomplete traversal of the frozen windows require cell-specific capacity/state mapping.
- **LEVEL 2 — required:** the large, SOC-shaped backbone residual requires a TRC-compatible voltage-reference/OCP basis before equilibrium and polarization errors can be separated.
- **LEVEL 3 — indicated but not identifiable here:** frozen dynamics do not remove the residual and surface support fails at low SOC, but the absent measured timebase prevents defensible dynamic calibration.
- **LEVEL 4 — applies to full identification from this dataset alone:** the own clean cycle is insufficient to identify electrode OCP/windows and dynamic/thermal parameters uniquely.

The minimum defensible *adaptation* before meaningful quantitative reduced-model validation is cumulative **LEVEL 1 + LEVEL 2**. A complete TRC dynamic parameterization would also need LEVEL 3 evidence from data that are not available here; it must not be inferred by fitting this provisional-time cycle.

## Identifiability matrix

| Parameter/component | Information required | Information available | Identifiable here? | Reason |
|---|---|---|---|---|
| Q_cell | complete charge/discharge throughput | complete clean discharge and recharge counters; Q_ref=2.8999 Ah | YES | Capacity is directly supported in capacity coordinates without measured time. |
| OCP / equilibrium-voltage representation | near-equilibrium low-rate or rested voltage across SOC, preferably both directions | loaded 1C discharge and CC-CV charge only | NO | Loaded terminal voltage mixes equilibrium voltage, resistance, kinetics, diffusion, and directional effects. |
| stoichiometric windows | cell-compatible electrode OCPs plus independent capacity/electrode constraints | cell capacity and loaded full-cell voltage only | NO | Many electrode windows can generate similar full-cell voltage; no electrode reference is available. |
| Q_n / Q_p grouped state capacities | Q_cell plus defensible electrode stoichiometric window widths | Q_cell is measured; electrode window widths are not identified | NO | The Panasonic Q_n/Q_p values cannot be rescaled uniquely from cell capacity alone. |
| R0 | time-resolved current step and immediate voltage response with fixture contribution characterized | adjacent-record R_app only; no measured sampling time | CONDITIONAL DIAGNOSTIC ONLY | R_app contains cell, fixture, wiring, timing, and fast polarization contributions and is not intrinsic R0. |
| alpha_n | measured-time relaxation/transient excitation with electrode separation | sustained 1C branch; provisional time only | NO | No measured timebase or electrode-specific excitation. |
| alpha_p | measured-time relaxation/transient excitation with electrode separation | sustained 1C branch; provisional time only | NO | No measured timebase or electrode-specific excitation. |
| d_n | time-resolved kinetic/polarization information over SOC | loaded curve with static/dynamic effects confounded | NO | Cannot separate negative kinetics from OCP, positive kinetics, diffusion, and setup effects. |
| d_p | time-resolved kinetic/polarization information over SOC | loaded curve with static/dynamic effects confounded | NO | Cannot separate positive kinetics from OCP, negative kinetics, diffusion, and setup effects. |
| thermal Cth | measured time, heat input, and calibrated temperature response | surface temperatures but no measured timestamps; heat generation not identified | NO | Thermal storage cannot be inferred in capacity coordinates alone. |
| thermal Gth | measured-time cooling/heating response and ambient boundary condition | ambient channel exists but no measured timestamps | NO | Heat-transfer rate requires a real timebase. |

## Accuracy-reporting policy

| Quantity | Status | Basis |
|---|---|---|
| Voltage-vs-SOC RMSE/MAE/bias | QUANTITATIVE / CONDITIONAL | Valid comparison of measured loaded voltage and prediction in capacity coordinates; not an OCP-validation metric. |
| Capacity error | QUANTITATIVE | Complete Chroma Ah counter supports 2.8999 Ah directly. |
| SOC trajectory error | NOT INDEPENDENTLY VALIDATED | SOC is a coulomb-counted reference coordinate, not independently measured true SOC. |
| Dynamic time-domain RMSE | PROVISIONAL / NOT PRIMARY | Chronology uses a capacity/current-derived uniform timebase, not measured time. |
| Temperature-vs-time RMSE | NOT SUPPORTED | Measured timestamps are absent and this audit is electrical only. |
| Apparent series-path response | QUANTITATIVE DIAGNOSTIC | Adjacent-sample DeltaV/DeltaI is reproducible but is not intrinsic cell R0. |
| Voltage-backbone residual | QUANTITATIVE DIAGNOSTIC | Numerically valid transfer diagnostic; TRC loaded 1C voltage is not equilibrium OCV. |
| OCP surface-domain violation count | QUANTITATIVE MODEL DIAGNOSTIC | Evaluated without clipping or extrapolation. |

## Explicit answers

1. **Is literal transfer adequate?** No.
2. **How much discrepancy exists at backbone level?** Overall backbone RMSE is `188.882 mV` with bias `179.277 mV`; the mismatch is already large and SOC-shaped before dynamics.
3. **Is the Panasonic state mapping coherent over the full TRC discharge?** No as a complete-window transfer. The TRC passes `96.825%` of Panasonic nominal capacity, finishes short of both nominal average-state endpoints, and produces low-SOC positive-surface support violations under load.
4. **Would dynamics also require adaptation?** Evidence indicates yes for quantitative terminal-voltage transfer, but the parameters cannot be identified defensibly with provisional time.
5. **What can be adapted defensibly?** Q_cell can be set from the complete 2.8999 Ah discharge. A capacity-coordinate empirical voltage reference can be constructed for diagnostics, but it must not be called OCP.
6. **What cannot be identified?** Electrode OCPs/windows, intrinsic R0, alpha_n/alpha_p, d_n/d_p, Cth, and Gth cannot be uniquely identified from this dataset.
7. **Minimum next adaptation level?** Cumulative LEVEL 1 + LEVEL 2. LEVEL 3 is indicated but unsupported for calibration; full identification remains LEVEL-4 data-limited.

## Stop point

No parameter optimization, OCP fit, voltage correction, affected-dataset calibration, thermal accuracy calculation, or climate/city matching was performed.
