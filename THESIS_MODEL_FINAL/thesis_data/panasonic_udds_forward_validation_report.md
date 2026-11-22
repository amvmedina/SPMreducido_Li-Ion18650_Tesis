# Panasonic 25 C UDDS forward-only validation

## UDDS TRANSFER LIMITED

The manifest-selected UDDS dataset was evaluated once with the frozen nominal conditional Panasonic parameterization. No optimization, baseline correction, clipping, parameter release, or model/OCP change was performed.

## Dataset audit

- Exact source: `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\Drive cycles\03-21-17_00.29 25degC_UDDS_Pan18650PF.mat`
- Manifest entry verified in: `C:\Users\amvel\Documents\Andes\THESIS_MODEL_FINAL\experimental_data_manifest.md`
- Columns: `TimeStamp, Voltage, Current, Ah, Wh, Power, Battery_Temp_degC, Time, Chamber_Temp_degC`
- Raw samples: 224,187; usable samples: 224,186. One exact duplicate final row was removed.
- Total recorded and usable duration: 22447.866000608 s.
- Positive sampling interval: median 0.098993629 s, mean 0.100130990 s, range 0.031995401-2.323000878 s. All actual positive time steps were preserved.
- Raw current convention: discharge negative, regenerative charge positive. Range -7.52923 to 4.12497 A. Sign conversion occurred only at the model-loader boundary.
- Measured voltage range: 2.50334-4.20779 V.
- Battery case temperature: 25.61949-27.11018 degC; chamber temperature: 25-25 degC.
- Exported `Ah` starts at -0.00000 Ah and ends at its minimum, -2.70149 Ah. It is non-monotonic within the drive (36,754 increasing and 124,282 decreasing increments) because regenerative charge is retained. Net exported discharge is 2.701490000 Ah.
- Current/time integration gives 2.697669057 Ah net, comprising 3.964335041 Ah gross discharge and 1.266665984 Ah regenerative charge.
- Source SOC is `1 + Ah/Q_cell`; this directly preserves cumulative-Ah continuity through regenerative intervals. Model states start at the frozen full-charge nominal effective endpoints and propagate continuously with the measured current and actual times.

## Frozen forward evaluation

- `R0=0.0259594900586 ohm`, `alpha_p=3204.22765511 s`, `d_n=0.000183029692927 1/s`
- Fixed `alpha_n=3000 s`, `d_p=0.0007 1/s`
- Frozen `Q_cell=2.994979138519 Ah`, `Q_n=23008.732269258333 C`, `Q_p=14780.898652565113 C`
- Frozen Guo-White MCMB negative OCP, 13-knot positive OCP, and nominal effective stoichiometric windows were used unchanged.

## Primary electrical metrics

Residual is predicted minus measured voltage. Primary metrics include the complete usable drive and its recorded 300 s post-cutoff recovery.

| Metric | Complete usable profile |
|---|---:|
| Samples | 224186 |
| RMSE [mV] | 30.988192 |
| MAE [mV] | 21.875765 |
| R2 | 0.983947098 |
| Mean signed error [mV] | 19.762413 |
| Maximum absolute error [mV] | 632.385889 |

Through the cutoff sample alone (221,186 samples), RMSE is 30.808410 mV, MAE is 21.808497 mV, and R2 is 0.983742293. The cutoff residual is 576.318635 mV and was not forced to match.

### By SOC region

| Region | Samples | RMSE [mV] | MAE [mV] | R2 | Mean signed [mV] | P95 absolute [mV] | Max absolute [mV] |
|---|---:|---:|---:|---:|---:|---:|---:|
| high [2/3, 1] | 85,164 | 14.360656 | 10.852647 | 0.981343199 | 5.638385 | 29.764529 | 117.712870 |
| middle [1/3, 2/3) | 81,634 | 21.522702 | 19.752408 | 0.940377255 | 19.677567 | 34.597696 | 127.032107 |
| low [0, 1/3) | 57,388 | 52.785444 | 41.254574 | 0.756271951 | 40.843214 | 104.888996 | 632.385889 |

## Residual versus operating condition

Near rest is defined as `|I_raw| <= 0.05 A`; all other points are loaded. Discharge and regenerative subsets use the same 0.05 A threshold.

| Subset | Samples | RMSE [mV] | MAE [mV] | Mean signed [mV] | Max absolute [mV] |
|---|---:|---:|---:|---:|---:|
| near rest | 11,466 | 30.637659 | 21.602310 | 19.568306 | 632.385889 |
| loaded | 212,720 | 31.006973 | 21.890505 | 19.772876 | 576.318635 |
| discharge | 168,083 | 33.237593 | 23.274812 | 21.795447 | 576.318635 |
| regenerative charge | 44,637 | 20.537993 | 16.677824 | 12.156774 | 224.682794 |
| post cutoff recovery | 3,000 | 42.182795 | 26.835342 | 26.835342 | 632.385889 |

### Current-magnitude bins

Loaded-sample empirical quartiles are 0.12249, 0.63778, and 1.32292 A. These data-derived edges provide four similarly populated loaded bins, plus the independently defined near-rest bin.

| Bin | |I| interval [A] | Samples | RMSE [mV] | MAE [mV] | Mean signed [mV] | P95 absolute [mV] |
|---|---:|---:|---:|---:|---:|---:|
| near_rest | [0.00000, 0.05000] | 11,466 | 30.637659 | 21.602310 | 19.568306 | 49.303103 |
| loaded_q1 | (0.05000, 0.12249] | 53,190 | 18.372057 | 14.893422 | 11.661159 | 36.142092 |
| loaded_q2 | (0.12249, 0.63778] | 53,232 | 25.695610 | 20.319548 | 18.653268 | 50.078319 |
| loaded_q3 | (0.63778, 1.32292] | 53,125 | 29.999833 | 22.719395 | 20.789476 | 59.553877 |
| loaded_q4 | (1.32292, infinity] | 53,173 | 44.140922 | 29.634384 | 27.992353 | 91.314328 |

Spearman correlation between absolute residual and decreasing SOC is +0.635873 overall and +0.642564 on loaded samples. The corresponding correlation with current magnitude is +0.234970 overall and +0.266361 on loaded samples.

Residual magnitude therefore grows with both decreasing SOC and increasing current magnitude, but the SOC association is stronger. Within each SOC region, the loaded-sample absolute-error/current-magnitude correlations remain positive, so the current trend is not solely an SOC-composition artifact.

## Largest current-step voltage changes

The events below are the ten largest adjacent-sample `|delta I_raw|` transitions. Voltage change is evaluated across exactly the same recorded sample pair; no baseline correction is applied.

| Rank | Time [s] | SOC | I before [A] | I after [A] | delta I [A] | Measured delta V [mV] | Predicted delta V [mV] | Delta-V error [mV] |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 22147.966994 | 0.097994 | -7.52923 | 0.00000 | +7.52923 | +180.250000 | +236.317254 | +56.067254 |
| 2 | 22066.957998 | 0.101626 | 0.10045 | 4.12497 | +4.02452 | +18.660000 | +127.813041 | +109.153041 |
| 3 | 20695.123999 | 0.166338 | 0.09718 | 3.98940 | +3.89222 | +23.810000 | +124.040303 | +100.230303 |
| 4 | 17951.105992 | 0.287107 | 0.09310 | 3.79585 | +3.70275 | +21.230000 | +110.730566 | +89.500566 |
| 5 | 16579.449994 | 0.345672 | 0.09228 | 3.78115 | +3.68887 | +17.370000 | +110.362723 | +92.992723 |
| 6 | 15207.329997 | 0.403592 | 0.09228 | 3.73705 | +3.64477 | +21.880000 | +108.627613 | +86.747613 |
| 7 | 19323.007995 | 0.227507 | 0.09555 | 3.73134 | +3.63579 | +20.590000 | +112.617228 | +92.027228 |
| 8 | 13835.588001 | 0.460851 | 0.09065 | 3.71664 | +3.62599 | +13.510000 | +109.840193 | +96.330193 |
| 9 | 12463.590997 | 0.517502 | 0.08983 | 3.66682 | +3.57699 | +12.870000 | +108.540774 | +95.670774 |
| 10 | 11091.585996 | 0.573426 | 0.08820 | 3.60884 | +3.52064 | +12.220000 | +111.062421 | +98.842421 |

Across these ten events, delta-voltage RMSE is 92.713765 mV and bias is +91.756212 mV. Most selected non-cutoff events are large regenerative-current onsets; the discharge-pulse-derived conditional resistance/kinetic representation substantially overpredicts their immediate upward voltage changes.

## State and OCP-domain audit

| State | Minimum | Maximum | Supported OCP interval where applicable |
|---|---:|---:|---:|
| `x_n_avg` | 0.298367533995 | 0.720451136448 | diagnostic state |
| `x_p_avg` | 0.224027933455 | 0.881065697692 | diagnostic state |
| `x_n_surf` | 0.284014880455 | 0.720935892646 | [0.01, 0.80] |
| `x_p_surf` | 0.223068293108 | 0.904693134286 | [0.18, 0.96] |

- Minimum source-Ah SOC: 0.097993717133; minimum predicted SOC from current/time integration: 0.099269500078.
- OCP-domain violations: **0**.
- Clipping count: **0**.

## Comparison with prior electrical validation

| Validation | Metric type | Result [mV] |
|---|---|---:|
| HPPC calibration, SOC >= 1/3 | Phase-balanced pulse-relative RMSE | 36.835198 |
| Withheld low-SOC HPPC | Phase-balanced pulse-relative RMSE | 85.729433 |
| Withheld low-SOC HPPC | Absolute-voltage RMSE | 59.582115 |
| Independent 1C, overall | Absolute-voltage RMSE | 114.362188 |
| Independent 1C, high / middle / low SOC | Absolute-voltage RMSE | 37.384 / 64.319 / 195.954 |
| UDDS, overall | Absolute-voltage RMSE | 30.988192 |
| UDDS, high / middle / low SOC | Absolute-voltage RMSE | 14.361 / 21.523 / 52.785 |

Metric definitions differ between pulse-relative HPPC and absolute-voltage continuous-profile validation, so their magnitudes are contextual rather than interchangeable. Nevertheless, UDDS shows the same direction of transfer loss: positive mean residual grows from high to middle to low SOC, so the frozen model increasingly overpredicts measured voltage as SOC falls.

## Interpretation and classification

Aggregate UDDS tracking is strong (`RMSE=30.988 mV`, `R2=0.983947`), but it is not uniform. Low-SOC RMSE is 52.785 mV versus 14.361 mV at high SOC, with a positive low-SOC bias of 40.843 mV. The cutoff residual reaches 576.319 mV. In addition, the ten largest current changes have 92.714 mV delta-voltage RMSE, exposing a transient/current-direction limitation that aggregate sample metrics dilute.

The dominant limitation is therefore a **combination**: a low-SOC, SOC-dependent positive/static voltage offset plus deficient transient transfer at the largest current changes, especially regenerative-current onsets; increasing current magnitude is a secondary aggravating factor. The limitation is not caused by OCP extrapolation or clipping.

**Final classification: UDDS TRANSFER LIMITED.**

No repair, optimization, baseline adjustment, thermal fit, or frozen-model change was attempted.
