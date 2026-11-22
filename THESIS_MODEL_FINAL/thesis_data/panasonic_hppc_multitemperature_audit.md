# Panasonic NCR18650PF multi-temperature HPPC audit

## Scope

Exploratory audit only. The frozen electrical and thermal models, OCPs, accepted parameters, and manifests were not modified. No Arrhenius dependence was implemented. UDDS, CC-CV, C/20, and 1C were not fitted.

## Source and channel audit

The five original MAT HPPC files were loaded directly from the newly local full dataset. Finite fractions are ordered as time/current/voltage/battery temperature/chamber temperature. The cold chamber channel is reported but never used as ambient or as a fit input.

| Nominal °C | File | Samples | Duration s | contiguous positive dt min/median/max s | Omitted gaps | finite fractions | Current sign |
|---:|---|---:|---:|---|---:|---|---|
| -20 | `06-15-17_11.31 n20degC_5Pulse_HPPC_Pan18650PF.mat` | 49655 | 58195.0 | 0.004003/0.998/1.009 | 9 | 1.000/1.000/1.000/1.000/0.000 | raw negative = discharge |
| -10 | `06-07-17_12.29 n10degC_5pulse_HPPC_Pan18650PF.mat` | 78494 | 78434.0 | 0.003993/0.998/1.008 | 10 | 1.000/1.000/1.000/1.000/0.000 | raw negative = discharge |
| 0 | `05-20-17_10.44 0degC_5pulse_HPPC_Pan18650PF.mat` | 85260 | 83098.8 | 0.006005/0.996/1.019 | 11 | 1.000/1.000/1.000/1.000/0.000 | raw negative = discharge |
| 10 | `03-27-17_09.06 10degC_5Pulse_HPPC_Pan18650PF.mat` | 91371 | 88792.3 | 0.008992/0.996/1.016 | 12 | 1.000/1.000/1.000/1.000/0.000 | raw negative = discharge |
| 25 | `03-11-17_08.47 25degC_5Pulse_HPPC_Pan18650PF.mat` | 102800 | 97599.4 | 0.007996/0.996/1.016 | 13 | 1.000/1.000/1.000/1.000/1.000 | raw negative = discharge |

Available channels by file:
- -20 °C: `TimeStamp; Voltage; Current; Ah; Wh; Power; Battery_Temp_degC; Time; Chamber_Temp_degC`
- -10 °C: `TimeStamp; Voltage; Current; Ah; Wh; Power; Battery_Temp_degC; Time; Chamber_Temp_degC`
- 0 °C: `TimeStamp; Voltage; Current; Ah; Wh; Power; Battery_Temp_degC; Time; Chamber_Temp_degC`
- 10 °C: `TimeStamp; Voltage; Current; Ah; Wh; Power; Battery_Temp_degC; Time; Chamber_Temp_degC`
- 25 °C: `TimeStamp; Voltage; Current; Ah; Wh; Power; Battery_Temp_degC; Time; Chamber_Temp_degC`

The chamber-temperature field exists in every MAT structure, but its finite fraction is zero at −20, −10, 0, and 10 °C; only the 25 °C channel is finite. It is therefore invalid for the cold campaigns and was not used. Large raw positive time increments are the documented omitted between-group discharges; the table reports only contiguous increments and counts the omitted gaps explicitly.

## Battery-temperature audit

| Nominal °C | Initial °C | Median °C | Min °C | Max °C | End-start K | Median ΔT in 4C/6C pulse K | Maximum ΔT K |
|---:|---:|---:|---:|---:|---:|---:|---:|
| -20 | -20.1260 | -19.9202 | -20.4003 | -18.4340 | +0.4230 | +0.00000 | +0.00000 |
| -10 | -10.1702 | -9.7103 | -10.1945 | -6.7688 | +0.8958 | +0.01210 | +1.13786 |
| 0 | 0.3358 | 0.5622 | 0.1207 | 4.3658 | +0.4415 | +0.83770 | +1.68673 |
| 10 | 10.5592 | 10.7677 | 10.1075 | 13.7905 | +0.4170 | +0.44011 | +1.29716 |
| 25 | 25.6307 | 25.8325 | 25.3393 | 27.9396 | +0.6052 | +0.43712 | +0.82941 |

Both nominal campaign temperature and median measured battery temperature are retained for the Arrhenius diagnostic. No choice between them is imposed here.

## Pulse reconstruction and common support

- -20 °C: 36 detected records; 26 individually usable; 26 in common support.
- -10 °C: 47 detected records; 36 individually usable; 26 in common support.
- 0 °C: 54 detected records; 48 individually usable; 26 in common support.
- 10 °C: 59 detected records; 56 individually usable; 26 in common support.
- 25 °C: 67 detected records; 64 individually usable; 26 in common support.
- Common discrete SOC levels: 100%, 95%, 90%, 80%, 70%, 60%, 50%, 40%, 30%, 25%.
- Common C-rates: 0.5C, 1C, 2C.
- Common matched combinations per temperature: 26; total matched pulse records: 130.
- Common SOC range: 25% to 100%.

| Nominal °C | Missing expected | Truncated | Abnormal V | Abnormal I | Insufficient rest | Invalid SOC | Nonfinite battery T | Support violation |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| -20 | 34 | 10 | 0 | 0 | 0 | 0 | 0 | 0 |
| -10 | 23 | 11 | 0 | 0 | 0 | 0 | 0 | 0 |
| 0 | 16 | 6 | 0 | 0 | 0 | 0 | 0 | 0 |
| 10 | 11 | 3 | 0 | 0 | 0 | 0 | 0 | 0 |
| 25 | 3 | 3 | 0 | 0 | 0 | 0 | 0 | 0 |

All detected and expected-but-missing records, flags, voltages, durations, temperatures, SOC assignments, and apparent transition resistances are in the pulse CSV. No pulse was silently discarded.

## Direct apparent-transition-resistance diagnostic

`R_app = abs((V_immediate - V_pre)/(I_pulse - I_pre))` is reported only as apparent transition resistance, never as intrinsic R0.

| Nominal °C | Median R_app [mΩ] | IQR [mΩ] |
|---:|---:|---:|
| -20 | 87.146614 | 85.107866–96.760793 |
| -10 | 62.442051 | 59.103139–68.762004 |
| 0 | 43.572663 | 42.069679–49.371787 |
| 10 | 31.118464 | 29.865800–34.158977 |
| 25 | 21.494079 | 20.891011–22.298231 |

Matched C-rate medians [mΩ]:

| C-rate | −20 °C | −10 °C | 0 °C | 10 °C | 25 °C |
|---:|---:|---:|---:|---:|---:|
| 0.5C | 85.995171 | 58.916868 | 44.033647 | 30.943615 | 22.073129 |
| 1C | 84.601587 | 62.179721 | 43.537185 | 31.569507 | 20.914108 |
| 2C | 99.756364 | 72.358064 | 48.180492 | 33.752153 | 22.104858 |

Matched nominal-SOC medians [mΩ]:

| SOC | −20 °C | −10 °C | 0 °C | 10 °C | 25 °C |
|---:|---:|---:|---:|---:|---:|
| 100% | 101.089207 | 77.138107 | 51.993791 | 39.774883 | 25.415661 |
| 95% | 99.756364 | 75.473976 | 54.429113 | 36.678165 | 24.211889 |
| 90% | 98.537320 | 72.358064 | 49.768886 | 34.294584 | 22.298527 |
| 80% | 99.978792 | 69.916529 | 43.608141 | 31.569507 | 21.992782 |
| 70% | 96.869973 | 68.473334 | 43.332184 | 29.572266 | 20.885284 |
| 60% | 96.995596 | 70.913144 | 42.662987 | 30.458779 | 20.992719 |
| 50% | 99.532212 | 69.802729 | 41.853343 | 30.010348 | 20.770291 |
| 40% | 87.963091 | 57.878567 | 42.237720 | 30.023572 | 21.382003 |
| 30% | 85.601582 | 58.767117 | 44.670074 | 32.255950 | 21.604001 |
| 25% | 81.766126 | 59.579165 | 45.348613 | 32.059331 | 22.290445 |

Matched-SOC and matched-C-rate values use the same pulse keys at every temperature; no different pulse population is substituted between temperatures. The raw transition metric shows a strong and consistent increase as temperature decreases, before any model fitting.

## Accepted 25 °C regression test

PASS: 40 pulses in 8 histories reproduced `R0=0.0259594900586 ohm`, `alpha_p=3204.22765511 s`, and `d_n=0.000183029692927 1/s`, with 12/12 successful starts. Accepted files were not overwritten.

## Independent common-support conditional fits

The exact accepted phase loss and the same 12 deterministic starts were used independently at every temperature. `alpha_n=3000 s` and `d_p=0.0007 1/s` remain fixed.

| °C | R0 Ω | alpha_p s | d_n 1/s | Balanced mV | Start mV | Pulse mV | Relaxation mV | Starts | Rank | Condition | Numerical bound | OCP edge |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| -20 | 0.0438789236 | 87494.8273 | 2.40588439e-06 | 175.1140 | 201.3560 | 223.1714 | 40.5584 | 12/12 | 3 | 15.149 | False | False |
| -10 | 0.0140763961 | 85907.2833 | 1.23142625e-05 | 104.4019 | 117.7829 | 135.3289 | 22.6398 | 12/12 | 3 | 22.671 | False | False |
| 0 | 1.99404719e-05 | 68754.9712 | 2.64855465e-05 | 66.7253 | 74.8087 | 86.9864 | 13.9212 | 12/12 | 3 | 14066 | False | False |
| 10 | 0.0222183521 | 12645.8211 | 4.02587113e-05 | 39.6857 | 42.8526 | 52.9351 | 9.2943 | 12/12 | 3 | 16.187 | False | False |
| 25 | 0.0187383449 | 3946.48114 | 7.52896455e-05 | 22.1214 | 25.8693 | 27.5510 | 6.3082 | 12/12 | 3 | 34.605 | False | False |

Practical-conditioning details:

- −20 °C: singular values 0.431729, 0.0593092, 0.0284981; sensitivity correlations `(R0,alpha_p)=+0.9543`, `(R0,d_n)=−0.8681`, `(alpha_p,d_n)=−0.8129`; log-bound fractions `0.1243/0.0099/0.3402`; near-best ranges `R0=[0.0438556,0.0439058]`, `alpha_p=[87494.8,87494.8]`, `d_n=[2.40225e−6,2.41443e−6]`.
- −10 °C: singular values 0.369403, 0.0347061, 0.0162940; sensitivity correlations `+0.9511/−0.8913/−0.8365`; log-bound fractions `0.2174/0.0115/0.4415`; near-best ranges `R0=[0.0140764,0.0187743]`, `alpha_p=[79629.3,85907.3]`, `d_n=[1.22869e−5,1.26767e−5]`.
- 0 °C: singular values 0.291999, 0.0239537, 2.07591e−5; sensitivity correlations `+0.9636/−0.9337/−0.8896`; log-bound fractions `0.2452/0.0308/0.4890`; near-best ranges `R0=[1.99405e−5,5.52488e−5]`, `alpha_p=[68755.0,68755.0]`, `d_n=[2.64855e−5,2.65688e−5]`.
- 10 °C: singular values 0.131052, 0.0132659, 0.00809625; sensitivity correlations `+0.9805/−0.9605/−0.9304`; log-bound fractions `0.1800/0.1782/0.4850`; near-best ranges `R0=[0.0222177,0.0253352]`, `alpha_p=[8873.00,12645.8]`, `d_n=[4.02533e−5,4.15603e−5]`.
- 25 °C: singular values 0.103019, 0.00532669, 0.00297702; sensitivity correlations `+0.9798/−0.9886/−0.9636`; log-bound fractions `0.1940/0.2795/0.4462`; near-best ranges `R0=[0.0171349,0.0187427]`, `alpha_p=[3946.48,5802.54]`, `d_n=[7.32346e−5,7.53339e−5]`.

These are practical local diagnostics, not structural-identifiability claims. At 0 °C the condition number exceeds 10⁴ and the near-best R0 spread is large. At 10 and 25 °C, distinct near-best `alpha_p` valleys remain. At −20 and −10 °C, `alpha_p` lies within about 1% of its log-range OCP-feasibility edge.

## Arrhenius diagnostics and leave-one-temperature-out robustness

| Parameter | Temperature definition | Slope [K] | Intercept | Signed Ea [J/mol] | R² | LOTO median abs. relative error | LOTO worst |
|---|---|---:|---:|---:|---:|---:|---:|
| R0 | nominal | 816.09323 | -8.1856904 | 6785.3767 | 0.003542 | 98.242% | 113024.876% |
| R0 | battery_median | 843.18064 | -8.2794464 | 7010.5939 | 0.003892 | 98.214% | 112882.754% |
| alpha_p | nominal | 5563.939 | -10.039865 | 46261.163 | 0.849120 | 63.607% | 298.610% |
| alpha_p | battery_median | 5476.3825 | -9.6845957 | 45533.178 | 0.846714 | 63.899% | 304.486% |
| d_n | nominal | -5505.3397 | 9.2671699 | 45773.942 | 0.922639 | 38.940% | 209.651% |
| d_n | battery_median | -5432.1313 | 8.9646873 | 45165.253 | 0.924589 | 39.083% | 206.719% |

A high R² is not treated as proof of Arrhenius physics. The five-point fits are judged together with local recoverability, boundaries, and leave-one-temperature-out errors.

## Campaign order and aging confounding

Chronological HPPC acquisition order recovered from filenames:
- 2017-03-11: 25 °C — `03-11-17_08.47 25degC_5Pulse_HPPC_Pan18650PF.mat`
- 2017-03-27: 10 °C — `03-27-17_09.06 10degC_5Pulse_HPPC_Pan18650PF.mat`
- 2017-05-20: 0 °C — `05-20-17_10.44 0degC_5pulse_HPPC_Pan18650PF.mat`
- 2017-06-07: -10 °C — `06-07-17_12.29 n10degC_5pulse_HPPC_Pan18650PF.mat`
- 2017-06-15: -20 °C — `06-15-17_11.31 n20degC_5Pulse_HPPC_Pan18650PF.mat`

Temperature decreases as the campaign progresses from 25 °C toward −20 °C. The dataset documentation explicitly warns that the cell aged during the campaign. Therefore temperature and campaign progression are confounded; no recovered trend is attributed exclusively to temperature, and no aging correction is attempted.

## Required decisions

- `R0(T)`: **PARAMETER NOT SUFFICIENTLY RECOVERABLE**
- `alpha_p(T)`: **PARAMETER NOT SUFFICIENTLY RECOVERABLE**
- `d_n(T)`: **PARAMETER NOT SUFFICIENTLY RECOVERABLE**

`R0` has essentially no Arrhenius-line support and becomes practically degenerate at 0 °C. `alpha_p` has moderate in-sample linearity but poor leave-one-temperature-out transfer and near-edge/multivalley solutions. `d_n` has the strongest in-sample linear appearance, but its leave-one-out worst error exceeds 200% and it remains embedded in the ill-conditioned, highly correlated three-parameter fit. Consequently none of the three independently recovered sequences is sufficiently robust to implement as a temperature law.

**NO — KEEP THE CURRENT TEMPERATURE-INDEPENDENT PARAMETERS**

No temperature-dependent law was implemented after this decision.
