# Panasonic joint 1C + HPPC one-node thermal calibration

## JOINT CALIBRATION IMPROVES THERMAL TRANSFER

Only `Cth` and `Gth` were fitted. The clean one-node architecture, frozen electrical/OCP/state mapping, capacity, entropy coefficient, and measured-voltage heat formulation were unchanged. UDDS was not used in fitting.

## HPPC thermal-continuity gate

- Raw MAT provenance: `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\5 pulse disch\03-11-17_08.47 25degC_5Pulse_HPPC_Pan18650PF.mat`
- Read-only electrical table: `C:\Users\amvel\Documents\Andes\Tesis\datasets\converted_csv\hppc_25C_converted.csv`
- Source rows: 102,800; audited raw continuous blocks: 14; retained entropy-supported blocks: **11**; retained samples: 81,637.
- Normal pulse median dt: 0.100995228 s; normal rest median dt: 0.996997207 s; retained duplicate-time increments: 120 (retained as zero-time observations, with no propagation).
- Thirteen omitted-data gaps exceed 10 s and span 1947.844997-3748.544997 s. No current, heat, ambient, or temperature was held across them.
- Global cumulative-Ah SOC is preserved without block resets. Retained SOC range is 0.117399528430-0.900000000000. The upper 0.90 boundary was interpolated inside block 3; the lower 0.10 boundary lies inside an omitted gap, so no interpolation across that gap was attempted.
- Battery and chamber temperature channels are finite, chamber temperature is available for every retained sample, and each retained block contains real pulse/rest thermal history. **Gate decision: HPPC SUITABLE FOR THERMAL CALIBRATION.**
- Each block is initialized with its measured starting battery temperature. This is a measured initial condition, not a fitted offset.

| Block | Source indices | SOC start | SOC end | Duration [s] | Samples | Median positive dt [s] | Max normal dt [s] | T start [degC] | T end [degC] | Chamber [degC] | Internal current [A] |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 3 | 17281-22904 | 0.900000000 | 0.866736968 | 3693.387328 | 5,625 | 0.996001 | 1.012004 | 25.630700 | 27.502470 | 25.000-25.000 | -0.00000 to 17.40053 |
| 4 | 22905-30539 | 0.806342558 | 0.769844808 | 4920.071996 | 7,635 | 0.996001 | 1.013000 | 26.235940 | 27.267100 | 25.000-25.000 | -0.00000 to 17.40053 |
| 5 | 30540-38174 | 0.709513836 | 0.673012747 | 4920.065005 | 7,635 | 0.996001 | 1.015002 | 26.034190 | 26.886020 | 25.000-25.000 | -0.00000 to 17.40053 |
| 6 | 38175-45809 | 0.612678437 | 0.576214077 | 4920.059000 | 7,635 | 0.996001 | 1.014006 | 25.630700 | 27.289510 | 25.000-25.000 | -0.00000 to 17.40053 |
| 7 | 45810-53444 | 0.515849716 | 0.479528929 | 4920.091006 | 7,635 | 0.996001 | 1.015002 | 25.630700 | 27.087770 | 25.000-25.000 | -0.00000 to 17.40298 |
| 8 | 53445-61079 | 0.419020995 | 0.382546617 | 4920.070004 | 7,635 | 0.996001 | 1.013000 | 25.630700 | 27.098970 | 25.000-25.000 | -0.00000 to 17.40053 |
| 9 | 61080-68714 | 0.322198952 | 0.285861470 | 4920.075003 | 7,635 | 0.996001 | 1.012004 | 25.630700 | 27.502470 | 25.000-25.000 | -0.00000 to 17.40135 |
| 10 | 68715-76349 | 0.273784591 | 0.237280163 | 4920.077005 | 7,635 | 0.996001 | 1.015997 | 25.630700 | 27.513680 | 25.000-25.000 | -0.00000 to 17.40053 |
| 11 | 76350-83984 | 0.225363553 | 0.188869141 | 4920.060006 | 7,635 | 0.996001 | 1.013000 | 25.641910 | 27.928380 | 25.000-25.000 | -0.00000 to 17.40053 |
| 12 | 83985-93209 | 0.176952531 | 0.155466572 | 6050.921997 | 9,225 | 0.996001 | 1.013996 | 25.630700 | 25.832450 | 25.000-25.000 | -0.00000 to 17.40053 |
| 13 | 93210-98916 | 0.128534831 | 0.117399528 | 3701.628997 | 5,707 | 0.996001 | 1.013996 | 25.630700 | 25.843660 | 25.000-25.000 | -0.00000 to 11.60008 |

Omitted-gap source-index pairs and durations are recorded individually in the metrics CSV. The physical heat input is `I*(Ueq_average-V_measured) - I*T_measured*dUeq/dT`; there is no absolute value, negative-heat clipping, automatic sign selection, empirical multiplier, regime weight, or offset.

## Joint calibration

The objective is `0.5*MSE_1C + 0.5*MSE_HPPC = 0.163745107442 K^2`. Forty deterministic log-space multistarts used only the specified numerical guards. 40 succeeded and 40 are within 0.01% joint-RMSE of the best.

- `Cth = 200.000000000000 J/K` (+70.179196% from 1C-only)
- `Gth = 0.078451554871 W/K` (-28.100035% from 1C-only)
- `tau_th = 2549.343991055122 s`

| Experiment | RMSE [K] | MAE [K] | Bias [K] | Max abs [K] | R2 |
|---|---:|---:|---:|---:|---:|
| 1C, joint fit | 0.254140709 | 0.202832997 | -0.032085637 | 0.601948511 | 0.889476100 |
| HPPC, joint fit | 0.512740397 | 0.420219534 | -0.417506303 | 1.874734903 | -2.172877629 |
| HPPC, 1C-only parameters | 0.617273179 | 0.547258199 | -0.545426806 | 1.610352765 | -3.598466911 |

### HPPC diagnostics (not reweighted)

| Subset | RMSE [K] | MAE [K] | Bias [K] | Max abs [K] | R2 |
|---|---:|---:|---:|---:|---:|
| load | 0.393602245 | 0.325547068 | -0.320845745 | 0.953480825 | -2.580046772 |
| rest | 0.519723918 | 0.426541176 | -0.423960698 | 1.874734903 | -2.169888572 |
| high | 0.458631340 | 0.366752461 | -0.359174452 | 1.657544304 | -1.929968332 |
| middle | 0.497547380 | 0.423820775 | -0.423381870 | 1.553304907 | -3.306764411 |
| low | 0.548835285 | 0.447565996 | -0.446162491 | 1.874734903 | -2.055687342 |

## Conditioning and multistart stability

**JOINT THERMAL PARAMETERIZATION BOUNDARY LIMITED** (practical local audit, not structural identifiability). The balanced numerical sensitivity matrix has rank **2/2**, singular values `1.91498945888, 0.390994091315`, condition number `4.89774526372`, sensitivity correlation `+0.949809434741`, and covariance-derived parameter correlation `-0.919550796871`.

Nearest log-bound fractions are 0.000000 (`Cth`) and 0.276350 (`Gth`). Near-optimal ranges are `Cth=199.999999999-200` and `Gth=0.078451548863-0.0784516201987`.

The previous 1C-only condition number was 6.120213126, with sensitivity correlation +0.304113503 and parameter correlation -0.866785425.

## Critical no-refit UDDS validation

Source: `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\Drive cycles\03-21-17_00.29 25degC_UDDS_Pan18650PF.mat`. The exact same entropy-supported SOC 0.90-to-0.10 measured-voltage heat trajectory was evaluated after freezing the joint parameters.

| Metric | 1C-only | Joint 1C+HPPC |
|---|---:|---:|
| Cth [J/K] | 117.523178000996 | 200.000000000000 |
| Gth [W/K] | 0.109112089746 | 0.078451554871 |
| tau_th [s] | 1077.086675492568 | 2549.343991055122 |
| 1C RMSE [K] | 0.143622065 | 0.254140709 |
| HPPC RMSE [K] | 0.617273179 | 0.512740397 |
| UDDS RMSE [K] | 0.526697652 | 0.378519120 |
| UDDS MAE [K] | 0.512960900 | 0.353605705 |
| UDDS bias [K] | -0.512714512 | -0.353251679 |
| UDDS max error [K] | 0.875239690 | 0.822319237 |
| UDDS R2 | -4.119410369 | -1.644069235 |
| Final measured [degC] | 26.469807973 | 26.469807973 |
| Final predicted [degC] | 26.112601243 | 26.319083223 |
| Measured rise [K] | +0.636968211 | +0.636968211 |
| Predicted rise [K] | +0.279761481 | +0.486243461 |

UDDS RMSE changes by +28.133509% relative to 0.526697652 K. Bias changes by +0.159462833 K from the prior -0.512714512 K.

## Decisions

1. **Does adding HPPC materially change Cth or Gth?** **Yes**. Changes are +70.179% and -28.100%, respectively.

2. **Does HPPC improve practical conditioning?** **No**. The condition number changes from 6.120 to 4.898, but the joint optimum is boundary-limited and the absolute parameter correlation worsens from 0.867 to 0.920; correlations and multistart spread are reported above.

3. **Does joint calibration improve independent UDDS transfer?** **Yes**. The no-refit RMSE change is +28.134%.

4. **Is the improvement large enough to replace the 1C-only calibration?** **No**. Replacement requires a clear independent improvement without pathological/boundary conditioning; the frozen baseline file was not overwritten.

5. **Does a persistent approximately constant UDDS underprediction remain?** **Yes**. Joint UDDS bias is -0.353252 K (93.3% of RMSE magnitude).

**Final classification: JOINT CALIBRATION IMPROVES THERMAL TRANSFER.**

No UDDS/CC-CV fitting, two-node state, empirical correction, ambient offset, electrical refit, or model-voltage heat calibration was performed.
