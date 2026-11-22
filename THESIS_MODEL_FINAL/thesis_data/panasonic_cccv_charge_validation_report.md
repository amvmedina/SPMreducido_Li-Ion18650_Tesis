# Panasonic 25 C CC-CV charge forward validation

## CC-CV CHARGE TRANSFER LIMITED

The manifest-selected associated charge was run forward once with the frozen discharge-side OCP representation and frozen conditional dynamics. Measured current was the input during both CC and CV; predicted voltage was never forced to the regulated CV voltage. No optimization, hysteresis, charge-specific OCP, baseline correction, clipping, or thermal fit was introduced.

## Dataset audit

- Exact raw source: `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\C20 OCV and 1C discharge tests_start_of_tests\03-09-17_19.03 3349_Charge2.mat`
- Exact converted source: `C:\Users\amvel\Documents\Andes\Tesis\datasets\converted_csv\charge_3349_2_converted.csv` (numerically verified against the raw channels)
- Raw columns: `TimeStamp, Voltage, Current, Ah, Wh, Power, Battery_Temp_degC, Time, Chamber_Temp_degC`
- Converted columns: `Time, Voltage, Current, Ah, Wh, Power, Battery_Temp_degC, Chamber_Temp_degC`
- Raw samples: 123; usable samples: 122 after removing one zero-time duplicate final electrical record. Its electrical fields and chamber temperature repeat exactly, while battery temperature changes by -0.20175 degC despite no elapsed time; the earlier record is retained.
- Duration: 7190.123995766 s. Positive sampling interval median 60.000006482 s, mean 59.422512362 s, range 0.010007992-60.008998588 s. The 0.010 s interval is retained because it is positive; all actual times are preserved.
- Raw current is positive on charge and zero at rest: 0.00000 to 2.89997 A. Conversion to model charge-negative current occurs only at the loader/simulation boundary.
- Voltage range: 3.21117-4.20007 V.
- Exported Ah is nondecreasing from 0.00000 to 2.78376 Ah and ends at 2.78376 Ah; no counter reversal occurs.
- Battery temperature: 25.60828-30.01311 degC; chamber temperature: 25-25 degC.

## Signal-derived experiment structure

- Initial rest: indices 0-10, through 540.006003 s (11 records). The first charged record is index 11 at 600.011994 s.
- CC phase: indices 11-58; current 2.89916-2.89997 A, mean 2.899565 A.
- CC-to-CV transition: the first current reduction occurs at index 59, time 3480.010002 s, where measured voltage is 4.20007 V and current is 2.56596 A. The transition metric uses the last CC sample and first two CV samples.
- CV phase: indices 59-111; regulated voltage range 4.19942-4.20007 V and current taper 2.56596 to 0.04982 A.
- Charge termination: last active record index 111 at 6590.110996 s and 0.04982 A; zero current first appears at index 112, 6650.118999 s.
- Final rest: indices 112-121, spanning 540.004997 s.

## Initial state and measured-current handling

The immediately preceding manifest-recorded 1C discharge removes 2.798236300099 Ah and is followed by recovery to 3.20796 V. The charge file begins at 3.21117 V with zero current. With frozen `Q_cell=2.994979138519 Ah`, the non-fitted initial SOC is therefore:

`SOC_initial = 1 - 2.798236300099/2.994979138519 = 0.065690887756`.

The frozen nominal effective trajectories give `x_n_avg=x_n_aux=0.282632548163` and `x_p_avg=x_p_aux=0.905559612762` at the first record. Auxiliary states equal average states because initialization is at the recorded zero-current recovery; the small measured rest drift is retained as validation error rather than used to adjust states.

The source records current at 60 s endpoints while exported Ah contains the charge passed over each interval. Accordingly, sample current drives the algebraic surface/kinetic/ohmic voltage, while `delta Ah/delta t` supplies the measured interval-average current for state propagation. This exactly preserves CC onset and CV taper charge instead of introducing a one-sample delay. It is a data-timing treatment, not a parameter or model-equation change.

Final SOC is `0.995166477151` after adding the measured 2.78376 Ah; neither initial nor final SOC was forced to zero or one.

## Primary voltage metrics

Residual is predicted minus measured voltage.

| Metric | Complete usable charge |
|---|---:|
| Samples | 122 |
| RMSE [mV] | 61.033423 |
| MAE [mV] | 59.550924 |
| R2 | 0.962510341 |
| Mean signed error [mV] | -47.526296 |
| Maximum absolute error [mV] | 96.022350 |
| P95 absolute error [mV] | 81.034958 |

### By experimental phase

| Phase | Samples | RMSE [mV] | MAE [mV] | R2 | Mean signed [mV] | P95 absolute [mV] | Max absolute [mV] |
|---|---:|---:|---:|---:|---:|---:|---:|
| Initial rest | 11 | 66.763466 | 66.682027 | -409.146201567 | 66.682027 | 72.121118 | 73.411118 |
| CC | 48 | 51.902289 | 50.184473 | 0.927905496 | -50.184473 | 69.739039 | 96.022350 |
| CV | 53 | 68.568852 | 68.001967 | -44655.095082357 | -68.001967 | 81.448077 | 82.120183 |
| CC-to-CV transition | 3 | 65.183296 | 65.116535 | -460.013252894 | -65.116535 | 68.270337 | 68.591900 |
| Final 5% taper | 38 | 64.489379 | 64.137839 | -39373.012421448 | -64.137839 | 77.133335 | 78.938169 |
| Final rest | 10 | 51.898469 | 51.875146 | -896.887375962 | -51.875146 | 54.483919 | 54.941880 |

CV voltage is externally regulated; its near-constant measured value is not evidence of model-predicted CV control. The very negative CV R2 reflects division by the tiny regulated-voltage variance and is not used as the transfer criterion.

### By cell SOC

| SOC region | Samples | RMSE [mV] | MAE [mV] | Mean signed [mV] | Max absolute [mV] |
|---|---:|---:|---:|---:|---:|
| low [0, 1/3) | 27 | 56.155984 | 53.445292 | 0.888211 | 96.022350 |
| middle [1/3, 2/3) | 21 | 48.316061 | 47.798658 | -47.798658 | 56.240615 |
| high [2/3, 1] | 74 | 65.789587 | 65.113757 | -65.113757 | 82.120183 |

Signed residual versus SOC has Spearman `rho=-0.479984`. The charge-side residual shifts toward increasing underprediction as SOC rises: mean residual changes from +0.888 mV at low SOC to -47.799 mV at middle SOC and -65.114 mV at high SOC. Low-SOC aggregation contains the positive initial-rest offset and the start of loaded CC, so its near-zero mean does not indicate uniformly small low-SOC error.

## CC current response

The measured CC current is 2.899565 A on average, range 2.89916-2.89997 A. With no correction, CC voltage RMSE is 51.902289 mV and predicted-minus-measured bias is -50.184473 mV. Its residual-versus-SOC trajectory is shown in the SOC figure.

## CV taper response

Using the measured taper current, CV RMSE is 68.568852 mV, MAE 68.001967 mV, bias -68.001967 mV, and maximum absolute error 82.120183 mV. At the final active sample (0.04982 A), residual is -57.417679 mV; the last five active taper records average -57.514070 mV.

## State and domain audit

| State | Minimum | Maximum | Frozen support where applicable |
|---|---:|---:|---:|
| `x_n_avg` | 0.282632548163 | 0.718186140717 | diagnostic state |
| `x_p_avg` | 0.227553746166 | 0.905559612762 | diagnostic state |
| `x_n_surf` | 0.282632548163 | 0.718677138248 | [0.01, 0.80] |
| `x_p_surf` | 0.226731237172 | 0.905559612762 | [0.18, 0.96] |

- Cell SOC range: 0.065690887756-0.995166477151.
- OCP-domain violations: **0**.
- Clipping count: **0**.

## Charge-direction transfer and prior validation

| Validation | Metric | Result [mV] |
|---|---|---:|
| HPPC calibration, SOC >= 1/3 | Phase-balanced pulse-relative RMSE | 36.835198 |
| Withheld low-SOC HPPC | Phase-balanced / absolute RMSE | 85.729 / 59.582 |
| Independent 1C discharge | Absolute RMSE | 114.362188 |
| UDDS | Absolute RMSE | 30.988192 |
| UDDS ten largest current steps | Delta-voltage RMSE | 92.713765 |
| Independent CC-CV charge | Absolute RMSE | 61.033423 |

The independent charge supports a current-direction limitation, but not as a single uniform sign offset. UDDS showed excessive predicted upward voltage changes at several large regenerative onsets. Here, after the initial rest, the sustained CC and CV predictions are predominantly below measurement. Thus the frozen discharge-calibrated representation can overstate an immediate regenerative step while understating sustained charge voltage; that contrast is evidence of direction/transient dependence rather than a correctable constant baseline.

## Interpretation and classification

The complete-charge RMSE is 61.033 mV. The model begins +66.682 mV high during the initial rest, then runs -50.184 mV low through CC and -68.002 mV low through CV. High-SOC bias reaches -65.114 mV, while the final 5% taper remains at 64.489 mV RMSE. These sign and SOC changes rule out a single static correction and show that discharge-side OCP/dynamic transfer is incomplete.

The dominant limitation is a **combination of high-SOC/static OCP mismatch, CC loaded-voltage underprediction, CV taper underprediction, and current-direction asymmetry**. The regulated CV trace is treated only as a measured reference, not as simulated charger control. No repair is attempted.

**Final classification: CC-CV CHARGE TRANSFER LIMITED.**

No optimization, parameter change, hysteresis, separate charge OCP, clipping, model-equation change, or thermal fit was performed.
