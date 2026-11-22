# Panasonic 25 C CC-CV charge thermal validation

## Strict no-refit scope

The frozen primary 1C-only one-node values `Cth=117.523178000996 J/K`, `Gth=0.109112089746 W/K`, and `tau=1077.086675492568 s` were used unchanged. No joint 1C+HPPC parameter, optimizer, offset, correction, or architecture change was used.
Validation B feeds its recursively predicted temperature into both the frozen Butler-Volmer voltage calculation and reversible heat after the shared measured initial condition; measured battery temperature is not used afterward, and no Arrhenius dependence is introduced.

## Dataset and supported interval

- Exact source: `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\C20 OCV and 1C discharge tests_start_of_tests\03-09-17_19.03 3349_Charge2.mat`
- Raw columns: `TimeStamp, Voltage, Current, Ah, Wh, Power, Battery_Temp_degC, Time, Chamber_Temp_degC`
- Raw samples: 123; duplicate rows removed using the previously audited rule: 1.
- Complete charge SOC: 0.065690887756-0.995166477151. Retained entropy-supported SOC: **0.100000000000-0.900000000000**, with both endpoints inserted by exact linear interpolation.
- Complete audited charge throughput: 2.783760000000 Ah; retained supported-interval throughput: 2.395983310815 Ah.
- Source start/end times: 667.590232452-3732.937565401 s; retained duration: 3065.347332949 s; samples: 53.
- Model-current range: -2.899970 to -1.404012 A; CC values are numerically negative after the single loader-boundary sign conversion.
- Measured voltage: 3.543739-4.200070 V; battery temperature: 26.900295-30.013110 degC; chamber temperature: 25.000000-25.000000 degC.
- CC-to-CV transition begins at source time 3480.010001734 s and SOC 0.855652984510. It is inside the supported interval; SOC=0.90 occurs during CV, so no later CV or final-rest data were admitted.

## Overall thermal validation

| Path | RMSE [K] | MAE [K] | Bias [K] | Max abs [K] | R2 | Final predicted [degC] |
|---|---:|---:|---:|---:|---:|---:|
| Validation A: measured voltage/measured T heat | 0.427780949 | 0.348334994 | +0.012743484 | 1.034230957 | 0.693423225 | 30.114973418 |
| Validation B: fully predictive electrothermal | 0.864386997 | 0.829400941 | -0.829400941 | 1.161652108 | -0.251736229 | 28.787310216 |
| Model-voltage / measured-T reversible diagnostic | 0.862878739 | 0.827836673 | -0.827836673 | 1.160671780 | -0.247371762 | 28.789629614 |

Measured initial/final temperatures are 26.900295185/29.080742461 degC. Validation A initial/final predictions are 26.900295185/30.114973418 degC; measured/predicted rises are +2.180447276/+3.214678233 K. Validation B initial/final predictions are 26.900295185/28.787310216 degC; its predicted rise is +1.887015032 K.

## CC/CV diagnostics

| Path/phase | Samples | Duration [s] | SOC range | Model current [A] | RMSE [K] | MAE [K] | Bias [K] | Max abs [K] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| validation_a / CC | 47 | 2752.425768 | 0.100000000-0.840207134 | -2.899970 to -2.899160 | 0.340973234 | 0.288234155 | -0.090198824 | 0.726524925 |
| validation_a / CV | 6 | 252.927564 | 0.855652985-0.900000000 | -2.565960 to -1.404012 | 0.840086065 | 0.819124899 | +0.819124899 | 1.034230957 |
| validation_b / CC | 47 | 2752.425768 | 0.100000000-0.840207134 | -2.899970 to -2.899160 | 0.896342598 | 0.868895583 | -0.868895583 | 1.161652108 |
| validation_b / CV | 6 | 252.927564 | 0.855652985-0.900000000 | -2.565960 to -1.404012 | 0.553553098 | 0.520026248 | -0.520026248 | 0.851510510 |
| model_voltage_measured_t / CC | 47 | 2752.425768 | 0.100000000-0.840207134 | -2.899970 to -2.899160 | 0.894888628 | 0.867447648 | -0.867447648 | 1.160671780 |
| model_voltage_measured_t / CV | 6 | 252.927564 | 0.855652985-0.900000000 | -2.565960 to -1.404012 | 0.551189501 | 0.517550699 | -0.517550699 | 0.848884737 |

## Charge-side heat audit

| Heat path | Eirr [J] | Erev [J] | Etotal [J] | Mean Qirr [W] | Mean Qrev [W] | Mean total [W] |
|---|---:|---:|---:|---:|---:|---:|
| Measured-voltage / measured-T | 1565.334960095 | 213.913534501 | 1779.248494596 | 0.510655006 | 0.069784436 | 0.580439442 |
| Fully predictive | 1128.899067677 | 213.286465915 | 1342.185533592 | 0.368277701 | 0.069579869 | 0.437857570 |

For the primary measured-voltage path, signed `Erev/Etotal=+0.120226902` and integral-absolute reversible fraction is `0.136892833`. Reversible heat is exothermic for 86.298682% and endothermic for 13.701318% of retained time.
The charge sign rule was verified numerically: **True**. With I<0, positive dUeq/dT gives positive/exothermic reversible heat and negative dUeq/dT gives negative/endothermic heat.
Relative to an irreversible-only measured-voltage propagation, reversible heat changes the final predicted temperature by `+0.654719124 K`, with maximum absolute trajectory change `0.774795986 K` and RMS change `0.523330022 K`.

## Electrical-error-to-heat diagnostic

`Qirr_measured-Qirr_model = I*(Vmodel-Vmeasured)` has mean `+0.142377305 W`, RMS `0.145930852 W`, and time integral `+436.435892418 J`. Numerical identity error is `8.327e-17 W`.
The frozen electrical model predominantly underpredicts charge voltage. Because model current is negative, `I*(Vmodel-Vmeasured)` is therefore positive: the measured-voltage path supplies more irreversible heat than the model-voltage path. This is consistent with the fully predictive trajectory being substantially cooler.

## Classifications and thesis interpretation

### Validation A: CHARGE-SIDE THERMAL TRANSFER LIMITED

### Validation B: COUPLED CHARGE ELECTROTHERMAL TRANSFER FAILED

1. **Does the discharge-calibrated thermal pair transfer to charge with measured-voltage heat?** **Only partially.** Validation A has RMSE `0.427781 K` and R2 `0.693423`, but its near-zero overall bias is cancellation between CC underprediction and CV overprediction. It overpredicts the final temperature by `+1.034231 K` and the temperature rise by `+1.034231 K`.
2. **Is failure primarily thermal, electrical, or mixed?** **Mixed overall, with electrical-voltage heat dominant in Validation B.** Validation A already shows charge/CV thermal-transfer limitation. Replacing measured voltage with model voltage changes RMSE to `0.862879 K`; adding predicted-temperature feedback changes it only to `0.864387 K`. Thus nearly all additional coupled failure precedes thermal feedback and follows the `+436.436 J` electrical heat deficit.
3. **Does error differ materially between CC and CV?** **Yes.** Validation A changes from CC RMSE/bias `0.340973/-0.090199 K` to CV `0.840086/+0.819125 K`, including a bias-sign reversal. Validation B remains negatively biased in both phases.
4. **Does reversible heat materially affect prediction?** **Yes, but it is not the dominant total heat source.** It contributes `12.023%` of signed total energy and shifts the final Validation-A prediction by `+0.654719 K` (maximum trajectory effect `0.774796 K`).
5. **Does this support transfer across current direction?** **It limits that claim.** The frozen pair retains useful CC shape information with measured-voltage heat, but the CV residual structure and the failed fully predictive path do not support general charge-direction transfer.

No thermal/electrical parameter, OCP, entropy coefficient, SOC mapping, or frozen model equation was modified.
