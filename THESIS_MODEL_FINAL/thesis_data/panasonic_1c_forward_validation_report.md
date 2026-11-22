# Panasonic 25 C 1C forward-only validation

## 1C TRANSFER LIMITED

The manifest-selected initial-series 1C discharge was evaluated once with the frozen conditional HPPC parameterization. No parameter was optimized, and no OCP, stoichiometric window, capacity, support interval, or model equation was changed.

## Data audit and preparation

- Exact raw source: `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\C20 OCV and 1C discharge tests_start_of_tests\03-09-17_17.59 3349_Dis1C_1.mat`
- Manifest entry verified: `C:\Users\amvel\Documents\Andes\THESIS_MODEL_FINAL\experimental_data_manifest.md`
- Raw columns: `TimeStamp, Voltage, Current, Ah, Wh, Power, Battery_Temp_degC, Time, Chamber_Temp_degC`
- Active branch: source indices 0-348 inclusive (349 samples, 3474.369004 s).
- Sample interval on the active branch: median 9.999996051 s, mean 9.983818978 s, range 4.367005080-10.010999814 s. Actual nonuniform times were preserved.
- Raw current is negative on discharge (-2.89982 to -2.89900 A). Its sign was changed only at the model boundary so model current is discharge-positive.
- Measured discharge capacity from zero-order-hold current/time integration: 2.798236300099 Ah. The relative exported `Ah` counter changes by 2.798180000000 Ah; its arbitrary absolute offset was not interpreted as prior discharge.
- Measured loaded voltage: 4.04420 V at the first active sample and 2.49948 V at cutoff.
- Battery case temperature: 24.98062-32.72549 degC. Chamber temperature: 25-26 degC.
- Relative cumulative exported Ah defined the SOC axis and initialized the discharge at SOC=1. The frozen grouped-state equations were then propagated with the measured current and actual time increments; measured battery temperature was prescribed only in the existing kinetic-voltage factor.

## Frozen forward evaluation

- `R0 = 0.0259594900586 ohm`
- `alpha_p = 3204.22765511 s`; `alpha_n = 3000 s`
- `d_n = 0.000183029692927 1/s`; `d_p = 0.0007 1/s`
- `Q_cell = 2.994979138519 Ah`; `Q_n = 23008.732269258333 C`; `Q_p = 14780.898652565113 C`
- Frozen nominal effective electrode windows and frozen Guo-White MCMB/13-knot positive OCP were used unchanged.

## Electrical metrics

Residual is defined as predicted minus measured voltage.

| Metric | Value |
|---|---:|
| Samples | 349 |
| RMSE [mV] | 114.362188 |
| MAE [mV] | 79.896393 |
| R2 | 0.855578032 |
| Maximum absolute error [mV] | 625.652809 |
| Mean signed error, predicted - measured [mV] | 79.896393 |

### By cell-SOC region

| SOC region | Samples | RMSE [mV] | MAE [mV] | R2 | Max abs [mV] | Mean signed [mV] |
|---|---:|---:|---:|---:|---:|---:|
| high [2/3, 1] | 124 | 37.383604 | 35.985907 | 0.872383944 | 52.194829 | 35.985907 |
| middle [1/3, 2/3) | 124 | 64.318954 | 62.553665 | 0.397326292 | 88.696767 | 62.553665 |
| low [0, 1/3) | 101 | 195.954040 | 155.098359 | 0.079637880 | 625.652809 | 155.098359 |

### Final cutoff region

The cutoff-region diagnostic is the last 1% of the recorded discharge capacity; it is reported separately and was not given special weight or forced to match.

| Metric | Value |
|---|---:|
| Samples | 5 |
| RMSE [mV] | 557.090944 |
| MAE [mV] | 554.313125 |
| R2 | -84.506916670 |
| Maximum absolute error [mV] | 625.652809 |
| Mean signed error, predicted - measured [mV] | 554.313125 |
| Final measured voltage [V] | 2.499480000 |
| Final predicted voltage [V] | 3.125132809 |
| Final-point error, predicted - measured [mV] | 625.652809 |

## State and domain audit

| Quantity | Minimum | Maximum | Frozen support |
|---|---:|---:|---:|
| `x_n_avg` | 0.282632548163 | 0.720451136448 | diagnostic state |
| `x_p_avg` | 0.224027933455 | 0.905559612762 | diagnostic state |
| `x_n_surf` | 0.257430419583 | 0.716850242708 | [0.01, 0.80] |
| `x_p_surf` | 0.230014864650 | 0.947461285261 | [0.18, 0.96] |

- OCP-domain violations: **0**.
- Clipping used: **NO**.
- Final source-Ah cell SOC: 0.065709685917.
- Predicted final cell SOC from measured current/time integration: 0.065690887756.
- Predicted final states: `x_n_avg=0.282632548163`, `x_p_avg=0.905559612762`, `x_n_surf=0.257430419583`, `x_p_surf=0.947461285261`.

## Interpretation

The high-SOC 1C RMSE is 37.384 mV, close in scale to the HPPC calibration-region phase-balanced RMSE of about 36.835 mV. Error then rises to 64.319 mV at middle SOC and 195.954 mV at low SOC. The residual is predominantly positive, meaning the model holds voltage above the measurement, and expands steeply near cutoff (final-point error 625.653 mV).

This independent continuous-discharge result therefore shows the same low-SOC transfer deterioration previously observed in withheld low-SOC HPPC pulses. The discrepancy is not an OCP-domain or clipping artifact: all surface stoichiometries remain supported. It is a systematic loaded-voltage/low-SOC transfer limitation of the frozen conditional representation, with the largest mismatch in the cutoff region.

**Final classification: 1C TRANSFER LIMITED.**

No repair, refit, parameter release, OCP/support/window change, model-equation change, UDDS run, or thermal fit was performed.
