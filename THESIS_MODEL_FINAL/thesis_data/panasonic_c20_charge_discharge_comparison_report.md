# Panasonic C/20 charge/discharge quasi-OCV comparison

## Scope and source

This is a data-only comparison. No electrochemical model or OCP function was run, no parameter was fitted, and no hysteresis term or branch correction was introduced.

- Source: `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\C20 OCV and 1C discharge tests_start_of_tests\05-08-17_13.26 C20 OCV Test_C20_25dC.mat`
- MAT fields: `TimeStamp, Voltage, Current, Ah, Wh, Power, Battery_Temp_degC, Time, Chamber_Temp_degC`
- Samples: 2453
- Frozen common capacity: 2.994979138519 Ah
- Sign observed from the signals: negative current is discharge; positive current is charge.
- Active/rest threshold used only for segmentation: |I| > 0.050 A.
- Source duration: 195824.477005 s (54.395688 h).
- Battery-temperature range: 11.41626 to 26.09024 °C.
- Chamber-temperature range: 10.00000 to 25.00000 °C.

## Reproducible signal segmentation

The script classifies every sample from measured current and run-length encodes the result. It verifies the five-part order rather than using the expected indices as inputs. Indices below are zero-based and inclusive.

| Segment | Indices | Samples | Time span (s) | Current range (A) | Voltage range (V) |
|---|---:|---:|---:|---:|---:|
| rest | 0–5 | 6 | 0.000000–240.010004 | 0.00000 to 0.00000 | 4.18398 to 4.18398 |
| C/20 discharge | 6–1246 | 1241 | 300.019002–74680.886005 | -0.14536 to -0.14454 | 2.49948 to 4.17030 |
| rest | 1247–1307 | 61 | 74740.900003–78280.903000 | 0.00000 to 0.00000 | 2.66300 to 2.86117 |
| C/20 charge | 1308–2390 | 1083 | 78340.916002–143255.048002 | 0.14455 to 0.14537 | 2.92679 to 4.20007 |
| rest | 2391–2452 | 62 | 143315.059998–195824.477005 | 0.00000 to 0.00000 | 4.15953 to 4.18591 |

The detected active branches are exactly discharge indices 6–1246 and charge indices 1308–2390.

## Common-capacity SOC comparison

The primary axes are

`SOC_dis = 1 - Q_discharge_since_start / Q_cell_Ah`

and

`SOC_chg = Q_charge_since_charge_start / Q_cell_Ah`.

The charge starts at nominal SOC = 0 because it follows the complete low-rate discharge. Its endpoint was **not** forced to SOC = 1. Both measured curves were linearly interpolated on 2001 points only over the exact common interval **[0, 0.872766248937]**. No extrapolation was used.

For `deltaV_qocv = V_C20_charge - V_C20_discharge`:

| Statistic | Gap (mV) |
|---|---:|
| Mean | 105.842 |
| Median | 109.944 |
| RMS | 111.723 |
| Minimum | 64.941 |
| Maximum | 427.310 |
| P5 | 66.231 |
| P95 | 156.666 |

### Gap by available SOC region

| Region | SOC interval | Grid points | Mean (mV) | Median (mV) | RMS (mV) | Min (mV) | Max (mV) | P5 (mV) | P95 (mV) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Low | [0, 1/3) | 764 | 90.543 | 77.180 | 99.472 | 64.941 | 427.310 | 65.604 | 170.131 |
| Middle | [1/3, 2/3) | 764 | 98.829 | 110.577 | 100.832 | 66.211 | 123.455 | 66.898 | 122.324 |
| High | [2/3, 0.872766249] | 473 | 141.882 | 149.319 | 142.831 | 111.887 | 174.807 | 113.384 | 163.409 |

The gap is positive at every comparison point and is strongly SOC-dependent. It has a pronounced endpoint excursion near nominal SOC = 0; that point compares two loaded measurements separated by the intermediate rest and must not be interpreted as an equilibrium material property.

## Endpoint and capacity diagnostic

| Quantity | Value |
|---|---:|
| Integrated discharge capacity | 2.994979138519 Ah |
| Integrated charge capacity | 2.613916708369 Ah |
| Charge/discharge capacity ratio | 0.872766249 (87.276625%) |
| Final charge SOC from frozen capacity | 0.872766248937 (87.276625%) |
| Discharge final loaded voltage | 2.49948 V |
| Discharge first post-current rest voltage | 2.66300 V |
| Discharge recovered voltage before charge | 2.86117 V |
| Charge initial loaded voltage | 2.92679 V |
| Charge final loaded voltage | 4.20007 V |
| Charge first post-current rest voltage | 4.18591 V |
| Charge final recorded rest voltage | 4.15953 V |

The 4.2 V loaded endpoint does not make the charge branch complete: coulomb counting reaches only 87.277% of the frozen capacity. The integrated discharge capacity differs from the frozen value by 0.000 nAh.

## SHAPE-ONLY NORMALIZED COMPARISON

As a secondary diagnostic only, each branch was mapped by its own transferred capacity. This is not the primary SOC comparison and is not used to redefine capacity. Its gap statistics are: mean 47.045 mV, median 39.203 mV, RMS 55.988 mV, minimum 29.770 mV, maximum 427.310 mV, P5 31.287 mV, and P95 109.597 mV. The lower central gap after this remapping shows that unequal effective capacity mapping contributes to the primary difference, but it does not remove the directional voltage separation.

## Relation to the existing CC–CV charge validation

1. **Systematic gap:** Yes. The measured low-rate charge voltage exceeds the discharge voltage over the entire shared SOC interval (minimum 64.941 mV; mean 105.842 mV).
2. **Comparable magnitude:** Yes. The mean measured gap is the same order as, and larger than, the prior CC bias magnitude (50.2 mV) and CV bias magnitude (68.0 mV: mean-gap ratios 2.11 and 1.56, respectively).
3. **SOC dependence:** Strong. Regional mean gaps are 90.543 mV (low), 98.829 mV (middle), and 141.882 mV (high), with the largest endpoint difference near nominal SOC = 0.
4. **OCP implication:** The positive, sustained low-rate charge/discharge quasi-OCV gap is consistent with a discharge-derived OCP producing a substantial part of the observed negative charge-direction voltage bias. It **could explain** a substantial part, but this diagnostic does not prove causation or uniquely allocate the error: both branches remain under approximately C/20 load, the endpoint protocols include rests, and the charge transfers less capacity.

This result is therefore described only as a **low-rate charge/discharge quasi-OCV gap**, not thermodynamic hysteresis.
