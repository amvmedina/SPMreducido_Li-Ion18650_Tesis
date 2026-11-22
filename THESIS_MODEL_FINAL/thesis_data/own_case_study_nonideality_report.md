# Own TRC18650 capacity-coordinate case-study reference

This reference is constructed directly from measured Chroma signals. No SPM or thermal model was run, no parameter was fitted, no voltage was corrected, and the loaded 1C curve is **not OCV**.

## Primary reference and coordinate

`cycle_cell1_ambient` is used as the clean relative reference because the completed audit placed it alone in the lowest apparent-resistance group. Its measured complete-discharge capacity is `2.899900000 Ah`; this is used as `Q_ref_Ah`. The coordinate is explicitly **coulomb-counted reference SOC**, with `SOC_ref_dis = 1 - discharged_Ah/Q_ref_Ah`.
The post-discharge recharge capacity is `2.899900000 Ah`. Raw discharge points and the comparison-only uniform `0.001` SOC interpolation are stored separately by `record_type`; interpolation does not replace or alter the raw observations.

Its median apparent rest-to-load series resistance is `0.035875537 ohm`. This is a relative setup diagnostic, not intrinsic cell resistance or model R0.

## Signal-derived ambient-cycle segmentation

The file does not start with a rest: it begins with a partial conditioning charge. CC/CV boundaries use the upper-voltage band together with declining positive current; discharge/rest identity uses step number and measured current sign.

| Segment | Step | Rows (zero-based, inclusive) | Samples | Current [A] | Voltage [V] | Step capacity [Ah] | Signed total-capacity change [Ah] | Inferred/provisional duration [s] | Temperature ranges [degC] |
|---|---:|---|---:|---|---|---|---:|---:|---|
| initial_conditioning_charge_cc | 1 | 0-290 | 291 | 2.8497 to 2.9034 | 3.9298 to 4.2038 | 0.0001 to 0.6994 | 0.6993 | 868.46 | Temp1 24.559-29.568; Temp2 24.197-28.991; Temp7 22.515-23.344 |
| initial_conditioning_charge_cv | 1 | 291-956 | 666 | 0.0475 to 2.8340 | 4.2017 to 4.2044 | 0.7017 to 1.1521 | 0.4504 | 1991.48 | Temp1 24.609-29.942; Temp2 24.302-29.398; Temp7 22.490-23.253 |
| pre_discharge_rest | 2 | 957-1157 | 201 | 0.0000 to 0.0000 | 4.1908 to 4.2010 | 0.0000 to 0.0000 | 0.0000 | 598.94 | Temp1 23.843-24.676; Temp2 23.544-24.395; Temp7 22.634-23.492 |
| main_1c_discharge | 3 | 1158-2361 | 1204 | -2.9017 to -2.8984 | 2.8342 to 4.0865 | -2.9000 to -0.0001 | -2.8999 | 3602.63 | Temp1 23.892-33.048; Temp2 23.638-32.198; Temp7 22.490-23.447 |
| post_discharge_rest | 4 | 2362-2562 | 201 | 0.0000 to 0.0000 | 2.9321 to 3.0754 | 0.0000 to 0.0000 | 0.0000 | 598.94 | Temp1 27.962-33.101; Temp2 27.427-32.198; Temp7 22.518-23.253 |
| recharge_cc | 5 | 2563-3603 | 1041 | 2.8497 to 2.9041 | 3.1779 to 4.2038 | 0.0001 to 2.5084 | 2.5082 | 3114.49 | Temp1 27.697-31.818; Temp2 27.243-31.071; Temp7 22.511-23.492 |
| recharge_cv | 5 | 3604-4039 | 436 | 0.1299 to 2.8438 | 4.2016 to 4.2044 | 2.5108 to 2.9000 | 0.3892 | 1302.70 | Temp1 26.359-31.884; Temp2 25.956-31.191; Temp7 22.684-23.447 |
| final_rest | 6 | 4040-4240 | 201 | 0.0000 to 0.0000 | 4.1805 to 4.1981 | 0.0000 to 0.0000 | 0.0000 | 598.94 | Temp1 24.725-26.442; Temp2 24.402-26.043; Temp7 22.624-23.015 |

All durations above are inferred/provisional. No raw time or timestamp exists, so no duration is presented as measured.

## Charge-discharge directional voltage separation

The normalized charge and discharge curves share the interval `SOC_ref = 0.000 to 1.000`. `DeltaV_dir = V_charge - V_discharge` has mean `353.182 mV`, median `332.893 mV`, RMS `358.403 mV`, minimum `116.800 mV`, maximum `495.833 mV`, P5 `294.800 mV`, and P95 `484.693 mV`.

| SOC region | Mean [mV] | Median [mV] | RMS [mV] | P5/P95 [mV] |
|---|---:|---:|---:|---:|
| high | 312.517 | 324.019 | 314.241 | 242.161 / 334.906 |
| middle | 323.628 | 315.510 | 323.977 | 308.683 / 350.291 |
| low | 423.312 | 425.890 | 426.120 | 354.707 / 489.885 |

This is **charge-discharge directional voltage separation**, not thermodynamic hysteresis. It can include ohmic, kinetic, diffusion, mapping, and true hysteretic contributions.

## Descriptive temperature evidence

During the clean discharge, Temp1 starts at `23.921 degC`, reaches `33.048 degC`, and ends at `33.046 degC`; its maximum rise above the discharge-start value is `9.128 degC`. Temp2 starts at `23.735 degC`, reaches `32.198 degC`, and ends at `32.198 degC`.
At the end of the post-discharge rest, the two surface sensors read `27.962` and `27.530 degC`. During discharge their mean signed difference (Temp1 - Temp2) is `0.567 degC`, and maximum absolute difference is `1.038 degC`.
The chamber/ambient range over the complete file is `Temp7: 22.490-23.492 degC`. Ten nonfinite ambient samples identified by the prior audit remain unmodified.

These are capacity-coordinate descriptions only. Thermal time constants, heat-transfer coefficients, and dynamic thermal errors cannot be inferred without measured timestamps.

## Complete-cycle comparison with clean reference

`DeltaV_reference = V_affected - V_clean_reference` at matched normalized coulomb-counted reference SOC. The deviation combines temperature, cell, setup, degradation/state, and electrochemical differences; it is not assigned solely to parasitic resistance.

| Dataset | Capacity [Ah] | Median R_app [ohm] | R_app ratio | Mean / median deviation [mV] | RMS / max abs [mV] | P5 / P95 [mV] | High/middle/low mean [mV] | Interpretation |
|---|---:|---:|---:|---:|---:|---:|---|---|
| cycle_cell1_ambient | 2.8999 | 0.035876 | 1.000x | 0.00 / 0.00 | 0.00 / 0.00 | 0.00 / 0.00 | 0.00 / 0.00 / 0.00 | CLEAN REFERENCE |
| cycle_cell1_30c | 2.4894 | 0.232494 | 6.481x | -428.60 / -437.17 | 433.98 / 601.50 | -522.45 / -320.89 | -498.79 / -434.77 / -352.28 | ELEVATED RESPONSE + TEMPERATURE DIFFERENCE |
| cycle_cell1_35c | 2.6171 | 0.222850 | 6.212x | -427.70 / -445.59 | 432.50 / 590.00 | -502.91 / -321.70 | -487.50 / -443.57 / -352.07 | ELEVATED RESPONSE + TEMPERATURE DIFFERENCE |
| cycle_cell2_undeclared | 2.7133 | 0.177540 | 4.949x | -370.59 / -387.97 | 374.96 / 452.46 | -441.68 / -270.44 | -426.00 / -384.14 / -301.67 | ELEVATED RESPONSE + CELL-TO-CELL DIFFERENCE |
| cycle_cell2_30c | 2.5994 | 0.215820 | 6.016x | -425.78 / -446.43 | 431.51 / 555.90 | -504.65 / -303.13 | -491.63 / -443.22 / -342.55 | ELEVATED RESPONSE + TEMPERATURE DIFFERENCE + CELL-TO-CELL DIFFERENCE |

## Repeated-transition pulse evidence

| Dataset | N | Median [ohm] | IQR [ohm] | Range [ohm] | Ratio to clean | Detector |
|---|---:|---:|---:|---|---:|---|
| pulse_cell1_part1 | 7 | 0.221348508 | 0.007900583 | 0.218786625-0.235910878 | 6.170x | EVIDENCE OF ELEVATED SERIES-PATH RESPONSE |
| pulse_cell1_part2 | 3 | 0.228938929 | 0.009566190 | 0.218013729-0.237146108 | 6.381x | EVIDENCE OF ELEVATED SERIES-PATH RESPONSE |

All ten pulse transitions lie in the elevated group. Part-1 and part-2 medians differ by only `7.590 milliohm`, so the elevation is **persistent across repeated transitions**, not isolated. Its physical source is not uniquely identified.

## Nonideality detector provenance

The detector reuses the audit's largest-adjacent-log-gap grouping of experiment-level R_app medians. It does not impose a universal resistance threshold. Each output row retains transition count, reference/experiment median and IQR, ratio, reference identity, and an evidence statement in the metrics CSV.

## Explicit answers

1. **Coherent clean reference?** Yes, as a relative experimental reference in capacity/SOC coordinates: it has a complete 2.8999 Ah discharge/recharge and the lowest transition-response group. It is not asserted to be ideal or OCV.
2. **Directional separation magnitude?** Mean `353.182 mV`, RMS `358.403 mV`, with the full distribution and SOC regions reported above.
3. **Affected-cycle deviation strength?** All four affected cycles show systematic capacity-matched voltage deviations; their numerical mean/RMS/maximum and regional values are tabulated above and exported without correction.
4. **Persistent pulse elevation?** Yes. Both pulse parts and every retained transition occupy the elevated relative-response group.
5. **Useful nonideality case studies?** Yes. They remain valid measurements of the complete experimental chain and are useful precisely because the audit detects a reproducible nonideality; they are not treated as calibration ground truth or discarded as failed measurements.
6. **Quantitative versus qualitative conclusions?** Capacity, normalized SOC coordinates, measured voltage/temperature values, voltage deviations, and adjacent-sample R_app are quantitative. Any time-dependent electrical or thermal inference, thermal constants, dynamic RMSE, or unique attribution of the excess resistance remains qualitative/not identifiable because measured time is absent.

## Stop point

No electrochemical or thermal model was run; no OCP, R0, diffusion, or thermal parameter was fitted; and no voltage correction was applied.
