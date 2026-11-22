# Universidad de los Andes TRC18650 final case-study assessment

## Methodological outcome

The case study is closed as an evidence-gated dataset assessment, not as a TRC parameter calibration. Data quality, model applicability, parameter transferability, and model accuracy eligibility are kept separate.

**Direct Panasonic-to-TRC transfer is not adequate. Q_cell is directly supported by the clean own experiment. Capacity/state mapping requires cell-specific adaptation. A cell-compatible equilibrium-voltage/OCP representation is required before rigorous TRC reduced-model validation. Available own data do not identify that OCP or the dynamic/thermal parameters uniquely. Therefore, NO FULL TRC PARAMETER REFIT IS PERFORMED.**

This is a methodological decision, not a failed optimization.

## Metric-role register

| Metric | Value | Unit | metric_role | Evidence basis |
|---|---:|---|---|---|
| clean_reference_Q_discharge | 2.899900000 | Ah | DATA-QUALITY DIAGNOSTIC | Complete Chroma discharge counter. |
| clean_reference_apparent_R | 0.035875537 | ohm | DATA-QUALITY DIAGNOSTIC | Adjacent rest-to-load DeltaV/DeltaI; not intrinsic R0. |
| Panasonic_backbone_RMSE | 0.188882141 | V | TRANSFER DIAGNOSTIC | Panasonic Ueq versus loaded TRC 1C curve; not OCP validation. |
| literal_full_transfer_RMSE | 0.082134289 | V | TRANSFER DIAGNOSTIC | Valid frozen-OCP-support samples using provisional time. |
| charge_transfer_RMSE | 0.040091853 | V | TRANSFER DIAGNOSTIC | Held-out directional transfer using provisional time. |

No metric in this final TRC case study is labelled MODEL ACCURACY because the direct model-transfer gate failed. Voltage comparisons remain transfer or data-quality diagnostics.

## A. Clean reference

`cycle_cell1_ambient` is the clean relative own-data reference. Its complete discharge capacity is `2.8999 Ah`; median apparent rest-to-load series resistance is `0.035875537 ohm`. The latter is a **DATA-QUALITY DIAGNOSTIC**, not intrinsic battery R0.

Electrical calibration: **NOT SUITABLE** — lowest relative apparent-resistance group and complete cycle, but no measured time channel and insufficient dynamic excitation. Electrical validation: **SUITABLE WITH CAVEAT** — cleanest current/voltage transitions and complete capacity coverage; timebase remains provisional.

## B. Affected complete-cycle datasets

Voltage deviations are affected-minus-clean at matched normalized reference SOC. They combine temperature, cell, state, electrochemical, and experimental/setup effects; no voltage correction is applied.

| Dataset | Cell / nominal T | Capacity [Ah] | Median R_app [ohm] | Ratio | Mean / RMS voltage deviation [mV] | Nonideality flag | Interpretation | metric_role |
|---|---|---:|---:|---:|---:|---|---|---|
| cycle_cell1_30c | cell 1 / 30 degC | 2.4894 | 0.232494022 | 6.481x | -428.60 / 433.98 | EVIDENCE OF ELEVATED SERIES-PATH RESPONSE | ELEVATED RESPONSE + TEMPERATURE DIFFERENCE | DATA-QUALITY DIAGNOSTIC |
| cycle_cell1_35c | cell 1 / 35 degC | 2.6171 | 0.222850133 | 6.212x | -427.70 / 432.50 | EVIDENCE OF ELEVATED SERIES-PATH RESPONSE | ELEVATED RESPONSE + TEMPERATURE DIFFERENCE | DATA-QUALITY DIAGNOSTIC |
| cycle_cell2_undeclared | cell 2 / undeclared | 2.7133 | 0.177540167 | 4.949x | -370.59 / 374.96 | EVIDENCE OF ELEVATED SERIES-PATH RESPONSE | ELEVATED RESPONSE + CELL-TO-CELL DIFFERENCE | DATA-QUALITY DIAGNOSTIC |
| cycle_cell2_30c | cell 2 / 30 degC | 2.5994 | 0.215820169 | 6.016x | -425.78 / 431.51 | EVIDENCE OF ELEVATED SERIES-PATH RESPONSE | ELEVATED RESPONSE + TEMPERATURE DIFFERENCE + CELL-TO-CELL DIFFERENCE | DATA-QUALITY DIAGNOSTIC |

These experiments provide evidence consistent with an elevated series-path response and experimental/setup nonideality. Fixture resistance was not uniquely identified. The datasets remain valuable because the methodology detects that they should not be treated as calibration ground truth.

## C. Repeated pulse evidence

| Dataset | Transitions | Median [ohm] | IQR [ohm] | Range [ohm] | Ratio | Finding | metric_role |
|---|---:|---:|---:|---|---:|---|---|
| pulse_cell1_part1 | 7 | 0.221348508 | 0.007900583 | 0.218786625-0.235910878 | 6.170x | Persistent elevated group | DATA-QUALITY DIAGNOSTIC |
| pulse_cell1_part2 | 3 | 0.228938929 | 0.009566190 | 0.218013729-0.237146108 | 6.381x | Persistent elevated group | DATA-QUALITY DIAGNOSTIC |

All ten retained pulse transitions remain in the elevated group. The elevation is persistent across both parts, not an isolated transition; its physical source remains non-unique.

## D. Direct Panasonic model transfer

| Diagnostic | Result | metric_role | Interpretation |
|---|---:|---|---|
| Voltage-backbone RMSE | 188.882 mV | TRANSFER DIAGNOSTIC | Panasonic Ueq versus loaded TRC 1C voltage; not OCP validation |
| Voltage-backbone bias | 179.277 mV | TRANSFER DIAGNOSTIC | Frozen backbone is systematically above loaded TRC voltage |
| Literal full-transfer RMSE | 82.134 mV | TRANSFER DIAGNOSTIC | Valid support only; provisional-time integration |
| Literal full-transfer bias | 56.982 mV | TRANSFER DIAGNOSTIC | Predicted minus measured |
| High/middle/low-SOC RMSE | 30.306 / 44.215 / 133.983 mV | TRANSFER DIAGNOSTIC | Deterioration is concentrated at low SOC |
| Surface-state support | 1182 valid; 22 invalid | TRANSFER DIAGNOSTIC | 22 low-SOC positive-surface violations; no clipping |
| Recharge RMSE | 40.092 mV | TRANSFER DIAGNOSTIC | Held-out directional check |
| Recharge CC / CV RMSE | 39.883 / 40.587 mV | TRANSFER DIAGNOSTIC | No direction-specific correction |

These results are not a calibrated TRC model and cannot be reported as model accuracy. The dynamic terms partially compensate the voltage-backbone mismatch, but capacity mapping remains incompatible and the discharge enters unsupported OCP surface states.

## Four separate assessment dimensions

| Dataset | Data quality | Model applicability | Parameter transferability | Model accuracy eligibility | Recommended use |
|---|---|---|---|---|---|
| pulse_cell1_part1 | E-cal NOT SUITABLE; E-val QUALITATIVE / DIAGNOSTIC ONLY; T-cal NOT SUITABLE; T-val QUALITATIVE / DIAGNOSTIC ONLY | DIAGNOSTIC CASE STUDY ONLY AFTER NONIDEALITY GATE | NOT ENTERED: experimental nonideality/data-quality gate | NOT ASSESSED | experimental-nonideality case study; qualitative temperature-response diagnostic |
| pulse_cell1_part2 | E-cal NOT SUITABLE; E-val QUALITATIVE / DIAGNOSTIC ONLY; T-cal NOT SUITABLE; T-val QUALITATIVE / DIAGNOSTIC ONLY | DIAGNOSTIC CASE STUDY ONLY AFTER NONIDEALITY GATE | NOT ENTERED: experimental nonideality/data-quality gate | NOT ASSESSED | experimental-nonideality case study; qualitative temperature-response diagnostic |
| cycle_cell1_ambient | E-cal NOT SUITABLE; E-val SUITABLE WITH CAVEAT; T-cal NOT SUITABLE; T-val QUALITATIVE / DIAGNOSTIC ONLY | CONDITIONAL ELECTRICAL TRANSFER DIAGNOSTICS IN CAPACITY/SOC COORDINATES | LEVEL 1: CAPACITY / STATE-MAPPING ADAPTATION REQUIRED; LEVEL 2: CELL-SPECIFIC VOLTAGE-REFERENCE / OCP ADAPTATION REQUIRED; LEVEL 3: DYNAMIC PARAMETER ADAPTATION ALSO REQUIRED; LEVEL 4: AVAILABLE DATA INSUFFICIENT FOR DEFENSIBLE QUANTITATIVE MODEL IDENTIFICATION | NO TRC MODEL ACCURACY CLAIM: INCOMPATIBLE DIRECT TRANSFER | strongest own quantitative voltage/capacity comparison with timing caveat; qualitative thermal comparison |
| cycle_cell1_30c | E-cal NOT SUITABLE; E-val QUALITATIVE / DIAGNOSTIC ONLY; T-cal NOT SUITABLE; T-val QUALITATIVE / DIAGNOSTIC ONLY | DIAGNOSTIC CASE STUDY ONLY AFTER NONIDEALITY GATE | NOT ENTERED: experimental nonideality/data-quality gate | NOT ASSESSED | experimental-nonideality case study; qualitative temperature-response diagnostic |
| cycle_cell1_35c | E-cal NOT SUITABLE; E-val QUALITATIVE / DIAGNOSTIC ONLY; T-cal NOT SUITABLE; T-val QUALITATIVE / DIAGNOSTIC ONLY | DIAGNOSTIC CASE STUDY ONLY AFTER NONIDEALITY GATE | NOT ENTERED: experimental nonideality/data-quality gate | NOT ASSESSED | experimental-nonideality case study; qualitative temperature-response diagnostic |
| cycle_cell2_undeclared | E-cal NOT SUITABLE; E-val QUALITATIVE / DIAGNOSTIC ONLY; T-cal NOT SUITABLE; T-val QUALITATIVE / DIAGNOSTIC ONLY | DIAGNOSTIC CASE STUDY ONLY AFTER NONIDEALITY GATE | NOT ENTERED: experimental nonideality/data-quality gate | NOT ASSESSED | experimental-nonideality case study; qualitative temperature-response diagnostic |
| cycle_cell2_30c | E-cal NOT SUITABLE; E-val QUALITATIVE / DIAGNOSTIC ONLY; T-cal NOT SUITABLE; T-val QUALITATIVE / DIAGNOSTIC ONLY | DIAGNOSTIC CASE STUDY ONLY AFTER NONIDEALITY GATE | NOT ENTERED: experimental nonideality/data-quality gate | NOT ASSESSED | experimental-nonideality case study; qualitative temperature-response diagnostic |

## Final transfer decision

The hierarchical gate returns: `LEVEL 1: CAPACITY / STATE-MAPPING ADAPTATION REQUIRED; LEVEL 2: CELL-SPECIFIC VOLTAGE-REFERENCE / OCP ADAPTATION REQUIRED; LEVEL 3: DYNAMIC PARAMETER ADAPTATION ALSO REQUIRED; LEVEL 4: AVAILABLE DATA INSUFFICIENT FOR DEFENSIBLE QUANTITATIVE MODEL IDENTIFICATION`.

- LEVEL 1 is required because the TRC 2.8999 Ah discharge does not traverse the frozen Panasonic capacity/state mapping completely.
- LEVEL 2 is required because the voltage-backbone residual is large and SOC-dependent.
- LEVEL 3 is indicated by residual dynamic mismatch but cannot be calibrated with provisional time.
- LEVEL 4 records that full OCP, dynamic, and thermal identification is unsupported by the available own data.

## Reporting policy

- Voltage-vs-reference-SOC errors are quantitative only as **TRANSFER DIAGNOSTICS** until a compatible TRC parameterization exists.
- Capacity and R_app are quantitative **DATA-QUALITY DIAGNOSTICS**; R_app is never intrinsic R0.
- SOC accuracy is not reported because no independent SOC ground truth exists.
- Dynamic time-domain results are provisional and not primary.
- **THERMAL DYNAMIC ACCURACY NOT ASSESSABLE** for every own dataset because measured time is unavailable.
- Temperature versus capacity/SOC remains **DESCRIPTIVE ONLY**.

## Provenance and stop point

This consolidation reads `own_data_quality_audit.csv`, `own_case_study_cross_experiment_metrics.csv`, `own_data_apparent_resistance_transitions.csv`, and `own_model_transfer_metrics.csv`. Detailed source analyses remain in `own_case_study_nonideality_report.md` and `own_model_transfer_audit.md`.

No model was rerun, no parameter was fitted, no OCP/window was rescaled, no voltage was corrected, and affected data were not used as calibration truth.
