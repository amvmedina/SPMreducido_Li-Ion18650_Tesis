# Panasonic 25 °C HPPC data audit and preparation

## Scope

Read-only segmentation and diagnostics only. `R_app` is an apparent instantaneous step-response quantity and is **not** automatically the model parameter `R0`. No dynamic parameters were optimized.

## Sources and signals

- Raw MAT path recorded in the manifest: `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\5 pulse disch\03-11-17_08.47 25degC_5Pulse_HPPC_Pan18650PF.mat`
- Read-only converted CSV loaded: `C:\Users\amvel\Documents\Andes\Tesis\datasets\converted_csv\hppc_25C_converted.csv`
- Columns: `Time, Voltage, Current, Ah, Wh, Power, Battery_Temp_degC, Chamber_Temp_degC`
- Rows: 102800
- Elapsed-time range: 0.000000 to 97599.399000 s
- Total logged duration: 97599.399000 s (27.110944 h)
- Voltage range: 2.49819 to 4.17497 V
- Current range: -17.40298 to 0.00000 A
- Current convention: raw negative current is discharge; raw positive current would be charge; zero is rest
- Battery case temperature: 25.33929 to 27.93959 °C
- Chamber temperature: 25.00 to 25.00 °C

## Sampling and gaps

- Median positive sample interval overall: 0.996001 s
- Median within-pulse sample interval: 0.100995 s
- Median contiguous-rest sample interval: 0.996997 s
- Duplicate elapsed-time increments: 153
- Long inter-group gaps greater than 10 s: 13, spanning 1947.845 to 3748.545 s
- The long gaps separate pulse groups. The cumulative `Ah` counter jumps consistently across them and is therefore used for SOC; direct integration of this pulse-only CSV would miss the omitted between-group discharge.

## Reproducible segmentation

- Rest: `abs(I) <= 0.05 A`
- Discharge pulse: `I <= -0.50 A`
- Charge pulse: `I >= +0.50 A`
- A time jump greater than 10 s starts a new segment
- Rest segments: 81; duration range 9.897 to 1199.941 s
- Discharge pulse segments: 67
- Charge/recovery pulse segments: 0 (recovery occurs by zero-current rest)
- Pulse counts by nominal magnitude: 0.5C: 14, 1C: 14, 2C: 14, 4C: 13, 6C: 12
- Pulse-duration range: 0.799 to 10.013 s
- Post-pulse rest-duration range: 58.998 to 1199.941 s

Pulse indices in `panasonic_hppc_pulse_summary.csv` are zero-based CSV data-row indices (the header is not counted). `V_pre` is the last rest sample before the pulse; `V_immediate` is the first loaded sample; `V_end` is the last loaded sample; `V_post_rest` is the final contiguous rest sample before the next pulse or long inter-group gap. `R_app` uses the median pulse current relative to the pre-pulse rest current.

## SOC and frozen nominal electrode mapping

- Fixed Panasonic capacity: `2.994979138519 Ah`
- `Q_discharge = -Ah` and `SOC = 1 - Q_discharge/Q_cell`
- Frozen nominal effective negative window: 0.720451136448 to 0.251849704066
- Frozen nominal effective positive window: 0.224027933455 to 0.953477823621
- Pulses outside frozen OCP supported domains: 0
- These mappings use nominal effective windows, not uniquely identified physical electrode stoichiometries.

## Usability assessment

- Detected discharge pulses: 67
- Usable complete pulses: **64**
- Low-voltage-cutoff-truncated pulses: 3 (IDs 60, 64, 67)
- Usable pulse-start SOC range: 7.878% to 100.000%
- Usable `R_app` range: 20.071 to 33.348 mΩ

**Calibration suitability:** The 64 complete pulses provide repeated 0.5C--6C excitation over a broad SOC range and appear suitable for a subsequent dynamic-parameter calibration after the three cutoff-truncated pulses are excluded or treated separately. The dataset also contains long omitted-data gaps, variable 0.1/1 s sampling, duplicate timestamps, and progressively missing high-current levels near low SOC; loaders must segment groups and must not integrate current across missing intervals.

No HPPC parameter fitting, UDDS, 1C validation, or model-equation change was performed in this audit.
