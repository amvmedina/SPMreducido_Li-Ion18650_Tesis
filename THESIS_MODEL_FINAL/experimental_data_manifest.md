# Experimental data manifest

Audit date: 2026-09-12  
Repository: `C:\Users\amvel\Documents\Andes\THESIS_MODEL_FINAL`  
Scope: read-only discovery and characterization of experimental data. No raw
file was copied, edited, renamed, moved, or deleted. No parameter was fitted.
The frozen model files were not modified.

## 1. Data sources and files found

### A. Own experimental data: TRC/Chroma measurements

These seven source CSVs are clearly the user's measurements. Cell identity is
available only as `cell1` or `cell2` in the filename; no Chroma channel number,
serial number, lot number, or manufacturer/model is populated in the files.

| Absolute source path | Cell | Experiment and condition |
| --- | --- | --- |
| `C:\Users\amvel\Documents\Andes\Tesis\TRC18650experimental data\amvTRC18650LiDisPulse1C cell1 part1.csv` | cell 1 | Approximately 1C pulse discharge/rest, ambient/chamber temperature recorded; interrupted part 1 |
| `C:\Users\amvel\Documents\Andes\Tesis\TRC18650experimental data\amvTRC18650LiDisPulse1C cell1 part2.csv` | cell 1 | Continuation of pulse discharge/rest; no header in this part |
| `C:\Users\amvel\Documents\Andes\Tesis\TRC18650experimental data\amvTRC18650LiIonCycle1C Tambient cell1.csv` | cell 1 | `Aire 22–24 °C`: approximately 1C discharge followed by approximately 1C CC-CV charge; initial charge excluded by confirmed notes; simultaneous with the previously undeclared cell-2 test and source of their shared measured-air channel `Temp7` |
| `C:\Users\amvel\Documents\Andes\Tesis\TRC18650experimental data\amvTRC18650LiIonCycle1C T30C cell1.csv` | cell 1 | 30 degC approximately 1C discharge followed by CC-CV charge |
| `C:\Users\amvel\Documents\Andes\Tesis\TRC18650experimental data\amvTRC18650LiIonCycle1C T35C cell1.csv` | cell 1 | 35 degC approximately 1C discharge followed by CC-CV charge |
| `C:\Users\amvel\Documents\Andes\Tesis\TRC18650experimental data\amvTRC18650LiIonCycle1C cell2 07232026.csv` | cell 2 | `Aire 22–24 °C`: approximately 1C cycle performed simultaneously with the first cell-1 laboratory-temperature test; both tests share the measured `Temp7` air record (22.49–23.49 °C), although sample-by-sample alignment is not established |
| `C:\Users\amvel\Documents\Andes\Tesis\TRC18650experimental data\amvTRC18650LiIonCycle1C T30C cell2 07242026.csv` | cell 2 | 30 degC approximately 1C cycle |

Supporting own chamber/controller logs were also found. They contain controller
temperature and setpoint, not battery voltage or current, and cannot be aligned
to Chroma data from the available CSV columns alone:

- `C:\Users\amvel\Documents\Andes\Tesis\Control Termico\datos\control_temperatura_30C_20260724_112526.csv`
- `C:\Users\amvel\Documents\Andes\Tesis\Control Termico\datos\control_temperatura_30C_20260724_131703.csv`
- `C:\Users\amvel\Documents\Andes\Tesis\Control Termico\datos\control_temperatura_30C_20260724_152319.csv`
- `C:\Users\amvel\Documents\Andes\Tesis\Control Termico\datos\control_temperatura_35C_20260725_095052.csv`
- `C:\Users\amvel\Documents\Andes\Tesis\Control Termico\datos\control_temperatura_35C_20260725_102916.csv`
- `C:\Users\amvel\Documents\Andes\Tesis\Control Termico\datos\control_temperatura_35C_20260725_104121.csv`
- `C:\Users\amvel\Documents\Andes\Tesis\Control Termico\datos\control_temperatura_35C_20260725_105316.csv`
- `C:\Users\amvel\Documents\Andes\Tesis\Control Termico\datos\control_temperatura_35C_20260725_113213.csv`
- `C:\Users\amvel\Documents\Andes\Tesis\Control Termico\datos\control_temperatura_35C_20260725_113859.csv`
- `C:\Users\amvel\Documents\Andes\Tesis\Control Termico\datos\control_temperatura_35C_20260725_134620.csv`
- `C:\Users\amvel\Documents\Andes\Tesis\Control Termico\datos\datos_temperatura.csv`
- `C:\Users\amvel\Documents\Andes\Tesis\Control Termico\datos\datos_temperatura_impulso.csv`

### B. External/reported dataset: Panasonic 18650PF

The source README identifies a new 2.9 Ah Panasonic 18650PF cell tested by Dr.
Phillip Kollmeyer at the University of Wisconsin-Madison using a Digatron tester
and thermal chamber. The selected 25 degC source files are:

| Profile | Absolute raw source path |
| --- | --- |
| Dataset description | `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Readme file - desc of tests performed.txt` |
| C/20 discharge/charge | `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\C20 OCV and 1C discharge tests_start_of_tests\05-08-17_13.26 C20 OCV Test_C20_25dC.mat` |
| Five-level HPPC, 0.5/1/2/4/6C | `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\5 pulse disch\03-11-17_08.47 25degC_5Pulse_HPPC_Pan18650PF.mat` |
| UDDS | `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\Drive cycles\03-21-17_00.29 25degC_UDDS_Pan18650PF.mat` |
| Initial-series 1C discharge | `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\C20 OCV and 1C discharge tests_start_of_tests\03-09-17_17.59 3349_Dis1C_1.mat` |
| Associated CC-CV charge | `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\C20 OCV and 1C discharge tests_start_of_tests\03-09-17_19.03 3349_Charge2.mat` |

Additional 1C capacity/aging candidates are present in
`C:\Users\amvel\Documents\Andes\Tesis\datasets\raw_mat`:

- `C:\Users\amvel\Documents\Andes\Tesis\datasets\raw_mat\03-10-17_23.36 3349_Dis1C_2.mat`
- `C:\Users\amvel\Documents\Andes\Tesis\datasets\raw_mat\07-22-17_22.44 4020_Dis1C_1.mat`
- `C:\Users\amvel\Documents\Andes\Tesis\datasets\raw_mat\07-24-17_07.00 4020_Dis1C_2.mat`

Existing working conversions, treated as derived data rather than raw sources:

- `C:\Users\amvel\Documents\Andes\Tesis\datasets\converted_csv\c20_ocv_25C_converted.csv`
- `C:\Users\amvel\Documents\Andes\Tesis\datasets\converted_csv\c20_25C_converted.csv` (same rows with fewer columns)
- `C:\Users\amvel\Documents\Andes\Tesis\datasets\converted_csv\hppc_25C_converted.csv`
- `C:\Users\amvel\Documents\Andes\Tesis\datasets\converted_csv\udds_25C_converted.csv`
- `C:\Users\amvel\Documents\Andes\Tesis\datasets\converted_csv\discharge_1C_25C_converted.csv`
- `C:\Users\amvel\Documents\Andes\Tesis\datasets\converted_csv\discharge_1C_25C_3349_2_converted.csv`
- `C:\Users\amvel\Documents\Andes\Tesis\datasets\converted_csv\discharge_1C_25C_4020_1_converted.csv`
- `C:\Users\amvel\Documents\Andes\Tesis\datasets\converted_csv\discharge_1C_25C_4020_2_converted.csv`
- `C:\Users\amvel\Documents\Andes\Tesis\datasets\converted_csv\charge_3349_2_converted.csv`

The full external source contains the same profile families at 25, 10, 0, -10,
and -20 degC plus temperature-rise experiments. Those additional files are not
recommended for the minimal first plan.

### C. Unknown/uncertain candidates

Seventeen ZIP archives were found under
`C:\Users\amvel\Documents\Andes\Tesis\Maryland Dataset`. Their filenames show
low-current OCV, initial-capacity, BJDST, DST, FUDS, and US06 tests on units
labelled SP20-1 and SP20-2 at 0, 25, and 45 degC. No local README or source
metadata was found that confidently establishes cell manufacturer/model,
capacity, sensor definitions, or provenance, so they are not selected now:

- `C:\Users\amvel\Documents\Andes\Tesis\Maryland Dataset\SP1_0C_LC_OCV_02_24_2016.zip`
- `C:\Users\amvel\Documents\Andes\Tesis\Maryland Dataset\SP1_25C_LC_OCV_11_5_2015.zip`
- `C:\Users\amvel\Documents\Andes\Tesis\Maryland Dataset\SP1_45C_LC_OCV_11_21_2015.zip`
- `C:\Users\amvel\Documents\Andes\Tesis\Maryland Dataset\SP1_Initial capacity_10_16_2015.zip`
- `C:\Users\amvel\Documents\Andes\Tesis\Maryland Dataset\SP2_0C_BJDST.zip`
- `C:\Users\amvel\Documents\Andes\Tesis\Maryland Dataset\SP2_0C_DST.zip`
- `C:\Users\amvel\Documents\Andes\Tesis\Maryland Dataset\SP2_0C_FUDS.zip`
- `C:\Users\amvel\Documents\Andes\Tesis\Maryland Dataset\SP2_0C_US06.zip`
- `C:\Users\amvel\Documents\Andes\Tesis\Maryland Dataset\SP2_25C_BJDST.zip`
- `C:\Users\amvel\Documents\Andes\Tesis\Maryland Dataset\SP2_25C_DST.zip`
- `C:\Users\amvel\Documents\Andes\Tesis\Maryland Dataset\SP2_25C_FUDS.zip`
- `C:\Users\amvel\Documents\Andes\Tesis\Maryland Dataset\SP2_25C_US06.zip`
- `C:\Users\amvel\Documents\Andes\Tesis\Maryland Dataset\SP2_45C_BJDST.zip`
- `C:\Users\amvel\Documents\Andes\Tesis\Maryland Dataset\SP2_45C_DST.zip`
- `C:\Users\amvel\Documents\Andes\Tesis\Maryland Dataset\SP2_45C_FUDS.zip`
- `C:\Users\amvel\Documents\Andes\Tesis\Maryland Dataset\SP2_45C_US06.zip`
- `C:\Users\amvel\Documents\Andes\Tesis\Maryland Dataset\SP2_Initial_capacity_10_16_2015.zip`

The A123/Guo material under `Experimental validation`, `third_party`, and the
closed synthetic-output directories is intentionally excluded from the usable
experimental-data candidates by the frozen-architecture instruction.

## 2. Own-data signal and test audit

### Column mapping and units

For the five cell-1 CSVs, the exported columns are `Description`, `Step No.`,
`Voltage(V)`, `Current(A)`, `Capacity(Ah)`, `Mode`, `Total Capacity(Ah)`,
`Temp1`, `Temp2`, and `Temp7`. Experimenter-confirmed metadata states:

- `Temp1`: cell surface temperature at the middle of the cell [degC].
- `Temp2`: cell surface temperature near one end [degC].
- `Temp7`: chamber-air or ambient temperature [degC].

The two cell-2 files instead contain `Temp3` and `Temp4`. Their locations are
not documented in the available metadata and no chamber-temperature column is
present inside those files. The previously undeclared cell-2 test was performed
simultaneously with the first cell-1 laboratory-temperature test and shares its
measured air thermocouple, `Temp7` (22.49–23.49 degC). This is a confirmed
shared measurement, not an inferred temperature. `Temp7` is not copied into
the cell-2 signal arrays, and different record counts plus the absent measured
timebase mean that sample-by-sample alignment is not established. `Temp3` and
`Temp4` remain separate measured signals of undocumented placement.

Presentation labels summarize measured-air ranges, not programmed references:
`Aire 22–24 °C`, `Aire 29–31 °C`, and `Aire 33–37 °C`. Numerical tables retain
the exact `Temp7` minima and maxima.

There is **no elapsed-time or date/time column in any own-data CSV**. Periods
shown below were inferred from exported capacity/current and are provisional,
not timestamps. `part2` has no header; its schema is positionally inferred from
`part1`. The pulse boundary has strong value continuity (delta V = -0.4 mV,
delta I = 0 A, delta total capacity = 0 Ah), but a lost-time gap cannot be ruled
out.

### Numeric characterization

| File/test | Samples and timing | Current [A] | Voltage [V] | Temperature evidence | Capacity evidence |
| --- | --- | ---: | ---: | --- | --- |
| Pulse cell 1, joined parts 1+2 | 74,265; about 0.10 s/sample and 7,383.9 s total, both provisional | -2.9018 to 0 | 2.5893 to 4.1031 | Temp1 26.15-30.44 degC; Temp2 26.15-30.85 degC; Temp7 23.58-24.86 degC | Partial pulse train removes 2.2497 Ah; not a complete capacity test |
| Ambient cell 1 cycle | 4,241; about 2.994 s/sample and 12,694.1 s total, provisional | -2.9017 to +2.9041 | 2.8342 to 4.2044 | Temp1 23.84-33.10; Temp2 23.54-32.20; Temp7 22.49-23.49 degC | Relevant discharge 2.8999 Ah; following CC-CV charge 2.8999 Ah |
| 30 degC cell 1 cycle | 4,372; about 2.992 s/sample and 13,077.6 s total, provisional | -2.9016 to +2.9042 | 2.5999 to 4.2047 | Temp1 29.82-37.87; Temp2 29.58-36.32; Temp7 29.38-30.92 degC | Discharge 2.4894 Ah; following charge 2.4653 Ah |
| 35 degC cell 1 cycle | 4,937; about 2.992 s/sample and 14,770.7 s total, provisional | -2.9017 to +2.9040 | 2.6000 to 4.2047 | Temp1 34.27-44.40; Temp2 34.06-43.01; Temp7 33.34-36.94 degC | Discharge 2.6171 Ah; following charge 2.5185 Ah |
| Cell 2 cycle, Aire 22–24 °C | 6,532; about 2.98 s/sample, provisional | -2.9016 to +2.9040 | 2.6000 to 4.2047 | Temp3 22.98-34.61; Temp4 23.08-32.64 degC, placement uncertain; shared measured air Temp7 22.49-23.49 degC from the simultaneous cell-1 test, without assumed sample alignment | Discharge 2.7133 Ah; following charge 2.6897 Ah |
| 30 degC cell 2 cycle | 4,501; about 2.98 s/sample, provisional | -2.9018 to +2.9039 | 2.5999 to 4.2048 | Temp3 29.72-39.36; Temp4 29.91-37.06 degC; placement uncertain | Discharge 2.5994 Ah; following charge 2.5054 Ah |

The relevant cycle sequence confirmed in existing notes is approximately 1C
discharge followed by approximately 1C CC-CV charge. Initial charge material in
each file is not part of the comparison. The measured discharge plateau is
approximately -2.90 A raw; charge is positive and tapers during CV.

The 30 and 35 degC cell-1 files show about 0.69 V instantaneous voltage changes
at the discharge transition, versus about 0.106 V for the ambient cell-1 file.
This is evidence of a setup/contact/sense or cell-condition inconsistency that
must be resolved before pooling these files.

## 3. External Panasonic signal and test audit

The raw dataset documentation defines `TimeStamp`, `Voltage`, `Current`, `Ah`,
`Wh`, `Power`, `Battery_Temp_degC`, `Time`, and `Chamber_Temp_degC`. `Time` is
elapsed seconds. Battery temperature is case temperature at the middle of the
cell, measured by thermocouple; chamber temperature is ambient. The Digatron
channel identifier is not stated.

| Selected converted profile | Sample interval and duration | Current [A] | Voltage [V] | Temperature [degC] | Notes |
| --- | --- | ---: | ---: | --- | --- |
| C/20 25 degC | Active samples about 60 s; active discharge 74,380.867 s | -0.14536 to +0.14537 | active discharge 4.17030 to 2.49948 | active chamber 25 | Full converted file includes long pauses/later chamber changes; use only the identified discharge branch, indices 6-1246 |
| HPPC 25 degC | Median positive increment 0.996 s; full logged span 97,599.399 s with long gaps between pulse groups | -17.40298 to 0 | 2.49819 to 4.17497 | battery 25.34-27.94; chamber 25 | Five discharge levels: 0.5, 1, 2, 4, and 6C across SOC; between-pulse discharges are in a separate source file |
| UDDS 25 degC | Median 0.098994 s; 22,447.866 s | -7.52923 to +4.12497 | 2.50334 to 4.20779 | battery 25.62-27.11; chamber 25 | Dynamic discharge plus regenerative charge |
| Initial 1C discharge 25 degC | About 10 s; active discharge 3,474.369 s | about -2.89982 to -2.89900 | 4.04420 loaded start; 2.49948 cutoff | battery up to 32.93; chamber 25-26 | Complete initial-series 1C discharge, 2.798235809 Ah by independent time integration |
| Associated CC-CV charge 25 degC | Median 60 s; 7,190.124 s | 0 to +2.89997 | 3.21117 to 4.20007 | battery 25.61-30.01; chamber 25 | Exported Ah reaches 2.78376 Ah |

The source README specifies charging at 1C to 4.2 V with 50 mA cutoff. At 25
and 10 degC, dynamic tests stop when voltage first reaches 2.5 V. The selected
C/20 and initial 1C discharges both reach 2.49948 V.

## 4. Cell capacity information

### Own TRC cells

- Exact manufacturer/model and manufacturer nominal capacity were **not found**.
  The label/recipe says TRC18650 and 1C operation is approximately 2.90 A,
  which suggests a 2.9 Ah rating but is not a manufacturer specification.
- Own CSVs do not contain explicit time, so
  `integral(abs(I) dt)/3600` cannot be independently reconstructed exactly.
- The Chroma-exported step-capacity values nevertheless provide measured
  delivered capacity. The strongest own baseline is the complete ambient cell-1
  discharge and recharge: 2.8999 Ah in each direction.
- Other exported discharge capacities are condition/setup dependent: 2.4894 Ah
  (30 degC cell 1), 2.6171 Ah (35 degC cell 1), 2.7133 Ah (undeclared-condition
  cell 2), and 2.5994 Ah (30 degC cell 2). These should not be averaged.
- The pulse test's 2.2497 Ah is partial and must not define full-cell capacity.

**Recommended own-cell source for initial `Q_cell_Ah`: 2.8999 Ah from the
ambient cell-1 complete discharge**, subject to confirming recipe endpoints,
sample period, and that later calibration/validation files use the same physical
cell and voltage-sense setup. This is a measured usable capacity, not a fitted
parameter and not a manufacturer nominal value.

### External Panasonic cell

- Manufacturer/model: Panasonic 18650PF, as identified by the dataset README,
  directory names, and profile filenames.
- Nominal manufacturer/rated capacity in the dataset metadata: 2.9 Ah; hence
  nominal 1C = 2.9 A.
- Independently integrated complete C/20 discharge:
  `Q_measured = 2.994979138519 Ah` from 4.17030 to 2.49948 V.
- Independently integrated initial 1C discharge:
  `Q_measured = 2.798235809 Ah` to the same 2.49948 V cutoff.
- The README reports that after about 110 cycles, 1C capacity had fallen from
  about 2.8 Ah to about 2.3 Ah; the 4020 files are therefore aging references,
  not interchangeable repetitions of the initial cell.

**Recommended Panasonic `Q_cell_Ah`: 2.994979138519 Ah from the complete C/20
discharge** for a Panasonic-specific analysis. Do not use this value for the own
TRC cell.

## 5. Current sign conventions and future loading boundary

Both usable source families use the opposite sign to the frozen model:

```text
Raw own TRC and Panasonic:  I_raw < 0 discharge, I_raw > 0 charge
Frozen model:               I_model > 0 discharge, I_model < 0 charge
Required loader mapping:    I_model = -I_raw
```

No raw file should be edited. Apply and document this sign inversion only in a
future loader/processed representation. Keep voltage unchanged. Convert own
temperatures from degC to kelvin at the loader boundary using
`T_K = T_degC + 273.15`.

## 6. Minimal experimental plan

The own and Panasonic cells must remain separate data streams; capacities and
parameters must not be pooled.

### A. Parameter/model information

1. For the own cell, use the ambient cell-1 full discharge's exported 2.8999 Ah
   as the initial measured usable capacity after confirming its recipe endpoints
   and time base.
2. For the external Panasonic case, use only the identified C/20 discharge
   branch for its measured capacity (2.994979138519 Ah) and low-current voltage
   information. It does not establish the own cell's capacity or OCP.

### B. First calibration

The recommended **first own-cell calibration experiment** is the joined cell-1
pulse discharge (`part1` then `part2`) because its repeated 1C load/rest
transients are more informative for resistance, kinetics, diffusion relaxation,
and thermal response than a constant-current cycle. Use Temp1 as the primary
surface-temperature observation, Temp2 as a spatial/checking signal, and Temp7
as ambient.

This recommendation is conditional: calibration must not start until the
actual Chroma sample period and part1/part2 time continuity are confirmed, and
the unusually large pulse voltage response is checked against contacts and
voltage-sense configuration. If the time base cannot be recovered, this dataset
is not calibration-ready.

For a separate reported-data demonstration, the smallest coherent Panasonic
calibration choice is the 25 degC HPPC profile, using its exact time, voltage,
current, case temperature, and chamber temperature together with the
Panasonic-specific C/20 capacity. Do not transfer its fitted cell parameters
silently to TRC.

### C. Independent validation

- Primary own-cell validation: the ambient cell-1 1C discharge segment, which
  is not used in pulse calibration and contains voltage, two surface
  temperatures, ambient temperature, and exported usable capacity.
- Minimal external reported validation: Panasonic UDDS at 25 degC after a
  Panasonic-only HPPC calibration. UDDS is dynamic and was not used for fitting.

### D. Minimal transfer experiment

Apply the own cell-1 pulse-calibrated parameters **without recalibration** to the
35 degC cell-1 1C discharge. This is the smallest temperature-transfer test and
uses measured surface and chamber temperatures. Before interpreting transfer
error, confirm that it is the same physical cell/channel and resolve the large
voltage-step discrepancy. Only after that test should a small subset such as
thermal boundary parameters be considered for readjustment.

## 7. Information missing before calibration

1. Actual elapsed-time/sample-period setting for every Chroma file, preferably
   a re-export containing timestamps.
2. Whether any time or samples were lost between pulse part1 and part2.
3. Confirmation that pulse, ambient, 30 degC, and 35 degC `cell1` refer to the
   same physical cell and Chroma channel.
4. Explanation of the much larger current-transition voltage steps in the
   pulse/30/35 degC runs than in the ambient run: contacts, two/four-wire sense,
   lead resistance, channel, recipe, or cell condition.
5. Exact voltage/current/capacity termination rules in the own recipes and the
   meaning of Chroma `Mode` codes 0 and 1.
6. Manufacturer, exact model, nominal capacity, age/cycle history, and serial
   identifiers for own cells 1 and 2.
7. Confirmation of Temp3 and Temp4 sensor locations. The shared measured
   `Temp7` air record for the simultaneous `Aire 22–24 °C` tests is confirmed,
   but sample-by-sample alignment remains unavailable without timestamps.
8. A defensible time alignment between chamber-controller logs and Chroma data,
   if those separate logs are to be used.
9. Chemistry compatibility of the frozen model's retained OCP functions with
   each selected cell. The data audit does not establish this.
10. Source/provenance and cell metadata for the Maryland ZIP archives before
    any use.

Until items 1-4 are resolved, the own pulse and temperature-transfer files are
valuable candidates but are not ready for quantitative dynamic calibration.

## 8. Corrected discharge-plot mapping (user-confirmed)

This correction applies **only** to the discharge voltage and recorded
discharged-capacity curves used in the Uniandes comparison figures. It does not
reassign charge records, current, capacity totals, temperature samples, or
other source-file metadata.

| Final plotted condition | Electrical source dataset | Electrical source fields | Thermal provenance retained separately |
|---|---|---|---|
| Rionegro · Aire 22–24 °C · reference setup | `cycle_cell1_ambient` | `Voltage(V)` and `Q_dis=abs(Q-Q_start)` from that file's `Capacity(Ah)` discharge branch | same cell-1 record; its `Temp7` is also the measured shared-air channel |
| Rionegro · Aire 22–24 °C · external-resistance setup | `cycle_cell2_undeclared` | `Voltage(V)` and `Q_dis=abs(Q-Q_start)` from that file's `Capacity(Ah)` discharge branch | cell-2 surface channels stay with cell 2; shared measured `Temp7` comes from the simultaneous cell-1 record, without claiming sample alignment |
| Cali · Aire 29–31 °C · external-resistance setup | `cycle_cell2_30c` | `Voltage(V)` and `Q_dis=abs(Q-Q_start)` from that file's `Capacity(Ah)` discharge branch | temperature provenance remains in its previously established records; no cross-file sample combination |
| Valledupar · Aire 33–37 °C · external-resistance setup | `cycle_cell1_30c` | `Voltage(V)` and `Q_dis=abs(Q-Q_start)` from that file's `Capacity(Ah)` discharge branch | `cycle_cell1_35c` remains the temperature record for this condition; its temperature columns are not transferred onto the electrical rows |

The `cycle_cell1_35c` voltage and discharged-capacity curve are excluded from the corrected
discharge comparisons, but the source file and its temperature measurements
are preserved. The prior axis forced the first and last records to SOC values
of 1 and 0 using `SOC_ref = 1 - abs(Q-Q_start)/abs(Q_end-Q_start)`. Because the
CSV files do not establish that these endpoints are fully charged and fully
discharged, the figures now use `Q_dis = abs(Q-Q_start)` and label the axis
**Discharged capacity (Ah)**. Each curve therefore ends at its actual recorded
capacity rather than at a normalized value of 1. The coordinate is not a
physical SOC or evidence of complete discharge. Detailed
mapping and endpoint audits are stored in
`thesis_data/uniandes_corrected_discharge_mapping_audit.csv` and
`thesis_data/uniandes_discharge_axis_audit.csv`.
