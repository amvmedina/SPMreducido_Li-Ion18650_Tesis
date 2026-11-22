# Own TRC18650 experimental-data quality audit

This is a read-only data-quality audit. No SPM/OCP/electrical/thermal parameter was fitted, no Panasonic parameter was applied, and no raw voltage was corrected.

## Discovery and standardized representation

Seven primary Chroma/TRC18650 CSVs were found under `C:\Users\amvel\Documents\Andes\Tesis\TRC18650experimental data` and audited individually. 12 chamber/controller CSVs were also found under `C:\Users\amvel\Documents\Andes\Tesis\Control Termico\datos`; they are supporting temperature logs, not distinct battery experiments, and the available columns do not provide a defensible absolute alignment to the Chroma files.

Every primary CSV contains voltage, current, step capacity, signed total capacity, step number/mode, and at least two temperature channels. The reusable loader exposes canonical voltage/current views while preserving the raw fields and their sign. None contains elapsed time or a timestamp. Therefore all duration, dt, current integration, and time-domain plots are explicitly **provisional**, using a uniform sampling interval inferred from active-step capacity/current. Raw timestamp duplicates, reversals, and long gaps are **not assessable**. The raw data remain unchanged.

Because the provisional timebase is inferred from the same capacity/current records, agreement between reconstructed current integration and the Ah counters is a consistency check, **not an independent capacity measurement**.

Signal evidence establishes `I_raw<0` for discharge and `I_raw>0` for charge; a future model loader must record `I_model=-I_raw` without overwriting source current.

The `cycle_cell1_ambient` and `cycle_cell2_undeclared` experiments were simultaneous and shared the measured air thermocouple `Temp7` (22.4896–23.4919 degC). This channel is associated with both tests as experimental context but is not copied into the cell-2 arrays. Different record counts and absent timestamps prevent sample-by-sample alignment. `Temp3` and `Temp4` remain separate measurements of undocumented placement.

## Suitability summary

| Dataset | Electrical calibration | Electrical validation | Thermal calibration | Thermal validation | Recommended thesis use |
|---|---|---|---|---|---|
| pulse_cell1_part1 | NOT SUITABLE | QUALITATIVE / DIAGNOSTIC ONLY | NOT SUITABLE | QUALITATIVE / DIAGNOSTIC ONLY | experimental-nonideality case study; qualitative temperature-response diagnostic |
| pulse_cell1_part2 | NOT SUITABLE | QUALITATIVE / DIAGNOSTIC ONLY | NOT SUITABLE | QUALITATIVE / DIAGNOSTIC ONLY | experimental-nonideality case study; qualitative temperature-response diagnostic |
| cycle_cell1_ambient | NOT SUITABLE | SUITABLE WITH CAVEAT | NOT SUITABLE | QUALITATIVE / DIAGNOSTIC ONLY | strongest own quantitative voltage/capacity comparison with timing caveat; qualitative thermal comparison |
| cycle_cell1_30c | NOT SUITABLE | QUALITATIVE / DIAGNOSTIC ONLY | NOT SUITABLE | QUALITATIVE / DIAGNOSTIC ONLY | experimental-nonideality case study; qualitative temperature-response diagnostic |
| cycle_cell1_35c | NOT SUITABLE | QUALITATIVE / DIAGNOSTIC ONLY | NOT SUITABLE | QUALITATIVE / DIAGNOSTIC ONLY | experimental-nonideality case study; qualitative temperature-response diagnostic |
| cycle_cell2_undeclared | NOT SUITABLE | QUALITATIVE / DIAGNOSTIC ONLY | NOT SUITABLE | QUALITATIVE / DIAGNOSTIC ONLY | experimental-nonideality case study; qualitative temperature-response diagnostic |
| cycle_cell2_30c | NOT SUITABLE | QUALITATIVE / DIAGNOSTIC ONLY | NOT SUITABLE | QUALITATIVE / DIAGNOSTIC ONLY | experimental-nonideality case study; qualitative temperature-response diagnostic |

## Apparent rest-to-load series resistance

Common rule: use adjacent records with `abs(I_pre)<=0.02 A`, `abs(I_loaded)>=2.5 A`, and `abs(delta I)>=2.5 A`; then calculate `R_app=abs(delta V)/abs(delta I)`. Both discharge and charge onsets are retained and labelled. This is the apparent response of the observed cell/fixture/wiring/measurement/electrochemical path, **not intrinsic cell resistance or model R0**.

The lowest group was selected from experiment medians by the largest adjacent gap in log(R_app), not by a universal resistance threshold. The data-driven split is `0.079808200 ohm`. Reference-group median/IQR are `0.035875537/0.000507537 ohm`.

| Dataset | N | Median [ohm] | Mean [ohm] | Std [ohm] | Range [ohm] | IQR [ohm] | Excess [ohm] | Ratio | Interpretation |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| pulse_cell1_part1 | 7 | 0.221348508 | 0.224567424 | 0.006841320 | 0.218786625-0.235910878 | 0.007900583 | 0.185472971 | 6.170x | substantially elevated relative to reference group; internally consistent elevated band |
| pulse_cell1_part2 | 3 | 0.228938929 | 0.228032922 | 0.009598313 | 0.218013729-0.237146108 | 0.009566190 | 0.193063392 | 6.381x | substantially elevated relative to reference group; internally consistent elevated band |
| cycle_cell1_ambient | 2 | 0.035875537 | 0.035875537 | 0.000717766 | 0.035368000-0.036383074 | 0.000507537 | 0.000000000 | 1.000x | lowest internally consistent reference group; insufficient transitions for distribution stability |
| cycle_cell1_30c | 2 | 0.232494022 | 0.232494022 | 0.011473781 | 0.224380834-0.240607211 | 0.008113188 | 0.196618486 | 6.481x | substantially elevated relative to reference group; insufficient transitions for distribution stability |
| cycle_cell1_35c | 2 | 0.222850133 | 0.222850133 | 0.020765590 | 0.208166644-0.237533623 | 0.014683490 | 0.186974596 | 6.212x | substantially elevated relative to reference group; insufficient transitions for distribution stability |
| cycle_cell2_undeclared | 2 | 0.177540167 | 0.177540167 | 0.010920110 | 0.169818483-0.185261851 | 0.007721684 | 0.141664630 | 4.949x | substantially elevated relative to reference group; insufficient transitions for distribution stability |
| cycle_cell2_30c | 2 | 0.215820169 | 0.215820169 | 0.015169459 | 0.205093741-0.226546596 | 0.010726427 | 0.179944632 | 6.016x | substantially elevated relative to reference group; insufficient transitions for distribution stability |

The elevated group is separated from the reference by a factor of roughly five to six. Within the repeated pulse records, every transition remains on the elevated side of the data-driven split; this is evidence of a persistent path condition rather than one isolated jump. Attribution is not unique, so the report uses **evidence consistent with parasitic series resistance** and does not subtract any inferred excess.

## Per-dataset report cards

### pulse_cell1_part1

- **Dataset:** `C:\Users\amvel\Documents\Andes\Tesis\TRC18650experimental data\amvTRC18650LiDisPulse1C cell1 part1.csv`
- **Cell / nominal test temperature / signal-verified type:** cell 1 / ambient/undeclared / pulse discharge/rest.
- **Time integrity:** no raw time; provisional dt `0.099787441 s`, duration `5986.149 s`; raw duplicate/non-monotonic timestamps and long gaps not assessable.
- **Current sign:** discharge_negative (high confidence); range `-2.90180` to `0.00000 A`.
- **Capacity consistency:** main discharge/charge counters `0.571900/not assessable Ah`; provisional integrations `0.571900/not assessable Ah`. Whole-file provisional signed/absolute integrals are `-2.030385/2.030385 Ah`, versus signed/absolute counter changes `-2.030200/2.030200 Ah`; differences `-0.000185/0.000185 Ah`. Complete sequence `False`; SOC status: `not justified`.
- **Voltage quality:** range `2.8109-4.1031 V`; main discharge start/end `3.3927/3.2076 V` with preceding/following rest `True/True`; largest adjacent jump `0.6840 V`; largest jump with `abs(delta I)<0.05 A` `0.0293 V`; near-extreme repeat count `2`; inconsistent transition signs `0`.
- **Thermal quality:** Temp1: initial/final 28.854800/28.180100, range 26.149900-30.435600 degC, apparent step 0.0063 degC, max adjacent jump 0.278100 degC, robust jump flags 0 above 0.500000 degC, missing 0 (cell surface middle), Temp2: initial/final 28.918300/28.349300, range 26.148800-30.849400 degC, apparent step 0.0045 degC, max adjacent jump 0.295200 degC, robust jump flags 0 above 0.500000 degC, missing 0 (cell surface end), Temp7: initial/final 24.457400/24.376900, range 23.581000-24.857600 degC, apparent step 0.0083 degC, max adjacent jump 0.416500 degC, robust jump flags 0 above 0.500000 degC, missing 0 (chamber/ambient); active runs with net heating `5`, rest runs with net cooling `8`. Sensor-pair offset/correlation/lag are `-0.140177 degC / 0.996786 / -12 samples`.
- **Apparent series resistance:** N `7`, median `0.221348508 ohm`; `substantially elevated relative to reference group`. evidence consistent with parasitic series resistance; source not uniquely identifiable.
- **Electrical calibration:** NOT SUITABLE — apparent rest-to-load resistance is 6.17x the clean reference median; no measured time channel.
- **Electrical validation:** QUALITATIVE / DIAGNOSTIC ONLY — voltage response contains a systematic elevated transition drop consistent with setup/path nonideality.
- **Thermal calibration:** NOT SUITABLE — unknown partition of extra electrical loss inside/outside cell plus absent measured timestamps.
- **Thermal validation:** QUALITATIVE / DIAGNOSTIC ONLY — temperature traces remain descriptive, but heat input is setup-confounded.
- **Recommended thesis use:** experimental-nonideality case study; qualitative temperature-response diagnostic.

### pulse_cell1_part2

- **Dataset:** `C:\Users\amvel\Documents\Andes\Tesis\TRC18650experimental data\amvTRC18650LiDisPulse1C cell1 part2.csv`
- **Cell / nominal test temperature / signal-verified type:** cell 1 / ambient/undeclared / pulse discharge/rest.
- **Time integrity:** no raw time; provisional dt `0.099767902 s`, duration `1424.087 s`; raw duplicate/non-monotonic timestamps and long gaps not assessable.
- **Current sign:** discharge_negative (high confidence); range `-2.90170` to `0.00000 A`.
- **Capacity consistency:** main discharge/charge counters `0.203400/not assessable Ah`; provisional integrations `0.203400/not assessable Ah`. Whole-file provisional signed/absolute integrals are `-0.219714/0.219714 Ah`, versus signed/absolute counter changes `-0.219500/0.219500 Ah`; differences `-0.000214/0.000214 Ah`. Complete sequence `False`; SOC status: `not justified`.
- **Voltage quality:** range `2.5893-3.4661 V`; main discharge start/end `2.7738/2.5893 V` with preceding/following rest `True/True`; largest adjacent jump `0.6877 V`; largest jump with `abs(delta I)<0.05 A` `0.0914 V`; near-extreme repeat count `3`; inconsistent transition signs `0`.
- **Thermal quality:** Temp1: initial/final 28.180100/26.578200, range 26.442400-28.773800 degC, apparent step 0.0005 degC, max adjacent jump 0.208500 degC, robust jump flags 0 above 0.500000 degC, missing 0 (cell surface middle), Temp2: initial/final 28.349300/26.667500, range 26.427700-29.085100 degC, apparent step 0.0149 degC, max adjacent jump 0.311500 degC, robust jump flags 0 above 0.500000 degC, missing 0 (cell surface end), Temp7: initial/final 24.376900/24.092600, range 23.898000-24.637500 degC, apparent step 0.0243 degC, max adjacent jump 0.249700 degC, robust jump flags 0 above 0.500000 degC, missing 0 (chamber/ambient); active runs with net heating `2`, rest runs with net cooling `4`. Sensor-pair offset/correlation/lag are `-0.108462 degC / 0.992180 / -12 samples`.
- **Apparent series resistance:** N `3`, median `0.228938929 ohm`; `substantially elevated relative to reference group`. evidence consistent with parasitic series resistance; source not uniquely identifiable.
- **Electrical calibration:** NOT SUITABLE — apparent rest-to-load resistance is 6.38x the clean reference median; no measured time channel.
- **Electrical validation:** QUALITATIVE / DIAGNOSTIC ONLY — voltage response contains a systematic elevated transition drop consistent with setup/path nonideality.
- **Thermal calibration:** NOT SUITABLE — unknown partition of extra electrical loss inside/outside cell plus absent measured timestamps.
- **Thermal validation:** QUALITATIVE / DIAGNOSTIC ONLY — temperature traces remain descriptive, but heat input is setup-confounded.
- **Recommended thesis use:** experimental-nonideality case study; qualitative temperature-response diagnostic.

### cycle_cell1_ambient

- **Dataset:** `C:\Users\amvel\Documents\Andes\Tesis\TRC18650experimental data\amvTRC18650LiIonCycle1C Tambient cell1.csv`
- **Cell / nominal test temperature / signal-verified type:** cell 1 / Aire 22–24 °C / combined cycle.
- **Time integrity:** no raw time; provisional dt `2.994702465 s`, duration `12697.538 s`; raw duplicate/non-monotonic timestamps and long gaps not assessable.
- **Current sign:** discharge_negative (high confidence); range `-2.90170` to `2.90410 A`.
- **Capacity consistency:** main discharge/charge counters `2.899900/2.899900 Ah`; provisional integrations `2.901963/2.900676 Ah`. Whole-file provisional signed/absolute integrals are `1.149646/6.958394 Ah`, versus signed/absolute counter changes `1.152000/6.952000 Ah`; differences `-0.002354/0.006394 Ah`. Complete sequence `True`; SOC status: `coulomb-counted reference SOC`.
- **Voltage quality:** range `2.8342-4.2044 V`; main discharge start/end `4.0865/2.8342 V` with preceding/following rest `True/True`; largest adjacent jump `0.1055 V`; largest jump with `abs(delta I)<0.05 A` `0.0265 V`; near-extreme repeat count `3`; inconsistent transition signs `0`.
- **Thermal quality:** Temp1: initial/final 24.631000/24.748000, range 23.842800-33.101000 degC, apparent step 0.0292 degC, max adjacent jump 0.545500 degC, robust jump flags 1 above 0.500000 degC, missing 0 (cell surface middle), Temp2: initial/final 24.197400/24.433800, range 23.544100-32.197900 degC, apparent step 0.0147 degC, max adjacent jump 0.218800 degC, robust jump flags 0 above 0.500000 degC, missing 0 (cell surface end), Temp7: initial/final 22.703000/22.818800, range 22.489600-23.491900 degC, apparent step 0.0192 degC, max adjacent jump 0.478200 degC, robust jump flags 0 above 0.500000 degC, missing 10 (chamber/ambient); active runs with net heating `3`, rest runs with net cooling `3`. Sensor-pair offset/correlation/lag are `0.546285 degC / 0.999523 / -1 samples`.
- **Apparent series resistance:** N `2`, median `0.035875537 ohm`; `lowest internally consistent reference group`. no relative elevation detected; still not an intrinsic-R0 certification.
- **Electrical calibration:** NOT SUITABLE — lowest relative apparent-resistance group and complete cycle, but no measured time channel and insufficient dynamic excitation.
- **Electrical validation:** SUITABLE WITH CAVEAT — cleanest current/voltage transitions and complete capacity coverage; timebase remains provisional.
- **Thermal calibration:** NOT SUITABLE — surface/ambient temperatures exist, but thermal dynamics cannot be calibrated defensibly without measured timestamps.
- **Thermal validation:** QUALITATIVE / DIAGNOSTIC ONLY — stable surface/ambient channels, but duration and gaps are reconstructed rather than measured.
- **Recommended thesis use:** strongest own quantitative voltage/capacity comparison with timing caveat; qualitative thermal comparison.

### cycle_cell1_30c

- **Dataset:** `C:\Users\amvel\Documents\Andes\Tesis\TRC18650experimental data\amvTRC18650LiIonCycle1C T30C cell1.csv`
- **Cell / nominal test temperature / signal-verified type:** cell 1 / Aire 29–31 °C / combined cycle.
- **Time integrity:** no raw time; provisional dt `2.993601328 s`, duration `13085.031 s`; raw duplicate/non-monotonic timestamps and long gaps not assessable.
- **Current sign:** discharge_negative (high confidence); range `-2.90160` to `2.90420 A`.
- **Capacity consistency:** main discharge/charge counters `2.489400/2.465300 Ah`; provisional integrations `2.490954/2.465105 Ah`. Whole-file provisional signed/absolute integrals are `-0.002559/4.984172 Ah`, versus signed/absolute counter changes `0.000300/4.979300 Ah`; differences `-0.002859/0.004872 Ah`. Complete sequence `True`; SOC status: `coulomb-counted reference SOC`.
- **Voltage quality:** range `2.5999-4.2047 V`; main discharge start/end `3.4850/2.5999 V` with preceding/following rest `True/True`; largest adjacent jump `0.6974 V`; largest jump with `abs(delta I)<0.05 A` `0.0439 V`; near-extreme repeat count `1`; inconsistent transition signs `0`.
- **Thermal quality:** Temp1: initial/final 29.883600/30.144700, range 29.824600-37.871000 degC, apparent step 0.0174 degC, max adjacent jump 0.201500 degC, robust jump flags 0 above 0.500000 degC, missing 0 (cell surface middle), Temp2: initial/final 29.743500/29.900500, range 29.581700-36.322900 degC, apparent step 0.0035 degC, max adjacent jump 0.265400 degC, robust jump flags 0 above 0.500000 degC, missing 0 (cell surface end), Temp7: initial/final 29.985600/29.718900, range 29.375300-30.924500 degC, apparent step 0.0134 degC, max adjacent jump 0.246900 degC, robust jump flags 0 above 0.500000 degC, missing 0 (chamber/ambient); active runs with net heating `2`, rest runs with net cooling `2`. Sensor-pair offset/correlation/lag are `0.653411 degC / 0.999450 / 0 samples`.
- **Apparent series resistance:** N `2`, median `0.232494022 ohm`; `substantially elevated relative to reference group`. evidence consistent with parasitic series resistance; source not uniquely identifiable.
- **Electrical calibration:** NOT SUITABLE — apparent rest-to-load resistance is 6.48x the clean reference median; no measured time channel.
- **Electrical validation:** QUALITATIVE / DIAGNOSTIC ONLY — voltage response contains a systematic elevated transition drop consistent with setup/path nonideality.
- **Thermal calibration:** NOT SUITABLE — unknown partition of extra electrical loss inside/outside cell plus absent measured timestamps.
- **Thermal validation:** QUALITATIVE / DIAGNOSTIC ONLY — temperature traces remain descriptive, but heat input is setup-confounded.
- **Recommended thesis use:** experimental-nonideality case study; qualitative temperature-response diagnostic.

### cycle_cell1_35c

- **Dataset:** `C:\Users\amvel\Documents\Andes\Tesis\TRC18650experimental data\amvTRC18650LiIonCycle1C T35C cell1.csv`
- **Cell / nominal test temperature / signal-verified type:** cell 1 / Aire 33–37 °C / combined cycle.
- **Time integrity:** no raw time; provisional dt `2.993682602 s`, duration `14776.817 s`; raw duplicate/non-monotonic timestamps and long gaps not assessable.
- **Current sign:** discharge_negative (high confidence); range `-2.90170` to `2.90400 A`.
- **Capacity consistency:** main discharge/charge counters `2.617100/2.518500 Ah`; provisional integrations `2.618842/2.518173 Ah`. Whole-file provisional signed/absolute integrals are `0.129680/5.372188 Ah`, versus signed/absolute counter changes `0.132700/5.367100 Ah`; differences `-0.003020/0.005088 Ah`. Complete sequence `True`; SOC status: `coulomb-counted reference SOC`.
- **Voltage quality:** range `2.6000-4.2047 V`; main discharge start/end `3.4965/2.6000 V` with preceding/following rest `True/True`; largest adjacent jump `0.6888 V`; largest jump with `abs(delta I)<0.05 A` `0.0641 V`; near-extreme repeat count `1`; inconsistent transition signs `0`.
- **Thermal quality:** Temp1: initial/final 34.316800/35.212200, range 34.267800-44.402200 degC, apparent step 0.0150 degC, max adjacent jump 0.299900 degC, robust jump flags 0 above 0.500000 degC, missing 0 (cell surface middle), Temp2: initial/final 34.128300/35.055100, range 34.059600-43.009900 degC, apparent step 0.0016 degC, max adjacent jump 0.248300 degC, robust jump flags 0 above 0.500000 degC, missing 0 (cell surface end), Temp7: initial/final 35.378000/35.214600, range 33.343200-36.937600 degC, apparent step 0.0023 degC, max adjacent jump 1.228300 degC, robust jump flags 6 above 0.500000 degC, missing 0 (chamber/ambient); active runs with net heating `4`, rest runs with net cooling `3`. Sensor-pair offset/correlation/lag are `0.584983 degC / 0.999213 / 0 samples`.
- **Apparent series resistance:** N `2`, median `0.222850133 ohm`; `substantially elevated relative to reference group`. evidence consistent with parasitic series resistance; source not uniquely identifiable.
- **Electrical calibration:** NOT SUITABLE — apparent rest-to-load resistance is 6.21x the clean reference median; no measured time channel.
- **Electrical validation:** QUALITATIVE / DIAGNOSTIC ONLY — voltage response contains a systematic elevated transition drop consistent with setup/path nonideality.
- **Thermal calibration:** NOT SUITABLE — unknown partition of extra electrical loss inside/outside cell plus absent measured timestamps.
- **Thermal validation:** QUALITATIVE / DIAGNOSTIC ONLY — temperature traces remain descriptive, but heat input is setup-confounded.
- **Recommended thesis use:** experimental-nonideality case study; qualitative temperature-response diagnostic.

### cycle_cell2_undeclared

- **Dataset:** `C:\Users\amvel\Documents\Andes\Tesis\TRC18650experimental data\amvTRC18650LiIonCycle1C cell2 07232026.csv`
- **Cell / nominal test temperature / signal-verified type:** cell 2 / Aire 22–24 °C / combined cycle.
- **Time integrity:** no raw time; provisional dt `2.992855899 s`, duration `19546.342 s`; raw duplicate/non-monotonic timestamps and long gaps not assessable.
- **Current sign:** discharge_negative (high confidence); range `-2.90160` to `2.90400 A`.
- **Capacity consistency:** main discharge/charge counters `2.713300/2.689600 Ah`; provisional integrations `2.714552/2.688465 Ah`. Whole-file provisional signed/absolute integrals are `1.356318/6.790243 Ah`, versus signed/absolute counter changes `1.360400/6.787200 Ah`; differences `-0.004082/0.003043 Ah`. Complete sequence `True`; SOC status: `coulomb-counted reference SOC`.
- **Voltage quality:** range `2.6000-4.2047 V`; main discharge start/end `3.6483/2.6000 V` with preceding/following rest `True/True`; largest adjacent jump `0.5370 V`; largest jump with `abs(delta I)<0.05 A` `0.0252 V`; near-extreme repeat count `1`; inconsistent transition signs `0`.
- **Thermal quality:** Temp3: initial/final 22.991200/23.361100, range 22.982900-34.612900 degC, apparent step 0.0044 degC, max adjacent jump 0.898700 degC, robust jump flags 2 above 0.500000 degC, missing 0 (temperature location unknown), Temp4: initial/final 23.139900/23.367800, range 23.084800-32.641100 degC, apparent step 0.0087 degC, max adjacent jump 0.469800 degC, robust jump flags 0 above 0.500000 degC, missing 0 (temperature location unknown); active runs with net heating `3`, rest runs with net cooling `3`. Sensor-pair offset/correlation/lag are `0.592631 degC / 0.999387 / 4 samples`.
- **Apparent series resistance:** N `2`, median `0.177540167 ohm`; `substantially elevated relative to reference group`. evidence consistent with parasitic series resistance; source not uniquely identifiable.
- **Electrical calibration:** NOT SUITABLE — apparent rest-to-load resistance is 4.95x the clean reference median; no measured time channel.
- **Electrical validation:** QUALITATIVE / DIAGNOSTIC ONLY — voltage response contains a systematic elevated transition drop consistent with setup/path nonideality.
- **Thermal calibration:** NOT SUITABLE — unknown partition of extra electrical loss inside/outside cell plus absent measured timestamps.
- **Thermal validation:** QUALITATIVE / DIAGNOSTIC ONLY — temperature traces remain descriptive, but heat input is setup-confounded and Temp3/Temp4 placement is undocumented; the shared measured Temp7 air channel is confirmed but cross-file sample alignment is not.
- **Recommended thesis use:** experimental-nonideality case study; qualitative temperature-response diagnostic.

### cycle_cell2_30c

- **Dataset:** `C:\Users\amvel\Documents\Andes\Tesis\TRC18650experimental data\amvTRC18650LiIonCycle1C T30C cell2 07242026.csv`
- **Cell / nominal test temperature / signal-verified type:** cell 2 / Aire 29–31 °C / combined cycle.
- **Time integrity:** no raw time; provisional dt `2.993712755 s`, duration `13471.707 s`; raw duplicate/non-monotonic timestamps and long gaps not assessable.
- **Current sign:** discharge_negative (high confidence); range `-2.90180` to `2.90390 A`.
- **Capacity consistency:** main discharge/charge counters `2.599400/2.505300 Ah`; provisional integrations `2.599496/2.505011 Ah`. Whole-file provisional signed/absolute integrals are `-0.053047/5.150768 Ah`, versus signed/absolute counter changes `-0.051600/5.147400 Ah`; differences `-0.001447/0.003368 Ah`. Complete sequence `True`; SOC status: `coulomb-counted reference SOC`.
- **Voltage quality:** range `2.5999-4.2048 V`; main discharge start/end `3.5306/2.5999 V` with preceding/following rest `True/True`; largest adjacent jump `0.6566 V`; largest jump with `abs(delta I)<0.05 A` `0.0191 V`; near-extreme repeat count `1`; inconsistent transition signs `0`.
- **Thermal quality:** Temp3: initial/final 29.730200/30.124400, range 29.724500-39.364300 degC, apparent step 0.0014 degC, max adjacent jump 0.242900 degC, robust jump flags 0 above 0.500000 degC, missing 0 (temperature location unknown), Temp4: initial/final 30.010000/30.220600, range 29.912200-37.060800 degC, apparent step 0.0014 degC, max adjacent jump 0.185300 degC, robust jump flags 0 above 0.500000 degC, missing 0 (temperature location unknown); active runs with net heating `2`, rest runs with net cooling `3`. Sensor-pair offset/correlation/lag are `0.656866 degC / 0.999356 / 3 samples`.
- **Apparent series resistance:** N `2`, median `0.215820169 ohm`; `substantially elevated relative to reference group`. evidence consistent with parasitic series resistance; source not uniquely identifiable.
- **Electrical calibration:** NOT SUITABLE — apparent rest-to-load resistance is 6.02x the clean reference median; no measured time channel.
- **Electrical validation:** QUALITATIVE / DIAGNOSTIC ONLY — voltage response contains a systematic elevated transition drop consistent with setup/path nonideality.
- **Thermal calibration:** NOT SUITABLE — unknown partition of extra electrical loss inside/outside cell plus absent measured timestamps.
- **Thermal validation:** QUALITATIVE / DIAGNOSTIC ONLY — temperature traces remain descriptive, but heat input is setup-confounded and sensor placement/ambient alignment are undocumented.
- **Recommended thesis use:** experimental-nonideality case study; qualitative temperature-response diagnostic.

## Methodological interpretation

1. **Strongest own experiment for quantitative model comparison:** `cycle_cell1_ambient`. It is the sole lowest-resistance-group experiment and has complete discharge/charge capacity plus documented cell-surface and ambient channels. Quantitative use should be restricted to voltage/capacity until an actual acquisition timebase is recovered; thermal use remains qualitative.
2. **Experiments showing setup-related voltage evidence:** pulse_cell1_part1, pulse_cell1_part2, cycle_cell1_30c, cycle_cell1_35c, cycle_cell2_undeclared, cycle_cell2_30c. Their relative separation is consistent with added series-path resistance, but the data cannot uniquely locate it in contacts, fixture, wiring, sensing, or cell state.
3. **Consistency of affected behavior:** repeated pulse transitions are internally consistent; heated cell-1 cycles occupy the same elevated band. Cell-2 values are also elevated but show more between-experiment spread, so attribution across physical cells remains cautious.
4. **Affected datasets retaining useful thermal information:** pulse_cell1_part1, pulse_cell1_part2, cycle_cell1_30c, cycle_cell1_35c retain two documented cell-surface channels plus ambient and show heating/cooling structure. They are useful descriptively, not as quantitative heat-ground-truth. Cell-2 temperature traces remain descriptive; Temp3/Temp4 placement is undocumented. The simultaneous Aire 22–24 °C cell-2 run shares the measured Temp7 air record, but cross-file sample alignment is not established.
5. **Best experimental-nonideality examples:** the two pulse parts provide repeated evidence across ten transitions; the 30/35 degC cell-1 cycles show that the elevated step also appears in nominally steady cycle tests.
6. **Framework reuse:** yes. Channel aliases/configuration, explicit current-sign metadata, time-source status, integrity/capacity/temperature checks, adjacent rest-to-load extraction, and data-driven cross-experiment grouping can be reused. New datasets with measured time improve the audit automatically without changing the methodology.

## Clean, affected, and ambiguous evidence

- **Relative clean reference:** `cycle_cell1_ambient` is the strongest clean candidate in this collection, not a certified intrinsic-resistance experiment. Only two comparable onsets are present, and the absence of raw time prevents a full dynamic-quality certification.
- **Clearly affected for electrical use:** the two cell-1 pulse parts and the 30/35 degC cell-1 cycles show persistent, approximately five-to-six-fold elevated apparent transition resistance relative to the reference. They are rejected for electrical identification.
- **Attribution remains ambiguous:** both cell-2 cycles are also elevated and therefore unsuitable for electrical identification, but comparison with cell 1 cannot separate fixture resistance from genuine cell-to-cell or temperature differences. Temp3/Temp4 placement is undocumented. The Aire 22–24 °C test has a confirmed shared measured Temp7 air record, without sample-level alignment.

The reported temperature `apparent step` is the modal nonzero adjacent digital increment, not a certified sensor resolution. The robust jump screen flags adjacent changes above `max(0.5 degC, median + 10 x 1.4826 x MAD)` for review; a flag alone is not an automatic rejection.

## Supporting controller logs

These files were discovered but are not independently classified as battery experiments because voltage/current/cell identity and a defensible alignment key are absent:

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

## Stop point

No SPM, OCP, R0, diffusion, or thermal parameter was fitted; no voltage correction or model-accuracy calculation was performed.
