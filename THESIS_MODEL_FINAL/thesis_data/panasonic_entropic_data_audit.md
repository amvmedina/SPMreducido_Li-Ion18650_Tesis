# Panasonic entropic-coefficient and thermal-data audit

## Decision

**PANASONIC ENTROPIC COEFFICIENT NOT IDENTIFIABLE FROM AVAILABLE DATA**

The local Panasonic tree contains one genuine C/20 experiment at one controlled temperature: the reference file has a 25 °C discharge and a subsequent 25 °C charge. Direction mixing is prohibited, so these two branches do not provide temperature replication. No second comparable low-rate discharge or charge curve exists at another controlled temperature. Consequently, no `dUeq/dT(SOC)` values were estimated, regressed, smoothed, plotted, or implemented.

No electrochemical model was run, no parameter was fitted, and no OCP or frozen electrical file was modified.

## Search coverage and reproducible candidate definition

- Dataset root: `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data`
- Dataset README inspected: `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Readme file - desc of tests performed.txt`
- MAT files inspected from their signals: **198 of 198**, with zero load failures.
- CSV files found: **59**; all are inside `EIS` folders and are impedance spectra, not low-rate capacity curves.
- Reference C/20 file: `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\C20 OCV and 1C discharge tests_start_of_tests\05-08-17_13.26 C20 OCV Test_C20_25dC.mat`
- Signals required from each MAT: `Time`, `Current`, `Voltage`, `Ah`, `Battery_Temp_degC`, and `Chamber_Temp_degC`.
- Active sign threshold: |I| > 0.020 A. Contiguous positive and negative current runs were examined independently.
- Deliberately broad candidate gate: at least 1.0 h, at least 0.10 Ah transferred, and median |I| no greater than 0.35 A.
- A curve is marked usable only when it is near-constant current, covers at least 75% of frozen `Q_cell_Ah = 2.994979138519 Ah`, has controlled temperature, a consistent Ah counter, and a broad voltage range.

The dataset README independently describes a C/20 discharge/charge only once, before HPPC/EIS/drive-cycle testing at 25, 10, 0, -10, and -20 °C. The conclusions below nevertheless come from the electrical and temperature signals rather than filenames alone.

## All signal-level low-current candidates

Indices are zero-based and inclusive. Current magnitude is shown as minimum/median/maximum over the active branch.

| ID | Exact path | Folder temperature | Direction | Indices | abs(I) min/median/max (A) | Duration (h) | Integrated Ah | Ah-counter Δ | Voltage range (V) | Most SOC? | Long rest in file? | Ah usable? | Disposition |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|---|
| C1 | `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\-10degC\Charges and Pauses\06-14-17_08.04 3787_Charge8.mat` | -10 °C | charge | 107–178 | 0.04982/0.27399/2.89997 | 1.175344 | 0.825877524 | 0.826250000 | 4.06431–4.20007 | no | yes (1.764 h max) | yes | variable-current/CC–CV segment; insufficient capacity coverage; missing or changing temperature; limited voltage/SOC interval |
| C2 | `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\-20degC\Charges and Pauses\06-23-17_23.36 3913_Charge2a.mat` | -20 °C | charge | 100–163 | 0.04982/0.17436/2.26543 | 1.041152 | 0.393884805 | 0.393090000 | 4.19942–4.20007 | no | yes (1.641 h max) | yes | variable-current/CC–CV segment; insufficient capacity coverage; missing or changing temperature; limited voltage/SOC interval |
| C3 | `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\C20 OCV and 1C discharge tests_start_of_tests\03-09-17_15.13 3349_Charge1.mat` | 25 °C | charge | 51–157 | 0.04982/0.22132/2.89997 | 1.758321 | 1.662763649 | 1.662920000 | 3.90668–4.20007 | no | no (0.825 h max) | yes | variable-current/CC–CV segment; insufficient capacity coverage; missing or changing temperature; limited voltage/SOC interval |
| C4 | `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\C20 OCV and 1C discharge tests_start_of_tests\05-08-17_13.26 C20 OCV Test_C20_25dC.mat` | 25 °C | discharge | 6–1246 | 0.14454/0.14536/0.14536 | 20.661352 | 2.994979139 | 2.994910000 | 2.49948–4.17030 | yes | yes (14.586 h max) | yes | usable low-rate reference |
| C5 | `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\C20 OCV and 1C discharge tests_start_of_tests\05-08-17_13.26 C20 OCV Test_C20_25dC.mat` | 25 °C | charge | 1308–2390 | 0.14455/0.14537/0.14537 | 18.031703 | 2.613916708 | 2.613900000 | 2.92679–4.20007 | yes | yes (14.586 h max) | yes | usable low-rate reference |

### Candidate temperature and current-stability detail

| ID | Chamber °C range (mean) | Battery °C range (mean) | Chamber samples finite | abs(I) P5–P95 (A); coefficient of variation |
|---|---:|---:|---:|---:|
| C1 | missing–missing (mean missing) | 11.150–16.779 (mean 15.383) | 0.0% | 0.05770–2.89997; CV=1.254 |
| C2 | missing–missing (mean missing) | 11.173–15.481 (mean 13.371) | 0.0% | 0.05892–1.45247; CV=1.234 |
| C3 | 12.000–24.000 (mean 19.168) | 11.393–24.151 (mean 20.009) | 100.0% | 0.06787–2.89997; CV=1.224 |
| C4 | 25.000–25.000 (mean 25.000) | 25.014–26.090 (mean 25.636) | 100.0% | 0.14454–0.14536; CV=0.003 |
| C5 | 24.000–25.000 (mean 24.999) | 24.801–25.463 (mean 25.099) | 100.0% | 0.14455–0.14537; CV=0.003 |

The three rejected candidates are not hidden low-rate OCV curves. Their current ranges extend into ordinary charge rates and then taper, their current coefficients of variation exceed 1.2, and they cover only restricted charge intervals. Two have no usable chamber-temperature channel; all three undergo substantial battery-temperature changes. Comparing them with the 25 °C C/20 branches would mix rate, direction/protocol, temperature history, and SOC coverage.

## Usable low-rate curves

| Exact path | Direction | Nominal temperature | Chamber °C | Battery °C | Current | Duration (h) | Transferred capacity | Voltage range (V) | Most SOC? | Rest context | Ah usable? |
|---|---|---|---:|---:|---:|---:|---:|---:|---|---|---|
| `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\C20 OCV and 1C discharge tests_start_of_tests\05-08-17_13.26 C20 OCV Test_C20_25dC.mat` | discharge | 25 °C | 25.000–25.000 | 25.014–26.090 | 0.14454–0.14536 A (discharge) | 20.661352 | 2.994979139 (100.000% Qcell) | 2.49948–4.17030 | yes | initial 0.067 h, intermediate 1.000 h, final recorded 14.586 h | yes |
| `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\C20 OCV and 1C discharge tests_start_of_tests\05-08-17_13.26 C20 OCV Test_C20_25dC.mat` | charge | 25 °C | 24.000–25.000 | 24.801–25.463 | 0.14455–0.14537 A (charge) | 18.031703 | 2.613916708 (87.277% Qcell) | 2.92679–4.20007 | yes | initial 0.067 h, intermediate 1.000 h, final recorded 14.586 h | yes |

The discharge is the previously used reference. Its chamber is exactly 25 °C and its current-integrated capacity is 2.994979139 Ah. The charge is also approximately C/20 at 25 °C, but it reaches only 2.613916708 Ah (87.277% of frozen cell capacity) before the 4.2 V cutoff. The branches are useful for the already documented direction/capacity comparison, not for a temperature derivative.

Long rest files at 0, -10, -20, and temperature-rise conditions were also encountered. They are isolated rest/conditioning records near limited states of charge, often with missing chamber values or large temperature transients. They do not supply full-SOC low-rate voltage curves and cannot replace the missing multi-temperature discharge series.

## Common-SOC and derivative gate

For any future discharge comparison, the required primary mapping remains

`SOC = 1 - Q_discharge_since_start / Q_cell_Ah`,

using actual transferred charge and `Q_cell_Ah = 2.994979138519 Ah`. Charge would use its own direction-consistent coulomb count but must not be mixed into a discharge temperature regression. Interpolation would be allowed only on the intersection of measured SOC supports, with no endpoint stretching and no extrapolation.

Here there is only one usable discharge temperature and one usable charge temperature, both 25 °C. Therefore:

- number of temperatures available for a discharge-only regression: **1**;
- temperature span: **0 K**;
- common multi-temperature SOC interval: **none**;
- per-SOC slope, R², and temperature count: **not computed**.

The acceptance conditions fail before regression. A slope from these files would necessarily mix charge/discharge direction or unrelated dynamic/variable-temperature experiments and would be dominated by protocol and capacity-axis differences. No coefficient is invented and no literature coefficient is substituted.

## Heat-generation implementation audit

The existing implementation already declares discharge-positive model current (`I > 0` discharge, `I < 0` charge) and contains the requested signs and SI-unit structure:

- `C:\Users\amvel\Documents\Andes\THESIS_MODEL_FINAL\models\thermal.py`, line 98: `Qdot_irr = I * (Ueq - V)`.
- `C:\Users\amvel\Documents\Andes\THESIS_MODEL_FINAL\models\thermal.py`, line 101: `Qdot_rev = -I * T * dUeq_dT`.
- `C:\Users\amvel\Documents\Andes\THESIS_MODEL_FINAL\models\thermal.py`, line 102: `Qdot_gen = Qdot_irr + Qdot_rev`.

The generic heat function can evaluate either irreversible quantity depending on the voltage supplied, but the caller must label them separately:

```text
Qdot_irr_measured = I * (Ueq_average - V_measured)
Qdot_irr_model    = I * (Ueq_average - V_model)
Qdot_rev          = -I * T * dUeq_dT(SOC)
Qdot_gen          = Qdot_irr + Qdot_rev
```

The coupled wrapper currently passes its predicted terminal voltage at `C:\Users\amvel\Documents\Andes\THESIS_MODEL_FINAL\models\electrothermal_spm.py`, line 205; its exported `heat_irreversible_w` is therefore the **fully predictive model value**, not the measured-voltage calibration diagnostic. There is no separate measured-voltage heat output, so a later thermal-calibration script must calculate and retain `Qdot_irr_measured` under that explicit name rather than silently replacing the model value.

There is also an important state-definition finding. The current wrapper forms `Ueq` at `C:\Users\amvel\Documents\Andes\THESIS_MODEL_FINAL\models\electrothermal_spm.py`, line 183 from OCP arrays that are evaluated using **surface stoichiometries** (`C:\Users\amvel\Documents\Andes\THESIS_MODEL_FINAL\models\grouped_spm.py`, line 342). The thermal specification in this task instead requires equilibrium voltage from average electrode state:

```text
Ueq_average = U_p(x_p_avg) - U_n(x_n_avg)
```

That distinction must be explicit in the later Panasonic thermal calibration/validation interface. It was not changed here because the electrical architecture is frozen and this task authorizes an audit, not a model-equation change.

Reversible heat remains disabled when no coefficient is supplied. Because the identifiability gate failed, no Panasonic entropic module, coefficient CSV, coefficient plot, or entropic tests were created.

## Final scope statement

This audit stops at the required gate. `Cth`, `Gth`, electrical parameters, OCPs, stoichiometric windows, Arrhenius terms, and validation profiles remain untouched.
