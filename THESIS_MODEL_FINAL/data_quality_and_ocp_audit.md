# Data quality and OCP provenance audit

Audit date: 2026-09-12  
Scope: evidence review only. No parameter was optimized, no OCP equation was
imported, and no frozen electrochemical or thermal model file was modified.

## Executive decision

- The ambient cell-1 cycle is the only own experiment classified **likely
  clean** by the common rest-to-load diagnostic.
- Pulse cell-1 parts 1 and 2, 30 degC cell 1, 35 degC cell 1, and both cell-2
  cycles are classified **likely affected by setup resistance**. Their
  apparent transition resistance is 5.1--6.6 times the ambient cell-1 value.
- The pulse experiment is **rejected for electrical calibration**. Its repeated
  transitions consistently contain the anomalous voltage drop; it is not an
  isolated outlier.
- `models/ocp.py` exactly reproduces the functions in Feng Guo's public
  `FrankSuperG/CPG-SPMT` software at the recorded Git blob. Project records
  describe the associated electrodes as NMC positive and graphite negative,
  but the functions do **not** match the OCP pair printed in Appendix A of the
  2025 Guo and Couto paper.
- No available evidence establishes that those electrode curves are compatible
  with either the own TRC18650 cell or the Panasonic 18650PF cell. The common
  18650 form factor is not chemistry evidence.
- The Panasonic C/20 trace is a full-cell quasi-OCV curve. It cannot by itself
  supply the two separate functions `U_p(x_p)` and `U_n(x_n)` required by the
  SPM.

This audit supersedes the provisional recommendations in
`experimental_data_manifest.md` to use the own pulse test for first calibration
and the 35 degC cycle for transfer. The resistance evidence below makes both
recommendations indefensible.

## A. Own-data quality audit

### Common transition definition

The same purely local definition was applied to every own CSV:

1. select an adjacent transition whose preceding sample has
   `abs(I_raw) <= 0.02 A`;
2. require the following sample to have `I_raw <= -2.5 A` (the approximately
   1C discharge level is about -2.9 A in the Chroma convention);
3. calculate

   `R_app = (V_rest - V_loaded) / abs(I_loaded)`.

This uses the last recorded rest sample and the first recorded loaded sample.
It neither reconstructs missing time nor extrapolates voltage to zero time.
`R_app` is the apparent resistance of the complete observed electrical path at
that transition. **It is not a measurement or estimate of true cell internal
resistance, and it must not be used as model `R0`.**

### Transition results

| Experiment | Number of qualifying transitions | `R_app` values [ohm] | Summary [ohm] | Classification |
| --- | ---: | --- | --- | --- |
| Ambient cell 1 cycle | 1 | 0.036383 | 0.036383 | **Likely clean** relative to this experimental set |
| Pulse cell 1, part 1 | 7 | 0.235911, 0.232522, 0.223341, 0.218928, 0.221135, 0.221349, 0.218787 | median 0.221349; range 0.218787--0.235911 | **Likely affected by setup resistance** |
| Pulse cell 1, part 2 | 3 | 0.218014, 0.237146, 0.228939 | median 0.228939; range 0.218014--0.237146 | **Likely affected by setup resistance** |
| 30 degC cell 1 cycle | 1 | 0.240607 | 0.240607 | **Likely affected by setup resistance** |
| 35 degC cell 1 cycle | 1 | 0.237534 | 0.237534 | **Likely affected by setup resistance** |
| Cell 2 cycle, condition undeclared | 1 | 0.185262 | 0.185262 | **Likely affected by setup resistance**; cell-to-cell differences prevent a cleaner attribution |
| 30 degC cell 2 cycle | 1 | 0.226547 | 0.226547 | **Likely affected by setup resistance** |

The directly comparable fully charged cycle starts are especially clear:

| Experiment | `V_rest` [V] | `V_loaded` [V] | `I_loaded` [A] | `R_app` [ohm] |
| --- | ---: | ---: | ---: | ---: |
| Ambient cell 1 | 4.1920 | 4.0865 | -2.8997 | 0.036383 |
| 30 degC cell 1 | 4.1824 | 3.4850 | -2.8985 | 0.240607 |
| 35 degC cell 1 | 4.1853 | 3.4965 | -2.8998 | 0.237534 |
| Cell 2, condition undeclared | 4.1853 | 3.6483 | -2.8986 | 0.185262 |
| 30 degC cell 2 | 4.1872 | 3.5306 | -2.8983 | 0.226547 |

These starts have nearly equal rest voltage and load current, while their
voltage steps differ by approximately a factor of five to seven. Temperature
or SOC therefore cannot plausibly explain the full separation. The cell-2
results cannot prove where the extra resistance lies because they use a
different physical cell, but they fall in the same anomalous band and are not
safe electrical-identification data.

The pulse result is stronger than a one-step comparison: all ten qualifying
steps, spanning both files and multiple states of charge, remain between
0.218014 and 0.237146 ohm. That repeatability is evidence of a persistent path
condition. Consequently:

> **Pulse cell-1 parts 1 and 2 are explicitly rejected for electrical
> calibration.**

They are also rejected for coupled thermal calibration. If the extra voltage
drop is in leads, contacts, or the tester path, an unknown portion of the
associated electrical loss occurs outside the cell. Assigning that loss to the
cell heat source would bias the thermal parameters. The same restriction
applies to the other setup-affected own cycles.

The ambient cell-1 value is only "likely clean" relative to the rest of this
set. It is not independently certified. The CSV also lacks elapsed time, so it
is not presently suitable for dynamic parameter identification. Its documented
temperature channels make it a candidate for thermal-only validation after the
sample time is recovered from existing acquisition records.

## B. OCP provenance audit

### Exact provenance of the frozen functions

`models/ocp.py` records both of the following:

- source repository/file:
  `https://github.com/FrankSuperG/CPG-SPMT/blob/main/cpg_spmt.py`;
- Git blob/commit identifier:
  `2b396b93713a9ab908fca27f4e4ff3b9a9bf9343`.

The coefficients and algebra in both `U_p` and `U_n` match the archived source
`C:\Users\amvel\Documents\Andes\Tesis\third_party\guo_cpg_spmt\cpg_spmt.py`
and the project's traceability copy
`C:\Users\amvel\Documents\Andes\Tesis\python\models\guo_grouped_spm_ocp.py`.
The public source identifies Feng Guo as copyright holder and asks users to cite
the related publications.

The traceability copy and
`C:\Users\amvel\Documents\Andes\Tesis\docs\GUO_SYNTHETIC_REPRODUCTION.md`
explicitly state that these functions **are not** the polynomial/exponential
OCP pair printed in Appendix A of Guo and Couto, *Journal of Power Sources* 649
(2025) 237309. The software functions were retained to reproduce the software,
not to claim exact reproduction of the paper equations.

Therefore the defensible provenance statement is:

> The exact equations currently frozen in the thesis model were copied from
> Guo's public CPG-SPMT software implementation at the recorded blob. Available
> project evidence does not establish the upstream experimental fit or a named
> earlier publication from which those exact coefficient sets originated.
> Because they differ from the 2025 paper appendix, they must not be described
> as equations transcribed directly from that appendix.

### Electrode chemistry represented

The source-code functions themselves are labelled only positive and negative
electrode OCPs. The repository's model documentation,
`C:\Users\amvel\Documents\Andes\Tesis\python\docs\Guo_Grouped_NMC_SPM.md`,
describes the grouped model as **NMC/graphite** and its literature initial values
as belonging to a 2.9 Ah NMC 18650 reference. On that documented basis:

- `U_p(x_p)`: NMC-family positive-electrode OCP;
- `U_n(x_n)`: graphite negative-electrode OCP.

This chemistry label is inherited from the Guo model documentation. It is not
a new chemical identification of either target cell.

### Compatibility with the two target cells

**Own TRC18650 cell.** No manufacturer model, cathode composition, anode
composition, electrode half-cell curve, or electrode stoichiometric operating
window was found. A 2.8999 Ah measured capacity and an 18650 can size do not
establish NMC/graphite compatibility. There is **no evidence** that the frozen
`U_p` and `U_n` are cell-compatible.

**Panasonic 18650PF cell.** The external dataset supplies high-quality
full-cell voltage/current data and a slow C/20 full-cell quasi-OCV curve, but no
separate positive/negative half-cell OCP measurements or target-electrode
stoichiometric windows. The project documentation explicitly says the Guo
values are not exact Panasonic parameters and that transfer to Panasonic was
not evaluated or supported. There is **no evidence** that the frozen pair is
Panasonic-compatible.

The C/20 trace can later test whether an assumed OCP pair and stoichiometric
windows reproduce the full-cell voltage shape, but a single full-cell curve
does not uniquely decompose into `U_p` and `U_n`.

### Previously considered local OCP sources

No available local source was found that is a verified TRC18650 or Panasonic
18650PF electrode-resolved OCP parameterization. The following candidates were
previously considered:

| Local source | Curves/chemistry | Relevance and limitation |
| --- | --- | --- |
| `python/docs/Khan_NMC_eSPM.md` and `python/data/chen2020_ocp_functions.py` | NMC positive and graphite negative; coefficients match the locally stored Chen2020 functions | Chemically plausible only at family level. The local provenance record identifies the source cell as an LG M50 21700 NMC811/graphite cell, not either target cell. |
| `python/data/panasonic_graphite_ocp.py` | Chen2020 graphite and Ecker2015 graphite alternatives | The file itself calls them sensitivity assumptions, not Panasonic fits. Chen2020 is LG M50 21700; Ecker2015 is a Kokam NMC/graphite pouch cell. Neither supplies a Panasonic-specific positive curve. |
| `python/docs/KHASIN_REDUCED_ELECTROCHEMICAL_MODEL_EVALUATION.md` | Khasin six-logistic OCP family associated with a LiCoO2/graphite cell | Required fitted coefficients are not published locally. The evaluation explicitly rejected a physical Panasonic interpretation and showed inferred diffusion depends on the assumed graphite OCP. |
| `python/docs/18650_lfp_parameters_from_saw2013.md` | LFP positive and graphite negative for a 1.3 Ah 18650 | Wrong documented cathode family and capacity basis unless a target cell were independently shown to be LFP; 18650 geometry alone is insufficient. |
| `python/docs/ferreira_tang_equations.md` | Positive and negative equilibrium-potential equations as functions of surface concentration | A prior equation source, but the available project notes do not establish complete target-cell parameterization or compatibility with TRC/Panasonic. |
| `datasets/converted_csv/c20_ocv_25C_converted.csv`, `python/data/panasonic_18650pf_ocv.csv`, and `python/data/panasonic_18650pf_quasi_ocv.csv` | Panasonic **full-cell** low-current OCV/quasi-OCV | Cell-specific voltage evidence, but not separable electrode functions and therefore not a direct replacement for `U_p` and `U_n`. |

No new OCP is selected or imported by this audit.

## C. Decision table

"Conditional" means the named unresolved prerequisite must be completed before
use; it does not authorize calibration now.

| Dataset / experiment | Electrical quality | Temperature quality | Capacity information | OCP compatibility | Usable for electrical calibration? | Usable for thermal calibration? | Usable for independent validation? | Reason |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Own ambient cell-1 cycle | **Likely clean**; `R_app=0.036383 ohm` | Good channel roles: center surface, end surface, ambient; no elapsed time | Complete discharge and charge each 2.8999 Ah | Not established for TRC | **No now** | **No now** | **Conditional: thermal-only** after timebase recovery | Best own run, but constant-current excitation, missing time, and unresolved OCP prevent electrical identification. |
| Own pulse cell 1, part 1 | **Likely setup-affected**; median 0.221349 ohm | Good sensor roles; no elapsed time | Partial test only; joined pulse removal 2.2497 Ah | Not established for TRC | **No -- rejected** | **No** | **No** for quantitative electrothermal validation | Seven repeated anomalous steps plus unknown partition of parasitic heat. |
| Own pulse cell 1, part 2 | **Likely setup-affected**; median 0.228939 ohm | Good sensor roles; no header/time; boundary continuity does not prove no time gap | Continuation of same partial test | Not established for TRC | **No -- rejected** | **No** | **No** for quantitative electrothermal validation | Three further anomalous steps confirm the condition persists across both parts. |
| Own 30 degC cell-1 cycle | **Likely setup-affected**; 0.240607 ohm | Center/end/ambient channels present; useful measurements but heat input is confounded | Discharge 2.4894 Ah; charge 2.4653 Ah | Not established for TRC | **No** | **No** | **No** | Near-equal start OCV/current to ambient run but about 6.6 times its apparent step resistance. |
| Own 35 degC cell-1 cycle | **Likely setup-affected**; 0.237534 ohm | Center/end/ambient channels present; heat input is confounded | Discharge 2.6171 Ah; charge 2.5185 Ah | Not established for TRC | **No** | **No** | **No** | About 6.5 times the ambient cell-1 diagnostic; unsuitable as the proposed temperature-transfer test. |
| Own cell-2 cycle, condition undeclared | **Likely setup-affected**; 0.185262 ohm | Temp3/Temp4 locations unknown; no ambient or elapsed time | Discharge 2.7133 Ah; charge 2.6897 Ah | Not established for cell 2 | **No** | **No** | **No** | Different cell prevents direct attribution, but result is about 5.1 times ambient cell 1 and metadata are incomplete. |
| Own 30 degC cell-2 cycle | **Likely setup-affected**; 0.226547 ohm | Temp3/Temp4 locations unknown; no ambient or elapsed time | Discharge 2.5994 Ah; charge 2.5054 Ah | Not established for cell 2 | **No** | **No** | **No** | Same anomalous electrical family and insufficient thermal metadata. |
| Panasonic 25 degC C/20 | Good reported current/voltage/time; slow full-cell trace | Case/chamber signals exist, but C/20 gives little thermal excitation | Measured C/20 capacity **2.994979138519 Ah** | Separate electrode OCPs absent | **No now**; use as OCV/capacity evidence | **No** | No if used to establish the OCP assumption | Cannot uniquely recover two electrode curves from one full-cell voltage curve. |
| Panasonic 25 degC HPPC | Good dynamic electrical excitation and exact time; long gaps must be segmented | Case/chamber measured; modest temperature span and long inter-group gaps | Use Panasonic-specific 2.994979138519 Ah SOC basis | Frozen Guo pair not yet justified | **Conditional** after explicit OCP-assumption gate | Not the minimal thermal choice | No if used for calibration | Best available electrical-identification profile once OCP limitation is accepted and documented. |
| Panasonic initial 1C discharge / associated CC-CV charge | Good current/voltage/time | Good case/chamber channels and useful temperature rise | 1C discharge integrates to 2.798235809 Ah; retain 2.994979138519 Ah C/20 SOC basis separately | Frozen Guo pair not yet justified | **Conditional**, but reserve from HPPC fit where possible | **Conditional** after the same OCP gate | Conditional if kept entirely outside fitting | Suitable steady thermal excitation, but not a source of electrode-resolved OCP. |
| Panasonic 25 degC UDDS | Good bidirectional dynamic current/voltage/time | Case/chamber measured | Use Panasonic-specific 2.994979138519 Ah SOC basis | Frozen Guo pair not yet justified | **No** (reserve it) | **No** (reserve it) | **Conditional -- preferred independent validation** | Dynamic profile should remain untouched by fitting; valid only after a transparent Panasonic-specific assumed-OCP workflow. |
| Maryland archives | Uncertain | Uncertain | Cell identity/capacity not established locally | Unknown | **No** | **No** | **No** | Missing source metadata, cell identity, capacity, and sensor definitions. |

For every own TRC cell-level SOC calculation retain:

`Q_cell_Ah = 2.8999 Ah`

For Panasonic-specific analysis retain separately:

`Q_cell_Ah = 2.994979138519 Ah`

The values must not be averaged, pooled, or transferred between cells.

## D. One minimum defensible path

The minimum path is a **Panasonic-specific quantitative branch plus a narrowly
labelled own-data thermal check**:

1. Exclude pulse cell 1, 30/35 degC cell 1, and both cell-2 cycles from all
   electrical identification. Do not use their voltage steps as cell `R0`, and
   do not use them to identify coupled thermal parameters.
2. Treat the frozen Guo-software NMC/graphite `U_p` and `U_n` as an **explicit
   literature approximation**, not as Panasonic electrode data. Before any
   optimization, make the future go/no-go check whether a physically admissible
   pair of electrode operating windows can reproduce the main Panasonic C/20
   full-cell voltage shape and endpoints. Failure of that check means the
   frozen SPM cannot support defensible Panasonic electrical identification;
   it must not be hidden by fitting kinetic or resistance parameters.
3. If that non-calibration OCP compatibility gate is accepted, use only the
   Panasonic branch for quantitative model work: C/20 capacity
   `2.994979138519 Ah`, HPPC for electrical calibration, the 1C run for thermal
   calibration, and untouched UDDS for independent electrothermal validation.
   Report the identified electrical quantities as effective parameters
   conditional on the assumed literature OCP, not electrode-physical truths.
4. Use the own ambient cell-1 run only as a separate, thermal-only check after
   its existing acquisition timebase is recovered. Keep `2.8999 Ah`, do not
   transfer Panasonic capacity or fitted electrical parameters to it, and call
   any cross-cell comparison qualitative unless TRC-specific electrical/OCP
   compatibility is independently established.

This path requires no new laboratory experiment and no model-architecture
change. It minimizes claims: one setup-clean own run is retained only where its
sensors are informative, setup-contaminated experiments are excluded, and the
external electrical result is explicitly conditional on an unverified
literature OCP approximation.

## Stop point

The audit ends here. No calibration or optimization has been performed.
