# Legacy two-node thermal-model audit

## Scope and provenance

The consolidated `THESIS_MODEL_FINAL` tree did not contain the legacy two-node source. The implementation was located read-only in earlier thesis trees. The canonical source used for this audit is:

- `C:\Users\amvel\Documents\Andes\Tesis\python\models\two_node_thermal.py`
- `C:\Users\amvel\Documents\Andes\Tesis\python\docs\two_node_18650_thermal_model.md`
- `C:\Users\amvel\Documents\Andes\Tesis\datasets\scripts\calibrate_two_node_thermal.py`
- `C:\Users\amvel\Documents\Andes\Tesis\datasets\scripts\calibrate_two_node_thermal_refit_discharge_soc_correction.py`
- `C:\Users\amvel\Documents\Andes\Tesis\datasets\scripts\plot_two_node_thermal_refit_discharge_soc_manual.py`
- `C:\Users\amvel\Documents\Andes\Tesis\datasets\scripts\report_fixed_r0_validation.py`

Byte-identical copies exist under `2026-1/Baterias/Tesis`, `2026-1/Control/Proyecto_Control_Baterias_MPC_Stackelberg`, `backup/Tesis`, and `backup/tesis2/Tesis`. The recorded SHA-256 hashes are:

| Legacy file | SHA-256 |
|---|---|
| `python/models/two_node_thermal.py` | `5EE4081782A7B288B3CC7AC2338422713B9DA738342126F50E07B66633CDCB11` |
| `datasets/scripts/calibrate_two_node_thermal.py` | `95F5B6DE5801291E840AF2EB9133AA2B24589D252310E69E3F827E99D74E9BD4` |
| `datasets/scripts/calibrate_two_node_thermal_refit_discharge_soc_correction.py` | `55BA575310BE21AF3A60EAF8C22040AE3859D4A3EC8E062DDEB39EDC140CD487` |
| `datasets/scripts/plot_two_node_thermal_refit_discharge_soc_manual.py` | `94486D0B1BFD2E8C7E412324E7151A01DDA5DF425A55FF3D7B7AF9DD2B064644` |

No legacy file was modified or copied over the clean implementation.

## Physical architecture

The states are core temperature `T_c` and surface temperature `T_s`. The measured output is `T_s`, compared against `Battery_Temp_degC`; `T_c` is unmeasured.

The exact continuous equations represented by the legacy Euler code are:

```text
C_c dT_c/dt = Q_gen - G_cs (T_c - T_s)
C_s dT_s/dt = G_cs (T_c - T_s) - hA (T_s - T_amb)
y = T_s
```

The four physical parameters are:

| Parameter | Meaning | Unit |
|---|---|---|
| `C_c` | core thermal capacitance | J/K |
| `C_s` | surface thermal capacitance | J/K |
| `G_cs` | core-to-surface conductance (`1/R_cs`) | W/K |
| `hA` | surface-to-ambient conductance (`1/R_sa`) | W/K |

The implementation used previous-sample zero-order-hold inputs with an explicit Euler update. It initialized `T_c(0) = T_s(0) = Battery_Temp_degC[0]`. Measured `Chamber_Temp_degC` was used as ambient when available; otherwise ambient was held at the initial battery temperature. No initial core/surface offset was fitted.

## Legacy fitting and datasets

`calibrate_two_node_thermal.py` fitted all four physical parameters. None was fixed from geometry or literature. Its stated guard bounds were `C_c=[20,100] J/K`, `C_s=[5,50] J/K`, `G_cs=[0.05,2] W/K`, and `hA=[0.005,0.5] W/K`. It used five manual seeds, 250 deterministic random candidates (`seed=42`), and coordinate search—not a least-squares sensitivity/identifiability analysis.

The objective gave equal weight to the temperature RMSE of:

- 1C discharge: `datasets/processed/discharge_1C_25C_trimmed.csv`;
- charge/rest: `datasets/converted_csv/charge_3349_2_converted.csv`.

It additionally penalized late 1C underprediction, core/surface separation above 2 degC, and violations of `T_c >= T_s` during heating. UDDS (`datasets/converted_csv/udds_25C_converted.csv`) was held out for validation. The later `report_fixed_r0_validation.py` applied the selected thermal model to C/20, 1C discharge, charge/rest, HPPC, and UDDS, but that was a reporting/validation stage rather than the original physical-parameter fit.

The baseline refit output records `C_c=78.0436825 J/K`, `C_s=23.59375 J/K`, `G_cs=0.466711436 W/K`, and `hA=0.105648286 W/K`. The subsequently hard-coded selected values were `C_c=79.266263 J/K`, `C_s=21.462628 J/K`, `G_cs=0.422606 W/K`, and `hA=0.105224 W/K`; the archived files do not document a separate reproducible fit that explains the small change from the baseline-refit CSV to these selected values.

## Legacy heat generation and nonphysical processing

The original documentation proposed `abs(I*(V-OCV))`, ignored entropic heat, and treated an entropy term as optional future work. The implemented scripts instead:

1. formed `I*(V_model-OCV_model)` using the then-current electrical model;
2. automatically selected a global sign multiplier `+1` or `-1` according to which choice made more active-current samples nonnegative;
3. applied `max(Q_signed, 0)` by default, clipping all negative heat to zero.

This differs from the now-verified Panasonic formulation. It used model voltage rather than measured voltage, omitted verified reversible heat, selected the heat sign empirically, and clipped negative heat.

## Every located empirical correction

The legacy module contains the global SOC multiplier

```text
lambda_SOC(SOC) = 1 + a * (1 - clip(SOC, 0, 1))^b
Q_corrected = lambda_SOC * Q_gen
```

with constraints `a >= 0` and `b >= 1`. It can be applied in two modes:

- globally to all samples; or
- only when model-convention current exceeds `0.05 A` (discharge-only/regime-dependent correction).

The correction-refit script varied `a` over `[0,7]` and `b` over `[1,12]` on the 1C discharge, using an objective containing RMSE plus late-underprediction, a 2750–3500 s underprediction penalty, and a core/surface-separation penalty. The officially selected plotting/reporting branch instead hard-coded `a=6`, `b=10` and applied it only during discharge. At SOC 0 the selected multiplier is 7; at SOC 1 it is 1.

No separate charge-specific multiplier was located. No additive empirical heat source, polynomial residual correction, fitted temperature bias, or initial-temperature offset was located. The nonphysical operations that did exist were the automatic heat-sign selection, negative-heat clipping, SOC clipping inside `lambda_SOC`, the SOC multiplier itself, and its discharge-only current-regime gate.

## Historical performance and interpretation limit

For the hard-coded selected physical parameters, the archived summary reports:

| Dataset | No SOC multiplier RMSE [degC] | Discharge-only `a=6,b=10` RMSE [degC] |
|---|---:|---:|
| 1C discharge | 0.727691 | 0.655786 |
| charge/rest | 0.473909 | 0.473909 |
| UDDS | 0.204643 | 0.222634 |

Thus the multiplier improved the archived 1C result but worsened the archived UDDS result. The low archived UDDS error did not depend on that multiplier. These historical numbers are not directly comparable to the current frozen Panasonic branch because they used different electrical/OCP parameters, model-voltage heat, automatic sign selection, negative-heat clipping, no verified entropy heat, different interval handling, and explicit Euler integration.

The legacy audit itself correctly noted that only surface temperature was measured and that the internal state and four thermal parameters were not independently observed. No structural-identifiability claim is supported by the legacy work.
