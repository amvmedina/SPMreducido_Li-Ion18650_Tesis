# Panasonic 25 °C 1C heat-generation audit

## Scope

This is a forward-only heat-component audit. No value was optimized; `Cth` and `Gth` were not fitted. The frozen Panasonic electrical parameters, OCP functions, nominal effective stoichiometric windows, capacities, and voltage equations were used unchanged.

- 1C source: `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\C20 OCV and 1C discharge tests_start_of_tests\03-09-17_17.59 3349_Dis1C_1.mat`
- Source columns: `TimeStamp, Voltage, Current, Ah, Wh, Power, Battery_Temp_degC, Time, Chamber_Temp_degC`
- Raw-current sign: negative on discharge; converted once at the model-loader boundary to the internal convention `I > 0` for discharge.
- Frozen cell capacity: `2.994979138519 Ah`.
- Full 1C forward voltage RMSE reproduced before the heat calculation: `114.362188001 mV`.

## Verified literature coefficient

The full-cell entropy coefficient is from He et al. (2024) Supplementary Information, Table S2, which explicitly cites Lu et al. (2019). The specified cell is Panasonic NCR18650PF, Li(NiCoAl)O2/carbon, 2.9 Ah, 2.5–4.2 V. No literature search or coefficient fitting was performed here.

`models/panasonic_entropic.py` stores the supplied table in mV/K and converts it to V/K. `dUeq_dT_panasonic(SOC)` uses linear interpolation only on the inclusive literature domain `0.10 <= SOC <= 0.90`; it raises outside that domain and never clips or extrapolates.

## Data and state preparation

- The complete active discharge was propagated once with actual nonuniform timestamps, measured discharge current, measured battery surface temperature, and frozen dynamics.
- Cell SOC was aligned exactly as in the closed 1C validation using relative raw Ah and frozen `Q_cell_Ah`.
- The audit interval was restricted to exact boundary-interpolated SOC endpoints `0.900000000000` and `0.100000000000`, preserving all 297 original interior timestamps.
- Evaluation points: `299`; duration: `2974.928411972 s` (`49.582140 min`).
- Measured battery temperature: `26.863600` to `30.986661 °C` (`300.013600` to `304.136661 K`).
- Discharge-positive current: `2.899000` to `2.899820 A`.
- Average negative stoichiometry: `0.298709853130` to `0.673590075729`.
- Average positive stoichiometry: `0.296974350671` to `0.880532825535`.
- Average-state OCP domain violations: **0**; clipping: **none**.

For all heat calculations in this audit,

```text
Ueq_average = U_p(x_p_avg) - U_n(x_n_avg)
```

was evaluated from the frozen OCPs. Surface stoichiometries were not used for `Ueq_average`. Its range over the supported interval is `3.332755745` to `4.051940205 V`.

## Heat definitions

```text
Qdot_irr_measured = I * (Ueq_average - V_measured)
Qdot_irr_model    = I * (Ueq_average - V_model)
Qdot_rev          = -I * T_K * dUeq_dT(SOC)
Qdot_total_measured = Qdot_irr_measured + Qdot_rev
Qdot_total_model    = Qdot_irr_model + Qdot_rev
```

The measured-voltage and model-voltage definitions remain separate throughout the CSV, plots, metrics, and energy accounting.

## Heat-component metrics

Means and RMS values are time-weighted using the actual timestamps. Energies are trapezoidal time integrals over exact SOC 0.90–0.10 support.

| Component | Minimum (W) | Maximum (W) | Mean (W) | RMS (W) | Integrated energy (J) |
|---|---:|---:|---:|---:|---:|
| `Qdot_irr_measured` | 0.473774835 | 1.226701312 | 0.593540539 | 0.606765687 | 1765.740613 |
| `Qdot_irr_model` | 0.328061753 | 0.464715815 | 0.389005478 | 0.391077003 | 1157.263449 |
| `Qdot_rev` | -0.147807256 | 0.070290654 | -0.071369783 | 0.090898128 | -212.319994 |
| `Qdot_total_measured` | 0.378540096 | 1.165847343 | 0.522170756 | 0.547158994 | 1553.420619 |
| `Qdot_total_model` | 0.196455607 | 0.503716613 | 0.317635695 | 0.328948882 | 944.943455 |
| `Qdot_irr_measured_minus_model` | 0.065684006 | 0.762217656 | 0.204535061 | 0.231634603 | 608.477164 |

### Integrated heat and reversible fractions

- Measured-voltage irreversible heat: `1765.740613 J`.
- Model-voltage irreversible heat: `1157.263449 J`.
- Net reversible heat: `-212.319994 J`.
- Measured-voltage total heat: `1553.420619 J`.
- Model-voltage total heat: `944.943455 J`.
- Signed reversible/total energy fraction, measured-voltage total: `-0.136679011` (`-13.667901%`).
- Signed reversible/total energy fraction, model-voltage total: `-0.224690687` (`-22.469069%`).
- Absolute reversible/total integral ratio, measured-voltage total: `0.156121307` (`15.612131%`).
- Absolute reversible/total integral ratio, model-voltage total: `0.256652455` (`25.665246%`).

The signed fractions retain the endothermic sign of net reversible energy; the absolute ratios describe reversible heat magnitude without cancellation.

## Physical sign and sign-change audit

For discharge, all evaluated current and absolute temperature values are positive. Numerically:

- at SOC 0.50, `dUeq/dT = 0.169000 mV/K` and `Qdot_rev < 0`: endothermic;
- at SOC 0.20, `dUeq/dT = -0.082000 mV/K` and `Qdot_rev > 0`: exothermic.

Under the specified piecewise-linear interpolation, the exact zero crossings are:

1. `SOC = 0.258571428571` between the 0.30 and 0.20 table points;
2. `SOC = 0.145695364238` between the 0.20 and 0.10 table points.

In discharge order from SOC 0.90 toward 0.10:

- endothermic reversible contribution: `0.90 >= SOC > 0.258571428571`;
- exothermic reversible contribution: `0.258571428571 > SOC > 0.145695364238`;
- endothermic reversible contribution: `0.145695364238 > SOC >= 0.10`.

`Qdot_rev` is zero at the two crossings. Direct evaluations immediately on both sides of each root verify the sign reversals; no sign or monotonicity constraint was imposed.

## Measured-versus-model irreversible heat

The explicitly reported difference is

```text
Qdot_irr_measured - Qdot_irr_model = I * (V_model - V_measured).
```

It ranges from `0.065684006` to `0.762217656 W`, has time-weighted mean `0.204535061 W`, RMS `0.231634603 W`, and integrates to `608.477164 J`. This difference is electrical-model voltage error expressed as heat-rate discrepancy; it was not pushed into any thermal parameter.

## Stop condition

The verified entropy coefficient and forward-only 1C heat audit are complete. No `Cth`, `Gth`, electrical, OCP, Arrhenius, UDDS, or CC–CV thermal calibration/validation was performed.
