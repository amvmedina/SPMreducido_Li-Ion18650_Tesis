# Panasonic HPPC low-SOC forward validation

## LOW-SOC TRANSFER LIMITED

Forward evaluation only. The nominal conditional parameters from the `SOC >= 1/3` diagnostic are used unchanged: `R0=0.0259594900586 ohm`, `alpha_p=3204.22765511 s`, and `d_n=0.000183029692927 1/s`. `alpha_n=3000 s` and `d_p=0.0007 1/s` remain fixed. No parameter was optimized.

- Withheld complete low-SOC pulses: **24** in 6 independently initialized histories.
- Eligibility: pulse-start cell `SOC < 1/3`; cutoff-truncated IDs 60, 64, and 67 remain excluded.
- Cumulative Ah initializes each history; omitted gaps are not integrated, raw negative current is converted to discharge-positive, and actual nonuniform time steps are retained.

## Dynamic metrics

The same pulse-relative, pulse-balanced, equal-phase metric used for calibration is applied without refitting.

- Phase-balanced dynamic RMSE: **85.729433 mV**.
- Calibration-region (`SOC >= 1/3`) RMSE: 36.835198 mV.
- Transfer/calibration RMSE ratio: 2.327378.

| Response phase | RMSE [mV] |
|---|---:|
| Immediate | 22.531228 |
| Within pulse | 145.725371 |
| Post relaxation | 17.466164 |

### By C-rate

| C-rate | Phase-balanced RMSE [mV] |
|---|---:|
| 0.5C | 39.547954 |
| 1C | 85.022724 |
| 2C | 99.664359 |
| 4C | 103.987157 |
| 6C | 101.181673 |

Voltage-step bias (`measured - predicted`): +15.100181 mV.

## Absolute-voltage diagnostic

- RMSE: 59.582115 mV.
- MAE: 40.733681 mV.
- R2: 0.801454959.

## Surface-state and OCP-domain audit

- `x_n,surf` range: 0.281877086592 to 0.402832594302.
- `x_p,surf` range: 0.718449833759 to 0.907127079394.
- Frozen OCP-domain violation: **NO**.
- No state or stoichiometry was clipped.

## Interpretation

The withheld-region error is 2.327 times the calibration-region phase-balanced RMSE. The largest diagnostic discrepancy is in the **loaded transient** (145.725 mV RMSE).

**Final classification: LOW-SOC TRANSFER LIMITED.**

No repair, refit, additional calibration, parameter release, OCP/support change, model change, UDDS run, or 1C validation was performed.
