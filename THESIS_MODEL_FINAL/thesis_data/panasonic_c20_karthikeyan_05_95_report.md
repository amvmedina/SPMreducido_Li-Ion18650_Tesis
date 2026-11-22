# Panasonic C/20 Karthikeyan 5--95% diagnostic

Only samples with `0.05 <= z <= 0.95` entered the voltage objective. The resulting linear stoichiometric trajectories were evaluated without modification over the complete discharge.

## Fit-region metrics (5--95%)

| Metric | 5--95% fit | Previous full-range fit evaluated on 5--95% |
| --- | ---: | ---: |
| RMSE | 37.587 mV | 40.875 mV |
| MAE | 30.293 mV | 34.063 mV |
| R² | 0.97166054 | 0.96648551 |

The fit region contains 1116 of 1241 samples, from `z=0.050011` to `z=0.949451`.

## Full-discharge evaluation

| Metric | 5--95% fitted windows, evaluated on 0--100% | Previous full-range fit |
| --- | ---: | ---: |
| RMSE | 58.262 mV | 53.882 mV |
| MAE | 38.169 mV | 40.303 mV |
| R² | 0.95509960 | 0.96159673 |

- Residual at `z=0`: 139.955 mV
- Residual near `z=0.50`: -57.399 mV
- Residual at `z=1`: -599.002 mV
- Interior 5%-bin mean residual range: -56.486 to 79.924 mV

## Stoichiometric-window comparison

| Endpoint | Previous 0--100% fit | New 5--95% fit | Change |
| --- | ---: | ---: | ---: |
| `x_n_start` | 0.238047229 | 0.224101124 | -0.013946105 |
| `x_n_end` | 0.028850400 | 0.047523802 | +0.018673402 |
| `x_p_start` | 0.401280996 | 0.411225489 | +0.009944493 |
| `x_p_end` | 0.883551177 | 0.885368379 | +0.001817202 |

Both trajectories obey `x_n_start > x_n_end` and `x_p_start < x_p_end`, remain inside the Karthikeyan literature domains over the full discharge, and use no clipping.

## Multistart stability

- Deterministic starts: 40
- Successful feasible solutions: 40
- Solutions within 1.0 mV fit-region RMSE of the best: 2

Every successful solution and its fit-region metrics is recorded in `panasonic_c20_karthikeyan_05_95_multistart.csv`.

## Diagnostic conclusion

Excluding the endpoint regions reduces their leverage, but the fit-region residual remains structured rather than noise-like. The endpoint/cutoff region amplifies the previous full-range error, but it is not the sole source of incompatibility: the Karthikeyan OCP pair retains an interior full-cell shape mismatch.

Panasonic capacity remained fixed and separate at `2.994979138519 Ah`.
