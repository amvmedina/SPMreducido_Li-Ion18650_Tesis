# Panasonic joint thermal boundary/profile audit

## JOINT OPTIMUM RESOLVED

This is a numerical-boundary diagnostic of the unchanged clean one-node model. Only the Cth numerical guard was expanded from 200 to 1000 J/K; Gth retained its 1e-4-to-1 W/K guard. The objective, datasets, heat generation, entropy support, HPPC block handling, and withheld UDDS procedure were unchanged.

## Expanded-bound joint optimum

- Effective lumped `Cth = 363.045556942567 J/K`
- `Gth = 0.018701208704 W/K`
- `tau_th = 19412.946119219945 s`
- Joint objective: `0.126164862688 K^2`
- 1C RMSE: `0.393235843 K`; HPPC RMSE: `0.312562469 K`
- Cth/Gth nearest log-bound fractions: `0.191235610`, `0.432032581`
- Successful multistarts: 40/40; near-optimal: 40. Successful ranges: Cth `363.044556709-363.047166679 J/K`, Gth `0.0187005863924-0.0187015831054 W/K`. Near-optimal ranges: Cth `363.044556709-363.047166679`, Gth `0.0187005863924-0.0187015831054`.
- HPPC continuity handling remained fixed at 11 blocks and 81,637 samples; no omitted gap was propagated.
- **Profile outcome A:** Cth converges to an interior optimum; it neither reaches the new 1000 J/K guard nor remains practically flat over the sampled high-Cth region.

Cth is an **effective lumped thermal capacitance**, not a measured material heat capacity. A very large effective value can compensate for incomplete heat-generation representation, boundary-condition mismatch, experiment-to-experiment differences, or lack of transferability of one common lumped parameter set. None of those effects was repaired here.

## Fixed-Cth objective profile

| Fixed Cth [J/K] | Optimal Gth [W/K] | tau [s] | Joint objective [K^2] | 1C RMSE [K] | HPPC RMSE [K] | Withheld UDDS RMSE [K] | UDDS bias [K] | UDDS max abs [K] |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 50.000000000 | 0.128006537251 | 390.605051 | 0.275243092092 | 0.382761216 | 0.635594238 | 0.586201058 | -0.578973789 | 0.848150315 |
| 75.000000000 | 0.119871471100 | 625.670139 | 0.229953423325 | 0.237520109 | 0.635209449 | 0.560612808 | -0.551727244 | 0.837550715 |
| 100.000000000 | 0.112306420376 | 890.421043 | 0.209172980369 | 0.163267011 | 0.625851295 | 0.536447329 | -0.524457002 | 0.854001991 |
| 117.523178001 | 0.106849604984 | 1099.893425 | 0.199642563356 | 0.151675735 | 0.613416333 | 0.517228608 | -0.502987659 | 0.863032553 |
| 150.000000000 | 0.096155850904 | 1559.967476 | 0.184990297878 | 0.182465275 | 0.580247378 | 0.473529985 | -0.455150846 | 0.857867995 |
| 200.000000000 | 0.078451554871 | 2549.343991 | 0.163745107442 | 0.254140709 | 0.512740397 | 0.378519120 | -0.353251679 | 0.822319237 |
| 300.000000000 | 0.041480233495 | 7232.360445 | 0.132021349057 | 0.354441445 | 0.372040267 | 0.134425607 | +0.017715255 | 0.482073181 |
| 500.000000000 | 0.000100000000 | 4999999.992452 | 0.189691229934 | 0.539037112 | 0.298029279 | 0.834985270 | +0.699750831 | 1.788294249 |
| 750.000000000 | 0.000100000000 | 7499999.999190 | 0.511165120792 | 0.961547046 | 0.312661989 | 0.494338611 | +0.401505267 | 1.054853754 |
| 1000.000000000 | 0.000100000000 | 9999999.999904 | 0.775971717067 | 1.203785424 | 0.320693135 | 0.334736875 | +0.252290101 | 0.763737956 |

Across the sampled high-Cth region (300-1000 J/K), the joint-objective span relative to its minimum is `487.762300%`. The stated practical-flatness diagnostic threshold is 1.0%.
At fixed Cth values of 500 J/K and above, the profile drives Gth to its unchanged 1e-4 W/K lower numerical guard while the joint objective rises. Those points do not constitute alternative joint optima.

UDDS values in this table were calculated only after Gth had been optimized against 1C+HPPC. They did not influence any optimization or parameter selection.

## Conditioning

| Point | Rank | Singular values | Condition number | Sensitivity correlation | Parameter correlation |
|---|---:|---:|---:|---:|---:|
| Previous 200 J/K boundary | 2/2 | 1.91498944976, 0.39099403228 | 4.89774597988 | +0.949809438034 | -0.919550819653 |
| Expanded-bound optimum | 2/2 | 1.84320063035, 0.0962833231804 | 19.143508652 | +0.963909123009 | -0.921857251733 |

## Withheld UDDS at expanded-bound optimum

RMSE `0.502850094 K`, bias `+0.420598577 K`, maximum absolute error `1.191183466 K`. The RMSE change relative to the frozen 1C-only value `0.526697652 K` is `+4.527751%`.
The improvement therefore persists relative to the frozen 1C-only baseline, but it shrinks from the previous 200 J/K profile-point RMSE `0.378519120 K` to `0.502850094 K`; the bias also reverses from `-0.353251679 K` to `+0.420598577 K`.

## Decisions

1. **Was Cth=200 J/K artificially truncated?** **Yes**. At fixed 200 J/K the profile objective is 0.163745107442 K^2; the expanded optimum is 0.126164862688 K^2.
2. **Is there a finite interior joint optimum?** **Yes**.
3. **Is the objective effectively flat toward large Cth?** **No**, using the reported 300-1000 J/K profile span and the explicit 1% diagnostic threshold.
4. **Does withheld UDDS improvement persist after relaxing the bound?** **Yes**.
5. **Is the joint parameterization numerically well-defined?** **Yes**.

**Final classification: JOINT OPTIMUM RESOLVED.**

The frozen thermal parameter module was not modified. No UDDS fitting, architecture change, empirical correction, electrical refit, OCP change, or entropy-coefficient change was performed.
