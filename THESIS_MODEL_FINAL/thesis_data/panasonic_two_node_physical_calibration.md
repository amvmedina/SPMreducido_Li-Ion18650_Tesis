# Panasonic clean physical two-node calibration and validation

## TWO-NODE PARAMETERIZATION NOT RECOVERABLE FROM THIS EXPERIMENT

Only the four legacy physical parameters were fitted, and only on the closed Panasonic 25 C 1C interval `0.10 <= SOC <= 0.90`. The independent UDDS interval was evaluated once after freezing the optimum. No electrical/OCP/entropy parameter, initial temperature, empirical heat term, or correction factor was fitted or modified.

## Clean architecture

```text
C_c dT_c/dt = Q_gen - G_cs (T_c - T_s)
C_s dT_s/dt = G_cs (T_c - T_s) - hA (T_s - T_amb)
measured output = T_s
```

The clean module preserves the legacy conductance parameterization (`G_cs=1/R_cs`, `hA=1/R_sa`) and uses an exact zero-order-hold state transition at every actual nonuniform interval. `T_c(t0)=T_s(t0)=T_measured(t0)` in calibration and validation, matching the legacy initialization. There is no interface for SOC multipliers, regime gates, heat clipping, additive heat, or temperature bias.

Heat used for both 1C fitting and UDDS validation is exactly `I*(Ueq_average-V_measured) - I*T_measured*dUeq/dT`, with discharge-positive internal current, frozen average electrode states/OCPs, and the verified Panasonic entropy table.

## Datasets

- 1C calibration source: `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\C20 OCV and 1C discharge tests_start_of_tests\03-09-17_17.59 3349_Dis1C_1.mat`. The exact entropy-supported discharge segment runs from SOC 0.90 to 0.10 (`0.10 <= SOC <= 0.90`), from source times 371.860459714 to 3346.788871686 s (299 samples; duration 2974.928411972 s).
- Independent UDDS validation source: `C:\Users\amvel\Documents\Andes\2026-1\wykht8y7tg-1\Panasonic 18650PF Data\Panasonic 18650PF Data\25degC\Drive cycles\03-21-17_00.29 25degC_UDDS_Pan18650PF.mat`. The exact continuous SOC 0.90-to-0.10 interval contains 194,744 samples over 19501.646379329 s.

## 1C calibration result

Forty deterministic log-space least-squares starts were run with broad positivity guard bounds. 40 succeeded; 27 are within 0.100% relative RMSE of the best.

| Parameter | Fitted value | Guard bounds | Nearest log-bound fraction |
|---|---:|---:|---:|
| `c_c_j_per_k` | 94.7844616277 | [1e-3, 5000] J/K | 0.257089 |
| `c_s_j_per_k` | 0.0001 | [1e-4, 2000] J/K | 0.000000 |
| `g_cs_w_per_k` | 0.889090880712 | [1e-5, 100] W/K | 0.293008 |
| `hA_w_per_k` | 0.110337426796 | [1e-5, 10] W/K | 0.326213 |

Modal eigenvalues are -9994.29141965 and -0.00103557093098 1/s, giving fast and slow time constants 0.000100057 and 965.650898537 s.

| Metric | Clean two-node 1C | Frozen one-node 1C |
|---|---:|---:|
| RMSE [K] | 0.138221456 | 0.143622065 |
| MAE [K] | 0.115748804 | 0.119569322 |
| Bias [K] | -0.017734800 | -0.008530653 |
| Maximum absolute error [K] | 0.338561724 | not supplied in closure prompt |
| R2 | 0.967306708 | 0.964702003 |

The corresponding resistances are `R_cs=1.124744412 K/W` and `R_sa=9.063107860 K/W`.

The clean two-node RMSE change relative to one-node is +3.760%. This is a small in-sample change, and the independent-transfer result below - not visual appearance - is used to decide whether the added state is worthwhile.

## Multistart stability

Near-optimal means RMSE no more than 0.1% above the best; the complete successful-start table is in the CSV output.

| Parameter | All-successful min | All-successful max | Near-optimal min | Near-optimal max | Near-optimal relative range |
|---|---:|---:|---:|---:|---:|
| `c_c_j_per_k` | 0.001 | 409.999758 | 94.7844616 | 94.784497 | 0.000037% |
| `c_s_j_per_k` | 0.0001 | 117.522167 | 0.0001 | 0.000100000006 | 0.000006% |
| `g_cs_w_per_k` | 0.889090881 | 100 | 0.889090881 | 0.889092127 | 0.000140% |
| `hA_w_per_k` | 2.0101711e-05 | 0.155750892 | 0.110337422 | 0.110337427 | 0.000004% |

## Local practical conditioning

The numerical sensitivity matrix is the `299 x 4` Jacobian of surface-temperature residual [K] with respect to dimensionless log-parameters. All 1,196 matrix entries are stored under `scope=sensitivity_matrix` in `panasonic_two_node_physical_metrics.csv`. Its numerical rank is **4/4**, singular values are `4.586881518e+01, 7.101541791e+00, 5.972396722e-01, 1.794971975e-06`, and its condition number is **2.555405645e+07**.

Sensitivity-column correlation matrix:

| | c_c_j_per_k | c_s_j_per_k | g_cs_w_per_k | hA_w_per_k |
|---|---:|---:|---:|---:|
| c_c_j_per_k | +1.000000 | +0.263669 | -0.731468 | +0.247306 |
| c_s_j_per_k | +0.263669 | +1.000000 | +0.227200 | +0.114882 |
| g_cs_w_per_k | -0.731468 | +0.227200 | +1.000000 | +0.231514 |
| hA_w_per_k | +0.247306 | +0.114882 | +0.231514 | +1.000000 |

Local parameter correlation matrix inferred from the pseudoinverse of `J^T J`:

| | c_c_j_per_k | c_s_j_per_k | g_cs_w_per_k | hA_w_per_k |
|---|---:|---:|---:|---:|
| c_c_j_per_k | +1.000000 | -0.613033 | +0.942437 | -0.528477 |
| c_s_j_per_k | -0.613033 | +1.000000 | -0.561445 | +0.128233 |
| g_cs_w_per_k | +0.942437 | -0.561445 | +1.000000 | -0.282755 |
| hA_w_per_k | -0.528477 | +0.128233 | -0.282755 | +1.000000 |

This is a local practical-conditioning audit based on one measured surface-temperature trace. It is not a structural-identifiability claim.

The fitted fast mode is 0.000100057118 s, whereas the smallest retained 1C interval is 6.79287615 s. The optimizer therefore pushes the surface node below the experiment's temporal resolution instead of identifying a resolved second thermal mode.

A distinct secondary multistart basin has RMSE 0.143622069 K (only 0.005400613 K above the best) at `C_c=0.001`, `C_s=117.522`, `G_cs=100`, `hA=0.109112`. Its surface trajectory differs from the best by only 0.039956089 K RMSE (maximum 0.093592697 K), despite radically different parameters. The four physical quantities are therefore not practically recoverable from this surface-only experiment.

## Critical no-refit UDDS validation

Both models use the same measured-voltage/reversible measured-temperature heat and exact SOC 0.90-to-0.10 interval. All fitted two-node parameters were frozen before this evaluation.

| Metric | Clean two-node UDDS | Existing one-node UDDS A |
|---|---:|---:|
| RMSE [K] | 0.532525609 | 0.526697652 |
| MAE [K] | 0.521553418 | 0.512960900 |
| Bias [K] | -0.521470836 | -0.512714512 |
| Maximum absolute error [K] | 0.858785101 | 0.875239690 |
| R2 | -4.233330651 | -4.119410369 |
| Final measured surface [degC] | 26.469807973 | 26.469807973 |
| Final predicted [degC] | 26.110904185 | 26.112601243 |
| Measured rise [K] | +0.636968211 | +0.636968211 |
| Predicted rise [K] | +0.278064423 | +0.279761481 |

The clean two-node UDDS RMSE change relative to one-node is -1.107%. Maximum predicted `T_c-T_s` is 0.162690851 K (UDDS) and 0.708001261 K (1C).

## Final decisions

1. **Did the old apparent thermal improvement depend on the empirical SOC/regime correction?** The archived SOC multiplier improved its legacy 1C RMSE (0.727691 to 0.655786 K), but worsened legacy UDDS RMSE (0.204643 to 0.222634 K). Therefore the old low UDDS error did not depend on the correction; the apparent late-1C improvement did.

2. **Does the clean physical two-node model outperform the current one-node model on independent UDDS?** **No**; its RMSE change is -1.107% with no refit or heat correction.

3. **Is the improvement large enough to justify the two-node model in the thesis?** **No**. The decision is based on independent UDDS transfer and recoverability: the clean model must actually improve held-out behavior and supply a defensible finite parameterization.

4. **Are the parameters practically recoverable from the available temperature data?** **TWO-NODE PARAMETERIZATION NOT RECOVERABLE FROM THIS EXPERIMENT.** The conclusion follows from rank, singular spectrum/condition number, correlations, multistart spread, and bound proximity; no structural-identifiability claim is made.

No model-voltage heat validation, CC-CV thermal validation, own-TRC analysis, or empirical repair was performed.
