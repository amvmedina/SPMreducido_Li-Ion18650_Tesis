# Panasonic HPPC local sensitivity and identifiability audit

## Scope and data product

No optimization or dynamic calibration was performed. The frozen OCP, nominal effective windows, grouped-SPM equations, and electrothermal architecture were not modified.

- Complete pulse/rest groups: **64** (IDs 60, 64, and 67 excluded).
- Independently initialized contiguous multi-rate histories: **14**, separated by the source file's long omitted-data gaps.
- Voltage-response observations in the Jacobian: 99461.
- Repeated source indices across independent groups: 0 (a boundary rest sample can close one group and initialize the next, but the initializing sample is excluded from the Jacobian).
- Each local history begins at the first usable pulse's `pre_pulse_index`, ends at the last usable pulse's `post_rest_index`, and uses local time zero at that block boundary. State is preserved across every pulse and rest inside the contiguous block.
- Raw negative current was converted at the loader boundary to model discharge-positive current.
- Cumulative exported `Ah` initializes each local average and auxiliary electrode state; current is never integrated through the 13 long omitted-data gaps.
- Actual nonuniform time intervals, including zero-length duplicate timestamps, were retained. No resampling was performed.
- Measured battery temperature is used only as the prescribed temperature in the existing kinetic factor; no thermal parameters are involved.

Phase definitions: `immediate` is the first loaded sample; `within_pulse` is every later loaded sample through pulse end; `post_relaxation` is every contiguous rest sample after pulse end. `all_combined` contains all three phases and excludes the pre-pulse initialization sample.

SOC bins are non-overlapping and use pulse-start cell SOC: low `[0,1/3)`, middle `[1/3,2/3)`, and high `[2/3,1]`.

## Local reference and finite differences

The model API has no built-in parameter defaults. The existing README/smoke-test fixture is therefore used only as a local numerical reference, not as a Panasonic parameter estimate:

| Parameter | Local reference |
|---|---:|
| `alpha_n` | 3000 |
| `alpha_p` | 1800 |
| `d_n` | 5e-05 |
| `d_p` | 0.0007 |
| `R0` | 0.02 |

Central scale-aware perturbations use `p*(1 +/- 0.0001)`. Each Jacobian entry is `dV/d(ln p) = p dV/dp` in volts, so columns correspond to equal fractional parameter changes.

## Sensitivity by response phase

Column RMS values below are mV per unit fractional parameter change.

### Immediate

Observations: 64

| Parameter | RMS normalized sensitivity [mV] |
|---|---:|
| `alpha_n` | 1.13589 |
| `alpha_p` | 12.8164 |
| `d_n` | 40.5784 |
| `d_p` | 15.2724 |
| `R0` | 185.93 |

### Within Pulse

Observations: 6400

| Parameter | RMS normalized sensitivity [mV] |
|---|---:|
| `alpha_n` | 1.19534 |
| `alpha_p` | 13.6248 |
| `d_n` | 40.6743 |
| `d_p` | 15.3528 |
| `R0` | 185.894 |

### Post Relaxation

Observations: 92997

| Parameter | RMS normalized sensitivity [mV] |
|---|---:|
| `alpha_n` | 0.0920395 |
| `alpha_p` | 1.90781 |
| `d_n` | 0 |
| `d_p` | 0 |
| `R0` | 0 |

### All Combined

Observations: 99461

| Parameter | RMS normalized sensitivity [mV] |
|---|---:|
| `alpha_n` | 0.31732 |
| `alpha_p` | 3.93114 |
| `d_n` | 10.3689 |
| `d_p` | 3.91371 |
| `R0` | 47.3905 |

## Combined identifiability

- Numerical rank: 5 of 5.
- Singular values: 15.3215, 1.40396, 0.638388, 0.110284, 0.0217429.
- Jacobian condition number: 704.665.
- Weakest right-singular direction (unit vector; sign is arbitrary): `alpha_n`=+0.99819, `alpha_p`=-0.04576, `d_n`=-0.00253, `d_p`=-0.03893, `R0`=-0.00025.

The reported parameter-correlation matrix is estimator-oriented: the pseudoinverse of `J.T @ J` is normalized by its diagonal. It is not the ordinary correlation of parameter values, because no parameters were estimated.

Pairs with `|correlation| >= 0.90` in the combined audit:

| Parameter pair | Correlation |
|---|---:|
| `d_p` / `R0` | -0.936179 |

## Stratified audit

| Stratum | Observations | Condition number | Weakest singular value |
|---|---:|---:|---:|
| current: `0.5C` | 25802 | 648.05 | 0.00207788 |
| current: `1C` | 25802 | 643.918 | 0.0039479 |
| current: `2C` | 23959 | 724.344 | 0.00625728 |
| current: `4C` | 22116 | 3546.11 | 0.00234266 |
| current: `6C` | 1782 | 6239.78 | 0.00188871 |
| soc: `high` | 37670 | 930.266 | 0.0108542 |
| soc: `middle` | 22602 | 1624.91 | 0.00480848 |
| soc: `low` | 39189 | 1202.59 | 0.00704513 |

High-correlation pairs in the stratified audits:

| Stratum | Parameter pair | Correlation |
|---|---|---:|
| current: `0.5C` | `d_n` / `R0` | +0.988675 |
| current: `1C` | `d_n` / `R0` | +0.996294 |
| current: `2C` | `d_n` / `R0` | +0.999105 |
| current: `4C` | `d_n` / `R0` | +0.999879 |
| current: `6C` | `d_n` / `R0` | +0.999863 |
| soc: `high` | `d_p` / `R0` | -0.992756 |
| soc: `middle` | `alpha_n` / `alpha_p` | -0.937463 |
| soc: `middle` | `d_p` / `R0` | -0.992751 |
| soc: `low` | `d_p` / `R0` | -0.963435 |

Detailed SVD, weakest-direction, and correlation results for every phase/current/SOC subset are in the machine-readable CSV files.

## Apparent resistance context

The 64 usable measured pulse steps have `R_app = 20.071-33.348 mOhm`. This observed step quantity contains the model's algebraic `R0 I` contribution together with instantaneous surface-state/OCP changes and the algebraic grouped kinetic contributions `eta_p-eta_n`. Consequently, `R_app` is contextual evidence for voltage-step scale only; setting `R0 = R_app` would assign non-ohmic voltage terms to `R0` and is not justified.

Consistently, the immediate-response normalized RMS sensitivities are 185.930 mV for `R0`, 40.578 mV for `d_n`, 15.272 mV for `d_p`, 12.806 mV for `alpha_p`, and 1.136 mV for `alpha_n`. The model therefore predicts that the measured instantaneous step is not uniquely ohmic at this local reference. The fixture `R0=20 mOhm` happens to sit at the lower edge of the measured `R_app` range, but that numerical proximity is not an identification result.

## Minimum defensible dynamic subset

| Parameter | Recommendation | Reason from implemented equations and this local audit |
|---|---|---|
| `R0` | **CALIBRATE** | It is the only purely ohmic candidate and acts directly on every loaded sample; it has no relaxation sensitivity. It must be separated jointly from algebraic kinetic voltage rather than assigned from `R_app`. |
| `alpha_n` | **WEAKLY IDENTIFIABLE / CONDITIONAL** | It controls negative diffusion, but has the smallest combined RMS sensitivity and contributes about 0.998 of the unit weakest singular direction. Do not release it in the smallest first fit. |
| `alpha_p` | **CALIBRATE** | It controls the positive surface displacement and positive diffusion-relaxation time scale; its post-rest RMS sensitivity is materially larger than that of `alpha_n`. |
| `d_n` | **CALIBRATE** | It has a strong loaded-sample signature and is not highly correlated with `R0` only when multiple current magnitudes are combined. The near-unity `d_n`/`R0` correlations within each fixed C-rate make the multi-rate design essential. |
| `d_p` | **WEAKLY IDENTIFIABLE / CONDITIONAL** | It enters only the algebraic positive kinetic overpotential, has no rest signature, and is strongly correlated with `R0` in the combined audit (`-0.936`). Do not release it in the smallest first fit. |

Thus the smallest defensible primary calibration subset is **`R0`, `alpha_p`, and `d_n`**, using all current magnitudes together. Hold `alpha_n` and `d_p` fixed at explicitly declared provisional values for that first fit, then release at most one conditionally in a later nested check only if stability improves. This is a subset recommendation, not a parameter estimate; the present generic fixture is not asserted to be the physically correct fixed value.

The earlier synthetic-identifiability phase is closed and its numerical artifacts are intentionally absent from this repository. Its general sensitivity conclusion cannot override this experimental audit: nonzero synthetic/local voltage sensitivity shows that a parameter affects the output, whereas practical recoverability additionally requires independent sensitivity directions under the actual HPPC currents, SOC coverage, sampling, noise, and parameter-reference point. The condition numbers, correlations, and weakest singular vectors here assess that distinction; they are still local diagnostics, not proof of global identifiability.

## Stop condition

No dynamic parameter optimization, HPPC calibration, UDDS/1C validation, thermal fitting, OCP change, window change, or model-equation change was performed.
