# Final implemented model equations

This document reports the equations as implemented in the accepted code. It distinguishes direct code transcription, conceptual interpretation, literature constitutive input, and project-specific calibration. The repository-wide internal sign convention is **positive current for discharge** (`I > 0` discharge, `I < 0` charge). Raw experimental signs are converted only at loader boundaries.

## Symbols and units

| Symbol | Meaning | Unit |
|---|---|---:|
| `I` | internal discharge-positive cell current | A |
| `t`, `Delta t` | time, interval length | s |
| `x_k,avg`, `x_k,aux`, `x_k,surf` | normalized solid concentration states for electrode `k in {n,p}` | 1 |
| `alpha_k` | grouped solid-diffusion time scale (`R_k^2/D_k` interpretation) | s |
| `Q_k` | grouped electrode capacity scale | C |
| `d_k` | grouped kinetic-rate parameter | 1/s |
| `U_n`, `U_p`, `U_eq`, `V` | electrode OCPs, cell equilibrium voltage, terminal voltage | V |
| `eta_n`, `eta_p` | grouped electrode overpotentials | V |
| `R0` | lumped series resistance | ohm |
| `R`, `F` | gas constant, Faraday constant | J/(mol K), C/mol |
| `SOC_cell` | independent cell-level coulomb-counted SOC | 1 |
| `Q_cell,Ah` | cell capacity used for SOC | Ah |
| `Qdot` | heat rate, positive into the cell | W |
| `T`, `T_amb` | predicted cell and prescribed ambient temperature | K |
| `Cth`, `Gth` | effective thermal capacitance and conductance | J/K, W/K |

## Grouped solid-state dynamics

State order in `models/grouped_spm.py::continuous_state_matrices` is `[x_n,avg, x_n,aux, x_p,avg, x_p,aux]`. With discharge-positive current:

```math
\dot{x}_{n,avg}=-I/Q_n
```

```math
\dot{x}_{n,aux}=\frac{30}{\alpha_n}(x_{n,avg}-x_{n,aux})-\frac{19I}{7Q_n}
```

```math
\dot{x}_{p,avg}=+I/Q_p
```

```math
\dot{x}_{p,aux}=\frac{30}{\alpha_p}(x_{p,avg}-x_{p,aux})+\frac{19I}{7Q_p}
```

Surface outputs are algebraic:

```math
x_{n,surf}=x_{n,aux}-\frac{\alpha_n I}{105Q_n},\qquad
x_{p,surf}=x_{p,aux}+\frac{\alpha_p I}{105Q_p}.
```

These formulas are direct code transcription. Conceptually, the code describes them as a two-state-per-electrode parabolic diffusion reduction, not a radial finite-difference discretization. **REDUCTION METHOD PROVENANCE REQUIRES MANUAL LITERATURE VERIFICATION**. A final thesis citation is required specifically for the reduction that produces the numerical factors `30`, `19/7`, and `1/105`, and for the interpretation `alpha_k=R_k^2/D_k`; the local repository establishes software provenance but not a complete literature derivation.

The linear subsystem is propagated exactly for left-sample zero-order-held current on every actual interval. `models/grouped_spm.py::exact_zoh_matrices` evaluates

```math
\begin{bmatrix}x_{k+1}\\1\end{bmatrix}
=\exp\!\left(\begin{bmatrix}A&B\\0&0\end{bmatrix}\Delta t_k\right)
\begin{bmatrix}x_k\\I_k\end{bmatrix}.
```

## Grouped kinetics and terminal voltage

`models/grouped_spm.py::grouped_overpotentials` implements

```math
\eta_k=\frac{2RT}{F}\operatorname{asinh}\!\left[
\frac{I}{6Q_kd_k\sqrt{x_{k,surf}(1-x_{k,surf})}}
\right],\quad k\in\{n,p\}.
```

`models/grouped_spm.py::simulate_grouped_spm` then evaluates

```math
U_{eq}=U_p(x_{p,surf})-U_n(x_{n,surf}),
```

```math
V=U_{eq}+\eta_p-\eta_n-R_0I.
```

These are direct implemented equations. `R0`, `alpha_p`, and `d_n` are the accepted conditional HPPC calibration quantities; `alpha_n` and `d_p` are fixed. Their accepted values are project-specific effective parameters, not universally identified material properties.

## Frozen Panasonic OCP representation

The Panasonic-facing functions are in `models/ocp_panasonic_qi_frozen.py` and enforce their supported domains without clipping or extrapolation.

```math
U_{eq}(z)=U_p[x_p(z)]-U_n[x_n(z)]
```

with nominal effective trajectories

```math
x_n(z)=x_{n,start}+z(x_{n,end}-x_{n,start}),\qquad
x_p(z)=x_{p,start}+z(x_{p,end}-x_{p,start}).
```

The negative OCP is the fixed Guo–White (2005) MCMB discharge spline implemented in `models/ocp_panasonic_qi.py::U_n`, supported on `0.01 <= x_n <= 0.80`. The positive OCP is a project-specific monotonic piecewise-linear 13-knot representation implemented by `models/ocp_panasonic_qi.py::U_p`, supported on `0.18 <= x_p <= 0.96`. Its knot voltages and nominal windows were identified from the Panasonic C/20 curve and frozen after selecting multistart 3. The endpoints are **nominal effective stoichiometric windows**, not uniquely identified physical stoichiometries.

## Independent cell SOC

`models/electrothermal_spm.py::_cell_soc_trajectory` implements

```math
SOC_{cell,k+1}=SOC_{cell,k}-\frac{I_k\Delta t_k}{3600Q_{cell,Ah}}.
```

`Q_cell,Ah` is independent of `Q_n` and `Q_p`. No SOC clipping is used. Panasonic uses `2.994979138519 Ah`; the own TRC cell uses `2.8999 Ah` and the values are never mixed.

## Heat generation

For discharge-positive current, `models/thermal.py::calculate_heat_generation` implements

```math
\dot Q_{irr}=I(U_{eq}-V),
```

```math
\dot Q_{ohmic}=I(R_0I)=R_0I^2,\qquad
\dot Q_{kinetic}=I(\eta_n-\eta_p),
```

```math
\dot Q_{irr}=\dot Q_{ohmic}+\dot Q_{kinetic},
```

```math
\dot Q_{rev}=-IT\frac{dU_{eq}}{dT},\qquad
\dot Q_{gen}=\dot Q_{irr}+\dot Q_{rev}.
```

The generic coupled wrapper uses surface-state `U_eq` and predicted terminal voltage. The accepted Panasonic thermal calibration/validation scripts deliberately use average-state equilibrium voltage,

```math
U_{eq,avg}=U_p(x_{p,avg})-U_n(x_{n,avg}),
```

and retain two explicitly labeled heat paths:

```math
\dot Q_{irr,measured}=I(U_{eq,avg}-V_{measured}),\qquad
\dot Q_{irr,model}=I(U_{eq,avg}-V_{model}).
```

The full-cell entropy coefficient is a literature constitutive input implemented in `models/panasonic_entropic.py::dUeq_dT_panasonic`, linearly interpolated only on `0.10 <= SOC <= 0.90`. Its code provenance is He et al. (2024) Supplementary Information Table S2, citing Lu et al. (2019). No entropy extrapolation or clipping is used.

## One-node thermal model

`models/thermal.py` implements

```math
C_{th}\frac{dT}{dt}=\dot Q_{gen}-G_{th}(T-T_{amb}).
```

For constant left-sample heat and ambient over one interval and `Gth > 0`, `models/thermal.py::lumped_thermal_step` uses the exact analytical update

```math
T_{k+1}=T_{amb,k}+\frac{\dot Q_{gen,k}}{G_{th}}+
\left(T_k-T_{amb,k}-\frac{\dot Q_{gen,k}}{G_{th}}\right)
\exp\!\left(-\frac{G_{th}\Delta t_k}{C_{th}}\right).
```

For `Gth=0`,

```math
T_{k+1}=T_k+\dot Q_{gen,k}\Delta t_k/C_{th}.
```

The accepted Panasonic `Cth` and `Gth` are project-specific effective one-node parameters calibrated on the 25 °C 1C discharge over SOC 0.10–0.90 using measured-voltage irreversible heat plus literature reversible heat. They are not independently measured material properties.

## Evaluation metrics

For residual `e_i = y_pred,i - y_meas,i`, the reporting scripts use

```math
RMSE=\sqrt{\frac{1}{N}\sum_i e_i^2},\qquad
MAE=\frac{1}{N}\sum_i |e_i|,
```

```math
bias=\frac{1}{N}\sum_i e_i,
```

```math
R^2=1-\frac{\sum_i e_i^2}{\sum_i(y_{meas,i}-\bar y_{meas})^2}.
```

RMSE is the primary error-scale summary, MAE describes typical absolute deviation, signed bias shows systematic over/underprediction, and `R²` is reported only with the variance context of the measured subset. A negative `R²` is possible for low-variance subsets such as CV plateaus and does not replace RMSE/bias interpretation.

The data-quality diagnostic is

```math
R_{app}=\left|\frac{\Delta V}{\Delta I}\right|.
```

`R_app` is an adjacent rest-to-load experimental-chain diagnostic. It must not be interpreted as true cell internal resistance or as the model parameter `R0`.

## Source-status summary

| Component | Status | Code source | Thesis citation action |
|---|---|---|---|
| Four grouped solid states and surface outputs | implemented directly | `models/grouped_spm.py` | manually verify/cite reduction derivation |
| Butler–Volmer-style grouped kinetics and voltage | implemented directly; inherited software architecture | `models/grouped_spm.py` | cite theoretical kinetics and software provenance |
| Panasonic MCMB negative OCP | literature constitutive input | `models/ocp_panasonic_qi.py` | Guo–White (2005) |
| Panasonic positive OCP | project-specific identified constitutive curve | `models/ocp_panasonic_qi_frozen.py` | cite this thesis identification workflow |
| Entropy coefficient | literature constitutive input | `models/panasonic_entropic.py` | He et al. (2024) Table S2; Lu et al. (2019) as cited there |
| One-node energy balance and exact update | implemented model equation | `models/thermal.py` | cite lumped thermal modeling convention if required by thesis style |
| Electrical/thermal parameter values | project-specific calibration | frozen Panasonic modules and reports | cite this thesis calibration chapters |
