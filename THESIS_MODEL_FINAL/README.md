# Reduced electrochemical-thermal cell model

## Aplicación local de auditoría

La interfaz para evaluar suficiencia de nuevos registros Chroma se ejecuta con:

```powershell
python -m streamlit run app.py
```

Consulte `docs/aplicacion_auditoria_baterias_streamlit.md` para instalación, alcance e interpretación.

This repository contains the definitive implementation for Thesis Objective 1:
a reduced grouped single-particle model (SPM), independent cell-level SOC by
coulomb counting, electrochemical heat generation, and a predicted lumped cell
temperature that is fed back into the electrochemical kinetics.

It contains no calibration, optimization, experimental data processing, or new
synthetic-identifiability work. Numerical values in tests and examples are
execution fixtures only and are not identified battery parameters.

## Public API and current convention

The coupled entry point is:

```python
from models import simulate_electrothermal_spm
```

Its default and recommended current convention is

```text
I > 0: discharge
I < 0: charge
```

The optional input convention `guo_charge_positive` is converted at the model
boundary. Returned current is always the internal discharge-positive current,
and all equations below use that convention.

The original prescribed-temperature electrical entry point remains available
and numerically unchanged:

```python
from models import simulate_grouped_spm
```

Time is in seconds, current in amperes, voltage in volts, heat rate in watts,
temperature in kelvin, thermal capacitance in J/K, thermal conductance in W/K,
cell capacity in Ah, and entropic coefficient in V/K.

## Electrochemical model

The state order is

```text
x = [x_n_avg, x_n_aux, x_p_avg, x_p_aux]
```

Here `avg` is normalized volume-average solid concentration and `aux` is the
transformed surface state from the grouped parabolic solid-concentration
approximation. For discharge-positive current, the exact implemented equations
are

```math
\dot{x}_{n,avg}=-I/Q_n,
```

```math
\dot{x}_{n,aux}=(30/\alpha_n)(x_{n,avg}-x_{n,aux})
                 -(19/(7Q_n))I,
```

```math
\dot{x}_{p,avg}=+I/Q_p,
```

```math
\dot{x}_{p,aux}=(30/\alpha_p)(x_{p,avg}-x_{p,aux})
                 +(19/(7Q_p))I.
```

The normalized surface concentrations are algebraic outputs:

```math
x_{n,surf}=x_{n,aux}-(\alpha_n/(105Q_n))I,
```

```math
x_{p,surf}=x_{p,aux}+(\alpha_p/(105Q_p))I.
```

This is a two-state-per-electrode parabolic diffusion approximation. It is not
a radial finite-difference model. The group `alpha_k = R_k^2/D_k` is a diffusion
time scale; particle radius and solid diffusivity are not separately resolved.

Current is zero-order held at the left sample of each time interval. The code
computes the matrix exponential of

```text
M = [[A, B],
     [0, 0]]
```

so propagation of the linear four-state subsystem is exact for each supplied
interval, including nonuniform time grids. A zero-length interval leaves the
state unchanged.

The grouped kinetic overpotentials are

```math
\eta_k=\frac{2RT}{F}\,\operatorname{asinh}\left(
\frac{I}{6Q_kd_k\sqrt{x_{k,surf}(1-x_{k,surf})}}
\right),\qquad k\in\{n,p\}.
```

Temperature enters the existing electrochemical model only through the factor
`2RT/F` in these kinetics. It does not enter the solid-state dynamics or the
current OCP functions. Equilibrium voltage and terminal voltage are

```math
U_{eq}=U_p(x_{p,surf})-U_n(x_{n,surf}),
```

```math
V=U_{eq}+\eta_p-\eta_n-R_0I.
```

The nine electrical parameters, in their mandatory order, are

```text
[alpha_n_s, alpha_p_s, q_n_c, q_p_c, d_n, d_p,
 soc0_n, soc0_p, r0_ohm]
```

`q_n_c` and `q_p_c` are grouped electrode capacity scales in coulombs. The
historical names `soc0_n` and `soc0_p` denote initial normalized electrode
concentrations; neither is silently reinterpreted as full-cell SOC.

## Independent cell-level SOC

Cell SOC is configured explicitly by `CellSOCParameters(q_cell_ah,
soc_cell_0)` and propagated by discharge-positive coulomb counting:

```math
SOC_{cell}(t)=SOC_{cell}(0)-
\frac{1}{3600Q_{cell,Ah}}\int_0^t I(\tau)\,d\tau.
```

With the discrete left-sample ZOH convention,

```math
SOC_{k+1}=SOC_k-\frac{I_k\Delta t_k}{3600Q_{cell,Ah}}.
```

`Q_cell_Ah` is independent of `Q_n` and `Q_p`; the implementation never infers
one from the other. SOC is not clipped. `cell_soc_valid`,
`cell_soc_domain_valid`, and `first_soc_violation` report values outside the
inclusive physical interval `[0,1]` while preserving the computed value.

## Heat generation

For the implemented voltage and discharge-positive current convention,
irreversible electrical loss is

```math
\dot Q_{irr}=I(U_{eq}-V).
```

Substitution of the terminal-voltage equation gives the exposed decomposition

```math
\dot Q_{ohmic}=I(R_0I)=R_0I^2,
```

```math
\dot Q_{kinetic}=I(\eta_n-\eta_p),
```

```math
\dot Q_{irr}=\dot Q_{ohmic}+\dot Q_{kinetic}.
```

No absolute value or clipping is applied. The diagnostic
`heat_irreversible_nonnegative` indicates whether the evaluated irreversible
heat is nonnegative to numerical tolerance. This is normally true for ordinary
loss-producing operation with a physically consistent parameter set.

Reversible heat is disabled by default because this repository contains no
chemistry-specific entropic data. If the caller explicitly supplies either a
constant or a function for `dU_eq/dT`, the discharge-positive convention uses

```math
\dot Q_{rev}=-IT\frac{dU_{eq}}{dT},
```

and

```math
\dot Q_{gen}=\dot Q_{irr}+\dot Q_{rev}.
```

The callable form receives `(x_n_surface, x_p_surface, temperature_k)` and
returns `dU_eq/dT` in V/K. With no supplied coefficient,
`heat_reversible = 0`, `heat_generation_total = heat_irreversible`, and
`reversible_heat_enabled` is false. The model does not claim chemistry-specific
reversible heat in that mode.

## Lumped thermal state and closed coupling

`LumpedThermalParameters` requires four explicit cell-level quantities:

```text
c_th_j_per_k   effective thermal capacitance C_th [J/K]
g_th_w_per_k   conductance to ambient G_th [W/K]
t_initial_k    initial predicted cell temperature [K]
t_ambient_k    externally specified ambient temperature [K]
```

The predicted cell-temperature state obeys

```math
C_{th}\frac{dT}{dt}=\dot Q_{gen}-G_{th}(T-T_{amb}).
```

For constant heat and ambient over an interval, the analytical ZOH update is

```math
T_{k+1}=T_{amb}+\frac{\dot Q_{gen,k}}{G_{th}}+
\left(T_k-T_{amb}-\frac{\dot Q_{gen,k}}{G_{th}}\right)
\exp\left(-\frac{G_{th}\Delta t_k}{C_{th}}\right)
```

when `G_th > 0`. For the adiabatic case `G_th = 0`, the exact update is

```math
T_{k+1}=T_k+\dot Q_{gen,k}\Delta t_k/C_{th}.
```

The wrapper evaluates electrical voltage and heat at predicted `T_k`, advances
the thermal state to `T_(k+1)`, and uses that predicted temperature in the next
kinetic evaluation. Measured cell temperature is not used as the predicted
output. No Arrhenius law or external temperature dependency has been added.

## Coupled outputs

`ElectrothermalSPMResult` explicitly returns arrays for:

- `time_s`, discharge-positive `current_a`, `terminal_voltage_v`;
- `cell_soc`, its validity array/flag, and its first violation;
- predicted `temperature_k`;
- `heat_generation_total_w`, `heat_irreversible_w`, `heat_ohmic_w`,
  `heat_kinetic_w`, and `heat_reversible_w`;
- equilibrium voltage `u_eq_v`, both electrode OCPs, overpotentials, and ohmic
  drop;
- average, auxiliary, and surface normalized concentrations for both
  electrodes;
- electrochemical domain validity and the original first-violation record.

All average, auxiliary, and surface normalized concentrations must remain
strictly inside `(0,1)`. The grouped SPM does not clip them; electrical outputs
from the first violation onward are `NaN`.

## Example configuration

```python
import numpy as np

from models import (
    CellSOCParameters,
    GroupedSPMParameters,
    LumpedThermalParameters,
    simulate_electrothermal_spm,
)

electrical = GroupedSPMParameters(
    alpha_n_s=3000.0,
    alpha_p_s=1800.0,
    q_n_c=11000.0,
    q_p_c=11000.0,
    d_n=5.0e-5,
    d_p=7.0e-4,
    soc0_n=0.60,
    soc0_p=0.40,
    r0_ohm=0.02,
)
soc = CellSOCParameters(q_cell_ah=3.0, soc_cell_0=0.80)
thermal = LumpedThermalParameters(
    c_th_j_per_k=500.0,
    g_th_w_per_k=0.5,
    t_initial_k=298.15,
    t_ambient_k=298.15,
)
time_s = np.arange(0.0, 21.0)
current_a = np.full(time_s.size, 0.2)
result = simulate_electrothermal_spm(
    time_s, current_a, electrical, soc, thermal
)
```

Every number above is an uncalibrated API example, not a target-cell value.

## Predicted quantities, inputs, and remaining limitations

Predicted/model-computed quantities are electrode states, surface states, OCP,
overpotentials, equilibrium and terminal voltage, independent cell SOC, heat
contributions, and cell temperature. Externally specified quantities are time,
current, ambient temperature, the nine grouped electrical parameters,
`Q_cell_Ah`, `SOC_cell_0`, `C_th`, `G_th`, and `T_initial`.

The OCP equations in `models/ocp.py` are the chemistry-specific functions used
by the public Guo grouped NMC/graphite implementation (documented there by URL
and commit). Their suitability for the target cell is not established. They
are the only retained chemistry-specific constitutive data; no code was imported
from that repository during this extension.

Experimental work is still required to identify or validate:

- the nine grouped electrical parameters and initial electrode states;
- target-cell capacity `Q_cell_Ah` and the initial cell SOC;
- effective thermal capacitance `C_th` and conductance `G_th` under the intended
  boundary conditions;
- the OCP functions for the target chemistry;
- chemistry-specific `dU_eq/dT` if reversible heat is to be enabled.

These are parameterization and validation tasks, not missing state equations in
the present Objective 1 implementation.

## Tests

Run the complete suite with:

```text
python -m pytest -q
```

The tests cover prescribed-temperature electrochemical regression, exact cell
SOC accounting and out-of-range reporting, thermal equilibrium, passive
cooling, heating, finite coupled execution, disabled reversible heat, and
temperature feedback into terminal voltage. They do not calibrate or optimize
parameters and do not access experimental or synthetic datasets.
