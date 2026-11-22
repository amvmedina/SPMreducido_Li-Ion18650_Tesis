# Panasonic grouped-SPM parameter mapping

## Scope

Code/equation audit and analytical mapping only. No parameter was optimized, and the frozen model and Panasonic OCP modules were not modified.

## Actual implemented parameterization

For discharge-positive model current `I`, the implemented state/output equations are:

```text
dx_n_avg/dt = -I/Q_n
dx_n_aux/dt = (30/alpha_n)(x_n_avg-x_n_aux) - 19 I/(7 Q_n)
x_n_surf = x_n_aux - alpha_n I/(105 Q_n)
dx_p_avg/dt = +I/Q_p
dx_p_aux/dt = (30/alpha_p)(x_p_avg-x_p_aux) + 19 I/(7 Q_p)
x_p_surf = x_p_aux + alpha_p I/(105 Q_p)
eta_k = (2RT/F) asinh[I/(6 Q_k d_k sqrt(x_k_surf(1-x_k_surf)))]
V = U_p(x_p_surf)-U_n(x_n_surf)+eta_p-eta_n-R0 I
```

| Parameter | Units | Exact role and equation locations | Category | Status for HPPC |
|---|---:|---|---|---|
| `alpha_n` | s | `30/alpha_n` negative auxiliary-state relaxation and `-alpha_n I/(105 Q_n)` surface correction. The documented group is particle radius squared over solid diffusivity. | solid diffusion/state output | dynamic candidate |
| `alpha_p` | s | `30/alpha_p` positive auxiliary-state relaxation and `+alpha_p I/(105 Q_p)` surface correction. | solid diffusion/state output | dynamic candidate |
| `Q_n` | C | Negative average and auxiliary state forcing, surface correction, and kinetic denominator. | electrode state evolution/kinetic scale | analytically frozen |
| `Q_p` | C | Positive average and auxiliary state forcing, surface correction, and kinetic denominator. | electrode state evolution/kinetic scale | analytically frozen |
| `d_n` | 1/s | Appears only in `6 Q_n d_n sqrt(x(1-x))` inside negative `asinh`. | grouped Butler-Volmer kinetics | dynamic candidate |
| `d_p` | 1/s | Appears only in `6 Q_p d_p sqrt(x(1-x))` inside positive `asinh`. | grouped Butler-Volmer kinetics | dynamic candidate |
| `SOC0_n` | dimensionless | Initializes both negative average and auxiliary normalized concentrations. It is not cell SOC. | thermodynamic alignment/initial state | analytically frozen |
| `SOC0_p` | dimensionless | Initializes both positive average and auxiliary normalized concentrations. It is not cell SOC. | thermodynamic alignment/initial state | analytically frozen |
| `R0` | ohm | Algebraic terminal drop `-R0 I`; also supplies `R0 I^2` to the coupled heat calculation. | ohmic resistance | dynamic candidate |

`models/electrothermal_spm.py` does not redefine these electrical parameters. It re-evaluates the grouped kinetic equations at predicted temperature, applies the same voltage equation, and keeps cell-level `Q_cell_Ah`/SOC separate.

## Exact frozen-window mapping

With total discharged charge `Delta q = 3600 Q_cell_Ah`, integration of the two average-state equations gives:

```text
Q_n = Delta q/(x_n_start-x_n_end)
Q_p = Delta q/(x_p_end-x_p_start)
SOC0_n = x_n_start
SOC0_p = x_p_start
```

- `Delta q = 10781.924898668400 C`
- `Q_n = 23008.732269258333 C` (6.391314519238 Ah-equivalent)
- `Q_p = 14780.898652565113 C` (4.105805181268 Ah-equivalent)
- `SOC0_n = 0.720451136448`
- `SOC0_p = 0.224027933455`

These are grouped effective electrode capacity/state quantities derived from the nominal effective windows, not separately measured physical electrode capacities or uniquely identified stoichiometries.

## Forward state check

Passing exactly `2.994979138519 Ah` through the average-state equations gives `x_n=0.251849704066` and `x_p=0.953477823621`.
The numerical endpoint errors are `0.000e+00` and `0.000e+00`. Thus the frozen mapping is reproduced to floating-point precision without optimization or clipping.

An independent forward check using the actual `exact_zoh_matrices` implementation applied 1 A for `10781.924898668400 s` (exactly the frozen cell capacity), then zero current for `300000 s`. The loaded average-state endpoint errors were `0.000e+00` and `0.000e+00`; after relaxation the maximum error across all four states was `0.000e+00`.

The exact loaded endpoint statement applies to the volume-average states. During discharge, auxiliary/surface states differ by the model's diffusion terms; after zero-current relaxation they converge to the unchanged average endpoints.
