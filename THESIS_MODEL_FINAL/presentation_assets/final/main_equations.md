# Main equations for the 20-minute presentation

Internal convention: **`I > 0` discharge; `I < 0` charge**.

```math
\dot x_{n,avg}=-I/Q_n,\qquad \dot x_{p,avg}=+I/Q_p
```

```math
\dot x_{k,aux}=\frac{30}{\alpha_k}(x_{k,avg}-x_{k,aux})+\text{signed current forcing}
```

```math
x_{n,surf}=x_{n,aux}-\frac{\alpha_n I}{105Q_n},\qquad
x_{p,surf}=x_{p,aux}+\frac{\alpha_p I}{105Q_p}
```

```math
\eta_k=\frac{2RT}{F}\operatorname{asinh}\left(\frac{I}{6Q_kd_k\sqrt{x_{k,surf}(1-x_{k,surf})}}\right)
```

```math
U_{eq}=U_p-U_n,\qquad V=U_{eq}+\eta_p-\eta_n-R_0I
```

```math
\dot Q_{irr}=I(U_{eq}-V),\qquad \dot Q_{rev}=-IT\frac{dU_{eq}}{dT}
```

```math
C_{th}\dot T=\dot Q_{irr}+\dot Q_{rev}-G_{th}(T-T_{amb})
```

The electrical state and thermal state are propagated with exact left-sample zero-order-hold updates on the actual nonuniform timestamps. No state/OCP clipping is used.
