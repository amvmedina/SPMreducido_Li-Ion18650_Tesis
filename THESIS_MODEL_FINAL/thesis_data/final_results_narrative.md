# Final results narrative

## Model contribution

The implemented framework couples a four-state grouped single-particle electrical model, independent cell-level coulomb counting, electrochemical heat generation, and an analytically propagated one-node thermal state. Current is internally discharge-positive, actual nonuniform sample intervals are retained, and invalid SOC/OCP states are reported rather than clipped. The final Panasonic parameterization combines a fixed MCMB negative-electrode OCP, a frozen 13-knot positive-electrode OCP, nominal effective electrode windows, conditional HPPC dynamic parameters, a literature full-cell entropy coefficient, and effective one-node thermal parameters.

## Calibration findings

The Panasonic C/20 representation is accepted on its declared 5–95% fit interval with 1.086 mV RMSE, while the electrode windows remain weakly identifiable. Dynamic calibration was deliberately conditional: only HPPC pulses beginning at SOC at or above one third entered the final three-parameter fit, giving 36.835 mV phase-balanced RMSE. The fitted `R0`, `alpha_p`, and `d_n` are therefore not claimed as globally identified Panasonic constants.

The primary thermal pair was calibrated on the 25 °C 1C discharge over SOC 0.10–0.90 using measured-voltage irreversible heat and literature reversible heat. It achieved 0.143622 K RMSE and `R²=0.964702`. These are effective lumped parameters under the tested boundary condition, not independently measured material properties.

## Validation findings

Withheld low-SOC HPPC transfer is limited: phase-balanced RMSE rises to 85.729 mV and the dominant error is the loaded transient. The independent 1C discharge shows the same direction of deterioration, with 114.362 mV overall RMSE and 195.954 mV at low SOC. UDDS voltage transfer is stronger overall at 30.988 mV RMSE, but its residual still grows from 14.361 mV at high SOC to 52.785 mV at low SOC and large regenerative current steps are overpredicted. CC–CV charge transfer remains limited at 61.034 mV RMSE with −47.526 mV bias, showing incomplete transfer across current direction.

Thermally, the one-node model underpredicts the independent UDDS temperature trajectory (0.526698 K RMSE, −0.512715 K bias). On CC–CV charge, the measured-voltage heat path gives 0.427781 K RMSE, whereas the fully predictive path degrades to 0.864391 K. The 436.438 J electrical heat deficit explains most of this added coupled error; the retained reversible contribution is 213.914 J. These results separate limitations of electrical heat prediction from limitations of the thermal state model.

## Why the one-node thermal model was retained

The clean two-node model marginally improves the fitted 1C RMSE (0.138221 K versus 0.143622 K) but does not improve independent UDDS transfer (0.532526 K versus 0.526698 K). Its surface capacitance reaches the numerical lower guard, its fast mode is about `1e-4 s`—far below the experiment’s temporal resolution—and its local Jacobian condition number is about `2.56e7`. The extra state is therefore not practically recoverable from the available surface-temperature record. The one-node architecture is retained because it is the more parsimonious identifiable representation, not because all thermal transfers are accurate.

## Uniandes data-quality lesson

The cleanest ambient cell-1 cycle has median `R_app≈0.03588 ohm`; affected cycles and pulse records lie around 0.178–0.232 ohm. This reproducible separation is evidence of elevated series-path response in the complete experimental chain, not a measurement of intrinsic cell resistance. The affected records are useful nonideality case studies but are rejected as electrical-identification truth. The cleanest cycle remains limited by an absent measured timebase and by its loaded 1C/CC–CV nature.

## Transfer findings

Literal Panasonic-to-TRC transfer is incompatible. The equilibrium-voltage backbone differs from the loaded TRC reference by 188.882 mV RMSE, and the provisional-time full transfer gives 82.134 mV RMSE only through compensating dynamic terms while losing positive-electrode OCP support near low SOC. The minimum defensible adaptation is cell-specific capacity/state mapping plus a TRC-compatible voltage-reference/OCP basis. Dynamic and thermal re-identification would also require measured-time evidence that is not available.

## Thesis contribution

The central contribution is not a claim of universal parameter transfer. It is a reproducible evidence hierarchy: audit the experiment, separate calibration from withheld validation, preserve model-domain checks, diagnose transfer before adaptation, and report unsupported identification as a limitation. This converts imperfect experimental data into defensible conclusions without silently correcting voltage, clipping states, or promoting diagnostics to physical parameters.

## Remaining limitations

- The grouped diffusion-reduction literature derivation still requires manual citation verification.
- Panasonic electrode windows and positive OCP knots are effective, weakly identifiable representations.
- Dynamic parameters are conditional on HPPC pulses with start SOC at or above one third.
- Low-SOC and current-direction transfer remain limited.
- Thermal calibration uses measured-voltage heat; fully predictive heat inherits electrical-model error.
- The Uniandes files have no measured timebase, preventing defensible dynamic or thermal identification.
