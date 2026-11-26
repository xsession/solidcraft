# Results of the last run

Every example in this folder, run by cfd-solvers-rs 0.18.0 through `real_model_examples_meet_their_checks`
(26 November 2026, laminar air at 20 °C). Each check's value, what it was held against, and
whether it passed.

| Example | Check | Value | Expected | Off by | Pass | Run (s) |
|---|---|---|---|---|---|---|
| Laminar flow in a square duct | pressure drop over 22 mm against the exact laminar value (Pa) | 0.008641 | 0.008803 | -1.8 % | yes | 175 |
|  | the fastest air, about twice the mean in a developed square duct (m/s) | 0.1002 | bound | – | yes |  |
| A round pipe on a box grid, with cut cells | pressure drop per unit of flow against Hagen–Poiseuille (Pa·s/m³) | 416.9 | 447.6 | -6.9 % | yes | 276 |
| Cells split around a block in a duct | air out through the outlet: all that comes in, 0.2 m/s over 32 × 32 mm (m³/s) | 0.0002048 | 0.0002048 | -0.0 % | yes | 29 |
| A powered chip on a board in a duct | power in the study (W) | 0.5 | 0.5 | +0.0 % | yes | 124 |
|  | heat out less heat released, as a share of it | 1.694e-09 | bound | – | yes |  |
|  | the chip is warmer than the inflow (K) | 328.3 | bound | – | yes |  |
| A heater on a stack, cooled by a cold plate | power in the study (W) | 5 | 5 | -0.0 % | yes | 4 |
|  | heat out less heat released, as a share of it | 2.097e-08 | bound | – | yes |  |
|  | the heater, above the 20 °C plate by the stack's drop (K) | 308 | bound | – | yes |  |
| A fan blowing into a duct | the fan's flow against its curve crossing the duct's laminar curve (m³/s) | 7.176e-06 | 7.998e-06 | -10.3 % | yes | 191 |
| A porous brick filling a duct | pressure drop across the brick against Darcy (Pa) | 0.2634 | 0.27 | -2.4 % | yes | 111 |
| A jet blown from a face inside the domain | air out: the jet's 0.2 m/s over 10 × 10 mm (m³/s) | 2e-05 | 2e-05 | -0.0 % | yes | 100 |
| Air rising off a heated block in a sealed box | the air moves, driven only by the heat (m/s) | 0.2455 | bound | – | yes | 373 |
