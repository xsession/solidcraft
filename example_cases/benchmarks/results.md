# Results of the last run

Every benchmark in this folder, run through `real_benchmarks_meet_their_references` (cfd-solvers-rs 0.18.0 with the phase-field fix, 26 November 2026).
Each check's value, what it was held against, and whether it passed.

| Benchmark | Domain | Check | Value | Expected | Off by | Pass | Run (s) |
|---|---|---|---|---|---|---|---|
| A thin dipole at its resonance | Antennas | length at resonance, in wavelengths | 0.4815 | 0.48 | +0.31 % | yes | 20 |
|  |  | input resistance at resonance (Ω) | 72.04 | 70 | +2.91 % | yes |  |
| Robertson's stiff reactions | Chemistry | A at 40 s | 0.7158 | 0.7158 | +0.00 % | yes | 0 |
|  |  | B at 40 s | 9.186e-06 | 9.186e-06 | +0.01 % | yes |  |
|  |  | C at 40 s | 0.2841 | 0.2842 | -0.01 % | yes |  |
| An RC low-pass filter | Circuits | gain at 1 kHz (dB) | -1.445 | -1.445 | -0.00 % | yes | 0 |
|  |  | gain at 10 kHz (dB) | -16.07 | -16.07 | +0.00 % | yes |  |
|  |  | output one time constant after the step (V) | 3.16 | 3.161 | -0.03 % | yes |  |
| Sod's shock tube | Compressible flow | pressure between the rarefaction and the shock, x = 0.6 | 0.3031 | 0.3031 | -0.00 % | yes | 1 |
|  |  | density left of the contact, x = 0.6 | 0.4263 | 0.4263 | -0.00 % | yes |  |
|  |  | density between the contact and the shock | 0.2656 | 0.2656 | -0.00 % | yes |  |
| The TE101 resonance of a rectangular cavity | EM | first resonance (Hz) | 3.342e+09 | 3.352e+09 | -0.29 % | yes | 126 |
| The strength of a phase-field (AT2) bar | Fracture | peak stress (Pa, on a unit section) | 0.07262 | 0.07262 | +0.00 % | yes | 0 |
| Natural convection in a square cavity, Ra = 1e3 | Heat | mean Nusselt number on the hot wall | 1.131 | 1.118 | +1.17 % | yes | 45 |
|  |  | heat out of the cold wall over heat into the hot one | -1 | -1 | -0.00 % | yes |  |
| Natural convection in a square cavity, Ra = 1e4 | Heat | mean Nusselt number on the hot wall | 2.277 | 2.243 | +1.54 % | yes | 31 |
|  |  | heat out of the cold wall over heat into the hot one | -1 | -1 | +0.00 % | yes |  |
| Natural convection in a square cavity, Ra = 1e5 | Heat | mean Nusselt number on the hot wall | 4.588 | 4.519 | +1.52 % | yes | 466 |
|  |  | heat out of the cold wall over heat into the hot one | -1 | -1 | -0.00 % | yes |  |
| A copper wall releasing heat, held at 300 K on both faces | Heat | mean rise above the faces (K) | 3.513 | 3.509 | +0.12 % | yes | 0 |
| Linear Landau damping at k = 0.5 | Kinetic | damping rate of the field (1/ω_p) | -0.1587 | -0.1533 | -3.53 % | yes | 437 |
| Decay of the Taylor–Green vortex | Flow | energy left, as a share of the start | 0.2134 | 0.2139 | -0.25 % | yes | 1 |
| Spherical aberration of a refracting sphere | Optics | rms spot radius at the paraxial focus (m) | 3.639e-06 | 3.639e-06 | +0.00 % | yes | 0 |
| The cold two-stream instability | Plasma | growth rate (1/s) | 0.3279 | 0.3536 | -7.26 % | yes | 4 |
| Uniaxial loading past yield with linear hardening | Plasticity | stress at 0.4 % strain (Pa) | 2.554e+08 | 2.554e+08 | +0.00 % | yes | 0 |
|  |  | equivalent plastic strain | 0.002723 | 0.002723 | +0.00 % | yes |  |
| A 75 Ω line between 50 Ω ports | RF | S21 real at 500 MHz | 0.8589 | 0.8589 | +0.00 % | yes | 0 |
|  |  | S21 imag at 500 MHz | -0.474 | -0.474 | -0.00 % | yes |  |
|  |  | S21 real at 2000 MHz | -0.2668 | -0.2668 | +0.00 % | yes |  |
|  |  | S21 imag at 2000 MHz | -0.8892 | -0.8892 | -0.00 % | yes |  |
| An abrupt silicon PN junction | Semiconductors | built-in potential (V) | 0.6132 | 0.6132 | -0.00 % | yes | 0 |
|  |  | current at 0.4 V over the current at 0.3 V | 42.55 | bound | bound | yes |  |
| Droplets evaporating by the d² law | Spray | evaporated share of the injected liquid | 0.6449 | 0.6447 | +0.03 % | yes | 1 |
| Axial modes of a fixed–free aluminium bar | Structures | mode 1 (Hz) | 2546 | 2546 | +0.00 % | yes | 0 |
|  |  | mode 2 (Hz) | 7638 | 7638 | +0.00 % | yes |  |
|  |  | mode 3 (Hz) | 1.273e+04 | 1.273e+04 | +0.01 % | yes |  |
|  |  | mode 4 (Hz) | 1.782e+04 | 1.782e+04 | +0.01 % | yes |  |
| A steel cantilever under its own weight | Structures | largest displacement, at the free end (m) | 5.755e-05 | 5.836e-05 | -1.38 % | yes | 0 |
