# Results of the last run

Every benchmark in this folder, run through `real_benchmarks_meet_their_references` (cfd-solvers-rs 0.18.0 with the phase-field fix, 27 November 2026).
Each check's value, what it was held against, and whether it passed.

| Benchmark | Domain | Check | Value | Expected | Off by | Pass | Run (s) |
|---|---|---|---|---|---|---|---|
| A Gaussian pressure pulse in water | Acoustics | pressure 3.05 mm from the centre (Pa) | 102.2 | 103 | -0.86 % | yes | 25 |
|  |  | pressure 3.25 mm from the centre (Pa) | 107.8 | 107.3 | +0.49 % | yes |  |
|  |  | pressure 3.45 mm from the centre (Pa) | 61.6 | 60.68 | +1.51 % | yes |  |
| A thin dipole at its resonance | Antennas | length at resonance, in wavelengths | 0.4815 | 0.48 | +0.31 % | yes | 10 |
|  |  | input resistance at resonance (Ω) | 72.04 | 70 | +2.91 % | yes |  |
| A cold Lennard-Jones chain | Atomistics | potential energy per atom less T/2 (ε) | -1.03 | -1.031 | +0.08 % | yes | 1 |
| A cell discharged, then left to rest | Batteries | charge delivered (Ah) | 0.025 | 0.025 | +0.00 % | yes | 0 |
|  |  | recovery in the rest's last half hour (V): it has settled | 0 | bound | bound | yes |  |
|  |  | recovery over the whole rest (V): the voltage comes back up | 0.02383 | bound | bound | yes |  |
| Tissue heated evenly against blood perfusion | Bioheat | rise after one time constant (K) | 1.263 | 1.265 | -0.16 % | yes | 0 |
|  |  | rise after ten time constants (K) | 2.002 | 2.002 | -0.00 % | yes |  |
| Robertson's stiff reactions | Chemistry | A at 40 s | 0.7158 | 0.7158 | +0.00 % | yes | 1 |
|  |  | B at 40 s | 9.186e-06 | 9.186e-06 | +0.01 % | yes |  |
|  |  | C at 40 s | 0.2841 | 0.2842 | -0.01 % | yes |  |
| An RC low-pass filter | Circuits | gain at 1 kHz (dB) | -1.445 | -1.445 | -0.00 % | yes | 1 |
|  |  | gain at 10 kHz (dB) | -16.07 | -16.07 | +0.00 % | yes |  |
|  |  | output one time constant after the step (V) | 3.16 | 3.161 | -0.03 % | yes |  |
| Sod's shock tube | Compressible flow | pressure between the rarefaction and the shock, x = 0.6 | 0.3031 | 0.3031 | -0.00 % | yes | 1 |
|  |  | density left of the contact, x = 0.6 | 0.4263 | 0.4263 | -0.00 % | yes |  |
|  |  | density between the contact and the shock | 0.2656 | 0.2656 | -0.00 % | yes |  |
| The TE101 resonance of a rectangular cavity | EM | first resonance (Hz) | 3.342e+09 | 3.352e+09 | -0.29 % | yes | 160 |
| The lid-driven cavity at Re = 100 | Flow | u at y = 0.8516 on the vertical centreline (m/s) | 0.2297 | 0.2315 | -0.78 % | yes | 228 |
|  |  | u at y = 0.6172 on the vertical centreline (m/s) | -0.1341 | -0.1364 | +1.72 % | yes |  |
|  |  | u at y = 0.5 on the vertical centreline (m/s) | -0.1974 | -0.2058 | +4.11 % | yes |  |
|  |  | u at y = 0.4531 on the vertical centreline (m/s) | -0.2002 | -0.2109 | +5.05 % | yes |  |
|  |  | u at y = 0.2813 on the vertical centreline (m/s) | -0.1514 | -0.1566 | +3.32 % | yes |  |
|  |  | u at y = 0.1719 on the vertical centreline (m/s) | -0.102 | -0.1015 | -0.48 % | yes |  |
| The strength of a phase-field (AT2) bar | Fracture | peak stress (Pa, on a unit section) | 0.07262 | 0.07262 | +0.00 % | yes | 0 |
| A sphere resting on a plane | Granular | how far it sinks, as a share of Hertz's | 1 | 1 | +0.00 % | yes | 0 |
| Natural convection in a square cavity, Ra = 1e3 | Heat | mean Nusselt number on the hot wall | 1.131 | 1.118 | +1.17 % | yes | 51 |
|  |  | heat out of the cold wall over heat into the hot one | -1 | -1 | -0.00 % | yes |  |
| Natural convection in a square cavity, Ra = 1e4 | Heat | mean Nusselt number on the hot wall | 2.277 | 2.243 | +1.54 % | yes | 34 |
|  |  | heat out of the cold wall over heat into the hot one | -1 | -1 | +0.00 % | yes |  |
| Natural convection in a square cavity, Ra = 1e5 | Heat | mean Nusselt number on the hot wall | 4.588 | 4.519 | +1.52 % | yes | 531 |
|  |  | heat out of the cold wall over heat into the hot one | -1 | -1 | -0.00 % | yes |  |
| A copper wall releasing heat, held at 300 K on both faces | Heat | mean rise above the faces (K) | 3.513 | 3.509 | +0.12 % | yes | 0 |
| A Gaussian carried once round by DG | High order | largest difference from the exact profile after one period | 1.417e-06 | bound | bound | yes | 0 |
|  |  | mass gained or lost, relative | 9.159e-16 | bound | bound | yes |  |
| Linear Landau damping at k = 0.5 | Kinetic | damping rate of the field (1/ω_p) | -0.1587 | -0.1533 | -3.53 % | yes | 338 |
| Decay of the Taylor–Green vortex | Flow | energy left, as a share of the start | 0.2134 | 0.2139 | -0.25 % | yes | 0 |
| A spinodal wave growing in a Cahn–Hilliard mixture | Microstructure | the wave's amplitude at the end | 0.001988 | 0.002008 | -1.02 % | yes | 0 |
| Spherical aberration of a refracting sphere | Optics | rms spot radius at the paraxial focus (m) | 3.639e-06 | 3.639e-06 | +0.00 % | yes | 0 |
| A silicon slab between silica, in a box | Photonics | fundamental effective index | 2.911 | 2.915 | -0.15 % | yes | 2 |
| The cold two-stream instability | Plasma | growth rate (1/s) | 0.3279 | 0.3536 | -7.26 % | yes | 3 |
| Uniaxial loading past yield with linear hardening | Plasticity | stress at 0.4 % strain (Pa) | 2.554e+08 | 2.554e+08 | +0.00 % | yes | 0 |
|  |  | equivalent plastic strain | 0.002723 | 0.002723 | +0.00 % | yes |  |
| A slug of solute spreading along a column | Porous media | solute at x = 0.5 | 0.7363 | 0.7364 | -0.02 % | yes | 0 |
|  |  | solute at x = 0.56 | 0.4045 | 0.4046 | -0.03 % | yes |  |
|  |  | solute at x = 0.6 | 0.1314 | 0.1314 | +0.05 % | yes |  |
| A five-bus load flow | Power networks | voltage at north (p.u.) | 1.011 | 1.011 | +0.00 % | yes | 0 |
|  |  | voltage at lake (p.u.) | 0.9958 | 0.9958 | +0.00 % | yes |  |
|  |  | voltage at main (p.u.) | 0.9917 | 0.9917 | +0.00 % | yes |  |
|  |  | voltage at elm (p.u.) | 0.9786 | 0.9786 | +0.00 % | yes |  |
|  |  | angle at north (°) | -4.11 | -4.11 | +0.00 % | yes |  |
|  |  | angle at elm (°) | -7.49 | -7.49 | +0.00 % | yes |  |
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
| Decaying turbulence in the k–ε model | Turbulence | turbulent kinetic energy after 2 s (m²/s²) | 0.009841 | 0.009842 | -0.01 % | yes | 0 |
