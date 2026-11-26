# Benchmarks

The repository's `example_cases/benchmarks` folder holds published problems with known answers,
one or more for each domain of the solver: a stiff chemistry problem, a shock tube, a cavity
resonance, a thin dipole, natural convection in a square cavity and more. Each is a complete
case with the checks its run must pass, and every expected value is worked out, from the closed
form or the published value named with it, by the folder's `make_benchmarks.py`. A test checks
that SolidCraft accepts every case; another runs them through the solver and judges them.

Where the case templates show that every solver runs, the benchmarks show that what it gives is
right, and how close.

## The last run

All 20 pass (solver companion 0.18.0, with the phase-field fix below).

| Benchmark | Reference | Result |
|---|---|---|
| Robertson's stiff reactions | Hairer and Wanner | A, B and C at 40 s within 0.01 % |
| RC low-pass filter | closed form | gain and step response within 0.03 % |
| Sod's shock tube | Toro's exact Riemann solution | pressure and densities between the waves within 0.01 % |
| TE101 cavity resonance | Pozar | −0.29 % |
| Thin dipole at resonance | Balanis; King | 0.4815 λ long (0.48), 72.0 Ω (about 70) |
| AT2 phase-field strength | Pham et al. (2011) | exact |
| Natural convection, Ra = 1e3, 1e4, 1e5 | de Vahl Davis (1983) | Nusselt number +1.2, +1.5, +1.5 % |
| Wall with heat generation | Incropera et al. | +0.12 % |
| Landau damping | Canosa (1973) | damping rate 3.5 % too fast, on 512 cells |
| Taylor–Green vortex | Taylor and Green (1937) | −0.25 % |
| Spherical aberration | Snell's law, traced exactly | exact |
| Two-stream instability | Birdsall and Langdon | growth rate −7.3 % |
| J2 plasticity, uniaxial | Simo and Hughes | exact |
| Mismatched transmission line | Pozar | exact |
| PN junction | Sze and Ng | built-in potential exact; ideality 1.03 |
| d² law | Godsave; Turns | +0.03 % |
| Fixed–free bar modes | Rao | within 0.01 % |
| Cantilever under its own weight | Timoshenko and Gere | −1.4 % |

## What they found

- **The phase-field bar broke at its ends first**, so its strength came out 4 to 6 % low: the
  damage update counted an end cell's own damage as the neighbour it lacks. Fixed; the bar now
  peaks at the AT2 strength exactly.
- **Landau damping needs fine cells**: the solver's first-order upwind steps add damping, +24 %
  on 64 cells across the wave and +3.5 % on 512.
- **The cavity's tetrahedra need near-cubic boxes**: boxes six times longer than wide put the
  resonance 6 % low.
- **The diode solver fails on devices longer than 2 µm**, so reverse-bias capacitance is not
  yet checked.
- **A thin dipole's 73 + j42.5 Ω** holds only for a sinusoidal current; the moment method's
  impedance at exactly half a wavelength is about 81 + j47 Ω, and its resonance matches the
  published one.

Domains with no benchmark yet: acoustics, atomistics, batteries, bioheat, the Cahn–Hilliard
equation, discontinuous Galerkin, granular flow, multiphase flow, photonics, power networks,
porous columns, turbulence and the lid-driven cavity (which needs velocity probes in the 3-D
flow solver's results).
