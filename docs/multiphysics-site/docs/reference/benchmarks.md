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

All 32 pass (solver companion 0.18.0, with the phase-field fix below).

| Benchmark | Reference | Result |
|---|---|---|
| Gaussian pulse in water | the 2-D wave's exact solution | within 1.5 % across the ring's crest |
| Cold Lennard-Jones chain | lattice sum and equipartition | −0.08 % |
| Discharge and rest | charge conservation; relaxation | charge exact; the voltage recovers 24 mV and settles |
| Uniform heating against perfusion | Pennes (1948) | −0.16 % after one time constant |
| Lid-driven cavity, Re = 100 | Ghia, Ghia and Shin (1982) | within 0.5 to 5 % down the centreline |
| Hertzian sphere at rest | Johnson | exact |
| DG advection | Hesthaven and Warburton | 1.4e-6 from the exact profile; error falls as h^4.9 |
| Spinodal growth | Cahn (1961) | −1.0 % |
| Slab waveguide | Okamoto; Marcuse | −0.15 % |
| Diffusing slug | Crank | within 0.05 % |
| Five-bus load flow | Newton–Raphson solved independently | within 1e-8 p.u. |
| k–ε decay | Launder and Spalding | −0.01 % |
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
- **The acoustic solver is first-order in time**: it starts its velocity at t = 0 rather than
  half a step before, so a pulse comes out 3 % low at a Courant number of 0.3, within 1 % at 0.075.
- **The 3-D flow solver's convection is first-order upwind**, the only scheme its momentum
  equation offers: the lid-driven cavity's vortex core turns 4 to 5 % slow on 65 × 65 cells.
- **The two-fluid (Euler–Euler) model has no pressure**: drag carries the liquid along with the
  bubbles, so a bubble column's published behaviour does not apply.
- **A thin dipole's 73 + j42.5 Ω** holds only for a sinusoidal current; the moment method's
  impedance at exactly half a wavelength is about 81 + j47 Ω, and its resonance matches the
  published one.

Every domain of the solver has a benchmark except multiphase flow (see above). Within some
domains, solvers still without one: the 1-D Maxwell and 3-D FDTD solvers, Joule heating, the
power-plane IR drop, galvanic corrosion, the Cosserat rod, the depth-dose transport, the 2-D heat
and Poisson solvers and scalar transport.
