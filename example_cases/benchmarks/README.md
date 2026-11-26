# Benchmarks

Published problems with known answers, one or more for each domain of the solver companion
(cfd-solvers-rs), each a complete case and the checks its run must pass. Where the
[case templates](../templates/README.md) show that every solver runs, these show that what it
gives is right.

Every expected value is worked out in [`make_benchmarks.py`](make_benchmarks.py), from the
closed form or the published value named beside it; the script writes the JSON files here.
Nothing is taken from another product: each case is written from the problem's statement.

| Benchmark | Domain | Reference | Checked |
|---|---|---|---|
| [Robertson's stiff reactions](chemistry-robertson.json) | Chemistry | Robertson (1966); Hairer and Wanner (1996) | the three concentrations at 40 s |
| [RC low-pass filter](circuits-rc-low-pass.json) | Circuits | Sedra and Smith | gain at 1 and 10 kHz; the step response at one time constant |
| [Sod's shock tube](compressible-sod-shock-tube.json) | Compressible flow | Sod (1978); Toro's exact Riemann solution | pressure and densities between the waves |
| [TE101 cavity resonance](em-cavity-te101.json) | EM | Pozar, *Microwave Engineering* | the first resonance |
| [Thin dipole at resonance](antennas-thin-dipole-resonance.json) | Antennas | Balanis; King | the resonant length and resistance |
| [AT2 phase-field strength](fracture-at2-strength.json) | Fracture | Pham, Amor, Marigo and Maurini (2011) | the peak stress |
| [Wall with heat generation](heat-wall-with-generation.json) | Heat | Incropera et al. | the steady mean temperature rise |
| [Natural convection, Ra = 1e3](heat-natural-convection-ra1e3.json), [1e4](heat-natural-convection-ra1e4.json), [1e5](heat-natural-convection-ra1e5.json) | Heat (3-D flow solver) | de Vahl Davis (1983) | the mean Nusselt number; the heat balance |
| [Landau damping](kinetic-landau-damping.json) | Kinetic | Canosa (1973) | the damping rate |
| [Taylor–Green vortex](lbm-taylor-green.json) | Flow (lattice Boltzmann) | Taylor and Green (1937) | the energy's viscous decay |
| [Spherical aberration](optics-spherical-aberration.json) | Optics | Hecht, *Optics*; Snell's law traced exactly | the rms spot at the paraxial focus |
| [Two-stream instability](plasma-two-stream.json) | Plasma | Birdsall and Langdon | the growth rate |
| [J2 plasticity, uniaxial](plasticity-j2-uniaxial.json) | Plasticity | Simo and Hughes | stress and plastic strain at 0.4 % |
| [Mismatched line](rf-mismatched-line.json) | RF | Pozar, ABCD parameters | S21 at two frequencies |
| [PN junction](semiconductors-pn-junction.json) | Semiconductors | Sze and Ng | the built-in potential; the ideal-diode growth |
| [d² law](spray-d2-law.json) | Spray | Godsave (1953); Turns | the evaporated share |
| [Bar modes](structures-bar-modes.json) | Structures | Rao, *Mechanical Vibrations* | four axial modes |
| [Cantilever under its own weight](structures-cantilever-self-weight.json) | Structures | Timoshenko and Gere | the tip deflection |

Each file is JSON: `title`, `domain`, `shows`, `source`, the `case` as the companion takes it,
and `checks`. A check's `value` is a run metric (`metric`), a value in a JSON result
(`json` and `path`), or a column of a CSV result (`file`, `column`) in a row (`row`: `first`,
`last`, `max`, `min`) or interpolated where another column has a value (`at`); then divided by
another such value (`over`), offset (`add`) and scaled (`divide_by`). It passes within
`tolerance` (a share) of `expect`, or between `at_least` and `at_most`.

Run them:

```sh
# Each benchmark is a case SolidCraft accepts (part of the tests):
cargo test -p solidcraft-multiphysics-host --lib benchmarks

# Each run by a companion and judged; SOLIDCRAFT_BENCHMARKS=chemistry,optics runs only those:
SOLIDCRAFT_TEST_CFD_SOLVE=/path/to/cfd-solve SOLIDCRAFT_TEST_OUTPUT_DIR=/tmp/benchmarks \
    cargo test --release -p solidcraft-multiphysics-host --lib real_benchmarks -- --ignored

# The same without building SolidCraft:
python3 example_cases/benchmarks/run_benchmarks.py /path/to/cfd-solve /tmp/benchmarks
```

[`results.md`](results.md) holds the last run's values. To add a benchmark, add it to
`make_benchmarks.py` with its source and its derivation, run the script, then both commands.

## What the benchmarks found

- **The phase-field bar broke at its ends.** Its damage update counted an end cell's own damage
  as the neighbour it lacks, so the ends' damage jumped to 1 at once and the bar's strength
  came out 4 to 6 % low. Fixed in the companion; the bar now peaks at the AT2 strength.
- **Landau damping needs fine cells.** The solver's first-order upwind steps add damping that
  falls with the cell size: 64 cells across the wave give +24 %, 512 give +3.5 %.
- **The cavity's tetrahedra need cubes.** Long, flat boxes split into tetrahedra converge far
  more slowly: boxes six times longer than wide put the TE101 resonance 6 % low; 10 mm cubes
  put it 0.3 % low.
- **The diode solver fails on devices longer than 2 µm** (its Gummel iteration diverges and a
  pivot vanishes), so the junction here is short, and its capacitance under reverse bias,
  whose depletion would reach the contacts, is not checked.
- **A thin dipole's 73 + j42.5 Ω** is the induced-EMF estimate for a sinusoidal current; the
  moment method converges to about 81 + j47 Ω at exactly half a wavelength, so the benchmark
  checks the resonance (0.48 λ, about 70 Ω) instead.
