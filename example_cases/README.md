# Example cases

Cases that show what SolidCraft Multiphysics and its solver companion (cfd-solvers-rs) can do,
each one checked by a test. Two kinds:

- [`model/`](model/): **studies of the model**, built from SolidCraft geometry with SolidCraft
  commands. One example for each capability of the 3-D flow solver a study of the model uses,
  and, where there is one, an exact or textbook answer its run must match.
- [`templates/`](templates/README.md): **case templates**, one small complete case for each of
  the companion's solvers (56 of them: flow, heat, structures, electromagnetics, plasma,
  chemistry, batteries, circuits and more), run as they come.

Nothing here comes from another product: every geometry is drawn by the commands in the
example, and every expected value is worked out from physics in the example's own words.

## Studies of the model

Each file in `model/` is JSON:

| Field | What it is |
|---|---|
| `title`, `shows` | What the example is and what it demonstrates |
| `capabilities` | The solver capabilities its case must ask for and its run must report |
| `steps` | SolidCraft commands, in order, with their parameters. `"$study"` is the study the steps made, `"$revision"` its revision when the step runs |
| `checks` | What its run must give: a metric of the run (`metric`, less `minus`, over `divide`, its size with `abs`) within `tolerance` (a share) of `expect`, or between `at_least` and `at_most` |

| Example | Capability | Check |
|---|---|---|
| [01 Square duct](model/01-square-duct.json) | air in a study of the model, scaled momentum, point and global goals | pressure drop against the exact laminar value (f·Re = 56.91) |
| [02 Round pipe, cut cells](model/02-round-pipe-cut-cells.json) | cut cells | pressure drop per unit of flow against Hagen–Poiseuille |
| [03 Local refinement](model/03-local-refinement.json) | cells split around a body | all the inflow leaves |
| [04 Board, conjugate heat](model/04-board-conjugate-heat.json) | conjugate heat, partial cells | heat out equals heat released |
| [05 Cold-plate stack](model/05-cold-plate-stack.json) | heat on a face, a fixed temperature on a solid-only side | energy balance; the stack's temperature drop |
| [06 Fan curve](model/06-fan-curve.json) | a fan on a face | the operating point against the duct's laminar curve |
| [07 Porous brick](model/07-porous-brick.json) | a porous body | pressure drop against Darcy |
| [08 Jet from a face](model/08-jet-from-a-face.json) | a condition on a face inside the domain | the jet's flow leaves |
| [09 Natural convection](model/09-natural-convection.json) | buoyancy | the heat sets the air moving |
| [10 Fan inside a duct](model/10-internal-fan.json) | a fan inside the domain | the operating point against the duct's laminar curve; all its flow leaves |
| [11 Turning cylinder](model/11-turning-cylinder.json) | a rotating region (the Magnus effect) | higher pressure where the surface moves against the stream |

Run them:

```sh
# Each example builds and asks for what it shows (no solver needed; part of the tests):
cargo test -p solidcraft-multiphysics-host --lib example_cases

# Each example run by a companion and judged by its checks:
SOLIDCRAFT_TEST_CFD_SOLVE=/path/to/cfd-solve SOLIDCRAFT_TEST_OUTPUT_DIR=/tmp/examples \
    cargo test --release -p solidcraft-multiphysics-host --lib real_model_examples -- --ignored
```

`SOLIDCRAFT_EXAMPLES=01,07` runs only those. The output folder gets `model-examples.json`, what
each example measured; [`model/results.md`](model/results.md) holds the last run's.

To add an example: write its file, run the first command (it must build and ask for its
capabilities), then the second (its checks must pass), and add it to the table above.

## Case templates

```sh
SOLIDCRAFT_TEST_CFD_SOLVE=/path/to/cfd-solve SOLIDCRAFT_TEST_OUTPUT_DIR=/tmp/templates \
    cargo test --release -p solidcraft-multiphysics-host --lib real_reference_runner_all_templates -- --ignored
python3 example_cases/templates/summarise_runs.py /tmp/templates /tmp/templates.json
python3 example_cases/templates/make_index.py /tmp/templates.json
```

## What is not here

The 3-D flow is laminar (the companion's turbulence models run on their own, in the `rans-*`
templates), single-phase in a study of the model (two-phase flow is in the `multiphase-*`
templates), and on a box grid of at most 32 768 cells. A rotating region is a turning frame for
steady turning, not a mesh that slides. Thermal devices (heat pipes, thermoelectric coolers) are
in the solver but not yet offered in a study of the model. See the Multiphysics documentation's *What it does not do* and *Training
cases* pages.
