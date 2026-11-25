# Training cases

A set of 33 public CFD training cases (the *holzmann-cfd.com* community cases, made for
OpenFOAM) was worked through to see what SolidCraft Multiphysics can do with each. This page
says how, what came out, and what is out of reach today and why.

## How the cases were built

The cases are published under the GNU GPL 3. SolidCraft takes nothing from them but shapes and
sizes:

- **Their setup files were not read.** Boundary conditions, fluids, solver settings and meshes
  are OpenFOAM dictionaries; SolidCraft's rules keep GPL code out of the project, read or
  copied, so none of them was opened.
- **Their geometry was measured, not copied.** Where a case comes with STL or STEP files, a
  short script of ours (`docs/multiphysics-site/tools/stl_parts.py`) measured each part's place
  and size, and SolidCraft drew the parts again with exact shapes: a battery cell as a cylinder,
  a channel as a box. Only the battery case's air channel, an irregular shape, was used as a
  surface, joined with its openings into one closed solid. No geometry file and no design that
  holds one is in the repository.
- **Settings are ours**: air or water at 20 °C, speeds that keep the flow laminar, powers and
  temperatures stated with each case. The detail pages give most of their own parameters as
  arbitrary; where they give one (the battery's 100 W per cell), it is said whether it was used.
- **Our own cases stand in where a classic shape has an exact answer**: a pipe, a square duct
  and a thin slot have exact laminar pressure gradients, which the downloaded cases' meshing
  tutorials (a pipe, a thin gap) do not report.

Everything was run in SolidCraft through its control channel by
`docs/multiphysics-site/tools/make_training_cases.py`, with the solver companion 0.18.0, and
the numbers below are that script's output.

## Results

All runs are laminar air at 20 °C (1.5e-5 m²/s, 1.2 kg/m³) on SolidCraft's box grid, at most
32 768 cells.

### Straight passages: the exact laminar pressure gradient

Where the flow is fully developed, its pressure falls at
dp/dx = (f·Re)·ν·U / (2·D<sub>h</sub>²), with f·Re = 64 for a round pipe, 56.91 for a square duct
and 82.34 for a 1 : 8 slot (Shah and London, *Laminar Flow Forced Convection in Ducts*, 1978).
Each case measures the gradient between two points well past the inlet.

| Case | Stands in for | Grid | Re | Measured / exact | Off by |
|---|---|---|---|---|---|
| Square duct, 8 × 8 mm | (exact shape on the grid) | 0.8 mm, 10 cells across | 27 | 0.327 / 0.333 m/s² | **−1.8 %** |
| Round pipe, 16 mm bore | Meshing a pipe (001) | 1.5 mm, 10.7 across | 53 | 0.096 / 0.092 | **+4.3 %** |
| Round pipe, cut cells | Meshing a pipe (001) | as above | 53 | 0.088 / 0.094 | **−6.6 %** |
| Thin gap, 2 × 16 mm slot | Thin gap meshing (006) | 0.4 mm, 5 across the gap | 12 | 2.21 / 2.44 | **−9.4 %** |

The square duct, which the grid holds exactly, lands within 2 %. The round pipe's staircase
wall costs a few per cent either way, and five cells across a gap leave about 10 %: the usual
rule of ten cells across the narrowest passage holds here too.

### A fan against a duct

A fan whose curve falls straight from 0.1 Pa at no flow to nothing at 1.2 L/min blows into the
square duct. Where its curve crosses the duct's laminar system curve, for fully developed flow,
it should deliver 8.0e-6 m³/s; the run settles at **7.2e-6 m³/s, 10 % less**. The system
curve leaves out the extra loss where the flow enters the duct and develops, which the run has
and which lowers the flow, so the difference is mostly the estimate's (stands in for *Fan
rotation and AMI*, 025, whose rotating fan is out of reach, see below).

### Heat: the energy balance and a conduction stack

| Case | Stands in for | What is checked | Result |
|---|---|---|---|
| Battery pack: six 46 × 80 mm cells, 5 W each, in a duct at 0.2 m/s | Battery cooling (029) | Heat carried out equals heat released | 30.000 W of 30 W, imbalance 1e-11; cells 63 to 78 °C at this low speed |
| TEG stack: the test device's parts, measured, 5 W in the heater, 20 °C under the cooler | Thermal-electric generator (023) | Drop across the 2 mm layer = P·L/(k·A) = 10.7 K | **12.8 K, +19 %**; energy imbalance 6e-10 |

The TEG layer is 2 mm thick in cells of 4.9 mm: it is held by partial cells, which mix the
layer and the copper beside it in series, exactly for heat across the layer. The 19 % comes
from the copper blocks' edges: 25 mm blocks are about five cells across and do not line up with
the cells, so the cells round their edges mix copper with air and pass less heat. Smaller cells
around the stack would bring it down.

### Spheres in a channel

| Case | Grid | Re | Wake behind the sphere |
|---|---|---|---|
| 10 mm sphere in a 30 mm channel (ours) | 1.7 mm, 6 cells across the sphere | 67 | none resolved |
| The case's 20 mm sphere in its 100 × 100 × 500 mm channel | 8 mm, 2.5 cells across | 100 | none resolved |

At Re 100 a sphere's closed wake reaches about 0.9 diameters behind it (Johnson and Patel,
*J. Fluid Mech.* 378, 1999). Neither grid resolves it: the case's long channel needs 8 mm cells
under the 32 768-cell limit, and a local refinement around the sphere, which would split whole
rows of the box grid, does not fit beside it. This is the box grid's limit at its clearest; a
refinement that splits only the cells near the sphere (an octree) is the remedy.

### The battery case's own channel

The case's air channel, joined with its inlet and outlet into one closed solid, converts and
cuts cleanly, but setting up the study on the resulting 10 000-facet wall took longer than an
hour: SolidCraft finds the exact face under each triangle of the wall, which grows with
triangles times faces. It was stopped; a faster lookup is a follow-up. The battery pack above
uses the case's cell size and pitch with our own duct.

## What the cases needed that SolidCraft lacked

Working the cases through found and fixed these in SolidCraft and its solver; each fix has its
own tests:

- **Studies had no real fluid.** A study of the model ran with its template's stand-in
  (0.02 m²/s, some 1 300 times as viscous as air) for 12 steps; a study now has a fluid (air by
  default) and an end time, and stops once it settles (see [the solver](solver.md)).
- **Small, slow flows froze.** The linear solves stopped at an absolute floor, so with
  millimetre cells the flow never moved: the square duct came out 42 % off before the fix
  (companion 0.18.0, `scaled_momentum_v1`) and 1.8 % after.
- **Curved and imported bodies were too finely tessellated** for a case (a 20 mm sphere in
  45 000 triangles); bodies are now tessellated for the study's cells.
- **A fixed temperature on a side only solid touches** (a cold plate under the model) was
  refused; it is accepted now.
- **Runs were stopped after five minutes**; they may take two hours, or
  `SOLIDCRAFT_CFD_TIMEOUT_S`.
- **Large heated runs passed the 16 MiB result limit**; a study writes its first and last flow
  frames only.
- **A fan's study took a needlessly short time step**; it now follows the fan's flow.
- **A design opened in a second tab could not run its study**: the solver was given to the
  first design only.

## Cases out of reach today

| Case | Needs | In SolidCraft today |
|---|---|---|
| Falling droplets, fluidic oscillator, sneezing simulation | Two-phase flow (a free surface, droplets) | Single-phase flow only |
| Kaplan turbine, vertical-axis wind turbine | A rotor turned by the flow (six degrees of freedom) | Rotating regions exist in the solver; flow-driven motion does not |
| Fan rotation and AMI, rotating rotor, 2-D AMI, rotating inlet, water pump, complex mesh motion | Sliding or moving meshes between rotating and fixed parts | A box grid; a rotating region is not yet offered in a study of the model |
| ACMI heat transfer | Heat across a moving, connecting and disconnecting interface | Fixed geometry |
| Tank with safety valve | A baffle that opens at a pressure | Fixed geometry |
| Combustion chamber (cold flow) | Mixing of a passive scalar | No scalar transport in the 3-D flow |
| Catalyst heat-up | Porous zones in a heated flow | Porous zones exist in the solver; a study of the model does not offer them yet |
| Solar chimney, gin and tonic | Natural convection (buoyancy); the drink also melting ice | Buoyancy exists in the solver; a study of the model does not offer it yet. No melting |
| Magnus effect | A cylinder whose wall turns | No moving walls in a study of the model |
| 2-D adaptive mesh refinement (blastFoam) | Compressible, shock-capturing flow | Incompressible flow |
| Pitot tube | Compressible flow up to 50 m/s, transient | Incompressible flow |
| Geometric variation, Tesla valve (DAKOTA) | An optimiser driving the geometry | Sweeps over one value; no optimiser |
| Fluent mesh conversion | Reading a Fluent mesh | Studies are built from the design's bodies |
| Suzanne's head, sphere and layer, feature edges, cell zones, helix | snappyHexMesh body-fitted meshes with boundary layers | A box grid with local refinement and cut cells; the pipe, the gap and the sphere are run above on it |
| Battery cooling at 100 W per cell | Turbulent air at about 10 m/s | No turbulence model in the 3-D flow: run laminar at 5 W per cell instead |

None of these is a small step; each is listed in the [limits](limits.md) or would be a feature
of its own.
