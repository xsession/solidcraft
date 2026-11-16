# Multiphysics workspace

Multiphysics adds saved numerical studies, actual companion execution and portable results under a separate workspace and command namespace.
The planned Simulation entry is unchanged. Existing Design and Render defaults are unchanged.

## Layout

The workspace is laid out like a CAD-embedded flow-analysis tool:

| Area | What it holds |
|---|---|
| Ribbon, **Analysis** tab | Project (Wizard, New project, General settings, Calculation control, Delete project), Input (Computational domain, Boundary condition, Sources, Engineering database, Goals, Coupling), Mesh (Mesh settings, Check project), Solve (Run, Solver monitor, Stop), Tools (Templates, Verification, Export case) |
| Ribbon, **Results** tab | Load results; Plots (Cut plots, Surface plots, Isosurfaces, Flow trajectories, Particle studies); Parameters (Point, Surface, Volume); Charts (XY plots, Goal plots); Output (Report, Animation, Compare, Export results) |
| Quick access (right of the tabs) | Save, Undo, Redo |
| Resource bar (far left) | Analysis tree, Engineering database, Templates |
| Analysis tree | Every project of the design; the active one (bold) expands into **General settings**, **Calculation control**, **Input data**, **Results** (with the loaded run) and **Verification**. See below |
| Graphics window (centre) | The CAD model, or the Results folder, a cut plot, an animation or the verification record when one is selected. A cue line at the bottom says what to do next |
| Floating windows | The feature dialog of the selected input node (its case settings, Apply and Cancel), result windows (goal plot, XY plot, point and volume parameters, report, compare), the Solver monitor, the template list and the Wizard |
| Status bar | Worker state, export progress, the last error or notice |

### Analysis tree

**Input data** lists flow-analysis folders. Computational domain, Fluid subdomains, Boundary
conditions, Goals and Mesh show in every project. Solid materials, Fans, Sources, Radiative
surfaces, Porous media, Perforated plates, Contact resistances, Initial conditions, Rotating
regions, Two-resistor components, Heat pipes, Printed circuit boards, Thermoelectric coolers,
LEDs, Electrical conditions, Batteries, Transferred boundary conditions and Coupling show once
the case has settings for them, or when **Customize tree…** turns them on.

Each folder lists its case entries as rows, and a list entry (velocity boundaries, solids, goals,
porous blocks…) lists one row per item, named from the item's `name`, `id`, `patch` or `type`.
Selecting a row edits only that entry in its feature dialog. Entries are filed by what their
names say; blocks that mix several kinds of settings (`config`, `thermal`, `turbulence`, …) are
filed child by child, so `thermal/solids` goes to Solid materials and `config/dt` to
Calculation control. Input data itself and the JSON editor always edit the whole case.

**Results** lists Cut plots (one row per VTU file), Point parameters, Volume parameters, XY plots
(one row per CSV file), Goal plots, Report, Animation (one row per PVD series) and Compare.
Surface plots, Isosurfaces, Flow trajectories, Particle studies and Surface parameters are not
built; Customize tree… can list them greyed, and the ParaView handoff covers them.

States: a project with structural-check errors shows an error icon; Results shows a warning icon
when the loaded run is older than the study; `•` marks unapplied edits.

| Row | Right-click |
|---|---|
| Project | Activate, Run, Clone project, General settings, Calculation control, Load results, Customize tree…, Delete project |
| General settings, Calculation control | Edit definition |
| Input data | Edit definition, Customize tree… |
| Input folder | Edit definition, Hide folder, Customize tree… |
| Case entry or item | Edit definition, Revert to applied, Copy as JSON |
| Results | Load run (each retained run), Unload results, Summary, Compare, Customize tree… |
| Results folder | Insert…, Hide folder |
| Result file | Show, Export results |

Clone project (`multiphysics.study.clone`) copies the committed case, question and setup into a
new study; retained runs stay with the source. Customize tree choices last for the session.

## Author and run

Choose **Multiphysics** in the workspace menu. Use **Wizard** (name, analysis type, template,
summary) or **New project** (the last template used) on the ribbon.
Write its engineering question in the project dialog, edit case fields or JSON in the Input data
dialogs, and choose Apply.
Each study retains its own unapplied draft, including invalid JSON, across study and workspace switches.
Cancel discards the selected draft. Apply, Remove and Keep result use document Undo and Redo.
Save the `.sc` design to retain applied studies and retained results.
Unapplied drafts remain only in the open document session; replacing or closing that document discards them.

The tree nodes are views of the same study (scripts may still address the stages Prepare,
Physics, Mesh, Couple, Run, Results and Verify; their handles are published on the matching
tree rows and on the Solver monitor button).
Mesh exposes actual case discretization settings. Coupling explains the unavailable general coupling adapter.
The structural check covers bounded inputs; the numerical companion checks its full case schema and prerequisites at execution.
An empty structural error list is not numerical or physical qualification.

The catalogue contains 34 design profiles. The 56 templates cover 45 bounded case types in the Rust port,
for all 34 profiles. The Wizard filters templates by analysis type, and **Use template** in the
template list opens the Wizard on its summary.
Catalogue presence is separate from executable support. Cases use generated or explicit numerical representations.
The current CAD model is not automatically converted into these cases.

## Native numerical companion

Use `cfd-solve` built from **cfd-solvers-rs 0.7.0**, based on upstream revision
`7cceaaaddafaf2515e01a86406e13d6037fdbf9a`.
The companion is a separate process, not a linked dependency or binary shipped by this change.
Versions 0.2.0 and 0.3.0 remain supported for cases that do not need newer capabilities.
Obstacle cases require `grid_obstacles_v1` (0.3.0 or 0.4.0). Fan cases require 0.4.0 with
`fan_curve_v1`. Porous-zone cases require 0.4.0 or later with `porous_blocks_v1`. Conjugate heat cases
require 0.5.0 or later with `conjugate_heat_v1`. Surface bodies, goals and natural
convection require 0.6.0 or later. The domain cases below require 0.7.0 with `domain_cases_v1`.
Other producer versions are refused.

The complete source archive includes the separate companion project. Run
`sh scripts/run-multiphysics.sh` on Linux or `scripts\run-multiphysics.bat` on
Windows to build it and launch SolidCraft with its executable path configured.
Both projects' Rust/toolchain, native linker and GUI prerequisites still apply.

```sh
export SOLIDCRAFT_CFD_SOLVE=/absolute/path/to/cfd-solve
cargo run --profile ci -p solidcraft -- --sample
```

```powershell
$env:SOLIDCRAFT_CFD_SOLVE = 'C:\tools\cfd-solve.exe'
cargo run --profile ci -p solidcraft -- --sample
```

An executable named `cfd-solve` on PATH is also discovered. Without it, authoring and case exports remain usable.
Browser builds prepare cases but cannot launch a local process or export a native directory.
Native CLI/MCP sessions inject the configured process host too.

Run freezes committed inputs, their SHA-256 and the study revision. It never silently applies a draft.
At most two workers run concurrently. Each uses a fresh private directory, bounded output and a five-minute timeout.
Cancel remains pending until the worker stops. Closing or replacing a document retires its work.
The native app waits briefly for cancellation acknowledgement at exit.
The transient history retains at most eight jobs; starting a later job can retire the oldest finished entry.
Keep important results in the design before they leave that session history.

## Post-processing and exports

Results displays actual metrics, manifests, artifacts, CSV plots and tables.
Keep result saves admitted current evidence in the design. Later study edits make that evidence historical.
Late or mismatched completions cannot become current results.
Every retained run keeps its original submitted case and hash independently of later edits.
Execution completion, linear convergence, steady convergence and physical validation remain distinct.

Export a case to a new JSON filename or a retained run to a new directory.
The handoff preserves original VTU/PVD/CSV/JSON/Touchstone artifacts, the producer manifest,
`case.json`, a metadata file and a short ParaView guide. Current-study metadata is labelled separately from the captured case.
Open `.pvd` collections or `.vtu` files in ParaView and save a state there.
The exporter does not fabricate `.pvsm` or a spatial field from a CSV plot.

Producer point/cell associations and values remain unchanged. Flow pressure is kinematic pressure in m²/s².
SIMPLE series use iteration, not physical seconds. Maxwell fields may use separate staggered grids.
Generic potential fields have no inferred physical units. Inspect conventions before comparisons.

### Result tools (FloEFD-style post-processing)

The Results tab reads retained evidence only; nothing here re-solves or edits the study.

| Tool | Where | What it shows |
|---|---|---|
| Result summary | Summary pane | Every numeric metric in the run manifest, with a unit hint where the name implies one |
| Goals | Results → Goals | `goals.csv` history with `goals.json` verdicts: final value, last relative change, per-goal and overall convergence, early stop step. One chart; each goal normalised to its own range, or a single goal in its own units when selected |
| Component temperatures | Results → Component temperatures | `thermal_report.json` rows: power, mean/max temperature, limit, margin and a pass/warning/fail/no-limit status with an icon |
| Compare studies | Results → Compare studies | Metric table of the latest run of every study, side by side (`multiphysics.study.compare`) |
| Colour map and range | Field slice (2D and 3D) | Viridis, Turbo, cool–warm, thermal and grayscale maps; gradient legend with nice ticks and unit; range can be locked so frames, fields and views share one scale |
| Frame slider | Artifact picker | Steps through the files of a run's PVD collection and re-inspects each VTU frame, keeping the field, plane and camera |
| Probe | Point, line and volume tools | Field value at a point: linear interpolation of point data, cell value for cell data; outside the mesh the nearest cell is reported and nothing is extrapolated |
| Min / max | Same | Extreme values and their locations |
| Line plot | Same | Up to 10 000 samples between two points; samples outside the mesh are gaps |
| Volume integral | Same | Integral, mean, min, max and measure (length, area or volume) over all cells |
| Difference | Same | `b − a` between two runs of one study on an identical mesh: largest absolute difference and its location, RMS, relative L2. Different meshes are refused rather than interpolated |

Scripts use `multiphysics.result.query` with `{"id", "run_id", "query": {"kind": ...}}` where
`kind` is `probe`, `min_max`, `line`, `integrate`, `compare`, `goals` or `components`. Every
answer carries the run's case hash, a `stale` flag when the study has changed since the run, and
the inferred unit. The logic lives in `crates/multiphysics/src/results.rs` and is unit tested
against hand-built meshes with known answers.

Still missing compared with FloEFD: isosurfaces, flow trajectories and particle studies, oblique
or arbitrary cut planes, surface plots on CAD faces, CAD-overlaid results in the main viewport,
generated reports and image/animation export. Use the ParaView handoff for those.

## Command interface

Every command uses `Session::execute` and typed input/output schemas. Programmatic calls never open dialogs.

| Command | Input | Result |
|---|---|---|
| `multiphysics.catalog` | `{}` | Catalogue, templates, host availability |
| `multiphysics.study.list` | `{}` | Project study/run summaries and case hashes, without case payloads |
| `multiphysics.study.create` | `template`, `name` | Created study authoring state and run summaries |
| `multiphysics.study.update` | `id`, `expected_revision`, `name`, `question`, `case` | Saved study authoring state and run summaries |
| `multiphysics.study.remove` | `id` | Removed identity |
| `multiphysics.preflight` | `id` | Structural issues, separate host status |
| `multiphysics.surface.inspect` | `mesh` with `vertices` and `triangles` | Bounded indexed-connectivity report; not a CAD or watertightness check |
| `multiphysics.run.start` | `id` | Asynchronous job |
| `multiphysics.run.status` | `{}` | Recent statuses, without large dataset bytes |
| `multiphysics.run.cancel` | `job_id` | Pending or acknowledged cancellation |
| `multiphysics.run.keep` | `job_id` | Study authoring state and retained-run summaries |
| `multiphysics.case.export` | `id`, `path` | New case file |
| `multiphysics.result.export` | `id`, `run_id`, `directory` | Original result handoff |
| `multiphysics.result.inspect` | `id`, `run_id`, `artifact` | Read-only metadata for a retained ASCII VTU or PVD artifact |
| `multiphysics.result.query` | `id`, `run_id`, `query` (`kind`: probe, min_max, line, integrate, compare, goals, components) | Read-only answer with case hash, stale flag and inferred unit |
| `multiphysics.study.clone` | `id`, optional `name` | New study with a copy of the committed case; runs stay with the source |
| `multiphysics.study.compare` | `ids` (empty for all) | Metric table of each study's latest retained run |

Exports respect existing agent path scopes. Scratch computations cannot launch external numerical work.
Projects are bounded to 64 MiB of serialized authoring and evidence; this is deliberately a small-case adapter.
Study command replies include the selected case for editing, but only run identity, stale state,
artifact counts and convergence metadata. Project listings omit numerical cases as well.
Original manifests and artifact text remain in the saved design and native result view;
use case export or result export when their complete bytes are needed.

Result inspection reports VTU point and cell counts, cell types, and field names, component
counts, declared numeric types and point/cell associations. PVD inspection reports producer
coordinates, file names, parts and whether each referenced file is retained with the run.
Coordinates may represent iterations or time, and units are not inferred. Inspection accepts
bounded ASCII unstructured-grid VTU and flat PVD references; binary/appended data and malformed
files are rejected. The command returns metadata; the native Results view renders supported cell slices separately. Inspection does not establish physical validation. Parsing happens
in the command caller; a UI must dispatch it off the window thread for large artifacts.

## Review boundaries

The portable `crates/multiphysics` contract crate sits at layer 1 and uses external serialization and hashing only.
The document owns `multiphysics` plus one pending format note. Engine commands and transient jobs live in their own modules.
IO writes handoffs through the safe writer. The UI module owns workspace presentation and drafts.
`apps/multiphysics-host` owns native processes and stays outside the portable library layers.
The integration leaves the CAD geometry kernel unchanged. The separate companion includes the bounded grid-obstacle adapter described below.

This is the first native implementation of the v1.4 design. The portable contract now has persistent CAD
selection and boundary-assignment records (see below). Remaining work includes:

- engine commands, an editor and a kernel snapshot provider for those records;
- arbitrary CAD volume meshing and the face-to-mesh patch mapping;
- conjugate heat transfer and coupling-graph execution;
- general scientific 3D field rendering and the full browser prototype filter pipeline;
- the unported solver, GPU and MPI methods.
Specialist RF, EMC, LED, detector and fan workflows retain the limits of actual connected case adapters.

Indexed surface inspection counts boundary, non-manifold and inconsistent-winding edges. It rejects
non-finite coordinates and invalid indices. It does not weld coincident coordinates or detect
collinearity or self-intersections. The optional companion does not expose a live residual stream;
the current job state must not be interpreted as convergence monitoring.

## Primitive electronics enclosure

Choose the `flow-electronics-enclosure` template. It models a 200 × 120 × 80 mm
fluid enclosure with a coarse PCB block and two component blocks. The inlet is
prescribed at 0.01 m/s. This is isothermal laminar airflow, without fan curves
or coupled solid/fluid temperatures. Obstacle min/max values are half-open
grid-cell indices. They are not CAD selections or millimetres.

After Run and Keep result, select a VTU artifact in Results and Inspect. The
preview supports XY/XZ/YZ slices of canonical axis-aligned HEX8 cells, scalar
fields and vector magnitudes. Cell values and holes are preserved without
point interpolation. Coordinates and field names retain producer conventions;
generic VTU units are not inferred. Frame coordinates come from retained PVD
collections and may represent physical time or iteration index.

Inspect mesh_report.json for excluded cells, boundary areas and signed fluxes.
The pressure field is kinematic pressure, in m²/s². Process completion and small
flux imbalance do not establish steady state, mesh independence or thermal validity.


## Orbitable cell slices

In native Results, choose **3D slice** after inspecting a supported VTU artifact.
The same field, XY/XZ/YZ plane and layer controls select a single planar cell slice.
The surrounding wireframe shows the numerical domain. Obstacle holes stay empty.
The initial view uses an orthographic camera; it is independent of the Design camera.
Use the configured SolidCraft navigation preset to orbit, pan and zoom.
**Fit result** restores the domain framing. **Vector glyphs** displays up to 128 sampled
original cell vectors when the selected field has three components.
Glyph lengths use a display scale relative to the largest vector magnitude in the dataset.
They are arrows, not integrated streamlines. Scalar fields disable the glyph control.

Color limits use the whole selected dataset, so moving a slice does not change
its color scale. The displayed coordinates and field arrays keep their original values.
The view normalizes coordinates internally for the shared camera's display scale;
this does not transform the case, retained artifacts or exported data.
The 2D view remains the initial choice. View controls do not edit a study or create
an Undo entry, and their camera settings are transient.

This is a single-slice presentation. It does not provide volume rendering,
CAD body overlays, general unstructured cell rendering or temperature coupling.


## Fan curve and exhaust filter

Choose the `flow-electronics-enclosure-fan` template. A shroud of four grid blocks
leaves a 40 × 40 mm opening in the left face. The fan blows inward through it.
A 10 mm porous zone across the right face represents an exhaust filter. The PCB and
component blocks are the same as in the prescribed-flow example.

The inlet flow is not prescribed. The case supplies a pressure-flow curve with an explicit
`static` or `total` convention, and the solver finds the operating point. After the
measured system pressure rise settles, it moves the flow towards the curve's crossing
with the measured resistance line. It never extrapolates the curve.

Results contains `fan_operating_point.json` and `.csv`, which hold the history,
corrections and final mismatch. They also contain `porous_zones.json`, with superficial
speed and dissipated power, and the usual VTU/PVD fields and `mesh_report.json`.
Fan pressures are in pascals. The field pressure is still kinematic, in m²/s².

A fresh 0.4.0 run converged after 48 corrections:

| Quantity | Value |
|---|---|
| Volume flow | 8.83 × 10⁻⁴ m³/s |
| System pressure rise | 0.879 Pa |
| Curve mismatch | below 10⁻¹³ |
| Inflow/outflow imbalance | 2.6 × 10⁻⁸ |

The curve and filter coefficients are illustrative, not product data. The flow is
isothermal and laminar on a coarse grid. Porous zones carry a first-order error at the
abrupt resistance jump; the companion tests show it shrinking under refinement.
None of this is thermal or physical validation.

## Geometry setup records

`solidcraft-multiphysics::setup` holds the records the setup research calls for:

- named face or body selections, by stable CAD identity;
- typed boundary assignments (normal-speed, vector, volume-flow or mass-flow inlets;
  static or total outlets and openings; walls with explicit slip, motion and thermal
  condition; symmetry);
- one-law thermal contacts and oriented orthotropic conduction;
- an optional explicit default wall.

`setup::audit` re-resolves references against an exact geometry snapshot. It reports
coverage, inherited defaults, conflicts, faces that do not bound the fluid, missing
pressure references and double-counted contacts. A deleted or split face is unresolved
and blocks Apply. It is never reattached by proximity.

The records are stored per study and omitted while empty. No command, editor or solver
adapter uses them yet, so they do not change any current case or run.


## Conjugate heat transfer (offline)

Choose `flow-electronics-enclosure-cht`. It is the fan enclosure with a `thermal` block:
- air properties and an inflow temperature in kelvin;
- one material entry per obstacle, with per-axis conductivity (the PCB conducts
  20 W/(m·K) in-plane and 0.4 through its thickness), power in watts and an optional limit;
- convection to ambient on four walls.

After the flow run, one steady energy equation covers every grid cell. Fluid cells conduct
and carry heat with the solver's own conservative face fluxes. Solid cells conduct and
generate heat. Interfaces use series conductance, so temperature and heat flux stay
continuous. The coupling is one-way: constant properties, no buoyancy, no radiation.

Results add the following files:
- `thermal.vtu`, with temperature in K and °C and a region index on all cells, solids included;
- `thermal_report.json`, with the energy balance, advected and conducted heat per patch,
  and the bulk rise;
- `thermal_components.csv`, with min, mean and maximum temperature and a status per component.

The status is pass, warn or fail. Warn means the margin is below 5% of the allowed rise
above the inflow temperature.

A fresh 0.5.0 run gives these results:

| Quantity | Value |
|---|---|
| Total power | 5.5 W |
| Advected with the air | 4.20 W |
| Conducted through the walls | 1.30 W |
| Relative energy imbalance | 4 × 10⁻¹³ |
| Processor | 83.9 °C, warn against an 85 °C limit |
| Power block | 70.1 °C, pass |
| PCB | 57.7 °C, pass |

The grid is coarse (10 mm cells) and every value is illustrative. None of this is
physical validation.

## Offline verification suite

`cfd-solve verify --out DIR` runs embedded checks locally in about a second:
- second-order accuracy of the energy operator on a manufactured solution;
- exact layered-wall conduction;
- the fan operating point against an independent fixed-flow run;
- refinement convergence of the Darcy plug pressure drop;
- the conjugate energy balance;
- bitwise determinism.

It writes `verification.json` and exits non-zero on any failure. No network access or
external data is used.


## Natural convection, goals, CAD bodies and the wizard (0.6.0, offline)

**Natural convection.** Add `buoyancy` (gravity vector, expansion coefficient, reference
temperature) to a `thermal` block. The energy equation then advances with every flow step,
using implicit Euler on all cells. The Boussinesq force `-β(T - T_ref)g` acts on the next
step. Flow frames carry `temperature_k`. Solids without density and specific heat use the
fluid's heat capacity, which is fine for steady answers but not for transient history; the
report says so. The verify suite checks the de Vahl Davis Ra = 10³ cavity:

| Rayleigh number | Grid | Nusselt number | Benchmark | Error |
|---|---|---|---|---|
| 10³ | 16 × 16 | 1.148 | 1.118 | +2.7% |
| 10³ | 20 × 20 | 1.141 | 1.118 | +2.0% |
| 10³ | 30 × 30 | 1.132 | 1.118 | +1.3% |
| 10⁴ | 20 × 20 | 2.318 | 2.243 | +3.3% |
| 10⁴ | 30 × 30 | 2.281 | 2.243 | +1.7% |

The error shrinks as the grid is refined.

**Goals.** `goals` monitors these quantities:
- patch volume flow;
- patch mean pressure;
- patch heat out;
- solid maximum and mean temperature.

`goal_control` sets the evaluation window, the relative tolerance (temperatures relative to
their rise above inflow) and an optional early stop. Results add `goals.json`, `goals.csv`
and `goal_<name>` metrics. Converged goals are not mesh independence.

**CAD bodies.** `bodies` takes closed, outward-wound triangle surfaces in metres, such as
SolidCraft's tessellation of a body. It checks them for watertightness and rasterizes them
by cell-centre winding numbers, with a tie rule so shared edges count once. Bodies become
solids like grid obstacles: no-slip walls, a material entry and power.
`mesh_report.json` records each body's represented volume against its true volume. The
representation is a staircase, not cut cells.

**Wizard.** `solidcraft_multiphysics::wizard::enclosure_case` builds a complete case from
mm, °C and W inputs and reports every automatic choice. The `multiphysics.wizard.enclosure`
command turns it into a saved study in one undo step and returns the report.

**Solver note.** The collocated PISO path became unstable above a viscous diffusion number
of about 0.25 in testing. Runs now report `max_viscous_diffusion_number`, and the wizard
keeps it at or below 0.05.

**Still missing compared with FloEFD:**
- immersed cut cells and local mesh refinement;
- turbulence wall treatment in the coupled solver;
- radiation;
- rotating regions;
- in-CAD interactive setup: the wizard is a command, without a form in the workspace;
- engineering-database breadth;
- parametric studies over CAD dimensions (sweeps change case values, not the CAD model).

## Domain cases and sweeps (0.7.0, offline)

cfd-solvers-rs 0.7.0 adds 27 case types for the electrical, mechanics, particle, kinetic
and coupled modules, each with a template:

| Area | Templates |
|---|---|
| Circuits and devices | SPICE RC filter (op, AC, transient, DC sweep); PN diode IV and CV; SPMe cell and four-cell pack |
| Electromagnetics | Half-wave dipole impedance (moment method); 3-D FDTD cavity; Nédélec FEM cavity port impedance; power-plane IR drop; silicon wire waveguide modes |
| Electrochemistry and power | Galvanic zinc/oxygen couple (mixed potential, polarization); five-bus power flow with PV curve, dispatch and swing |
| Multiphase | Euler–Euler bubble column; stirred compressible VOF; evaporating spray in a cross-flow |
| Structures | Plate cantilever under gravity (FEM stresses); axial bar modes; J2 cyclic plasticity; phase-field fracture; Cosserat rod; DEM settling |
| Particles and kinetics | Two-stream PIC; slab depth dose; Landau damping; Lennard-Jones MD |
| Coupled and other | Pennes bioheat SAR hot spot; busbar Joule heating; Cahn–Hilliard; DG advection; THMC column |

Each case writes CSV, JSON and, where the method has a field, VTU/PVD results, so the
Results view plots and slices them like the flow cases. The companion's
`tests/domain_cases.rs` checks each template against a closed form or a conservation
law, for example:
- the dipole resonates at 0.465 λ with 71.9 Ω;
- the cavity's first resonance is within 3.6 % of TE101 on 144 tetrahedra;
- two-stream growth is within 8 % of cold-beam theory;
- the bar's first mode is within 0.1 % of `c/4L`.

Details and limits are in the companion's `docs/DOMAIN_CASES.md`.

**Sweeps.** `multiphysics.study.sweep` takes a study, a JSON pointer to an existing scalar
in its case (for example `/protocol/0/current_a`) and 1 to 16 values of the same type. It
creates one checked study per value in one undo step, or none if any value is rejected.
Run each study as usual and compare their metrics.
