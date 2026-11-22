# Multiphysics workspace

Multiphysics adds saved numerical studies, actual companion execution and portable results under a separate workspace and command namespace.
The planned Simulation entry is unchanged. Existing Design and Render defaults are unchanged.

The workspace's interface is its own crate, `crates/multiphysics-ui`. The apps call
`solidcraft_multiphysics_ui::register()` at start-up, which installs its hooks (draw, tick, reset,
close, the per-document view) into `ui-egui`'s `multiphysics_workspace` module; that module keeps
the tab, the routing and the state slot. Editing the workspace rebuilds only its crate and relinks
the app, and its tests (`cargo test -p solidcraft-multiphysics-ui`) do not build `ui-egui`'s.

## Layout

The workspace is laid out like a CAD-embedded flow-analysis tool:

| Area | What it holds |
|---|---|
| Ribbon, **Flow analysis** tab | One row: Project (Wizard, New project, Model study, Settings, Domain), Insert (Inlets & walls, Goals), Mesh (Global mesh, Local mesh), Solve (Run, Monitor, Stop), Results (Results, Cut plots, Surface plots, Flow paths), Tools (Check, Library), Help (Validate). Short labels; the full name is the button's accessible name. Each group's drop-down (and the menu bar) also holds the rest: Clone, Component control, Add from components, Batch run, Calculation control, Isosurfaces, Particle studies, Point and Surface parameters, Create openings, Templates |
| Quick access (right of the tabs) | Save, Undo, Redo |
| Resource bar (far left) | Analysis tree, Engineering database, Templates |
| Analysis tree | Every project of the design; the active one (bold) expands into **General settings**, **Calculation control**, **Input data**, **Results** (with the loaded run) and **Verification**. See below |
| Graphics window (centre) | The CAD model, or the Results folder, a cut plot, an animation or the verification record when one is selected. A cue line at the bottom says what to do next |
| Floating windows | The feature dialog of the selected input node (its case settings, Apply and Cancel), result windows (goal plot, XY plot, point and volume parameters, report, compare), the Solver monitor, the template list and the Wizard |
| Status bar | Worker state, export progress, the last error or notice |

### Engineering database and Wizard

The **Engineering database** fills the graphics window. On the left, the database tree lists the
internal reference library, by family (metals, polymers, glass, wood, other) with a count each,
and the design's own materials. On the right, **Items** is a table of the chosen folder (name,
family, number of properties) under the search row, and **Item properties** is the chosen
reference's ruled Property | Value | Unit | Basis table. Hovering a row shows its conditions,
uncertainty, original value and source. Double-clicking an item opens its properties. While a
reference is chosen, a bar under both pages names its copy in the design, asks for the review
and copies it (`material.create`). The design folder lists the design's materials in a Name |
Origin | Density table (built-in first) above the material editor; choosing a row opens it there,
and the table steps aside while a material is being edited.

The **Wizard** shows its seven steps in a list on the left: passed steps are ticked, and a
passed step, or any step once a template is chosen, can be opened by clicking it. The page on
the right shows the step's title and hint, then its settings as a ruled Parameter | Value
table: a band per case group, a heading per entry or list item (named from its `name`, `id`,
`patch` or `type`), and one row per value. Short tuples are one row of fields. The template
step shows every category as a filter chip above a table of templates and their solver
profiles; double-clicking a template moves on. The
units step lists the SI units every case uses. The summary repeats the name, template and
profile, and shows the whole case read only. Back, Next, Finish and Cancel draft sit at the
bottom right.

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

Use `cfd-solve` built from **cfd-solvers-rs 0.14.0**, based on upstream revision
`7cceaaaddafaf2515e01a86406e13d6037fdbf9a`.
The companion is a separate process, not a linked dependency or binary shipped by this change.
Versions 0.2.0 and 0.3.0 remain supported for cases that do not need newer capabilities.
Obstacle cases require `grid_obstacles_v1` (0.3.0 or 0.4.0). Fan cases require 0.4.0 with
`fan_curve_v1`. Porous-zone cases require 0.4.0 or later with `porous_blocks_v1`. Conjugate heat cases
require 0.5.0 or later with `conjugate_heat_v1`. Surface bodies, goals and natural
convection require 0.6.0 or later. The domain cases below require 0.7.0 or later with `domain_cases_v1`. Thermal devices need 0.8.0 with
`thermal_devices_v1`, global, point and equation goals need 0.8.0 or later with `goals_v2`, a local mesh needs
0.9.0 or later with `local_mesh_v1`, transferred boundary conditions need 0.10.0 or later with
`transferred_bc_v1`, rotating regions need 0.11.0 or later with `rotating_regions_v1` and
`piso_neighbour_update_v1`, internal fans need 0.12.0 or later with `internal_fans_v1` and
`balanced_body_force_v1`, conditions on faces of bodies inside the domain need 0.13.0 or later
with `body_face_conditions_v1`, and cut cells need 0.14.0 with `cut_cells_v1`.
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

After a run, **Results** (ribbon or menu) keeps the finished run if it is not kept yet and
opens it on its last flow field, which is read at once (no separate Inspect step). The
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

## Thermal devices and goals v2 (0.8.0, offline)

`collocated_flow3d` conjugate-heat cases can now carry these devices on existing thermal solids:

| Analysis tree / Insert | Case field | Companion model |
|---|---|---|
| Contact resistances, Thermal joints | `thermal.solids[].contact_resistance_m2_k_per_w` | Series resistance on faces shared with a different solid |
| Two-resistor components | `thermal.solids[].junction_to_case_k_per_w` | Junction = case max + R·P; the verdict uses the junction |
| Heat pipes | `thermal.solids[].heat_pipe {axis, effective_resistance_k_per_w}` | Axial conductivity L/(R·A) |
| Thermoelectric coolers | `thermal.thermoelectric_coolers[] {name, cold_solid, hot_solid, heat_pumped_w, electrical_power_w}` | Fixed-duty cooler; COP is reported |

The Goals folder's Insert global, point, surface, volume and equation goal entries open the Goals panel
with that goal type picked. *Add* drafts a goal. Global goals are maximum fluid temperature. Point
goals sit at the domain centre: temperature for thermal cases, pressure otherwise. Volume goals
are the maximum temperature of the first solid. Equation goals refer to the first non-equation goal.
Edit the draft and apply it as usual. The admission check refuses points outside
the domain, equations that refer to other equations, unknown solids, or non-positive resistances.
Any of these features requires producer 0.8.0.

**Batch run** (Solve) starts every project that is ready and not already queued or running. It
skips the rest and reports how many it started.

Still not available, with the reason shown on hover: rotating regions. They need a
rotating-frame model in the solver.

## Working in the graphics window (3-D airflow cases)

For a `collocated_flow3d` case the graphics window draws the case itself: the computational
domain, every obstacle, porous block and surface body, and each domain face labelled with its
role (inlet, outlet, wall, symmetry, fan or default). The Analysis tree node decides what the
mouse does there. Every change goes into the project's draft. *Apply* commits it, and *Show CAD
model* switches back to the CAD window.

| Node | In the graphics window | In the panel |
|---|---|---|
| Computational Domain | An arrow on each face; drag it to move that face out or in by whole cells (a label shows the cells and metres) | Size in cells and metres; − cell / + cell per face |
| Solid Materials, Contact resistances, Thermal joints, Two-resistor, Heat pipes, other input nodes | Click picks a part, Ctrl+click adds or removes, a second click on the same spot takes the part behind | The picked parts. Solid Materials: a conductivity, *Apply to picked parts* (creates or updates each part's solid record). Devices: one picked part fills the solid name |
| Boundary Conditions, Fans | Click picks a domain face (same Ctrl and repeat-click rules) | Inlet (speed along the inward normal), Outlet (pressure), Wall, Symmetry for every picked face. *Lids*: with one face picked, an opening in cells and *Create lids* |

Moving a minimum face (left, bottom or front) shifts every grid part, surface body and point goal
by the same number of cells, so the model stays where it is in the domain. A move that would cut a
part is refused. The fan's own face cannot be given another role; edit the fan instead.

*Create lids* (Tools, or the Boundary panel) closes the picked face with one-cell-thick obstacles
`<face>_lid_1..4` everywhere except the opening. The face's condition then acts only through the
opening, the grid version of a lid with an inlet face in it.

### Studies of the loaded model

*Study of the model* (Project) makes a study whose geometry is the open design itself, not a
copy: the case is built from the design's bodies whenever its settings or conditions change, and
again before a run when the design has changed since (`multiphysics.model.*` commands).
- *Internal*: the fluid is inside the bodies. The domain is the bodies' exact box; domain faces
  that touch fluid are walls unless a condition opens them, and faces the bodies cover entirely
  take no condition. *External*: the fluid is around them, in their box plus a margin.
- The cell size is automatic (at most half the thinnest body) or typed. The build reports how
  much of the bodies' volume the grid holds; walls thinner than a cell are lost.
- Conditions are picked on faces in the 3-D view. A flat face on the domain's boundary opens that
  side of the domain (a duct's end, a vent in an outer wall). Any other face becomes a condition
  on its own body (companion 0.13.0, `body_face_conditions_v1`): the body's grid faces nearest the
  face's triangles. Not yet with heat transfer.
- *Mesh* (the Mesh node):
  - *Level of initial mesh*, 1 to 7: the automatic cell size (level 3) made about 1.4 times
    smaller per level up. Moving it drops a typed cell size. The grid's cell counts are shown,
    and the cell limit says when it made the cells bigger than asked.
  - *Show the grid on the model*: the grid's planes drawn on the domain's three low faces.
  - *Refine where the run changes most* (solution-adaptive): reads the shown or last kept run's
    velocity, flags the tenth of cells with the largest jump to a face neighbour, and splits the
    cells around each cluster of them (up to four boxes, `adapt_1`, …, ×2). Adapting again
    replaces the earlier boxes; run again to use the finer grid. By command
    (`multiphysics.model.adapt`) the field, share, number of boxes and split can be chosen.
  - Smaller cells around a selected body or face, split 2 to 4 times, plus one cell around it.
  - After any refinement the time step is fitted again to the smallest cells (local mesh, 0.9.0).
- *Goals* (the Goals node): volume flow through or mean pressure on a selected face, pressure at
  a point, or an equation of other goals. A face goal uses the face's condition patch, the
  domain side it lies on, or a wall patch made on its body for the goal. Goal plots show them.
- *Find a leak*: with a face inside an enclosure and one outside it selected (External, with a
  margin), the shortest path through the fluid is drawn; none means the inside is sealed.
- Every body is its triangle surface, voxelised at the cell centres (a staircase), unless *Cut
  cells* is on (companion 0.14.0, `cut_cells_v1`): then each connected piece of fluid in a cell
  the surface crosses is a cell of its own (so a wall thinner than a cell still separates the
  fluid on its sides), with the open part of its faces and a wall face, measured on a lattice up to 8 times
  finer per axis (volumes and open areas within half a sub-cell), with cells under 30 % fluid
  merged into a neighbour. Thin walls and slanted or curved faces keep their volume and open
  area. Not yet with heat transfer, porous media, fans, rotating regions or transferred
  conditions. The cut-cell report compares the fluid volume with the exact one.
- *Fan* (General): a fan curve (`flow pressure; …`, m³/s and Pa) on a flat face of the model on
  the domain boundary; it blows in against an outlet there (the companion's fan curve). One fan
  per study, not with cut cells (`multiphysics.model.fan`).
- Global goals (Goals node): the fluid's largest and mean speed and its largest, smallest and
  mean pressure (companion 0.16.0, `goals_v3`), and with heat on its largest and mean fluid
  temperature and largest solid temperature.
- *Heat transfer* (the Solid Materials node, always listed for a study of the model): the bodies
  conduct heat and the temperature is solved after the flow (conjugate heat transfer).
  - The inflow temperature (°C); the fluid is air unless set by command.
  - A material on the selected bodies: conductivity and power spread evenly through the body.
    A body without one conducts like the starter material (1 W/(m·K)), with a note.
  - Heat sources on faces: power on a picked face goes into the body's cells behind its grid
    faces (companion 0.15.0, `face_heat_sources_v1`).
  - Thermal conditions (temperature, heat flux, convection, adiabatic) on a flat face of the
    model on the domain boundary, such as a duct's open end. A face inside the domain is refused
    with the reason: give it a heat source instead.
  - Goals: a body's highest or mean temperature (the selected body) and the temperature at a
    point.
  - Not yet with cut cells, or with inlets and outlets on faces inside the domain (a wall goal on
    such a face works).
- Results on the model (once a run's field file is shown):
  - Cut Plots offers *On the model*: the slice is drawn on its plane in the design's 3-D view,
    with the design's camera. *Vectors* adds arrows of the in-plane velocity, scaled to the
    largest (at most 600; cells are skipped evenly beyond that).
  - Flow Trajectories, Particle Studies and Isosurfaces draw on the model, without a cut plot:
    in the graphics window while the tool is open in the property panel beside it, else in a
    view below the tool. The seeds and *Trace* (or *Compute*) come first; the lines are
    coloured by speed on the panel's colour map and range.
  - The view on the model is framed on the model the first time it is shown for a study, and
    only its part on screen is drawn, so a scrolled view still lines up with the model.
  - Surface Plots colours each wetted face of the model by the fluid cell next to it (looked up
    half a cell to one and a half cells off the face). Faces with no fluid beside them (outer
    faces of an internal study, solid-side faces) are drawn see-through, so wetted faces behind
    them show.
  - A click on the slice or on a coloured face probes the result there: a marker and the value
    on the model, and the point goes into the Point Parameters probe. A click where there is
    no fluid says so.

### Bringing CAD bodies into the case

*Add from components* (Insert, or the From Components tree node) turns bodies of the open model
into case parts. Use the body check list or *Use the model selection*. The bodies are tessellated
through the OBJ export (millimetres, welded closed).
- *Automatic* makes a body that fills its bounding box (within 2 %) a grid block snapped to whole
  cells, and any other body a surface body with its exact triangles. You can also force *Grid
  blocks* or *Surface bodies*.
- The first import fits the domain around the bodies plus a margin. The cell size is either
  typed, or chosen from the case's current cell count but no larger than half the thinnest
  body and no larger than the margin, within the solver's cell limit. A note says when a body
  is still thinner than two cells, or when a surface body holds no cell centre (the solver
  refuses it). The import records which CAD point is the domain origin.
- After fitting, a flow case's time step is divided by a whole factor until the diffusion and
  Courant numbers on the new cells are at most one, with that many more steps for the same
  simulated time where the step limits allow; a note gives both values.
- Grid blocks that fill every cell leave no fluid: the structural check refuses that case.
- Later imports can keep the domain and place bodies by that origin. A body whose case name
  already exists replaces that part's geometry and keeps its material.
- *Replace the case's other parts* starts from an empty domain.
- In a thermal case, new parts get the starter solid material. Materials and goals of removed
  parts are dropped, with a note.

With *Show CAD model* on, the CAD selection works as a picker too:
- Solid Materials and the device nodes offer *Pick the bodies selected on the model*.
- Boundary Conditions offers *Pick the domain faces of the model faces selected*: a model face
  lying on a domain face (within half a cell) picks that domain face.

The CAD origin is kept for the session only, because the solver case has no field for it. After
reopening, fit the domain again, or keep using the case parts.

## Flow result tools

The Results folder's Surface Plots, Isosurfaces, Flow Trajectories, Particle Studies and
Surface Parameters nodes work on a run's VTU field file. They draw in the 3-D field view, which
opens automatically and inspects the file if needed, and list their numbers in the panel.

| Node | Inputs | Shows |
|---|---|---|
| Flow Trajectories | Vector field, a seed line (start, end, count), steps, *Also trace upstream* | Lines coloured by speed. Counts of lines that reached the boundary, stopped in still flow, hit the step limit or started outside. Mean and maximum residence time |
| Particle Studies | As above, plus particle diameter, density, fluid viscosity, gravity (−Y) | Particle paths with Stokes drag: response time ρd²/(18μ), released at the local flow velocity, exact drag update per step |
| Isosurfaces | Field, iso level, optional colour field | The surface (marching tetrahedra), its triangle count and area |
| Surface Plots | Field, whole boundary or one side (X−…Z+) | The boundary faces coloured by the field |
| Surface Parameters | As Surface Plots | Area, area-weighted mean, minimum, maximum and integral. For a vector field, the outward flux ∫v·n dA (the volume flow for velocity) |

Streamlines use fourth-order Runge–Kutta with steps of a quarter of the local cell. Point data
are interpolated linearly in the cells' tetrahedral split; cell data are constant per cell.
The same tools are available to scripts as the `trajectories`, `isosurface` and `surface` kinds
of `multiphysics.result.query`.

## PCB, LEDs, batteries and perforated plates

These device nodes compile their inputs into ordinary case fields, so any producer with conjugate
heat transfer or porous zones runs them. Pick the part, or for plates the domain face, in the
graphics window. Fill in the panel and choose *Apply to the draft*.

| Node | Inputs | Written to the case |
|---|---|---|
| PCB | Board thickness, copper layers as `µm:coverage` | Solid conductivity: in-plane Σkᵢtᵢ/t, through-plane t/Σ(tᵢ/kᵢ). Each copper layer mixes copper (390 W/m·K) and FR-4 (0.3 W/m·K) by its coverage. The board normal is the part's thinnest direction |
| LEDs | Electrical power, optical efficiency, optional junction-to-case resistance | `power_w = P (1 − η)`, plus `junction_to_case_k_per_w` (needs producer 0.8.0) |
| Batteries | Current, internal resistance, entropic coefficient dU/dT, cells | `power_w = n (I²R − I·T·dU/dT)`, with T the inflow temperature. A net heat sink is refused |
| Perforated Plates | Name, free-area ratio φ | A one-cell porous layer on the face, with inertial resistance C = K/h and K = (0.707(1−φ)^0.375 + 1 − φ)²/φ² (Idelchik, thin sharp-edged plate) |

*Show Basic Mesh* (the Mesh folder's menu, or the graphics window's check box) draws the
grid's cell lines on the domain faces, thinned to 64 lines per axis.

## Local mesh (cfd-solvers-rs 0.9.0)

The grid can be refined locally. Open *Global Mesh* (or Insert Local Mesh, or Mesh → Local Mesh
in the ribbon). The graphics window then shows the cell lines. Pick parts there, choose a split
(×2, ×3 or ×4) and a margin in cells, and click *Refine around picked parts*. The result is a
draft, as with every other edit. Repeat for finer levels or other parts; *Remove local mesh* returns
to a uniform grid at the coarsest cell size.

- The grid stays a tensor product of node lines. Refining around a part splits every cell row
  that crosses the part's box along x, y and z. The refinement therefore runs through the
  whole domain on the other two axes, unlike FloEFD's octree cells.
- Existing grid nodes never move, so obstacles, porous blocks and lids keep their exact
  positions; their indices are remapped. Surface bodies and point goals are in metres and are
  unchanged.
- The case stores `mesh.coordinates` (`x`, `y`, `z` node lists) and needs producer 0.9.0 with
  `local_mesh_v1`. Admission checks that each list starts at 0, rises strictly, ends at the
  length and has one more entry than the cells. The cell limit (32 768) still applies.
- A domain-face drag adds or removes cells of the edge cell's width. Adding components from CAD
  snaps parts to the graded cells; fitting the domain to components starts a new uniform grid.

In the solver, the mesh builder, porous zones, surface-body voxelization, the body volume
report, heat-pipe lengths, the diffusion-number metric and `thermal.vtu` all use the actual
cell sizes. A uniform grid gives the same results as before. In a heated channel refined around
a chip, the chip's temperature rise matched the uniform grid to within 0.8 %.

## Transferred boundary conditions: zooming in (cfd-solvers-rs 0.10.0)

Run a project, open *Transferred Boundary Conditions* (Insert, or the tree), pick the parts to
study more closely in the graphics window, and choose a margin in parent cells, a refinement
(×1–×4) and a name. *Create zoom-in project* copies the project and drafts the child case for
the copy. Review it and *Apply*, then run it.

- The child box is the picked parts' box plus the margin, snapped outward to the parent's grid
  nodes. Each child cell is then split by the refinement.
- Parts inside are kept: grid parts are clipped, surface bodies and point goals are shifted.
  Parts outside, or surface bodies crossing the box, are left out with their materials and
  goals, and the panel lists them.
- A child face on the parent's own boundary keeps the parent's condition. A face inside the
  parent takes the parent run's results: where the parent flow enters on average, a
  fixed-value velocity with the velocity sampled at every child face centre. Elsewhere, a fixed
  pressure with the sampled pressure profile and a zero-gradient velocity. Faces lying in a
  parent solid become walls in the velocity profile.
- The parent is the run's own input case, so the child matches the results it is fed from.
  Temperatures are not transferred: the child keeps the parent's inflow temperature.
- The case stores the profiles under `transferred`. The solver applies them face by face in the
  momentum equation, the Rhie–Chow fluxes and the pressure equation.

In the solver's test, a child channel fed this way matched the parent's flow to within 1.7 % of
the inflow speed. The same child with a uniform inlet was off by 41 %.

## Rotating regions (cfd-solvers-rs 0.11.0)

Open *Rotating Regions* (Insert, or the tree) and pick the parts that turn in the graphics
window: a rotor, an impeller or fan blades. Give the region a name, a speed in rpm (negative
turns the other way), the axis (X, Y or Z) and a margin in cells, then choose *Add rotating
region around picked parts*. Review the draft and *Apply*. The graphics window draws each
region as a wire box with its axis and speed, and the panel lists the regions with a
*Remove* button.

- The region is the cell block covering the parts plus the margin. It stays at least one cell
  off the domain boundary, so the outer walls stand still. The axis runs through the block's
  centre.
- The model is the multiple reference frame (MRF) model. In the region the flow is solved
  relative to the turning frame, and the walls of parts inside it turn with it. It gives the
  steady flow a rotor drives. Blade passing is not resolved.
- A speed the time step cannot follow is refused by the solver above 0.2 rad per step. Adding a
  region therefore shortens the time step to 0.1 rad per step and raises the step count, so the
  simulated time stays the same.
- Adding a region turns on the PISO neighbour update (`config.piso_neighbour_update`), which
  keeps the solver stable at larger time steps. Removing the region leaves it on. With
  companion 0.12.0, also turn on balanced body forces (`config.balanced_body_force`, set by
  adding an internal fan, or in the JSON editor). The frame forces then balance exactly: a
  box turning with its frame reproduces solid-body rotation to round-off, against 0.5 %
  without.
- A run writes `rotating_regions.json`. For each region it gives the cell count, Ω and the
  swirl ratio: the mean tangential speed over the frame speed. A solid-body rotation gives 1.

## Internal fans (cfd-solvers-rs 0.12.0)

Open *Fans* and choose *Inside the domain*. Pick the part that stands for the fan, such as a
fan box brought in from the CAD model. Give it a name, the direction it blows (+x … −z) and
its curve as `flow L/s : pressure Pa` pairs, for example `0:40, 10:20, 16:0`. Then choose
*Make picked part an internal fan*. *On a domain face* keeps the external fan on a domain
face, as before.

- The part's cell block becomes open fluid that the fan drives. The part is removed, with its
  solid material and the goals on it. The block must not overlap solids, porous blocks or
  another fan.
- The fan is an actuator block. Each step it measures its flow through its cross-section and
  takes the pressure rise its curve gives at that flow. It applies that rise along its
  direction, so air crossing the block gains it. The curve is not extrapolated: a flow beyond
  it takes the end value, and the report flags it. Swirl and blade effects are not modelled.
- Adding a fan turns on the PISO neighbour update and balanced body forces. Without balanced
  forces, the cells beside the fan block would show a checkerboard of spurious velocity and
  the fan would read the wrong flow.
- A run writes `internal_fans.json` with each fan's operating point (flow, pressure rise,
  whether it is on the curve) and its history. The metrics
  `internal_fan_<name>_volume_flow_m3_per_s` and `_pressure_rise_pa` can be plotted like goals.
- In the solver's test, a fan in a channel against a fixed back-pressure delivered the curve's
  flow at that pressure to five digits.

## Fluid volume

*Tools → Fluid Volume* opens the Computational Domain panel and checks the case. *Check
fluid volume* runs the check again. It is redone on its own whenever the draft changes.

- It counts the fluid and solid cells and their volumes in litres. A cell is solid when a grid
  part covers it, or when a surface body encloses its centre. This is the same rule the
  solver uses, so the counts match its mesh report. Porous blocks and perforated plates stay
  fluid.
- The fluid is split into regions joined through cell faces, largest first. Each region lists
  the inlets, outlets and fan it reaches. A region none of them reaches is *sealed*: air
  trapped inside a closed part or between parts. The graphics window outlines sealed regions
  in the warning colour.
- A sealed region still conducts and stores heat, but no flow can reach it. *Fill sealed
  pockets with solid* turns each one into grid blocks named `pocket_<n>_<m>`, merged into as
  few boxes as it can. Review the draft and *Apply*. If a pocket is really meant to hold
  flowing air, open it to the rest of the fluid instead.
- The panel also warns about parts that overlap an earlier part, and about surface bodies too
  thin to enclose any cell centre, which the solver leaves out.

## Results mesh

The Results folder's *Mesh* node draws the computational grid of a run's field file in the
3-D field view. The cell lines lie on the faces of its bounding box, thinned to 64 per axis. The
panel lists the cell count and, per axis, the number of cells and the smallest and largest width.

## Checking the Multiphysics changes, online or offline

`scripts/multiphysics-check.ps1` (Windows) or `scripts/multiphysics-check.sh` runs format,
`cargo check --all-targets`, clippy with `-D warnings`, and the multiphysics tests of the
multiphysics crate, the engine and the UI, then the companion's tests. It writes every step's
output and exit code to `multiphysics-check.log`, and keeps going after a failure, so one log
shows every problem.

To build on a machine with no network, run `scripts/offline-vendor.ps1` (or `.sh`) once where
crates.io is reachable. It packs every dependency of the workspace and the companion into
`solidcraft-vendor.zip`. On the offline machine, unpack it at the repository root, copy
`vendor-crates/cargo-config.toml` to `.cargo/config.toml`, and run the check with `-Offline`
(or `--offline`). The workspace needs Rust 1.95 or newer (`rust-version`); the repository pins
1.99.1 for CI, and `RUSTUP_TOOLCHAIN=stable` overrides that pin for a local offline build.
