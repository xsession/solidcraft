# Multiphysics authoring contracts

This independent, safe Rust crate stores studies and portable evidence for the
Multiphysics workspace. It does not invoke a solver, compute geometry, or change
the existing Simulation workspace.

The original 23 case templates in `data/templates/` come from the MIT-licensed
`cfd-solvers-rs` 0.2.0 source package. Their original license is retained beside
them. They select 18 bounded numerical case types, not all upstream solver
methods. `data/solver-catalog.json` is the original MIT workspace design v1.4
inventory of 34 broad profile families; `data/DESIGN-LICENSE` retains its notice.
Catalogue presence does not establish an installed or qualified execution path.

The additional original `flow-electronics-enclosure` template is a discretized,
isothermal primitive example, with prescribed inlet velocity, a coarse PCB block
and component blocks. Coordinates use SI metres and grid indices are half-open
`min`/`max` bounds. This example is neither CAD fluid extraction nor a thermal or
fan-curve calculation. In particular its 10 mm PCB block is deliberately coarse,
not a representation of a thin physical PCB. Obstacles must be grid-aligned,
non-overlapping, uniquely named and within the declared mesh. Contact between
blocks is permitted. Fluid connectivity and numerical mesh checks remain the
solver's responsibility.

Study names, questions and case objects are authored document state. Incomplete
cases remain savable; `preflight` reports bounded structural checks separately
from the companion's detailed schema, numerical and physical validation.
Updates use an expected revision and commit atomically. No-op Apply preserves
the revision. Removal never reuses a study identifier.

Case limits are 2 MiB, depth 32, 200,000 JSON nodes, 32,768 generated mesh
cells/nodes, 20,000 requested steps/iterations and 50 million mesh-step work
units. These bounds protect the authoring interface; the execution host must
also limit runtime, output, memory and detailed method-specific work.

Runs retain the raw submitted case, its hash, study revision, native manifest and
bounded text artifacts. Each run has at most 128 flat ordinary filenames and
16 MiB of artifact text; no individual artifact exceeds 8 MiB. A project holds
at most 64 studies and 12 runs per study. The complete compact serialized project,
including original cases, normalized manifests and escaped artifact text, is
limited to 64 MiB. Every mutation validates that final budget before committing.
These are deliberate portable example limits, not large-dataset post-processing
support. The host validates fresh directory ownership and symlinks before admission.

Admission understands `cfd-solvers-rs/run/1`, producer versions 0.2.0 to 0.7.0, the pinned upstream source-baseline revision and `native_rust_serial_cpu`.
Each authored extension needs its own capability: nonempty obstacles need
`grid_obstacles_v1` (0.3.0 or 0.4.0); a `fan` needs `fan_curve_v1` (0.4.0);
nonempty `porous_blocks` need `porous_blocks_v1` (0.4.0); a `thermal` block needs
`conjugate_heat_v1` (0.5.0) and a material entry for every obstacle and surface body
(no implied solid). Surface `bodies` (closed, outward-wound triangle surfaces in metres,
checked for watertightness here), `goals` and thermal `buoyancy` need `surface_bodies_v1`,
`goals_v1` and `natural_convection_v1` (0.6.0 or later), and the 27 domain solvers in
`airflow::DOMAIN_SOLVERS` require `domain_cases_v1` (0.7.0). `airflow::requirements`
lists them for a case. Legacy 0.2.0 runs remain admissible for other cases.
Fan and porous admission checks the curve's monotonicity and range, the explicit
pressure convention, single ownership of the fan patch and its fixed-pressure
reference. It also checks porous bounds, coefficients and overlaps with obstacles. The normalized
manifest case must preserve every explicitly authored key/value; additional
solver defaults are permitted. Declared D2Q9 single-precision scalars use the
exact f32-to-f64 encoding of the pinned runner. This does not add a numerical
tolerance to other settings. The host independently binds the original
submitted case bytes to the job. This is structural consistency, not a
cryptographic signature or independent proof that a solver executed.

Editing keeps earlier immutable runs as historical evidence. Attaching a new
run requires the current case hash and revision. Execution completion, linear
convergence, steady convergence and physical validation remain distinct. A
manifest cannot establish physical validity, convergence quality, field units,
or a CAD-to-volume-mesh connection.

## Geometry setup records

`setup` holds the portable records for geometry-based setup. A `Study` keeps them
under `setup`. They are omitted from the saved study while they are empty.

- **Named selections** of faces or bodies. Each member uses a stable CAD identity:
  occurrence path, body, and a persistent face lineage or exact split-piece name.
- **Typed boundary assignments.** An inlet is one of normal speed, a fixed vector in a
  named frame, volume flow, or mass flow with explicit density. An outlet or opening
  states static or total gauge pressure and its backflow policy. A wall states slip,
  motion and thermal condition. Symmetry is also available.
- **Thermal contacts**, each with exactly one law: total K/W, area-specific m²K/W,
  conductance W/(m²K) or a thin layer.
- **Oriented conduction.** Isotropic, or orthotropic with right-handed principal axes;
  `world_tensor` returns R diag(k) Rᵀ.
- **An optional explicit default wall.**

`setup::audit` re-resolves every reference against a `TopologySnapshot` of the exact
current geometry. The document or kernel layer supplies that snapshot; this crate never
evaluates geometry. A missing reference, or a lineage now split into several pieces, is
unresolved. It is never reattached to a nearby face. The report shows:

- face counts and areas, and the volume flow implied by each inlet;
- uncovered fluid boundary faces, and faces using the default wall (reported as
  inherited, not reviewed);
- faces owned by two assignments, as conflicts;
- faces that do not bound the fluid, a missing pressure reference, contacts defined twice
  or overlapping, area mismatches between contact sides, and bodies with two conductivities.

`Project::apply_setup` refuses unresolved or invalid setups atomically. Coverage gaps do
not block Apply; they stay visible until Run.

No numerical adapter consumes these records yet. Engine commands, the workspace editor
and the kernel snapshot provider are the next integration steps. A ready audit is
preparation evidence only. It does not establish mesh validity, solver admission or
physical validity.

## Enclosure wizard

`wizard::enclosure_case` turns a short description into a complete case. The description
gives:
- the enclosure size in mm;
- fan cooling (inlet face, opening, curve in L/s and Pa, exhaust face, optional filter) or a
  sealed box with natural convection;
- the ambient temperature in °C and the external wall coefficient;
- components (box in mm, power in W, per-axis conductivity, optional limit in °C);
- a mesh level from 1 to 5.

It chooses near-cubic cells within the level's budget and snaps components to the grid.
It builds the fan shroud, exhaust filter, walls, materials and buoyancy. The time step keeps
the Courant number at or below 0.5 at the expected peak speed and the viscous diffusion
number at or below 0.05. It adds a maximum-temperature goal per component, with early stop.

The `WizardReport` records each choice: the grid, cell size, time step, and the requested
versus represented box and largest shift for every component, with warnings. A
wizard-generated fan case and natural-convection case were each run through the real
companion and admitted as current evidence. The result is an ordinary case, so it can be
edited, saved and exported like any other.
