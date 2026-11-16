# Electronics airflow and cooling architecture

Research date: 2026-11-13. Status: proposed architecture.

Implementation status, 2026-11-13:

- **Delivery 1 (obstacle airflow) and 2 (actual screenshot):** implemented with the
  grid-obstacle companion and the 2D/3D cell slices.
- **Delivery 4, partial:** companion 0.4.0 adds an external fan curve (`fan_curve_v1`)
  with an operating point that follows the measured system resistance. It also adds
  grid-aligned Darcy–Forchheimer porous zones (`porous_blocks_v1`). Both are checked
  against independent fixed-flow runs and analytic plug drops.
- **Contracts for delivery 3:** the portable DomainDefinition-side records exist
  (named selections, boundary assignments, contacts, conduction and a coverage audit).
  Exact CAD fluid extraction, meshing and their editor do not.
- **Still planned:** internal momentum-source fans, coupled heat balance and CHT,
  credibility studies and scale.

## Decision

Extend SolidCraft's existing command-driven Multiphysics workspace with separate,
versioned geometry, mesh, run and result contracts. Keep the Rust study/model layer
portable to wasm and keep native processes in the existing host. Reuse the current
viewport, camera, selection and theme; add a scientific rendering layer rather
than building a second editor or encoding CFD values as CAD face colors.

## Current source baseline

- `crates/multiphysics/src/lib.rs`: revisioned studies, case hashes and retained run
  evidence; current adapter is pinned to solver 0.2.0 and an upstream revision.
- `crates/doc/src/config_styles.rs`: certified `ConfigFaceTarget` references exist.
  They are a starting point for selections, not an implemented CFD boundary model.
- `crates/multiphysics/src/mesh_diagnostics.rs`: bounded triangle-surface inspection;
  it does not certify a CAD fluid volume or finite-volume mesh.
- `crates/engine/src/multiphysics_jobs.rs` and `apps/multiphysics-host/`: asynchronous
  jobs and native process execution already provide the orchestration seam.
- `crates/multiphysics/src/vtk.rs`: bounded ASCII VTU/PVD ingestion preserves cell
  and point associations; the workspace draws 2D and orbitable 3D cell slices, and
  `crates/multiphysics/src/results.rs` adds probes, extrema, line samples, integrals,
  run differences, goals and component tables. CAD-overlaid rendering is not done.
- `crates/viewport/src/scene.rs`: existing viewport meshes, lines and clipping.
- Companion solver `src/fvm/poly_mesh.rs`: cells, owner/neighbour faces and patches.
  Construction checks do not yet establish closure, overlap freedom or connectivity.
- Companion solver `src/case/poly.rs`: flow case generation is box/sheared Cartesian;
  solver mesh and VTU box topology are currently constructed separately.

The separate flow, solid heat and thermal transport case types do not establish
conjugate heat transfer. New capabilities require a new truthful producer identity.

## Contracts and ownership

1. **DomainDefinition**: document-owned body/occurrence references, named regions,
   certified face groups, openings, intended fluid seeds, simplification choices,
   units and assembly transforms. CAD edits leave unresolved groups visible and
   block execution rather than silently selecting a nearby face.
2. **GeometrySnapshot**: immutable evaluated geometry revision, SI conversion,
   region/interface identities, source hash and simplification provenance. Exact
   CAD subtraction and classification belong in kernel/solidcore operations.
   Display triangles cannot decide which CAD body or face a boundary means.
3. **MeshArtifact**: coordinates, cell/face connectivity, oriented area vectors,
   cell volumes, region IDs, boundary patch IDs, source-face mapping, mesher identity,
   input hash and quality report. Solver addressing and exported result topology
   derive from the same mesh, with explicit cell ordering or a checked permutation.
4. **PhysicsDefinition**: typed materials, boundaries, sources, initial conditions,
   reference pressures and model choices. Separate preparation from admission to
   a backend; unsupported turbulence, fan or heat-transfer models fail explicitly.
5. **FrozenRun**: geometry/mesh/physics hashes, backend identity, schema version,
   capability set, numerical controls, execution limits and reproducible input files.
6. **ResultDataset**: mesh identity, region, field name, unit, component count,
   cell/point association, time, original arrays, diagnostics and provenance.
   Derived interpolation, glyphs and streamlines identify their source and method.

Use modules within the existing crates initially. Introduce an independent crate
only when a dependency or reuse boundary requires it; do not add a service per stage.

## Mesh strategy

First deliver a discretized primitive enclosure: grid-aligned PCB/component blocks,
real fluid-cell exclusion, named no-slip obstacle patches and named openings.
Require obstacle extents to coincide with grid planes; reject or explicitly preview
snapping. Report represented dimensions and unresolved gaps. Reject sealed or
disconnected fluid components lacking the required pressure reference/flow path.
This is a bounded numerical demonstration, not arbitrary CAD meshing.

Next add exact CAD fluid-domain extraction, leak diagnostics and a replaceable
body-fitted mesher interface. Plan for surface/gap refinement and boundary layers.
Do not implement general cut-cell, adaptive, turbulent CHT and rotating blades as
the first milestone. Mesher dependencies require a separate licence review; public
OpenFOAM/BARAM documentation is behavioral reference, not source for copying.

Before solve, check finite values, index bounds, face orientation, positive volume,
cell closure, patch coverage, duplicate faces, connectivity and supported quality
limits. Geometric intersection/overlap checks need suitable geometric validation;
local topology checks alone cannot prove non-overlap. Thresholds are backend and
model specific, not one universal quality score.

## Electronics models

- Start with prescribed flow rate/velocity and pressure outlet. Label this honestly.
- Add an external fan pressure-flow curve and later internal momentum-source fans.
  Preserve whether the manufacturer curve uses static or total pressure, flow units,
  direction, reference conditions and the allowed interpolation/extrapolation range.
  The operating flow must respond to system resistance.
- Permit vents/filters as explicit geometry or calibrated porous resistance. Store
  coefficient convention, units and local axes; open-area percentage is insufficient.
- Add solid conduction and fluid energy only with explicit interface coupling.
  Preserve matching interface temperature and balanced normal heat flux, or a
  declared contact-resistance law, with convergence and integrated energy checks.
- Support PCB conductivity tensors in board-local coordinates and power sources in
  watts. Do not claim junction temperature without an appropriate package model.
- Add buoyancy, turbulence and radiation under separate verified capabilities.

## Execution and storage

Keep the UI thread free of meshing, solving, parsing large arrays and spatial-index
construction. Use cancellable jobs with bounded progress messages and stale-result
checks against immutable input hashes. Distinguish process completion, numerical
convergence, conservation assessment, refinement evidence and experimental validation.

Keep small portable example evidence. For larger datasets use a content-addressed
artifact directory and compact document references, with explicit portable bundle
export, missing-artifact handling and retention policy. Current limits (8 MiB per
artifact, 16 MiB per run) should not simply be removed. Bound decoded bytes, cells,
timesteps and memory independently of compressed input sizes.

Use VTU/PVD interchange first. Plan controlled binary/compressed and partitioned
ingestion later; do not feed appended binary data through the current ASCII parser.
Treat artifact paths and external references as untrusted. Preserve original arrays.

## Workspace and viewing

Use the existing SolidCraft shell with stages Geometry, Regions, Boundaries, Mesh,
Physics, Run and Results. Each stage exposes prerequisites and actionable failures.
The left tree groups regions, boundaries, mesh and runs; the right inspector edits
the selection; the bottom area contains job logs and physical monitors.

Implement cell-valued slices and a legend first, then probes, vector glyphs,
isosurfaces and streamlines. Keep CAD selection and result-cell selection distinct
but linked through provenance. Display field units, time, region, source run and
stale status. Do not interpolate across solid/fluid boundaries. Do not confuse
steady streamlines with transient particle pathlines. Share camera/clipping and
coordinate transforms with the existing CAD viewport.

## Delivery and acceptance

1. **Obstacle airflow**: generate a single connected fluid-only mesh around explicit
   grid-aligned blocks; verify no fluid cells occupy obstacles, each boundary face
   has one assignment and exported cell order matches the solver. Use an analytical
   channel regression and flux/continuity diagnostics before the enclosure example.
2. **Actual screenshot**: run that case, load the retained field and display a real
   velocity slice with geometry, units, time and provenance. Verify save/reopen and
   stale results after geometry edits. Label prescribed flow and isothermal physics.
3. **CAD enclosure**: certified CAD groups, exact fluid extraction, realistic mesh
   generation, quality visualization and remapping failures after topology edits.
4. **Fan and cooling**: curve operating point, porous models and coupled heat balance;
   compare against suitable published benchmarks and later measurements.
5. **Credibility and scale**: coarse/medium/fine studies with consistent geometry,
   temporal refinement when relevant, partitioned artifacts and remote execution.

Acceptance limits must be chosen and recorded per benchmark. A zero exit code,
small linear residual or attractive screenshot is insufficient evidence of physical
accuracy. A refinement study estimates discretization sensitivity; measurements
are needed for physical validation of a particular real device.

## Public research sources

- BARAM region seeds and fluid/solid assignment:
  https://baramcfd.org/en/manual-en/manual-barammesh-en/region-barammesh-en/
- BARAM boundary-layer controls:
  https://baramcfd.org/en/manual-en/manual-barammesh-en/boundary-layer-barammesh-en/
- BARAM workspace and mesh quality:
  https://baramcfd.org/en/manual-en/manual-barammesh-en/installation-barammesh-en/
- BARAM physical monitors:
  https://baramcfd.org/en/manual-en/manual-baramflow-en/monitor-en/
- OpenFOAM meshing overview (public documentation):
  https://doc.openfoam.com/2312/tools/pre-processing/mesh/generation/snappyhexmesh/
- OpenFOAM CHT solver overview (public documentation):
  https://doc.openfoam.com/2212/tools/processing/solvers/rtm/heat-transfer/chtMultiRegionFoam/
- SimScale fan boundary and internal momentum-source behavior:
  https://www.simscale.com/docs/simulation-setup/boundary-conditions/fan/
  https://www.simscale.com/docs/simulation-setup/advanced-concepts/momentum-sources/
- SimScale porous resistance:
  https://www.simscale.com/docs/simulation-setup/advanced-concepts/porous-media/
- SimScale orthotropic PCB materials:
  https://www.simscale.com/knowledge-base/how-to-set-anisotropic-thermal-properties/
- ParaView dataset associations and hierarchy:
  https://docs.paraview.org/en/latest/UsersGuide/understandingData.html
- VTK interchange specification:
  https://docs.vtk.org/en/v9.7.2/vtk_file_formats/vtkxml_file_format.html
- NASA verification/validation and spatial refinement:
  https://www.grc.nasa.gov/www/wind/valid/tutorial/tutorial.html
  https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html
