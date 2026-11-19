# Multiphysics reference workflow implementation

Base: `50145b4482f531e2b4e99f10deb0c449ac8ef2dc` (SolidCraft 0.44.0).

The user-supplied Screenpresso archive contains 79 FLOEFD-for-Creo captures. The implementation uses original SolidCraft widgets, theme, commands and illustrative data. Vendor artwork, proprietary material tables and private model identifiers are excluded. The second archive mixes product documentation and proposed workflows; its two chat/report illustrations are not evidence of running commercial software. Creo Flow Analysis and FLOEFD-for-Creo are distinct products.

## Functional mapping

| Observed workflow | SolidCraft entry | Current behavior |
|---|---|---|
| Project / Insert / Mesh / Solve / Results / Tools / Help | Seven grouped menus above the existing ribbon | Shared dispatcher and native UI kit; disabled actions explain missing capability |
| New project wizard | Project → wizard | Seven pages: identity, analysis/template, SI units, general/control/domain/mesh, materials/boundaries, initial conditions, summary |
| Preserve setup while moving between wizard pages | Wizard candidate | Edits remain transient; Finish creates a template-compatible study in one undo step; Cancel creates none |
| Clone a study | Project → Clone active project | Copies committed case and setup through the existing engine command; numerical evidence stays with source |
| Create a porous region | Insert → Porous | Prepares one unoccupied cell on a bounded primitive grid; requires review and Apply; surface-body bounds require explicit authoring |
| Fan pressure/flow curve | Insert → Fans | Prepares an external singleton fan with an original illustrative curve; replacing prescribed patch conditions requires an explicit checkbox |
| Component material properties | Insert → Solid materials | Prepares a thermal record for an existing named numerical body; defaults are illustrative; existing internal material catalogue remains available under Tools |
| Goals | Insert → Goals | Prepares a named outlet volume-flow goal; duplicate names are refused |
| Vector and tuple property editing | Structured case fields | Scalar array entries use real number/text/bool widgets, with stable RFC 6901 field targets |
| Material database | Tools → Engineering database | Existing reference catalogue and custom design materials; no copied vendor catalogue |
| Solve / monitor / stop | Solve menu | Existing asynchronous companion job commands; preparing a feature never launches a job |
| Result tools | Results menu | Existing supported plots/queries; unsupported surface, isosurface, trajectory and particle tools are disabled |

The units page explicitly states the case's existing SI convention. It does not implement arbitrary display-unit conversion. Wizard pages edit only fields present in the selected numerical template. Bounded authoring validation is separate from full companion schema admission and physical validation.

## Remaining capability gaps

Arbitrary curved CAD boundary adapters and electrical contacts remain unavailable in this primitive-flow authoring path. Internal fans (actuator blocks on a fan curve), the fluid volume check (fluid and solid cells, connected fluid regions, sealed pockets), lids, local mesh refinement, TECs, heat pipes, two-resistor, PCB and LED devices, perforated plates, transferred boundary conditions (zooming) and rotating regions (MRF) are implemented on 3-D airflow cases; see the workspace guide. A tree folder or menu label does not establish a numerical implementation. Unsupported models are not silently mapped to generic fluid conditions.

Radiation and electrical case fields are editable when a selected adapter actually exposes them; their presence in a generic editor does not add coupled radiation or Joule heating to collocated_flow3d. Solid thermal records describe the numerical case, while design material catalogue entries belong to the CAD document; this change does not invent an automatic coupling between them.

## Code and acceptance

- `crates/multiphysics/src/authoring.rs`: bounded, read-only candidate preparation.
- `crates/engine/src/cmd/multiphysics.rs`: typed preparation command and atomic optional case override for study creation.
- `crates/ui-egui/src/multiphysics_menus.rs`: grouped native menus and capability reasons.
- `crates/ui-egui/src/multiphysics_shell.rs`: existing dispatcher and seven-page wizard.
- `crates/ui-egui/src/multiphysics_workspace.rs`: feature form, explicit Apply and structured scalar-array widgets.
- Existing test files contain `reference_feature_*`, `reference_feature_prepare_*`, `reference_menus_*` and `reference_seven_step_*` acceptance cases.

See the delivered verification report for checks actually completed. A test's existence is not a passing result. Screenshots produced by the UI acceptance helper are headless egui widget captures, not native OS-window recordings.

## Checkpoint (18 November 2026)

The previous scratch checkout was removed between turns. This reconstruction starts from the pinned current base, not from an assumed surviving patch. Source and tests have a separately saved checkpoint.

Last completed checks: 59 core tests passed (one real-VTU fixture test ignored), 15 locale runtime tests passed, 16 workspace GUI tests passed (one real-enclosure fixture test ignored), 14 host tests passed (four companion tests ignored), strict core clippy passed, and layers / mesh-decisions / assets / UI-container gates passed. The native UI/host build completed successfully. The focused engine transaction check passed (1 passed; 0 failed; 0 ignored) against the real compiled engine library. The full engine unit-test binary was killed during compilation twice and did not execute. `xtask lists` reports pre-existing ordering problems in 13 untouched root locale catalogs. Full CI has not passed or been claimed.

The four new workflow captures at 1440×900 and 1280×800 were inspected: menus are visible, structured numeric tuple editing is reflected, and wizard navigation/Finish controls are on screen. The feature sheet floats over part of the ribbon at the smaller height while keeping its footer on screen. These are actual headless egui captures, with no numerical worker installed in that test session; they show setup, not a solved flow field.

The focused engine transaction check is complete; the full engine unit-test build remains blocked by compilation-time SIGKILL. The source ZIP and changed-source checkpoint include the final report; the checkpoint also includes verification logs and captures. The new GUI tests are verified; numerical qualification of the illustrative new feature values is outside this change.

## Reference layout pass (18 November 2026, second round)

Read from the same FLOEFD-for-Creo captures: the fully customised Analysis tree, the ribbon
groups and their dropdowns, the tree right-click menus and the docked property panel.

| Observed | SolidCraft now |
|---|---|
| One Flow Analysis tab: Project, Insert, Mesh, Solve, Results, Tools, Help | One ribbon tab with those groups. Project: Wizard, New, Clone, General settings, Component control, Computational domain. Insert: Boundary condition, Goals. Mesh: Global mesh, Local mesh. Solve: Run, Batch run, Calculation control, Solver monitor, Stop. Results: Load results, Cut plots, Point parameters, XY plots, Goal plots, Report, Export results. Tools: Check, Engineering database, Create lids, Templates. Help: Validation. Batch run, Local mesh and Create lids are disabled with the reason on hover |
| Input Data order: Computational Domain, Component Control, Fluid Subdomains, Solid Materials, Boundary Conditions, Transferred Boundary Conditions, Fans, Heat Sources, Radiative Surfaces, Radiation Sources, Contact Resistances, Thermoelectric Coolers, Heat Sink Simulations, Porous Media, Perforated Plates, Two-Resistor Components, Heat Pipes, Thermal Joints, Printed Circuit Boards, Electrical Conditions, Initial Conditions, Goals, Mesh (Global Mesh), From Components | Same order, then SolidCraft's Rotating regions, LEDs, Batteries and Coupling. Component control, Radiation sources, Heat sink simulations, Thermal joints and From components were added as folders. Mesh always lists a Global mesh row. General settings and Calculation control left the tree (ribbon and project menu, as in the reference) |
| Results (Not loaded): Mesh, Cut Plots, Surface Plots, Isosurfaces, Flow Trajectories, Particle Studies, Point/Surface/Volume Parameters, XY Plots, Goal Plots, Report, Export Results | Same order, then Animation and Compare. "Results (not loaded)" when the project has no run. Mesh is listed greyed (no mesh view yet); Export results opens the report's export |
| Goals folder: Insert Global / Point / Surface / Volume Goals…, Insert Equation Goal… | Same list; global and surface goals open the goal feature, the others are disabled with the reason |
| Mesh folder: Global Mesh, Insert Local Mesh…, Show Basic Mesh, Basic Mesh Color…, Create Mesh…, Delete All Local Meshes | Same list; only Global mesh is active, the rest explain what is missing |
| Item: Edit Definition…, Clone, Delete…, Save As…, Copy to Project…, Properties… | Edit definition, Clone and Delete (list items only), Copy to project (submenu of the other projects), then SolidCraft's Revert to applied and Copy as JSON. Save As and Properties are not listed |
| Feature dialogs docked at the right with ✓ / ✕ | Right-hand resizable property panel. ✓ applies pending case edits and closes; ✕ closes and keeps the draft. The dialog bodies (and their Apply / Cancel handles) are unchanged |

This round was edited without a compiler (the pinned toolchain could not be downloaded in that
session). Formatting, the raw-widget guard and locale-key coverage were checked; the native
build and the GUI suite were not run. Three raw `selectable_label` calls in the wizard were
replaced by the kit's `list_picker`, which the widget guard requires.

## Backend wiring (19 November 2026)

- Insert → Contact resistance, Thermal joint, Two-resistor, Heat pipe and TEC are enabled for
  `collocated_flow3d` cases with a thermal block. They are disabled for other cases, with the
  reason on hover. Each opens its tree node. Its panel drafts the device on a named thermal
  solid, or for a TEC, between the first two solids, through
  `multiphysics.case.feature.prepare`.
- The Goals folder menu enables all five Insert goal entries. The panel has a Goal type picker.
- Batch run is enabled in the ribbon and in the Solve menu.
- Engine `multiphysics.case.feature.prepare` accepts kinds `contact_resistance`,
  `two_resistor`, `heat_pipe`, `thermoelectric_cooler`, `global_goal`, `point_goal`,
  `volume_goal` and `equation_goal`.
- Verified here: the companion (414 tests) and the multiphysics crate (drafting and admission
  tests). The engine and egui changes were checked only with rustfmt, the widget guard and the
  locale-key check; they were not compiled.

## Picking in the graphics window (19 November 2026, second round)

- `multiphysics_case_view.rs` draws a 3-D airflow case in the graphics window. It shows
  domain-face drag arrows for Computational Domain, part picking for the solid-material and
  device nodes, and face picking for Boundary Conditions and Fans. Ctrl+click adds; a repeated
  click on the same spot steps to the next hit, like Select Other.
- `solidcraft_multiphysics::case_scene` holds the geometry and edits: `scene`, `pick`,
  `move_domain_face`, `assign_solid_material`, `set_patch_role` and `create_lids`. It has 7 unit
  tests.
- Create lids (ribbon Tools, Tools menu) is enabled for 3-D airflow cases and opens Boundary
  Conditions.
- Uncompiled here: the egui view, the panel tools and the GUI test
  `case_view_drags_domain_faces_and_sets_picked_faces_and_parts`.

## Add from components (19 November 2026, third round)

- Engine commands `multiphysics.cad.bodies` (model body names) and
  `multiphysics.case.component.import` (tessellate the bodies and return an unapplied case
  candidate plus the CAD origin used). Neither changes the document.
- `case_scene::parse_obj`, `import_components` and `patch_at` have 3 new unit tests. The engine
  test `cad_bodies_become_case_parts_without_touching_the_document` is uncompiled here.
- Insert → Add from components, a ribbon button, and the From Components panel. The Parts and
  Faces panels offer the CAD selection as a picker.

## Flow result tools (19 November 2026, fourth round)

- `solidcraft_multiphysics::flowpost` adds a binned cell `Locator`, `trajectories` (RK4 tracers and
  Stokes-drag particles), `isosurface` (marching tetrahedra) and `surface` (boundary faces,
  area-weighted statistics, outward flux). Its 5 unit tests check a linear field, the circle
  radius under solid-body rotation, the Stokes terminal velocity, plane and sphere areas, and
  the inlet and outlet flux.
- `results::Query` gains `Trajectories`, `Isosurface` and `Surface`.
  `flow_result_queries_trace_cut_and_integrate_retained_fields` needs the real XML reader, so
  it runs in the full build only.
- Surface Plots, Isosurfaces, Flow Trajectories, Particle Studies and Surface Parameters are
  enabled in the tree, the Results menu and the ribbon. The 3-D field view gains an overlay.

## Devices and basic mesh (19 November 2026, fifth round)

- `solidcraft_multiphysics::devices` has `pcb_conductivity`/`apply_pcb`, `led_heat`/`apply_led`,
  `battery_heat`/`apply_battery`, and `perforated_loss`/`apply_perforated_plate`. Its 3 unit
  tests check the hand-computed laminate values, the LED and battery heat, the Idelchik
  coefficient (about 4.4 at φ = 0.5 and 52.6 at φ = 0.2) and that ½·C·u²·h equals K·½u².
- Insert → PCB, LEDs, Batteries and Perforated Plates are enabled for 3-D airflow cases; the
  first three need a thermal block.
- Show Basic Mesh is enabled in the Mesh folder and the Global Mesh menus.

## Local mesh (19 November 2026, sixth round)

- Companion 0.9.0: `MeshCoordinates` and `PolyMeshSpec::coordinates`,
  `make_rectilinear_hexa_mesh`, `VtuMesh::rectilinear_cells`, and the capability
  `local_mesh_v1`. All 418 tests pass, 4 of them new in `tests/local_mesh.rs`: exact agreement
  with the uniform builder, closed graded cells, mass and energy balance on a graded channel and
  agreement with the uniform result, porous and surface-body cell counts, and refusals.
- `case_scene::refine_region`, `refine_around`, `uniform_mesh` and `nodes`. Scene, domain drag,
  component import, face matching and the perforated plate now use node coordinates. Admission
  gains `LOCAL_MESH`. Two new unit tests.
- Local Mesh in the ribbon, the Mesh menu and the Mesh folder's Insert/Delete Local Mesh items is
  enabled for 3-D airflow cases. The Global Mesh panel has the refinement tools.

## Transferred boundary conditions and results mesh (19 November 2026, seventh round)

- Companion 0.10.0: `transferred_bc_v1`, per-face profiles on fixed-value patches. All 420
  tests pass, 2 of them new in `tests/transferred.rs`.
- `solidcraft_multiphysics::zoom::zoom_case` and `results::zoom_from_run` have 2 unit tests,
  covering exact transfer of linear fields, part clipping and the refusals. Admission
  validates `transferred` (requirement `TRANSFERRED`).
- Insert → Transferred Boundary Conditions is enabled for 3-D airflow projects with a run. The
  Results Mesh node is enabled. The only remaining unavailable item is Rotating Regions.

## Rotating regions and PISO neighbour updates (19 November 2026, eighth round)

- Companion 0.11.0: `rotating_regions_v1`, the multiple reference frame model, and
  `piso_neighbour_update_v1`. All companion tests pass, 6 of them new in `tests/rotating.rs`
  and `tests/piso_neighbour_update.rs`. A box turning with its frame settles into solid-body
  rotation to 0.49 % of ΩR (second order in the cell size), with the pressure curvature ½Ω²
  to 0.01 %. A rotor drags the fluid the same amount either way.
- Found while testing: the C++-parity PISO correctors never rebuilt H(u)/a_P. Each step was then
  one unrelaxed SIMPLE iteration, which diverged above a diffusion number ν·Δt/h² of about 0.3.
  The opt-in update makes a 20² cavity stable at 0.4 and a 16² rotating box at 1.3. The default
  keeps parity with the C++ code, so old cases give identical results.
- `case_scene::add_rotating_region` and `remove_rotating_region`, and `Scene::rotating`.
  Admission validates `rotating_regions` (requirement `ROTATING`) and the flag (`PISO_UPDATE`).
  Two new unit tests.
- Insert → Rotating Regions is enabled for 3-D airflow cases.

## Fluid volume check (19 November 2026, ninth round)

- `solidcraft_multiphysics::fluid_volume`: `fluid_volume` classifies the cells as the
  companion meshes them, splits the fluid into face-connected regions and lists the inlets,
  outlets and fan each region reaches. `fill_sealed_pockets` turns the sealed regions into
  grid blocks. Five unit tests, one of them against the exact generalized winding number.
- Found while testing: the companion's surface-body voxelizer counted a cell centre lying on
  an edge shared by two triangles twice or not at all, because each triangle rounded the edge
  test with its endpoints in the other order. On a prism split along a diagonal through cell
  centres, the companion gave 619 cells against an exact 540. It now evaluates every edge in
  a fixed order and gives 540 (`tests/body_voxels.rs`). The fix is in companion 0.11.0.
- Tools → Fluid Volume is enabled for 3-D airflow cases. It opens the Computational Domain
  panel and runs the check. Every item of the reference menus now opens a working tool.

## Internal fans and balanced body forces (19 November 2026, tenth round)

- Companion 0.12.0: `internal_fans_v1` adds internal fans as actuator blocks following a
  pressure-flow curve, with the flow measured from face fluxes. `balanced_body_force_v1`
  adds balanced-force pressure and body-force coupling.
- Found while testing: where a body force jumps (the edge of a fan block), the cell Gauss
  pressure gradient put errors of about d·f/4 into the neighbouring cells. In the test
  channel that was 2.1 m/s of checkerboard against a 0.5 m/s flow, and the fan read a flow 4 %
  off. The balanced form rebuilds ∇p − a from the face values. The channel is then uniform to
  1e-9 and the fan sits on its curve to five digits. It also made the rotating-region
  solid-body test exact to round-off, where it was 0.49 % off.
- A pre-existing issue left as it is: porous blocks (implicit Darcy drag) still checkerboard
  at their edges, because the mobility jumps there. Their pressure drop overall is right,
  but the cell velocities beside a plug are not.
- `case_scene::add_internal_fan` and `remove_internal_fan`, and `Scene::fans`. Admission
  validates `internal_fans` (`INTERNAL_FANS`) and the flag (`BALANCED_FORCES`). One new unit
  test.
- Fans → *Inside the domain* picks a part and makes it an internal fan. The 3-D view draws
  each fan with an arrow.
