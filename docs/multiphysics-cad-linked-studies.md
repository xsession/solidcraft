# Multiphysics studies on the loaded model

Research date: 2026-11-20. Status: phases 1, 2 and 3 implemented (2026-11-20). Partial cells are
measured on a finer lattice, not clipped exactly; the remaining work is listed under phase 3.

Implemented:
- Phase 1: `multiphysics.model.create`, `.set`, `.condition`, `.build`, `.status`
  (`crates/engine/src/cmd/multiphysics_model.rs`, `crates/multiphysics/src/model.rs`). The case is
  built from the design's bodies at each change; Run refuses a case older than the model and the
  workspace's Run builds it again first. The workspace keeps the 3-D view of the model, with the
  domain and the conditions drawn over it (`crates/ui-egui/src/multiphysics_model_ui.rs`).
- Phase 2: conditions on any face of a body (companion 0.13.0, `body_face_conditions_v1`) and leak
  tracking (`multiphysics.model.leak`).
- Phase 3, part: `multiphysics.model.refine` (local mesh around a body or face) and
  `multiphysics.model.goal` (face flow and pressure, point pressure, equations); they have the
  Mesh and Goals nodes of the analysis tree.

## The problem

A Multiphysics study today works on its own numerical case, not on the design that is open:

- A new study copies a template case (`crates/multiphysics/data/templates/*.json`). Its geometry is
  the template's grid blocks, so a study of an enclosure shows the template's enclosure while
  the user's part sits unused in the window.
- The graphics window of the workspace (`crates/ui-egui/src/multiphysics_case_view.rs`) draws the
  case's own box in the content panel. The main 3-D view, its camera, selection and the CAD model
  play no part; *Show CAD model* only overlays tessellated bodies on that drawing.
- The CAD model enters only through *Add from components*: a one-time tessellation of chosen
  bodies into grid blocks or surface bodies (`case_scene::import_components`). After that the
  case is a copy. Editing the part does not change the study, and nothing says the study is out
  of date.
- Boundary conditions go on the six faces of the domain box (`left` … `back`), not on model
  faces. A model face counts only when it lies on a domain face (`case_scene::patch_at`).
- `crates/multiphysics/src/setup.rs` already defines what a CAD-linked setup needs: named
  selections of real faces and bodies by stable identity, typed boundary conditions on them, and
  an audit against the model's current topology. No adapter turns it into a case, and no UI uses it.

## How CAD-embedded flow tools work

Read from public material on FloEFD and on SOLIDWORKS Flow Simulation, which shares its technology
(sources at the end). Only the workflow is taken; no text, icons, templates or data.

1. **The project belongs to the model.** A project is attached to the open model (and its
   configuration). There is no separate geometry: the solid bodies of the model are the solids
   of the analysis, and the fluid is the space they leave.
2. **Internal or external.** The user says whether the fluid is inside the model (an enclosure,
   a duct) or around it. For internal flow the cavity must be closed: openings are covered by
   *lids*, which are real bodies added to the model, and the inlet or outlet is the lid's inner
   face. Cavities without flow conditions can be left out.
3. **Check geometry.** A check reports the fluid and solid volumes and whether the cavity is
   closed. Zero fluid volume means a leak; a leak tracker draws a path from a chosen inside face
   to an outside face through the first gap. It can also create the fluid region as a body to
   inspect.
4. **Computational domain.** A box in the CAD window, drawn around the model. For an external
   analysis it is the model's box plus margins; for an internal one it shrinks to the fluid
   region. Its faces are dragged in the main view or typed. With heat conduction in solids, all
   solids are inside the domain.
5. **Conditions on model faces.** Inlets (velocity, volume or mass flow), outlets (environment
   or static pressure), walls and heat sources are picked on faces of the model in the main view.
   The face must lie between fluid and solid ("face is not laying on the boundary between solid
   and fluid region" is the refusal). One condition per face. The symbol is drawn on the face,
   with arrows for the flow direction, and the condition follows the face when the model changes.
6. **Materials and sources on bodies.** Solid materials, volume heat sources, fans and porous
   media are assigned to bodies or components of the model. Component control turns bodies off
   for the analysis without deleting them.
7. **Mesh from the model.** The mesh is a Cartesian grid over the domain that the solver builds
   from the CAD geometry itself (immersed boundary, with partial cells where the surface cuts a
   cell). The user sets the basic mesh and local refinement on bodies, faces or regions; there is
   no manual meshing. The basic mesh is shown as lines in the CAD window.
8. **Goals and results on the model.** Goals sit on model faces and bodies (average temperature
   of a body, flow through a face). Cut plots, surface plots on chosen model faces and flow
   trajectories are drawn in the CAD window, with the model shaded or transparent around them.

## Gaps

| Flow-tool behaviour | SolidCraft today | Gap |
|---|---|---|
| Project attached to the model | Study copies a template case | No link to the document's bodies |
| Model edits update the study | One-time import, no staleness | Study geometry is a copy |
| Internal / external choice | Not a study setting | Missing |
| Lids as model bodies | `create_lids` blocks cells on a domain face | Not model geometry; only on domain faces |
| Check geometry, leak tracking | `fluid_volume.rs` finds fluid regions and sealed pockets on the grid | Not run on the model; no leak path |
| Domain in the main view | Separate drawing in the content panel | Not in the 3-D view, no model around it |
| Conditions on model faces | Only the six domain faces | Records exist (`setup.rs`), no adapter, no UI |
| Materials and sources on bodies | On case parts after import | Not on model bodies |
| Mesh from the model | Grid blocks snapped to cells, surface bodies voxelised by cell centre | No partial cells; no refinement on faces |
| Results on the model | 2-D slice and a separate 3-D field view | Not in the main view over the model |

## Proposal

### Phase 1: the study uses the model (no solver change)

The solver (`collocated_flow3d`) already takes surface bodies and voxelises them, takes openings
as parts of the six domain faces, and reports fluid regions. Phase 1 uses only that.

1. **Model-based study.** A new kind of study stores a `setup::Setup` (bodies in, analysis type,
   named face selections, conditions, materials) instead of fixed geometry. Templates stay as
   they are for the solver demonstrations.
2. **Case built at run time.** Run (and Check) builds the solver case from the document as it is:
   tessellate the bodies in the study (the existing OBJ path), fit or keep the domain, place
   surface bodies, map conditions. The built case is frozen with the run, as today, so results
   keep their exact input. When the model changes, the study shows that its last run is older
   than the model, through the existing `TopologySnapshot` audit, instead of keeping a copy.
3. **Internal / external.** External: the domain is the bodies' box plus margins. Internal:
   the domain is the fluid region's box; *Check* runs the existing fluid-volume flood fill on
   the built grid and reports fluid volume, sealed pockets and openings without a condition.
4. **Domain and conditions in the main 3-D view.** The Multiphysics workspace keeps the normal
   viewport with the model (one camera path, `navigation::handle`). The domain box is a viewport
   overlay with drag handles; conditions are picked on model faces with the normal selection,
   which already resolves faces by exact geometry. The tree lists the conditions with the faces
   they hold, and the faces show the condition's colour and flow arrows.
5. **Which faces a condition may use in phase 1.** A picked model face that lies on a domain
   face becomes an opening there (the lids path, sized from the face's exact extent). Any other
   face is refused with a reason (the solver has no conditions on faces inside the domain yet),
   rather than being moved to the nearest domain face.
6. **Bodies.** Materials and heat sources go on bodies of the model by identity; *component
   control* excludes bodies. Lids become real bodies of the design, added by a command, so the
   user sees and edits them like any other body.
7. **Results over the model.** Cut planes and the field view draw in the main viewport with the
   model shaded or transparent, reusing the existing slice data.

### Phase 2: conditions on any model face (companion change, done)

- Boundary conditions on the surface of immersed bodies: the solver marks the fluid cells next
  to a tagged triangle set and applies inlet, outlet, wall-temperature or heat-flux conditions
  there. This is what makes a pipe inlet or a vent inside the domain work. It needs a new
  capability name in `cfd-solvers-rs` and a version bump, like the earlier ones.
- Leak tracking: a path through fluid cells from an inside face to an outside face.

### Phase 3: partial cells (done as a sampled version, companion 0.14.0)

Done: `cut_cells_v1` (`companion/cfd-solvers-rs/src/case/cut_cells.rs`), opt-in per model study.
Fractions and open areas come from a lattice up to 8 times finer per axis (within half a
sub-cell), not exact polygon clipping; each cut cell's wall face closes it exactly; cells under
30 % fluid merge. Checked against a sphere's volume, a constriction's open gap and its flow
balance, and a condition on a cut face (`tests/cut_cells.rs`). Left: exact clipping, polyhedral
VTU output (cells show as their grid box now), and cut cells with heat transfer, porous media,
fans, rotating regions and transferred conditions (refused for now). The plan as written:

The companion's flow mesh is already a general polyhedral mesh (`fvm::PolyMesh`: cells, and
faces with owner, neighbour, centre and area vector), so cut cells need no new solver core.
What they need:

1. **Cutting** (`case/cut_cells.rs`, new): for each hex the surface crosses, the fluid part's
   volume and centroid, each grid face's open area and centroid, and the cut polygon's area
   vector and centroid as a new wall (or condition) face. Exact polygon clipping of the body's
   triangles against the hex, with the winding classification the voxeliser already uses for
   cells the surface misses. Degenerate cases (a vertex or edge on a grid plane) need a fixed
   tie rule, as voxelisation has.
2. **Small cells**: a fluid part below about a third of its hex makes the time step and the
   pressure solve stiff. Merge it into its largest fluid neighbour (cell agglomeration) and
   report the merged volume.
3. **Results**: VTU has no hex for a cut cell. Export the cut cells as VTK polyhedra (cell type
   42) or keep the hex with a volume-fraction field; the workspace's slice view reads hexes today.
4. **Conditions**: the cut polygon of a condition's triangles becomes that patch's face, so a
   curved vent or inlet gets its true area instead of a staircase.
5. **Checks**: closure of every cut cell (Σ area vectors = 0), the sum of fluid volumes against
   the exact fluid volume, a sphere and a tilted channel against analytic drag and flow.

A new capability (`cut_cells_v1`) and companion version, opt-in per case at first. Estimate:
the largest single piece of this plan, roughly as much as everything above together.

### Phase 3: accuracy (remaining)

- Partial cells where the surface cuts a cell, instead of whole-cell voxels.
- Local refinement on bodies and faces, using the existing `mesh.coordinates` grading.
- Goals and surface plots on model faces.

### Rules this has to respect

- Which face or body a pick means is decided by exact geometry (`cargo xtask mesh-decisions`);
  the tessellation is only the solver's input, made at run time.
- Programmatic calls never open dialogs; every step is a command with schemas.
- A changed default keeps the old behaviour behind a preference: template studies stay.

## Open questions

1. Phase 1 only allows conditions on faces that lie on the domain boundary. Is that useful
   enough for your models (enclosures with vents on their outer walls), or is phase 2 needed
   before the feature is worth shipping?
2. Should lids be ordinary bodies of the design (visible in Design too), or bodies that exist
   only for the analysis?
3. Do template studies stay in the same tree, or move to a separate "examples" list?

## Sources

- [CFD Online: computational domain and lid tool (FloEFD)](https://www.cfd-online.com/Forums/floefd-floworks-flotherm/172351-computational-domain-lid-tool-print.html)
- [CFD Online: checking geometry (FloEFD)](https://www.cfd-online.com/Forums/floefd-floworks-flotherm/98510-checking-geometry-print.html)
- [Siemens support: how to access and inspect the fluid region body (FloEFD)](https://support.sw.siemens.com/en-US/okba/KB000130372_EN_US/FLOEFD-How-to-access-and-inspect-the-fluid-region-body/index.html)
- [SOLIDWORKS FloXpress help: preparing the model](https://help.solidworks.com/2018/english/solidworks/floxpress/t_preparing_model_floxpress.htm)
- [SOLIDWORKS FloXpress help: flow inlet conditions](https://help.solidworks.com/2018/english/SolidWorks/floxpress/r_setting_flow_inlet_conditions.htm)
- [GoEngineer: lids in Flow Simulation](https://www.goengineer.com/blog/creating-lids-solidworks-flow-simulation)
- [CATI: first surface plot in Flow Simulation](https://www.cati.com/?p=131985)
- [SOLIDWORKS white paper on Cartesian SmartCells meshing](https://www.solidworks.com/sites/default/files/2018-10/DS-18028-Enabling-Fast-Accurate-CFD-Whitepaper-02.pdf)
- [Immersed boundary method (overview)](https://en.wikipedia.org/wiki/Immersed_boundary_method)

## The workflow in three tutorial videos (2026-11-20)

Read from the auto-generated transcripts of three public videos (summarised, not quoted):
a FloEFD training session on conjugate heat transfer in an electronics enclosure (project
setup); FloEFD for Solid Edge on an audio amplifier; and FloEFD meshing tips (an aerofoil and a
three-way valve). What a user does, in order:

1. **Wizard:** name, unit system (they pick °F because their inputs are in °F), internal
   analysis, *heat conduction in solids* on, fluid (air), default solid material, the outer
   walls' heat transfer coefficient and ambient temperature, initial conditions.
2. **Check the domain:** it is found from the model automatically and then hidden.
3. **Materials on bodies:** pick parts in the model and give them a material from the library
   (aluminium heat sink, an anisotropic PCB, a user-defined material kept in a shared library).
   *Select all identical parts* picks every instance at once. Lids get an *insulator*.
4. **Openings:** lids close the openings; the lid's inner face gets a fan (a curve from the fan
   library, optionally with swirl) or an environment-pressure opening.
5. **Heat sources and thermal walls:** a total power on picked components (20 W over four
   transistors), a temperature on the mounting faces or on several faces at once.
6. **Goals:** global (mean static pressure, mean fluid temperature), surface goals on the fan's
   face (the previous selection reused), mass or volume flow through faces, and volume goals
   (max and mean temperature of chips). Goals also decide when the solve has converged.
7. **Mesh:** an automatic level slider with the basic mesh drawn over the model; or manual cell
   counts with control planes that grade the grid; local meshes on a helper body (a transparent
   envelope) to a refinement level, levels splitting a cell into 8; solution-adaptive refinement
   up to a maximum level, which doubles as a mesh-sensitivity study.
8. **Run:** the solver monitor shows goal convergence and a live preview on planes.
9. **Results on the model:** cut plots with velocity vectors and colours, switched to
   temperature; surface plots on the enclosure and on the components; a probe; animated flow
   trajectories, all in the CAD view.
10. **Design study:** change the design (another configuration), clone the project to it, run,
    and compare: the same plot side by side and a table of the goals.

### Against SolidCraft now

| Step | SolidCraft | Gap |
|---|---|---|
| 1 Wizard | Study of the model (internal/external, cell size); SI units only | No unit choice; no heat transfer in model studies |
| 2 Domain | Built from the model and drawn over it | — |
| 3 Materials | Only in template/thermal cases, not on model bodies | Materials on bodies, identical parts, insulator |
| 4 Openings | Inlet/outlet on any face; no lid helper | Fan (curve, swirl) on a model face; create lids |
| 5 Heat | Not in model studies (refused with face conditions and cut cells) | Heat sources and thermal walls on bodies and faces: needs CHT with face conditions and cut cells in the companion |
| 6 Goals | Face flow, face pressure, point pressure, equations | Global goals, volume (body) goals, mass flow, goals that stop the run |
| 7 Mesh | Cell size; local refinement ×2–4 around a body or face | Level slider, basic mesh drawn on the model, control planes, solution-adaptive refinement, 2-D |
| 8 Run | Solver monitor (status); goal plots after | Live goal convergence and preview |
| 9 Results | Cut plots in a panel; "On the model" cut plot being added | Vectors, surface plots on model faces, probe, trajectories in the 3-D view |
| 10 Design study | Clone, compare metrics | Study per configuration; side-by-side plots |

### Order proposed

1. Finish results on the model: cut plots (in progress), then surface plots on model faces and
   a probe. This is what every video ends with.
2. Heat transfer in model studies: materials and heat sources on bodies, temperature and heat
   transfer coefficient on faces, volume goals on bodies. Needs the companion's conjugate heat
   to accept body face conditions and cut cells.
3. Mesh: a level slider with the grid drawn on the model; solution-adaptive refinement.
4. Fans on model faces, lids, global goals, configurations.
