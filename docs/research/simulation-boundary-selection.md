# Boundary selection for complex simulation geometry

Research date: 2026-11-13.

This review uses public official documentation. It describes documented workflows,
then proposes an original SolidCraft interaction model. It does not establish that
SolidCraft currently implements the proposed selection, fluid extraction or solver
coupling features. No commercial application was operated, and no product source,
screenshots, help text, templates or material catalogs were copied.

## Documented workflows

The cited documentation spans several releases. Ansys references identify their
release in the URL; COMSOL references identify 6.3 or 6.4. The Siemens articles are
older official workflow descriptions, not verification of the newest release's
precise screen layout. SimScale documentation is maintained online.

| Environment | Documented workflow | Implication for a new interface |
| --- | --- | --- |
| Ansys Fluent | Geometry and mesh tools filter selection by entity type. Hide/show, transparency and clipping expose interior regions. Boundary conditions operate on boundary zones. [S1, S2, S18] | Keep viewport selection, CAD surface identity and generated solver-zone identity distinct. |
| Ansys CFX-Pre | A boundary has a type and a Location selected from regions or picked in the viewer. Multiple regions can share a condition. Selection modes include single, box, polygon and flood selection. Overlapping picking initially selects the nearest displayed object. [S3, S4] | Offer a persistent location list and an overlap candidate chooser alongside direct mouse selection. |
| Simcenter STAR-CCM+ | Geometry parts are assigned to regions. Their surfaces can become separate boundaries. Surface/volume filters and tags help organize complicated assemblies; hiding classified surfaces supports incremental setup. [S7, S8] | Preserve the mapping from component to named surface group to region boundary. |
| COMSOL Multiphysics | Entity-level filters, graphics context menus, named selections and visibility controls support physics assignment. Reusable selections can serve several model nodes. A selection list helps resolve overlapping entities. [S9, S10] | Treat reusable selections as model objects, separate from temporary highlighting and visibility. |
| SimScale | Create a boundary condition, then click the permitted entities. Saved selections group reusable faces or volumes. Fluid-region extraction is a separate geometry preparation step. [S12, S13, S14] | Separate domain preparation from assignment, and show whether a control expects a face or a volume. |

## Curved and hidden surfaces

COMSOL documents adjacency-based selection with optional tangent-continuity
grouping and an angular tolerance. This can extend a seed across a continuous
sheet. Its geometry adjacency documentation also warns that entities in different
assembly parts are not automatically adjacent. [S11, S16]

Fluent's display tools include clipping, cut-edge display, filters, transparency
and local display aids. These help expose interiors without changing the mesh.
Its meshing interface also supports hiding selected surfaces. [S1, S2]

For SolidCraft, a curved wall should be selectable as a CAD face or an explicit
face set. Render triangles should only help locate that topology; they should
not become the durable identity of the assigned wall. A cylindrical or filleted
wall must not be reduced to one of six box sides.

Recommended selection tools:

- **Face / component filter:** prevent accidental component selection while
  assigning a wall, and accidental face selection while assigning a bulk material.
- **Add / remove:** maintain a visible list of selected references with individual
  removal and a clear-selection action.
- **Hide / isolate / restore:** expose an enclosure's inner walls without deleting
  the enclosure or changing its material and boundary assignments.
- **Section view:** inspect an interior with an adjustable viewing plane. A display
  section must not create an inlet, wall or physical cap.
- **Pick through:** list overlapping selectable entities with owner, face label,
  depth and a temporary highlight. Let the engineer confirm the intended one.
- **Grow selection:** optionally include adjacent or tangent-connected faces with
  an explicit tolerance and preview. Stop at component or region boundaries unless
  the user requests otherwise.
- **Review selection:** show member count, total area, owner components, bounding
  extent and unresolved references. Support zoom to the selected set.

These are design recommendations, not claims about current SolidCraft behavior.

## Internal and external fluid domains

SimScale's flow-volume documentation distinguishes the space occupied by fluid
from the surrounding solid CAD. Internal extraction uses a face touching the
desired cavity as a seed. Openings can require boundary faces; closed cavities can
be extracted without such openings. External extraction creates a surrounding
fluid volume with specified extents. Domain size must account for the relevant
flow development while controlling computational cost. [S14]

Its CAD preparation documentation separately describes internal caps and their
analysis-type restriction. A feature named "cap" therefore cannot be assumed to
work for every solver or domain representation. [S17]

Recommended SolidCraft preparation flow:

1. Choose **Internal cavity** or **External environment**.
2. Identify participating solid components and the intended fluid region.
3. For an internal cavity, select a seed and identify its openings. Preview any
   proposed physical caps before accepting them.
4. For an external case, size the surrounding domain and preview the remaining
   fluid after the component volumes are excluded.
5. Inspect disconnected regions, leaks, tiny gaps and excluded components.
6. Preserve solids when a supported conjugate heat-transfer workflow requires
   their conduction domains. A fluid-only workflow may need a different handoff.
7. Create named inlet, outlet/opening, wall, symmetry and interface groups on the
   resulting simulation topology.

A connected preview or a closed visual shell is insufficient evidence of valid
fluid extraction. Exact topology checks, meshing and solver-specific validation
remain separate requirements.

## Persistent selection and boundary identity

COMSOL named selections can be reused by material and physics nodes. Editing a
selection updates the nodes that reference it. Named selections may also control
visibility. [S10]

SimScale saved selections group multiple entities for later setup assignments.
Their geometric entity dimension matters: volumes, faces, edges and nodes serve
different purposes. [S12]

Recommended SolidCraft data relationships:

- A **selection set** stores topology references and an optional selection rule.
- A **boundary assignment** references a selection set and contains physical data.
- A **mesh mapping** records which generated patches implement each selection set.
- A **view state** stores hiding, transparency, clipping and hover feedback.

Visibility must not remove a physical assignment. CAD regeneration must re-resolve
references against the current exact topology. Missing or ambiguous references
must produce an unresolved state requiring review. Reusing a face number, nearest
triangle or screen location after regeneration is not a sufficient repair policy.

One logical group may produce several solver patches. Ansys documents a case
where faces grouped across different bodies fragment into separate Fluent face
zones. The application should preserve the logical assignment and display its
one-to-many solver mapping. It should not force repeated manual condition entry
or merge distinct material domains merely to simplify the interface. [S15]

## Normals, velocity and frames

CFX-Pre separates the coordinate frame from the boundary's stationary/rotating
interpretation. Boundary values can be evaluated relative to the rotating domain
or an absolute stationary frame. [S3]

Recommended SolidCraft controls must distinguish:

| Quantity | Required meaning |
| --- | --- |
| Normal inlet speed | A scalar interpreted against the boundary's local oriented normal. |
| Inlet or freestream velocity | A vector with explicit coordinate frame and units. |
| Initial fluid velocity | An initial condition, separate from the maintained boundary condition. |
| Moving-wall velocity | Motion of a selected wall, with its reference frame. |
| Translational frame velocity | Motion of the chosen simulation reference frame. |
| Rotational frame velocity | Angular velocity plus axis, origin and frame definition. |

For a curved face, the normal varies across the surface. A constant normal speed
and a constant Cartesian velocity vector therefore describe different fields.
Display oriented arrows at several sample locations and make any sign convention
visible. Preserve the distinction between the CAD face's orientation and the
simulation region's outward normal, especially on two-sided interfaces.

These control requirements are recommendations. Their availability must depend
on the capabilities of the selected solver adapter.

## Practical mouse workflow

The following sequence is proposed for SolidCraft:

1. Open the prepared fluid region and choose **Add boundary**.
2. Choose its physical role. Activate a face-only selection filter.
3. Hover to preview an exact face and its owner. Click to add it; use an explicit
   remove gesture or list action to remove it.
4. For an enclosed wall, hide or section the obstructing geometry. Use the overlap
   list if several faces remain under the pointer.
5. Optionally grow through adjacent tangent faces. Review the highlighted extent
   before accepting the expanded set.
6. Name the group. Review face count, area, owners, interface side and normals.
7. Enter physical values with units and frame. Show conditional parameters only
   when the chosen model supports them.
8. Apply one undoable assignment. Escape cancels the draft and releases pointer
   capture without changing the study.
9. Re-select the assignment in the tree to highlight all its faces. Edit, undo,
   save and reopen must preserve the same model references and physical values.
10. Run a coverage review before meshing or solving: show uncovered boundaries,
    explicit assignments, inherited defaults, conflicts and unresolved references.

The engineering acceptance test must include a curved internal wall and a hidden
enclosure face. A flat external box alone cannot verify this interaction.

## Physical and interaction pitfalls

CFX distinguishes an **Opening** from a predominantly inflowing inlet or outflowing
outlet. An opening permits both inflow and outflow. Choosing the wrong type can
alter reverse-flow behavior. [S6]

CFX Quick Setup supplies a no-slip adiabatic wall for otherwise unassigned domain
boundaries. This is a documented default, not evidence that every such wall is a
correct physical assumption. [S5]

For a new interface, show inherited defaults explicitly. A coverage audit should
not treat "assigned by fallback" as equivalent to "reviewed by the engineer".
A fluid-solid contact also needs its own interface and thermal semantics; an
exterior insulating-wall assumption must not silently replace that interface.

Keep selection ownership clear. Navigation gestures must not change the boundary
set. Clicking ribbon controls must not pick the geometry underneath. Invalid
inputs, unresolved topology and unsupported solver features must reject Apply
without partially changing the study.

## Official sources

- **S1 — Ansys Fluent 2024 R2, Examining the Mesh:** selection filters, clipping and
  display aids. <https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_ug/tgd_user_display.html>
- **S2 — Ansys Fluent 2025 R2, Meshing-Based Graphical Interface:** surface visibility
  and clipping tools. <https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/flu_ug/ch65s01.html>
- **S3 — Ansys CFX-Pre 2025 R1, Boundary Basic Settings:** location, coordinate
  frames and stationary/rotating interpretation.
  <https://ansyshelp.ansys.com/public/Views/Secured/corp/v251/en/cfx_pre/i1310482.html>
- **S4 — Ansys CFX-Pre 2026 R1, Picking Mode:** viewer selection and overlaps.
  <https://ansyshelp.ansys.com/public/views/secured/corp/v261/en/cfx_pre/i1302734.html>
- **S5 — Ansys CFX-Pre 2026 R1, Boundary Definition:** creation and default wall.
  <https://ansyshelp.ansys.com/public/Views/Secured/corp/v261/en/cfx_pre/i1300386.html>
- **S6 — Ansys CFX 2024 R2, Fluid Boundaries:** inlet, outlet and opening semantics.
  <https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/cfx_mod/CHDIJEGH.html>
- **S7 — Siemens, A new user's guide to STAR-CCM+ simulation, Part 4/5:** parts,
  regions and boundary setup; older official workflow article.
  <https://community.sw.siemens.com/s/article/A-new-user-s-guide-to-STAR-CCM-simulation-Part-4-5-Physics-setup-and-Solvers>
- **S8 — Siemens, Fully setup CFD in 1 hour:** selection filters, tags and hiding;
  older official workflow article.
  <https://community.sw.siemens.com/s/article/Fully-setup-CFD-in-1-hour-STAR-CCM-sim-template-with-surface-wrapper>
- **S9 — COMSOL, How to Select Geometry:** selection and hide/show workflows.
  <https://www.comsol.com/support/learning-center/article/how-to-select-geometry-34341>
- **S10 — COMSOL 6.4, Named Selections:** reuse and overlap disambiguation.
  <https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_ref_visualizationselection.22.14.html>
- **S11 — COMSOL 6.3, Adjacent Selection (Mesh Sequences):** adjacency and tangent
  grouping. <https://doc.comsol.com/6.3/doc/com.comsol.help.comsol/comsol_ref_visualizationselection.22.36.html>
- **S12 — SimScale, Saved Selections:** reusable geometric sets; page updated
  2025-05-15. <https://www.simscale.com/docs/topological-entity-sets/>
- **S13 — SimScale, Boundary Conditions:** creation, assignments and entity types;
  page updated 2026-02-12. <https://www.simscale.com/docs/simulation-setup/boundary-conditions/>
- **S14 — SimScale, How to Create Flow Volume Extraction:** internal and external
  fluid-domain preparation. <https://www.simscale.com/knowledge-base/flow-volume-extraction/>
- **S15 — Ansys Knowledge, How to avoid segmented faces zones in ANSYS Fluent?:**
  multi-body named-selection segmentation.
  <https://innovationspace.ansys.com/knowledge/forums/topic/how-to-avoid-segmented-faces-zones-in-ansys-fluent/>
- **S16 — COMSOL 6.4, Adjacent Selection (Geometry Sequences):** assembly adjacency
  and selection availability.
  <https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_ref_visualizationselection.22.26.html>
- **S17 — SimScale, CAD Edit:** fluid volumes and analysis-specific internal caps.
  <https://www.simscale.com/docs/cad-preparation/cad-editing/>
- **S18 — Ansys Fluent 2024 R2, Cell Zone and Boundary Conditions Overview:**
  zone-based assignment and internal versus external boundary topology.
  <https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/flu_ug/flu_ug_bcs_sec_overview.html>

OnScale is not included in this comparison because this review did not establish
an official, relevant CFD boundary-selection workflow for it. Absence from this
table is not a claim that the product lacks a capability.
