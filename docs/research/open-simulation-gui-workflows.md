# Boundary and material workflows in open simulation tools

Research date: 2026-11-13. These findings describe documented user behavior,
not an audit of every release. Only public user documentation was consulted.
No simulator source, material catalog, screenshots, templates or help text was
copied. The SolidCraft interaction sequences below are original proposals;
they do not claim that CAD boundary assignment or thermal adapters are implemented.

## Documented workflows

| Environment | User interaction described in its documentation | Consequence for SolidCraft design |
| --- | --- | --- |
| ElmerGUI | Create a condition, then Alt/AltGr-double-click a surface or edge to assign it. A toolbar mode provides an alternative. Shift-double-click opens body properties for equations, material, body force and initialization. Boundary names in the object browser also open assignment. [E1] | Make body and boundary picking distinct. The viewport and browser must resolve to the same authoritative assignment. |
| ElmerGUI, interior boundaries | Select an exterior face and use View → Hide/show selected; the tutorial then selects the revealed internal boundary. [E2] | Hidden internal geometry needs reversible visibility tools and a browser fallback. |
| BARAM boundary editor | A region boundary list highlights the chosen patch. A context action chooses the boundary type; double-click or Edit opens parameters. A name filter reduces large lists. [B1] | Provide searchable patches, two-way highlighting and a clear type-versus-values editor. |
| BARAM geometry preparation | STL import can split surfaces by feature angle and merge small regions using an area threshold. Named STL SOLIDs stay separate. A Hex6 domain exposes six external faces. [B2] | Retain imported patch names; preview grouping before acceptance. Do not confuse an external domain with a component surface. |
| BARAM visibility | Display controls include hiding, opacity and cut settings. Clip/Slice can expose an interior; individual objects can be exempted from cutting. [B3] | Section views and per-component visibility belong beside selection tools. |
| SALOME geometry groups | A group selects a main object and subshape type. Multiple faces can be added or removed, isolated or hidden. Groups support union, intersection and subtraction. [S1] | Store named selections with explicit membership and visible editing feedback. |
| SALOME mesh clipping | A clipping plane has position and orientation, can be inverted and previewed, and can be manipulated with the mouse. [S2] | Keep display clipping independent of geometry modification and boundary identity. |
| code_saturne | Define boundary zones before their physical conditions. Existing mesh groups can be imported from preprocessing. Semantic regions can combine groups and geometric criteria. [C1, C2] | Treat a patch or volume zone as a durable setup object, rather than a temporary screen selection. |
| OpenFOAM | Boundary conditions refer to mesh patches; patch groups can collect several patches. Initial internal values and boundary values are distinct. [O1, O2] | Validate the generated condition set for every active field and every resolved patch. |

A curved surface does not require a planar-boundary approximation in this
selection model. One semantic patch can contain many mesh faces with different
local normals. Conversely, coloring a few display triangles does not create
a valid solver patch. The CAD-to-mesh adapter must preserve or explicitly
resolve the named selection.

## Proposed mouse workflows for SolidCraft

The following sequences are design proposals, not implemented capabilities.
Their names and interaction states are original.

### Curved inlet or wall

1. Activate boundary selection and hover a CAD face or imported boundary patch.
2. Preview its name, owning occurrence, adjacent region, area and local normal arrows.
3. Click to add it. Use additive selection to include further faces, or remove an unwanted face.
4. Optionally preview connected or tangent-continuous growth. Show the resulting membership before acceptance.
5. Save the selection with a meaningful name, choose its physical role, edit values and Apply.
6. Validate the resolved solver boundary and retain the mapping with the study.

Normal speed on a curved inlet must use each boundary location's normal.
A specified global velocity vector is a different input. Display the direction
and sign convention instead of silently switching between the two.

### Wall inside an electronics enclosure

1. Select the cover and hide it, isolate the internal component, or enable section clipping.
2. Activate boundary selection. Hover and inspect the owning component and region before clicking.
3. If several surfaces overlap under the cursor, show a pick-through candidate list with preview highlighting.
4. Choose the actual internal face, then assign its boundary or interface condition.
5. Restore visibility and review the assignment in the study tree.

A section-plane cap is a display artifact, not a new physical boundary.
Internal fluid-solid interfaces need paired region ownership and an interface
model. They must not silently become exterior pressure outlets.

### Internal and external domains

For an internal study, show the retained fluid region, opening caps and any
leak/connectivity diagnostics before meshing. For an external study, show the
far-field enclosure independently from the solid objects. A seed or explicit
volume selection must identify the intended fluid region. Region choices and
boundary mappings need revalidation after geometry or mesh changes.

## Materials, thermal conditions and motion

BARAM adds materials to a case, offers property editing and supports user
database import/export as CSV. Available properties depend on phase and active
physics. This is useful behavioral evidence for a catalog and project-material
workflow; it is not permission to copy that database or its numbers. [B4]

BARAM separates wall motion and slip from thermal conditions. Its thermal
choices include prescribed temperature, heat flux, adiabatic behavior and
exterior convection/radiation. Coupled walls represent paired interfaces;
thin-wall resistance differs from resolving a solid region. [B5]

Its initialization editor separately sets internal velocity and temperature,
and can use a specified boundary as a starting-value source. [B6]

code_saturne separates zones, physical models, body forces, properties,
initialization and boundary conditions in its setup sequence. Enabling thermal
and radiation models exposes their additional settings. Its radiation
documentation distinguishes medium absorption from wall conditions. [C3, C4]

For SolidCraft, use a read-only starter catalog, editable project-local material
snapshots and explicit component/region assignments. Store units, source,
reference state, valid range and revision with the snapshot. Surface finish
and radiative properties should be assignable independently of bulk density,
heat capacity and conductivity. Show whether each property is actually consumed
by the selected solver adapter.

Gravity, initial fluid velocity, inlet velocity, moving-wall velocity and
reference-frame motion need separate typed settings. A single ambiguous
system-velocity field cannot represent these different physical meanings.
Thermal and radiation authoring remains preparation-only until an executable
adapter consumes the configuration and produces verified numerical evidence.

## Source limitations

The retrieved BARAM material documentation contains a suspicious Sutherland-law
temperature exponent. Do not use its printed equation or material numbers as
an authoritative property source. The library needs independently checked
primary data and explicit provenance.

Direct opening of the ElmerGUI manual failed in the research browser, but search
retrieval exposed the cited user-manual sections. The interior-boundary workflow
was corroborated by CSC's separate GUI tutorial. Software documentation describes
intended behavior; it does not establish that every version implements it identically.

## Sources

- **E1:** [CSC ElmerGUI Manual, Model definitions](https://bio.nic.funet.fi/index/elmer/doc/ElmerGUIManual.pdf).
- **E2:** [CSC Elmer GUI Tutorials, interior boundaries example, page 58](https://www.ftp.funet.fi/index/elmer/doc/GetStartedElmer.pdf).
- **B1:** [BARAM: Types of Boundary Conditions](https://baramcfd.org/en/manual-en/manual-baramflow-en/boundary-conditions-en/types-of-boundary-conditions-en/).
- **B2:** [BARAM: Geometry](https://baramcfd.org/en/manual-en/manual-barammesh-en/geometry-barammesh-en/).
- **B3:** [BARAM: Progress, display and mouse controls](https://baramcfd.org/en/manual-en/manual-barammesh-en/installation-barammesh-en/).
- **B4:** [BARAM: Materials](https://baramcfd.org/en/manual-en/manual-baramflow-en/materials-en/).
- **B5:** [BARAM: Wall Boundary Conditions](https://baramcfd.org/en/manual-en/manual-baramflow-en/boundary-conditions-en/wall-boundary-conditions-en/).
- **B6:** [BARAM: Initialization](https://baramcfd.org/en/manual-en/manual-baramflow-en/initialization-en/).
- **S1:** [SALOME 9.16: Working with groups](https://docs.salome-platform.org/9/gui/GEOM/work_with_groups_page.html).
- **S2:** [SALOME Mesh 9.16: Clipping](https://docs.salome-platform.org/9/gui/SMESH/clipping.html).
- **C1:** [code_saturne 9.0: Manage boundary conditions](https://code-saturne.org/doc/code_saturne-9.0/base_setup_boundary_conditions.html).
- **C2:** [code_saturne 9.1: Mesh properties and selection](https://code-saturne.org/doc/code_saturne-9.1/cs_ug_mesh_select_c.html).
- **C3:** [code_saturne 9.0: Setup of the main parameters](https://code-saturne.org/doc/code_saturne-9.0/base_setup_param.html).
- **C4:** [code_saturne 9.1: Radiative thermal transfers](https://code-saturne.org/doc/code_saturne-9.1/advanced_radiative_thermal.html).
- **O1:** [OpenFOAM 13 User Guide: Mesh boundary](https://doc.cfd.direct/openfoam/user-guide-v13/boundaries).
- **O2:** [OpenFOAM 13 User Guide: Boundary conditions](https://doc.cfd.direct/openfoam/user-guide-v13/boundary-conditions).
