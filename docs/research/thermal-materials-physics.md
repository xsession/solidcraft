# Thermal materials and physics setup research

Research date: 2026-11-13. Scope: public vendor documentation and public property-reference documentation, followed by an original design proposal for SolidCraft. No vendor source code, material catalog, screenshot, help text or simulation template is incorporated here. Documentation describes a workflow; it does not establish capability in SolidCraft.

## Existing SolidCraft boundary

The workspace already has a 59-profile reference catalog in `assets/render/physical/profiles.json`, with `sources.json`, `quantities.json` and `exports/native-cards.json` in the same catalog tree. It includes copper, aluminium, FR-4 and polymers. A reference profile is not an automatically assigned, fully characterized native CFD material. The native engineering card in `crates/doc/src/physical_material.rs` supports density, three principal thermal conductivities in material-local axes, specific heat, expansion and mechanical properties. It carries property sources and selected conditions. Its validation establishes numerical consistency, not source authenticity, applicability, certification or solver execution.

The existing card has no native fluid-viscosity property, surface-emissivity property, general temperature-function evaluator or automatic CFD material binding. A conductivity tuple with a material-axis label does not itself create an oriented body-to-solver assignment. Existing `material.create` and `material.edit` commands with a `physical_card` are the authoritative project-storage route. Reuse those commands and design-owned cards instead of creating an unrelated second native material store. CAD mass density, engineering candidate data and rendering appearance must remain distinguishable.

This document proposes contracts and interaction behavior. The current library-navigation increment is scoped to reusing `render.library.list`, `render.library.describe` and `render.library.physical.prepare`, followed by explicit candidate review, `material.create` and the existing material editor. Its scope adds no new physics-card variants. This document does not implement conjugate heat transfer, radiation, buoyancy, contact heat transfer or moving-frame flow in the current numerical companion.

## Comparison of documented workflows

| Environment | Material and physics setup pattern | Useful lesson for SolidCraft |
|---|---|---|
| Simcenter FLOEFD | A reusable Engineering Database defines solids; component assignments refer to those definitions. Its tutorial treats PCB directional conductivity, global gravity and regional fluid initialization separately. [1] | Separate reusable definitions from assignments and distinguish global physics from regional initial conditions. |
| SOLIDWORKS Flow Simulation | The official student workbook describes predefined and user-defined fluid/solid properties, including tabular dependencies, contact resistance, radiative surfaces and device definitions. This is historical workflow evidence rather than a current-version interface audit. [2] | A library should accommodate different property families without presenting them as one interchangeable bulk-material model. |
| Ansys Icepak | Conductivity can use object-relative orthotropic axes. Material properties can have constant, linear or piecewise temperature dependence. Surface radiation and contacts are configured independently. [3–7] | Expose orientation, property evaluation and assignment scope explicitly. |
| SimScale | Material definitions are assigned to volumes. Default definitions can be customized and reusable user/company libraries are available. Its electronics tutorial calls for exactly one material per part and uses visibility controls during assignment. [8–10] | Show coverage, missing assignments and duplicates; keep assigned bodies easy to reveal and inspect. |

## Three assignment types

### Volume material

A volume assignment identifies the solid or fluid region, the exact material revision and an orientation where needed. Show selected components in the viewport, the number of assigned regions and their names before Apply. Provide a filter for unassigned regions and a command to reveal the assigned objects. The material editor should show only properties relevant to the selected analysis, while preserving the full definition.

Bulk conduction uses conductivity. Transient heat storage uses density and specific heat. Fluid momentum may require dynamic viscosity, density and a density model. These are separate requirements; a card with density alone is not sufficient to configure a thermal-fluid case.

For PCB and laminate work, support isotropic conductivity, orthotropic principal values and a through-thickness specialization. An orientation preview should identify the material axes on the selected component. SimScale documents custom orthogonal axes and a visible cross-plane direction; Icepak relates anisotropy to each object's orientation. [3, 11] An original SolidCraft implementation should let several occurrences share a material definition while retaining their own assignment orientations.

Do not present a generic FR-4 name as an automatically characterized populated board. Provide a place for the stackup, homogenization method and supporting measurement or supplier record.

### Surface finish and radiation

A surface assignment identifies faces, the participating side and a finish revision. Store emissivity independently from bulk conductivity and visual roughness. Different faces of the same solid may have different coatings or finishes.

The radiation model must distinguish exchange among modeled surfaces from exchange with a prescribed environment. The user also needs to know whether a surface is opaque, transparent or partially transmitting, where the connected model supports those distinctions. SimScale documents these alternatives; Icepak exposes exchange with objects versus a reference-temperature environment. [6, 12] A material checkbox alone does not supply view factors, occlusion handling or a radiation solver.

Proposed validation includes finite emissivity within its model's bounds, explicit temperature and wavelength applicability when known, compatible radiation-model support and face-side coverage. Missing data should be reported as missing, rather than silently replaced by a universal coating value.

### Contact between regions

Model a contact as an interface with two identifiable sides and an overlap area. Show both sides highlighted, along with the interpretation and units of its resistance. Offer explicit alternatives:

| Input | SI unit | Meaning |
|---|---|---|
| Total resistance | K/W | Resistance of the whole specified contact |
| Area-specific resistance | m²·K/W | Resistance normalized by contact area |
| Area conductance | W/(m²·K) | Reciprocal of area-specific resistance |
| Thin layer | thickness in m and conductivity in W/(m·K) | Layer contribution before any additional interface resistance |

Icepak distinguishes total and area-based quantities and warns about selected areas that differ from actual contact areas. SimScale documents that assigning its thin-layer definition on both contact faces adds two layers. [5, 13] SolidCraft should represent one physical interface once, detect duplicate assignments and avoid requiring users to remember a product-specific double-counting rule.

## Property evaluation contracts

Extend the existing engineering card through the established format/version mechanism. Proposed additions should include:

- An explicit phase and intended physical model; unknown properties remain absent.
- SI storage with display-unit conversion at input/output boundaries.
- Constant, temperature-table and subsequently formula-backed property models, each explicitly typed.
- A reference temperature and pressure where applicable, a validity interval and a defined out-of-range policy.
- Strictly ordered temperature coordinates, finite values and property-specific positivity checks.
- An interpolation rule and a decision to reject extrapolation unless a reviewed model explicitly supports it.
- Material orientation separate from world placement and from the shared library definition.
- Property-level source, condition, uncertainty and quality information.

Icepak's public documentation demonstrates constant, linear and piecewise temperature dependence. [4] The proposed SolidCraft contract should additionally make evaluation failures observable: an unsupported table or an out-of-range case must not quietly fall back to an unrelated constant.

The selected solver adapter should declare which properties and models it consumes. Preflight should distinguish a valid saved definition, complete study assignments, adapter support and actual numerical execution. Unsupported physics remains preparation-only.

## Gravity, buoyancy and velocity

Gravity is a named-frame vector with acceleration units, not a setting inferred from the camera. Show its direction in the viewport and resolve it through the study coordinate frame. Natural convection also needs a temperature-dependent density or supported buoyancy approximation and an energy equation. Icepak's setup documentation separates gravity and solution initialization, and requires gravity for natural-convection configurations. [7]

If a Boussinesq model is added, its reference density, reference temperature and thermal expansion coefficient must be explicit. SimScale documents these as part of its buoyancy approximation. [14] Applicability checks belong to the selected physical model; do not turn one product's example air values or broad temperature guidance into universal limits.

The phrase “system velocity” should resolve to a specific assignment:

| Setting | Scope | Required distinction |
|---|---|---|
| Initial velocity | All fluid regions or selected subregions | Starting state; does not continuously impose inflow |
| Inlet or freestream velocity | Boundary selection | Imposed flow condition with an explicit direction convention |
| Moving wall velocity | Wall faces | Surface motion, separate from a fluid-region frame |
| Rotating wall | Wall faces, axis and origin | Surface rotational speed, not a rotating fluid zone |
| Moving frame | Fluid region | Translation, angular velocity, origin and absolute/relative convention |

Fluent documentation explicitly distinguishes absolute from cell-zone-relative initialization. SimScale documents moving-wall and rotating-wall inputs and states that a rotating wall does not represent a rotating zone. [15, 16] The SolidCraft UI should display the frame and velocity interpretation next to every vector field. It must not use one vector to set all these concepts at once.

## Reusable internal library

Build on the existing native material/card identity and command system. The proposed library should distinguish shipped references, editable user definitions, project-local revisions and immutable run snapshots. Duplicating a reference creates a user definition; editing that definition should not rewrite retained evidence from an earlier run.

Every imported property should retain its source title, publisher or author, URL/document/page, retrieval date, units, measured or calculated status and conditions. Record a redistribution claim separately from a numerical validation result. The current `Source.license` field is already documented as a caller's claim rather than permission verification.

NIST's WebBook offers fluid thermophysical properties, but availability is not blanket redistribution permission. NIST's SRD notice states restrictions on reproduction and storage; check the applicable dataset terms before bundling any records. [17–19] IAPWS provides recommended water-property formulations, including viscosity and conductivity releases, but this research did not establish blanket permission to redistribute its documents or an implementation. [20–22]

A practical first library release can expose native cards and reusable user records without claiming a complete certified property database. Add curated property records only with reviewed provenance and usable conditions. Never copy another simulator's catalog or populate missing values with zero. A material-library action should show its stored data and limitations, not imply that a disconnected solver now uses them.

## Official-source bibliography

Accessed 2026-11-13. These are behavior and reference-data sources, not licenses to copy a product's assets or catalog.

1. [Siemens, Simcenter FLOEFD for Solid Edge tutorial](https://static.sw.cdn.siemens.com/siemens-disw-assets/public/1qqiCDdtFFrqUsQ2cvnAUH/en-US/simcenter-floefd-solid-edge-README.pdf).
2. [SOLIDWORKS, Flow Simulation Student Workbook](https://www.solidworks.com/sw/images/content/Training/SolidWorks_Flow_Simulation_Student_Workbook.pdf).
3. [Ansys Icepak, Defining Anisotropic Thermal Conductivity Tensors](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v242/en/Subsystems/Icepak/Content/Materials/DefiningAnisotropicThermalConductivityTensors.htm).
4. [Ansys Icepak, Defining Properties Using Temperature-Dependent Functions](https://ansyshelp.ansys.com/public/Views/Secured/corp/v242/en/ice_ug/ice_ug_sec_model_prop_temp_depend.html).
5. [Ansys IcepakFEA, Thermal Resistance Types](https://ansyshelp.ansys.com/public/views/secured/electronics/v261/en/subsystems/IcepakFEA/Content/IcepakFEA/ThermalContact_ResistanceTypes.htm).
6. [Ansys Icepak, Stationary Wall Properties](https://ansyshelp.ansys.com/public/Views/Secured/Electronics/v252/en/Subsystems/Icepak/Content/Icepak/StationaryWallProperties.htm).
7. [Ansys Icepak, Solve Setup Dialog](https://ansyshelp.ansys.com/public/views/secured/electronics/v252/en/subsystems/icepak/content/Icepak/IcepakSolveSetupDialog.htm).
8. [SimScale, Materials](https://www.simscale.com/docs/simulation-setup/materials/).
9. [SimScale, How to Define a Custom Material](https://www.simscale.com/knowledge-base/custom-materials-in-simscale/).
10. [SimScale, Electronics Box Cooling Tutorial](https://www.simscale.com/docs/tutorials/thermal-management-tutorial-cht-analysis-of-a-raspberry-pi-using-ibm/).
11. [SimScale, How to Set Orthotropic Thermal Parameters](https://www.simscale.com/knowledge-base/how-to-set-anisotropic-thermal-properties/).
12. [SimScale, Radiation Behavior](https://www.simscale.com/docs/analysis-types/convective-heat-transfer-analysis/radiation/).
13. [SimScale, Contact Resistance](https://www.simscale.com/docs/simulation-setup/advanced-concepts/contact-resistance/).
14. [SimScale, Boussinesq Approximation](https://www.simscale.com/docs/simwiki/cfd-computational-fluid-dynamics/what-is-boussinesq-approximation/).
15. [Ansys Fluent, Initializing the Solution](https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/flu_ug/flu_ug_sec_solve_initialize.html).
16. [SimScale, Wall Boundary Conditions](https://www.simscale.com/docs/simulation-setup/boundary-conditions/wall/).
17. [NIST Chemistry WebBook, Thermophysical Properties of Fluid Systems](https://webbook.nist.gov/chemistry//fluid/).
18. [NIST, Standard Reference Data Public Law](https://www.nist.gov/srd/public-law).
19. [NIST, Copyright and Licensing Statements](https://www.nist.gov/open/copyright-fair-use-and-licensing-statements-srd-data-software-and-technical-series-publications).
20. [IAPWS, Releases and Guidelines](https://iapws.org/documents/release).
21. [IAPWS, Ordinary-Water Viscosity Release](https://iapws.org/technical-guidance/release/viscosity).
22. [IAPWS, Ordinary-Water Thermal Conductivity Release](https://iapws.org/technical-guidance/release/ThCond).
