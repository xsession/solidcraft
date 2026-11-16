# Rust case files

The top-level `solver` selects a typed JSON variant. Unknown solver names and unknown fields fail explicitly. Each `examples/*.json` is runnable with `cfd-solve run FILE --out DIRECTORY`; the same examples are embedded in the executable under the names printed by `--help`.

The output directory must not already contain a completed run's `manifest.json`. Failed numerical calls exit with status 2 and do not write a success manifest. Partial field files from an interrupted run may remain. A run's `execution_completed` means the requested steps completed. `linear_solver_converged` refers only to a reported linear solve. `steady_state_converged` is set only for SIMPLE with explicit velocity-change and continuity tolerances; a completed fixed number of SIMPLE iterations can have `steady_state_converged:false`.

## Case variants

The 24 domain case types added in 0.7.0 are described in [DOMAIN_CASES.md](DOMAIN_CASES.md); their
field names are in the embedded examples (`cfd-solve example <name>`).

| `solver` | Main fields | Spatial/data output |
| --- | --- | --- |
| `poisson1d` | `elements`, `length`, `source`, `left`, `right` | Nodal potential VTU/CSV |
| `poisson2d` | `mesh`, `diffusivity`, `source`, `boundary` | Triangle nodal potential VTU |
| `poisson3d` | `mesh`, `diffusivity`, `source`, `boundary` | Tetrahedron nodal potential VTU |
| `scalar_diffusion2d` | `mesh`, constant `diffusivity`, `reaction`, `source`, patch `boundaries` | Triangle nodal potential VTU |
| `heat2d` | `mesh`, `config`, `steps`, `output_interval`, `initial`, `boundary`, `source` | Nodal temperature VTU/PVD series |
| `lbm_d2q9` | `config`, `steps`, `initial`, `walls_y`, `top_wall_speed` | Cell-centred density/velocity/solid VTU |
| `maxwell1d` | `config`, `steps`, `center_fraction`, `width_fraction` | Separate E/H staggered VTU/CSV |
| `euler1d` | `config`, `duration`, primitive `left`/`right`, `interface_fraction` | Cell density/momentum/energy/pressure VTU/CSV |
| `diffusion2d` | `config`, `iterations`, `source`, `boundary:[left,right,bottom,top]` | Nodal field VTU |
| `rf_tem` | `line`, `frequencies_hz`, `reference_impedance` | Touchstone S2P and complex S21 CSV |
| `optics_sequential` | `object_space_index`, `surfaces`, `rays`, optional `image_z` | Ray trace JSON and optional RMS spot metric |
| `chemistry_batch` | `species`, `reactions`, `concentrations`, `temperature`, `duration`, `config` | Initial/final concentrations JSON |
| `acoustic_kspace2d` | `config`, `sound_speed`, `density`, `steps`, Gaussian parameters, `sensor_cells` | Cell pressure/velocity VTU and sensor CSV |
| `poly_scalar_transport3d` | `mesh`, `config`, `output`, `initial`, `initial_gradient`, prescribed `velocity`, `source`, `boundaries` | Hexa scalar CellData VTU/PVD and diagnostics |
| `collocated_flow3d` | `mesh`, `config`, `output`, `method:simple/piso/pimple`, initial fields, acceleration, velocity/pressure patches and steady tolerances | Hexa velocity/kinematic-pressure/derived CellData VTU/PVD and solver diagnostics |
| `solid_heat3d` | `mesh`, `config`, `output`, `initial`, `initial_gradient`, volumetric `source`, `boundaries` | Hexa temperature CellData VTU/PVD and energy diagnostics |
| `thermal_transport3d` | `mesh`, `config`, `output`, `initial`, gradient, prescribed `velocity`, heat `source`, `absorption`, `environment_temperature`, `boundaries` | Hexa temperature CellData VTU/PVD and sensible-energy diagnostics |
| `rans_transport3d` | `mesh`, `output`, prescribed `velocity` and row-major `velocity_gradient`, tagged `turbulence` model/initial fields/config/patches | Hexa turbulence/eddy-viscosity/stress CellData VTU/PVD and step diagnostics |

Configuration structures implement `serde(default)` where documented in code: omitted configuration members use upstream-derived defaults. Examples deliberately override small grid sizes. `LinearConfig` also accepts partial controls. Library callbacks support spatially varying functions and richer source/material models than these first JSON forms.

The new FVM cases use `mesh:{"dimensions":[8,6,2],"lengths":{"x":1,"y":0.6,"z":0.1},"xy_shear":0.2}` and `output:{"steps":20,"output_interval":5}`. The shear transforms x to x + shear*y. Patches are `left`, `right`, `bottom`, `top`, `front`, `back`. Scalar conditions use `{"patch":"left","kind":"fixed_value","value":300}` or `zero_gradient`; thermal zero-gradient means adiabatic. Velocity kinds are `fixed_value`, `zero_gradient`, `slip`. Pressure kinds are `fixed_value`, `zero_gradient`. Duplicate patch assignments fail. An omitted patch keeps its solver default, so supply boundary conditions appropriate to the problem.

`collocated_flow3d` also accepts optional `obstacles` (grid_obstacles_v1), `fan` (fan_curve_v1) `porous_blocks` (porous_blocks_v1) `thermal` (conjugate_heat_v1, with optional `buoyancy` under natural_convection_v1), `bodies` (surface_bodies_v1) and `goals`/`goal_control` (goals_v1); see the README. A thermal case writes `thermal.vtu`, `thermal_report.json` and `thermal_components.csv`; temperatures are kelvin (°C also exported). A fan case writes `fan_operating_point.json/.csv`; a porous case writes `porous_zones.json`. Their pressures are in pascals (kinematic pressure times `config.density`). Fan and porous fields use SI units and reject unknown members.

RANS model tags are `spalart_allmaras`, `k_epsilon`, `k_omega_sst`, `reynolds_stress`; see the four `examples/rans-*.json` files. Their velocity is prescribed throughout the run. Each field's boundary list is explicit, e.g. `k_boundaries` and `epsilon_boundaries`. SA/SST wall distance and SST grid scale are user inputs. Configured model constants and temporal schemes are used by the numerical equations.

New FVM fields are cell-centred on real hexahedra. Flow pressure is **kinematic pressure** in m²/s²; multiply by constant density when physical pressure in Pa is needed. SIMPLE PVD coordinates are iteration counts; PISO/PIMPLE and thermal/RANS PVD coordinates are physical seconds. Reynolds stress exports six components in the order stated in its array name. No arbitrary-polyhedron connectivity is fabricated from face centres.

## Meshes and patches

For 2D FEM:

```json
{"type":"rectangle","nx":20,"ny":12,"width":0.1,"height":0.06}
```

An explicit mesh uses `{"type":"explicit","mesh":{...}}` with `nodes:[{"x":...,"y":...}]`, `triangles:[{"node":[0,1,2]}]`, `boundary_edges:[{"node":[0,1],"patch":0}]` and `boundary_node:[true,...]`. Indices are zero based. The mask has one Boolean per node. Validation rejects nonfinite coordinates, invalid indices and degenerate triangles. Explicit meshes must already be generated and consistently tagged; a general CAD mesher/importer is not included.

For 3D FEM use `{"type":"box","nx":4,"ny":4,"nz":4,"width":1,"height":1,"depth":1}` or an explicit mesh with 3D `nodes`, `tetrahedra:[{"node":[0,1,2,3]}]` and the boundary mask.

Rectangular boundary patches are 0=left, 1=right, 2=bottom, 3=top. Scalar mixed conditions use, for example:

```json
{"patch":1,"type":"robin","value":5.0,"alpha":2.0}
```

The scalar operator is `-div(k grad u) + c u = f`; Robin assembly adds `alpha*u` to the boundary matrix and `value` to its load. Omitted patches have natural zero-flux boundary terms. Supply conditions that make the problem well posed; a zero residual alone does not remove a pure-Neumann nullspace.

## Units and representation

Physical API fields use SI units when their names or upstream contracts specify SI. Optical geometric lengths use metres, while optical wavelength/material/thin-film interfaces explicitly use nanometres. Euler examples are nondimensional reference problems. LBM uses lattice units; its exported quantities are named accordingly and no physical-unit conversion is inferred. Maxwell energy is field energy per transverse area, and 2D acoustic energy is per out-of-plane depth.

The `manifest.json` records the exact resolved input, port version, source revision, backend, metrics and artifact list. There is no HTTP service, CAD study adapter, arbitrary commercial-project importer or general restart loader in this release.
