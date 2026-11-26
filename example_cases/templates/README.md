# Case templates

Every case template SolidCraft offers for a new study (**New project** in the Multiphysics
workspace, or `multiphysics.study.create {"template": …}`), each a small, complete case of
one solver of the companion. `real_reference_runner_all_templates` in
`apps/multiphysics-host` runs all of them through a real companion; this page is written by
`make_index.py` from the templates themselves and that run's results.

56 templates.

| Family | Template | Shows | Solver | Ran |
|---|---|---|---|---|
| Acoustics | `acoustic-kspace2d` | 2D k-space | `acoustic_kspace2d` | yes, 0.0 s |
| Antennas | `wire-dipole` | half-wave dipole impedance sweep | `wire_antenna` | yes, 0.4 s |
| Atomistics | `lennard-jones-md` | Lennard-Jones chain MD | `lennard_jones_md1d` | yes, 0.3 s |
| Battery | `battery-discharge` | SPMe cell 1C discharge and rest | `battery_cell` | yes, 0.0 s |
| Battery | `battery-pack` | four-cell series pack cycle | `battery_cell` | yes, 0.0 s |
| Bioheat | `bioheat-hotspot` | SAR hot spot (Pennes) | `pennes_bioheat2d` | yes, 0.1 s |
| Chemistry | `chemistry-batch` | isothermal batch | `chemistry_batch` | yes, 0.0 s |
| Circuits | `spice-rc-filter` | SPICE RC filter with diode clamp | `spice_circuit` | yes, 0.0 s |
| Coupled | `joule-busbar` | busbar Joule heating | `joule_heating2d` | yes, 0.1 s |
| Diffusion | `diffusion2d` | structured 2D | `diffusion2d` | yes, 0.0 s |
| EM | `edge-cavity` | Nédélec FEM cavity port impedance | `edge_cavity3d` | yes, 4.5 s |
| EM | `fdtd3d-cavity` | 3D FDTD cavity pulse | `maxwell3d` | yes, 0.0 s |
| EM | `maxwell1d` | Yee Maxwell 1D | `maxwell1d` | yes, 0.0 s |
| Electrochemistry | `galvanic-corrosion` | galvanic zinc/oxygen couple | `galvanic_corrosion` | yes, 0.0 s |
| Flow | `flow-cavity-pimple` | transient cavity (PIMPLE) | `collocated_flow3d` | yes, 0.0 s |
| Flow | `flow-cavity-simple` | steady cavity (SIMPLE) | `collocated_flow3d` | yes, 0.0 s |
| Flow | `flow-channel` | channel (PISO) | `collocated_flow3d` | yes, 0.0 s |
| Flow | `flow-electronics-enclosure` | primitive electronics enclosure (isothermal) | `collocated_flow3d` | yes, 4.1 s |
| Flow | `flow-electronics-enclosure-fan` | fan-curve enclosure with exhaust filter (isothermal) | `collocated_flow3d` | yes, 34.6 s |
| Flow | `lbm-taylor-green` | D2Q9 Taylor–Green | `lbm_d2q9` | yes, 0.0 s |
| Flow | `sod-weno` | Sod shock tube (WENO) | `euler1d` | yes, 0.1 s |
| Fracture | `phase-field-bar` | phase-field notched bar | `phase_field_fracture1d` | yes, 0.0 s |
| Granular | `dem-settling` | DEM spheres settling in a box | `dem_granular` | yes, 0.5 s |
| Heat | `flow-electronics-enclosure-cht` | fan-cooled enclosure with conducting, powered components | `collocated_flow3d` | yes, 7.2 s |
| Heat | `flow-electronics-enclosure-natural` | sealed enclosure, natural convection | `collocated_flow3d` | yes, 42.5 s |
| Heat | `heat2d` | transient FEM 2D | `heat2d` | yes, 0.0 s |
| Heat | `solid-heat` | solid conduction FVM 3D | `solid_heat3d` | yes, 0.0 s |
| Heat | `thermal-radiation` | prescribed flow and radiation | `thermal_transport3d` | yes, 0.2 s |
| High order | `dg-advection` | DG advection | `dg_advection1d` | yes, 0.0 s |
| Kinetic | `landau-damping` | Landau damping (Vlasov–Poisson) | `landau_damping1d` | yes, 0.8 s |
| Microstructure | `cahn-hilliard` | Cahn–Hilliard spinodal decomposition | `cahn_hilliard1d` | yes, 0.0 s |
| Multiphase | `multiphase-bubble-column` | Euler–Euler bubble column | `multiphase_box` | yes, 0.1 s |
| Multiphase | `multiphase-vof-stir` | stirred compressible VOF | `multiphase_box` | yes, 2.8 s |
| Optics | `optics-lens` | sequential lens | `optics_sequential` | yes, 0.0 s |
| Photonics | `photonic-waveguide` | silicon wire waveguide modes | `photonic_waveguide` | yes, 0.7 s |
| Plasma | `pic-two-stream` | two-stream instability (PIC) | `pic_two_stream1d` | yes, 0.8 s |
| Plasticity | `j2-cyclic` | J2 cyclic uniaxial stress | `j2_uniaxial` | yes, 0.0 s |
| Porous | `thmc-column` | coupled THMC column | `thmc_column1d` | yes, 0.0 s |
| Power | `power-network` | five-bus flow, PV curve, dispatch and swing | `power_network` | yes, 0.4 s |
| RF | `rf-tem` | TEM two-port sweep | `rf_tem` | yes, 0.0 s |
| Radiation | `transport-depth-dose` | slab depth dose (condensed history) | `particle_transport_slab` | yes, 0.2 s |
| Rods | `cosserat-cantilever` | Cosserat cantilever under gravity | `cosserat_cantilever` | yes, 0.2 s |
| SI/PI | `pdn-ir-drop` | power-plane IR drop | `pdn_ir_drop` | yes, 0.1 s |
| Scalar | `poisson1d` | Poisson 1D | `poisson1d` | yes, 0.0 s |
| Scalar | `poisson2d` | Poisson FEM 2D | `poisson2d` | yes, 0.1 s |
| Scalar | `poisson3d` | Poisson FEM 3D | `poisson3d` | yes, 0.1 s |
| Scalar | `scalar-mixed` | mixed-boundary FEM 2D | `scalar_diffusion2d` | yes, 0.0 s |
| Semiconductors | `tcad-pn-diode` | PN diode IV and CV | `tcad_diode1d` | yes, 0.2 s |
| Spray | `spray-evaporation` | evaporating injection in a cross-flow | `spray_injection` | yes, 1.9 s |
| Structures | `bar-modes` | axial bar modes | `bar_modes1d` | yes, 0.0 s |
| Structures | `elasticity-cantilever` | plate cantilever under gravity (FEM) | `elasticity2d` | yes, 0.0 s |
| Transport | `poly-scalar` | scalar FVM 3D | `poly_scalar_transport3d` | yes, 0.2 s |
| Turbulence | `rans-kepsilon` | k–epsilon transport | `rans_transport3d` | yes, 0.3 s |
| Turbulence | `rans-reynolds-stress` | Reynolds-stress transport | `rans_transport3d` | yes, 0.3 s |
| Turbulence | `rans-sa` | Spalart–Allmaras transport | `rans_transport3d` | yes, 0.2 s |
| Turbulence | `rans-sst-des` | SST/DES transport | `rans_transport3d` | yes, 0.0 s |
