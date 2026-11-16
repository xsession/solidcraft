# Changes

## Unreleased

- 0.7.0: added 27 `cfd-solve` case types with embedded examples and the
  `domain_cases_v1` capability: `spice_circuit`, `tcad_diode1d`, `battery_cell`,
  `wire_antenna`, `maxwell3d`, `edge_cavity3d`, `pdn_ir_drop`, `galvanic_corrosion`,
  `elasticity2d`, `bar_modes1d`, `j2_uniaxial`, `phase_field_fracture1d`,
  `cosserat_cantilever`, `dem_granular`, `pic_two_stream1d`, `particle_transport_slab`,
  `pennes_bioheat2d`, `joule_heating2d`, `cahn_hilliard1d`, `dg_advection1d`,
  `landau_damping1d`, `lennard_jones_md1d`, `photonic_waveguide`, `thmc_column1d`,
  `power_network`, `multiphase_box` and `spray_injection`.
  `tests/domain_cases.rs` checks them against closed forms and conservation laws. The
  antenna case uses a Harrington pulse-basis moment method written for it, because the
  ported upstream Pocklington reference gives an unreliable feed reactance. See
  `docs/DOMAIN_CASES.md`.
- Ported the upstream `cfd/multiphysics` group (13 units): field registry with SI unit signatures,
  conservative 1-D remap and face-flux-to-cell-rate transfer, FEM/FVM/structured projections,
  partitioned fixed-point coupling with Aitken relaxation, block-coupled GMRES systems, staged
  pipelines, conjugate interface flux, electrochemical heat, axial optical deformation, lumped
  porous electro-thermal flow, EM/circuit port iteration, Pennes bioheat with SAR averaging and
  feedback, and one-way Joule heating. Also ported `fem/electromagnetics2d`
  (Electrostatics2D, DcConduction2D).
- Ported six FVM reacting/multiphase units: `species_transport` (MultiSpeciesTransport),
  `reacting` (reacting_cell_sources over a kinetics network), `premixed_flame1d` (Fisher-KPP
  progress-variable flame with stability-bounded explicit stepping), `flame_models` (premixed and
  non-premixed interfaces), `soot_radiation` (one-equation soot with optically thin radiation) and
  `vof` (bounded compressive VOF, CSF surface tension, contact-angle wall normal).
- Ported the remaining nine FVM solver units: `particles` (Stokes-drag cloud with two-way source,
  Saffman lift, virtual mass, collision, coalescence, Weber breakup, D^2 evaporation, volume
  fraction), `particle_interactions` (cell-linked deterministic collision/coalescence/breakup),
  `dense_particle_rheology` (hindered drag with regularised mu(I)), `spray` (injected parcels with
  drag, evaporation and latent heat), `euler_euler` (two-fluid transport with interphase drag),
  `projection2d` and `incompressible2d` (periodic staggered projection and Navier-Stokes),
  `thin_film` (capillary lubrication film) and `compressible_vof` (barotropic two-phase VOF).
  The upstream matrix-free CG (`core/conjugate_gradient`) is now `matrix_free_conjugate_gradient`.
- Ported the LBM CPU units: lattice descriptors (D2Q9 in-place, D3Q19, D3Q27), momentum-exchange
  and Darcy-Forchheimer diagnostics, surface voxelization with rigid transforms, regularized /
  cumulant / MRT D2Q9 collisions, D2Q5 passive scalar, 3-D VOF free surface, immersed-boundary
  point particles, visualization writers (legacy VTK, live JSON, PGM, offline HTML viewer,
  vorticity, Q-criterion, streamlines), Esoteric-Pull in-place streaming with solids, local
  forcing and Smagorinsky LES, two-lattice BGK/TRT pull, 16-bit compressed populations (half,
  bfloat16, shifted half, packed int16) and f32/f64 precision-generic pull. `io` gains OBJ and
  ASCII STL triangle-surface readers. The SYCL and distributed LBM variants remain unported.
- Ported the FEM solver group (Brinkman, penalty and mixed Stokes, thermoelastic block, SUPG/DG
  stabilization, Newmark dynamic bar, time-harmonic eddy currents, Darcy, axisymmetric
  elasticity, nonlinear magnetostatics with B-H curves, coils and Maxwell stress, cubic nonlinear
  Poisson), Green-Lagrange kinematics, phase-change enthalpy, the shared material record
  (`physics`), nodal DG (`dg`), Vlasov and Vlasov-Poisson (`kinetic`), 1-D molecular dynamics
  (`atomistics`) and acoustic post-processing.
- Ported the hybrid DAE graph, power flow (two-bus, multi-bus with continuation, swing
  transients, dispatch), the workflow runner, acoustic coupling, porous THMC / Allen-Cahn /
  fractured networks / Cahn-Hilliard with CALPHAD tables, aqueous acid speciation and the
  electrochemistry group (kinetics, corrosion models, corrosion cell, mixed potential, implicit
  electrode cell, 1-D and PolyMesh Nernst-Planck with electroneutral and Poisson potentials).
- Ported the remaining chemistry units (reaction paths, well-stirred reactor networks, the
  monolithic speciation block, NASA-7/Shomate thermodynamics, Pitzer/Davies/SIT activity models,
  binary Gibbs minimization, sorption, and PHREEQC/Cantera-YAML readers) and `fem/adaptivity`
  (residual estimator, Dorfler marking, conforming longest-edge refinement).
- Ported the twelve remaining optics units: spot diagrams and ray fans, 1-D Fraunhofer/PSF/MTF,
  Zernike fitting, birefringence, Kerr polarization, rigid transforms, field/wavelength/aperture
  models, polynomial freeform surfaces, ZMX/SEQ/LEN interchange, JSON serialization, a
  non-sequential plane scene with ghost ranking, and optimization (coordinate descent,
  differential evolution, sensitivities, Monte-Carlo tolerancing). `core::Mt19937_64` and
  `core::StdNormal` reproduce GCC libstdc++ random sequences bit for bit, so seeded runs match
  upstream.
- Library code no longer contains `unwrap`, `expect`, `panic!` or `unreachable!`: dual-number
  division follows IEEE semantics (use `checked_div` for an error), boundary faces without a patch
  and zero-corrector PIMPLE/SIMPLE calls return errors, and `Vec3::component` yields NaN for an
  out-of-range axis.
- Ported the CPU FDTD group: 3-D Cartesian Yee (PEC/PMC, anisotropy) with a byte-identical VTK
  writer, axisymmetric TM and multigrid levels, time/DFT monitors, 1-D wave ports and
  S-parameters, NF2FF, 1 g / 10 g SAR, multiconductor lines and EMC helpers. `physics::apply_to_fdtd`
  completes the material bridge. HDF5/XDMF output reports `Unsupported` because no HDF5 library
  is linked.
- Ported `rf/nport` (N-port Z/S and power-wave conversions, Touchstone SnP I/O, stability,
  gain, noise and load-pull, passivity, causality, adaptive sweeps, rational fitting, mixed mode,
  L-match), `rf/peec` with SPICE export, `rf/conductor_models` and `rf/modal_elements`. An
  upstream defect in the finite-ground Fresnel formula is kept for parity and documented; the
  new `fresnel_ground_reflection` gives the physical coefficients.
- Ported the `cfd/em` group: rectangular waveguide ports, multiconductor cables and shield
  transfer, GSTC metasurfaces, 1-D driven and cavity solves, signal/power integrity (eye, decap
  optimisation, PDN IR drop, pulses), Nedelec edge elements in 2-D and 3-D with ports, Thevenin
  and resonator Q, face-jump adaptivity with conforming Tet4 refinement, and photonic band and
  waveguide modes. `core::erf`/`erfc` add pure-Rust error functions.
- Ported the remaining core runtime units (backend selection, rayon-backed parallel helpers,
  slab decomposition, transfer stats, distributed CSR partitions and owner-routed assembly,
  mixed-precision refinement, profiler, scratch arena, SIMD packs, SoA fields, task graphs, Tri3
  mesh partitioning, distributed reverse AD, AMG adapter, kernel-tuning database) and
  `distributed/repartition`. NUMA placement, aligned allocation and the SYCL CSR facade are
  partial: thread binding and custom alignment need `unsafe`, and SYCL is replaced by a CPU
  reference.
- Ported the remaining CPU fvm units: octree hexahedral AMR with balanced refinement and
  coarsening, conservative remap and an adapt step, ideal-gas mixture thermo, dimensioned
  fields and time histories, rigid/ALE mesh motion, and Morton reordering.
- Ported `io/checkpoint` (binary restart files), the OpenFOAM polyMesh and field reader/writers
  and the XDMF mesh contract. The FVM HDF5 checkpoint validates its inputs and reports
  `Unsupported` because no HDF5 library is linked.
- `battery`: lithium-ion particles, SPM/SPMe with degradation and thermal, drive cycles and EIS; pseudo-2D Doyle-Fuller-Newman; balanced series packs; 3-D voxel conduction with electrothermal pack coupling. Matches upstream values; steps are transactional.
- `rf` MoM and antennas: thin-wire MoM with junctions, ground, loads and multiport S; RWG surface EFIE with far field; dielectric SIE and Helmholtz treecode; physical optics and SBR; dipole, array, polarization, match and optimizer utilities. Matches upstream values.
- `tcad`: 1-D Scharfetter-Gummel drift-diffusion device (equilibrium, DC, IV, transient, small-signal, C-V, Schottky and gate contacts), Fermi-Dirac statistics, MOS capacitor, long-channel MOSFET, diffusion-derived BJT, IGBT, Gummel-Poon, MOS level 2/3, VDMOS, hysteretic core and SPICE-ready static evaluators. Matches upstream values including Gummel iteration counts.
- `mechanics`: Neo-Hookean and J2 constitutive updates, objective finite-strain J2 with calibration audit, total-Lagrangian Tet4/Tri3 elements, nonlinear and peridynamic bars, penalty impact, meshfree neighbor graphs with halos, phase-field fracture and its coupled operator, and stable Cosserat rods (with the exact upstream multibody vector/quaternion helpers). Matches upstream values.
- `multibody`: exact vector/quaternion helpers, rigid bodies, joints and motors, three broad phases, sphere/plane/mesh and GJK/EPA contact, penalty/Hertz-Mindlin/impulse laws, bonded particles, rigid-body and DEM systems, implicit generalized-alpha integration, articulated trees, flexible-body interfaces, resolved/unresolved CFD-DEM coupling and transport-agnostic distributed DEM. Matches upstream values; the BVH also reports pairs upstream drops. MPI/SYCL resident units map to these CPU references.
- `distributed` and LBM decomposition: Cartesian process grids, bricks and neighbors, halo SoA grids, lattice-selective halos, device assignment, an in-process 26-neighbor halo exchange, and brick-decomposed pull LBM with virtual full/selective exchange and checkpoints. MPI and SYCL resident units map to these CPU references (native_partial).
- `workflow`: case configs with unit validation, material database, provenance and the validation catalogue; factorial and Latin-hypercube campaigns with templated case folders, a TSV registry, native/Docker/Apptainer/Slurm command plans, a runtime doctor, local process launching and log parsing; file-based control directives and Slurm parsing; and multi-server deploy, supervision, controller and Docker scaffolds generated byte-identically to upstream. Generation only; nothing runs remotely.
- `circuit`: SPICE-class MNA simulator (linear, controlled, diode/MOS/BJT/JFET/switch, N-port lines, static/dynamic/electrothermal compact devices, behavioral sources) with DC homotopy, sparse GMRES, AC, sweeps, fixed and adaptive transient and noise; a SPICE reader with expressions, .FUNC, subcircuits, .INCLUDE/.LIB and .IC/.NODESET; S-parameter, sensitivity, distortion, Fourier/THD, Monte Carlo, pole/zero and PSS analyses; RAW files; OSDI/OpenVAF/ngspice seams. OSDI modules are inspected without loading them (native_partial).
- Added the `particle` module (pushers, Monte-Carlo collisions, 1-D/2-D/3-D PIC, 3-D domain exchange, nuclear cross sections and depletion, slab Monte Carlo, plasma chemistry, wakefields and condensed-history transport with dose scoring), ported from upstream `cfd/particle` and checked against `parity/goldens-particle.json`.
- Added `python` (safe counterparts of upstream `python/c_api.h`: version string and the triangle-mesh inline XDMF writer) and `io::write_xdmf_inline`, byte-compared with `parity/goldens-python.json`; the C ABI symbols themselves stay unported because they need `unsafe`.
- `tests/multiphysics.rs`, `tests/fvm_reacting.rs`, `tests/fvm_multiphase.rs`,
  `tests/lbm_extended.rs`, `tests/fem_batch.rs`, `tests/small_physics.rs`,
  `tests/power_workflow.rs`, `tests/porous_echem.rs`, `tests/chemistry_advanced.rs`,
  `tests/fem_adaptivity.rs`, `tests/optics_design.rs`, `tests/fdtd_extended.rs`, `tests/rf_network.rs`,
  `tests/em_extended.rs`, `tests/special_functions.rs`, `tests/core_runtime.rs`, `tests/fvm_amr.rs`, `tests/io_formats.rs`, `tests/battery.rs`, `tests/rf_mom.rs`, `tests/tcad.rs`, `tests/mechanics.rs`, `tests/multibody.rs`, `tests/distributed_halo.rs`, `tests/workflow.rs`, `tests/circuit.rs`, `tests/particle.rs` and `tests/python_api.rs` check closed forms, conservation and bounds.
  The `tests/upstream_*_parity.rs` files compare value by value against goldens printed by the
  `parity/oracle_*.cpp` drivers (multiphysics, fvm_reacting, fvm_particles, fvm_flow2d, lbm,
  lbm_pull, fem_batch, small, core_power, porous_echem, chemistry, fem_adaptivity, optics, stdrand, fdtd, rf, em, core_runtime, fvm_amr, io, battery, rf_mom, tcad, mechanics, multibody, distributed, lbm-dist, workflow, circuit, particle, python; the LBM file writers byte for byte), built
  locally against the pinned upstream revision.
- Coverage is now 252 native reference, 33 partial and 0 not ported, out of 285 public units.
  Every upstream `cfd/solvers/fvm` unit now has a native reference. The ported units are recorded
  in `parity/port-map.json`, and `scripts/update_coverage.py` regenerates `parity/status.json` and
  `docs/COVERAGE.md` from it.

## 0.6.0

- `natural_convection_v1`: an optional `thermal.buoyancy` makes the coupling two-way. The energy
  equation takes an implicit Euler step with every flow step, and the Boussinesq acceleration
  `-beta (T - T_ref) g` drives the next one. Solids may declare `density_kg_per_m3` and
  `specific_heat_j_per_kg_k`; without them, the fluid capacity stands in and is reported as pseudo.
- `goals_v1`: monitored goals with a windowed relative-change test and optional early stop.
  Kinds are patch volume flow, patch mean pressure, patch heat out, and solid maximum and mean
  temperature. Outputs are `goals.json` and `goals.csv`.
- `surface_bodies_v1`: closed triangle surfaces become solids by cell-centre winding numbers,
  with a tie rule, after watertightness and orientation checks. The mesh report gives
  represented versus true volume.
- The energy system is row-scaled, so its tolerance is in kelvin.
- The verify suite adds the de Vahl Davis Ra = 1e3 cavity (Nusselt number within 5% at 16x16).
- Flow runs report `max_viscous_diffusion_number`. The collocated PISO path is not reliable above
  about 0.25; keep it near 0.05.
- Example: `flow-electronics-enclosure-natural`. The CHT example now uses goals with early stop.

## 0.5.0

- Added `conjugate_heat_v1`: an optional `thermal` block on `collocated_flow3d` solves one steady energy
  equation on all background cells after the flow run. Fluid cells conduct and advect heat with the
  conservative solver fluxes (first-order upwind). Grid obstacles become solids with per-axis
  conductivity and uniformly distributed power. Interfaces use series conductance. The temperature
  rise above the inflow is the unknown, so a small flow continuity error is not a heat source.
  Outputs are `thermal.vtu`, `thermal_report.json` (energy balance per patch) and
  `thermal_components.csv` (pass/warn/fail against optional limits).
- Added `cfd-solve verify --out DIR`, an offline self-verification suite with six embedded checks.
- Added the `flow-electronics-enclosure-cht` example. Its values are illustrative.

## 0.4.0

- Added `fan_curve_v1`: an optional `fan` on `collocated_flow3d` drives one planar exterior patch from a
  manufacturer pressure-flow curve. The volume flow is not prescribed; after the measured system pressure
  settles, the flow moves to the curve's crossing with the measured resistance line. Static and total
  pressure conventions are explicit. The curve is never extrapolated. `fan_operating_point.json/.csv`
  record the curve, convention, history and final mismatch.
- Added `porous_blocks_v1`: grid-aligned isotropic Darcy-Forchheimer zones (`viscous_resistance_per_m2`,
  `inertial_resistance_per_m`) with an implicit, Rhie-Chow-consistent drag on the momentum diagonal.
  `porous_zones.json` reports cells, superficial speed and dissipated power. The abrupt resistance jump has a
  first-order interface error that the tests show converging under refinement.
- Added `flow-electronics-enclosure-fan`: a shrouded inlet fan, three component blocks and an exhaust filter zone.
  Its fan curve and filter coefficients are illustrative, not characterized products. The flow is isothermal and laminar.

## 0.2.0

- Added native block CSR/Jacobi/GMRES, complex sparse solves, structured stencil and Poisson multigrid methods.
- Added general face-based PolyMesh validation, Cartesian/sheared hexa geometry, Gauss/least-squares/corrected operators, five face schemes and reusable workspaces.
- Added implicit scalar transport and SIMPLE/PISO/PIMPLE with Rhie–Chow flux, non-orthogonal correction, pressure gauge and Euler/BDF2/CN history paths.
- Added explicit fluid heat, implicit solid conduction, material/source callbacks, Boussinesq force, gray/optically thin radiation and enclosure radiosity.
- Added SA, k–ε, SST/DES and six-stress LRR/epsilon transport, plus local Smagorinsky/WALE/mixing-length closures.
- Added flow-derived fields, forces/probes/y+, residual callbacks and step diagnostics.
- Added five JSON case types and ten examples, with actual hexahedral cell VTU/PVD time series. SIMPLE steady convergence is distinct from iteration completion.
- Hardened RMS against overflow and NaN masking, rejected the source's ineffective collocated high-order selector, and made RANS multi-equation commitment atomic.
- Added C++ FVM/core/RANS oracle fixtures, bringing the suite to 78 scenarios / 14,185 numeric values; the Rust suite contains 93 integration tests and 23 executable examples.
- Expanded public-unit coverage from 30 reference / 10 partial / 245 unported to 53 / 11 / 221. Full Rust clone and full backend/application qualification remain unfinished.

## 0.1.0

Native Rust migration baseline with shared numerics, scalar/structural FEM, D2Q9 LBM, 1D Euler/FDTD, RF networks, sequential optics, chemistry batch, k-space acoustics and basic exports. Pinned source inventory, vendored dependencies and isolated C++ numerical fixtures included.
