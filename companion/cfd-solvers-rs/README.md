# cfd-solvers-rs 0.7.0

## Domain cases

`cfd-solve` now runs 45 case types. The 27 new ones expose the circuit, semiconductor,
battery, antenna, 3-D FDTD, Nédélec FEM, PDN, power-network, corrosion, elasticity,
plasticity, fracture, rod, DEM, multiphase, spray, PIC, particle transport, bioheat, Joule
heating, Cahn–Hilliard, DG, Vlasov, molecular dynamics, photonic mode and THMC modules. Each has an embedded example and
reference checks; see [docs/DOMAIN_CASES.md](docs/DOMAIN_CASES.md).

```sh
cargo run --release -- example wire-dipole --out results/dipole
```


## Natural convection, goals and surface bodies

New `collocated_flow3d` fields:

- `thermal.buoyancy`: `{gravity_m_per_s2, expansion_per_k, reference_temperature_k}`.
  Two-way transient coupling.
- `goals`: a list of `{name, kind, target}`.
- `goal_control`: `{window_steps, relative_tolerance, stop_when_converged, min_steps}`.
- `bodies`: a list of `{name, vertices, triangles}`, each a closed, outward-wound surface
  in metres.

```sh
cargo run --release -- example flow-electronics-enclosure-natural --out results/natural
```

`cfd-solve verify` now includes the de Vahl Davis natural-convection benchmark. Keep the
reported `max_viscous_diffusion_number` small: the collocated PISO path became unstable
above about 0.25 in testing.


## Conjugate heat transfer and offline verification

Version 0.5.0 adds `conjugate_heat_v1`. An optional `thermal` block on `collocated_flow3d`
holds:

- `fluid`: `specific_heat_j_per_kg_k` and `conductivity_w_per_m_k`;
- `inflow_temperature_k`;
- `solids`: one entry per obstacle, each with a per-axis `conductivity_w_per_m_k`
  (`[x, y, z]`), `power_w` and an optional `max_temperature_k`;
- `boundaries`: `{patch, condition}` pairs with `adiabatic`, `temperature`, `heat_flux`
  or `convection`. Unlisted exterior patches are adiabatic.

The energy solve is steady and one-way, on the settled flow. It assumes constant
properties, no buoyancy and no radiation.

`cfd-solve verify --out DIR` runs the embedded verification checks offline. It writes
`verification.json` and exits with status 2 on failure.

```sh
cargo run --release -- example flow-electronics-enclosure-cht --out results/enclosure-cht
cargo run --release -- verify --out results/verify
```

This program never opens a network connection. Vendored dependencies allow
`cargo build --offline --locked`.


## Fan-curve and porous-zone extension

Version 0.4.0 adds `fan_curve_v1` and `porous_blocks_v1` to `collocated_flow3d`
and keeps `grid_obstacles_v1`. Existing cases are unchanged.

`fan` drives one planar exterior patch from a pressure-flow curve. Supply
`patch`, an ambient `reference_patch` with a fixed-value pressure, an explicit
`pressure_convention` (`static` or `total`) and 2–64 `curve` points with
strictly increasing `volume_flow_m3_per_s` and strictly decreasing
`pressure_rise_pa`. The fan owns its patch velocity, so that patch must not
appear in `velocity_boundaries`. The flow is not prescribed. Once the measured
system pressure rise settles (`settle_tolerance`, default 0.002 per step), at
most every `update_interval` steps (default 25) the flow moves by
`relaxation` (default 0.7) towards the curve's crossing with the measured
resistance line. The curve is never extrapolated; `within_curve_range` reports
clamping. `fan_operating_point.json` and `.csv` record the curve, convention,
history, corrections and final mismatch. This is a lumped external fan with a
uniform normal inlet, not a rotor, swirl or internal momentum-source model.

`porous_blocks` are grid-aligned isotropic Darcy–Forchheimer zones. Each has
half-open `min`/`max` cell indices, `viscous_resistance_per_m2` (D = 1/K) and
`inertial_resistance_per_m` (C). The momentum sink per unit mass is
`-(nu*D + 0.5*C*|u|) u` on the superficial velocity, assembled implicitly on the
momentum diagonal. Zones must not overlap each other or obstacles.
`porous_zones.json` reports cells, mean superficial speed and dissipated power.
Porosity storage/convection, anisotropy and thermal effects are not modelled.
The collocated interpolation at an abrupt resistance jump has a first-order
error. On a 4-cell plug the pressure drop is about 6% low. Refinement reduces it
(2% at 8 cells, 0.9% at 16).

```sh
cargo run --release -- example flow-electronics-enclosure-fan --out results/enclosure-fan
```

That example has a shrouded 40 × 40 mm inlet opening, the PCB/component blocks
and a 10 mm exhaust filter zone. Its fan curve and filter coefficients are
illustrative values, not characterized products. The flow is isothermal and
laminar, with no conjugate heat transfer.

## Grid-obstacle airflow extension

Version 0.3.0 adds `grid_obstacles_v1`; the upstream revision below identifies
the original port baseline, not these local changes. Existing cases remain
supported. `flow-electronics-enclosure` is a prescribed-flow, isothermal laminar
example with a coarse PCB block and two component blocks. It does not model a
fan curve, conjugate heat transfer, or arbitrary CAD geometry.

```sh
cargo run --release -- example flow-electronics-enclosure --out results/enclosure
```

`obstacles` is an optional array of `{name,min,max}` objects on
`collocated_flow3d`. Bounds are half-open grid-cell indices, not physical lengths.
The mesh requires zero shear, at most 32,768 background cells and 128 disjoint
blocks. Each exposed obstacle patch is a stationary no-slip wall. Solver and VTU
cells derive from the same retained fluid cells. Disconnected fluid domains,
invalid bounds, and incompatible obstacle boundary overrides are rejected.

`mesh_report.json` records represented geometry, fluid/excluded cell counts,
boundary areas, signed volume/mass flux and relative boundary flux imbalance.
It is numerical evidence for this discretized geometry, not experimental
validation. The exported pressure is kinematic pressure, in m²/s².

A native Rust port of selected numerical methods from [xsession/cfd_solvers](https://github.com/xsession/cfd_solvers), pinned to revision `7cceaaaddafaf2515e01a86406e13d6037fdbf9a`.

**Every upstream public header unit now has a native Rust implementation; 33 remain partial.** Of 285 upstream public header units, 252 have native reference implementations, 33 have partial native coverage, and 0 remain unported. Header counts are an inventory of public units, not independent solver counts or a percentage of implementation effort. Every unit and its method selectors are retained in [parity/status.json](parity/status.json).

The Rust library and executable contain no C++ solver linkage or FFI. They use safe Rust, with `unsafe_code = "forbid"` in this crate. The C++ drivers under `parity/` generate test fixtures only; neither a C++ compiler nor an upstream checkout is required for normal Rust builds or tests. Third-party Rust dependencies retain their own licenses and implementation policies.

## Build and run

Install Rust 1.85.1 with the platform's normal native linker. The toolchain version is pinned in `rust-toolchain.toml`. Windows MSVC builds need the Windows SDK and C++ Build Tools linker; the solver implementation itself is Rust.

```sh
cargo test --locked --all-features
cargo build --locked --release --all-features
cargo run --release -- example sod-weno --out results/sod
cargo run --release -- example heat2d --out results/heat
cargo run --release -- example flow-channel --out results/channel
cargo run --release -- example flow-cavity-pimple --out results/cavity
cargo run --release -- example rans-sst-des --out results/turbulence
cargo run --release -- run examples/maxwell1d.json --out results/em
cargo run --release -- coverage
```

`build.sh` and `build.bat` run tests and create the native executable. The `parallel` feature enables Rayon for D2Q9 cell updates on grids of at least 8,192 cells. Other methods remain serial. `step_serial` is also available for comparisons. The feature does not provide MPI, GPU, CUDA or SYCL execution.

The source archive includes vendored Rust dependency sources. With the pinned toolchain and linker already installed, `cargo build --offline --locked --release --all-features` builds without fetching crates. Any included Linux executable is a convenience build for the documented environment; build from source on other targets.

The supplied archive also includes `bin/linux-x86_64/cfd-solve`, built for Linux x86_64 with glibc 2.35 or newer, and actual outputs under `sample-results/`. For example, open `sample-results/heat2d/heat.pvd` for a ready-made post-processing dataset, or run:

```sh
./bin/linux-x86_64/cfd-solve example lbm-taylor-green --out results/lbm
```

## Implemented methods

| Area | Native methods in this release | Access |
| --- | --- | --- |
| Shared numerical core | CSR/block CSR, pivoted block Jacobi, complex sparse solves, Jacobi/ILU(0), PCG, BiCGStab, restarted GMRES, Poisson multigrid, structured stencil, dense real/complex and tridiagonal solves | Library |
| Nonlinear methods | Newton with line search, explicit finite differences, actual matrix-free JFNK, Jacobian actions, adjoint actions, implicit DAE steps and events | Library |
| Differentiation and spectral methods | Forward/reverse AD, generalized inverse-iteration eigenmodes, radix-2 FFT in 1D/2D, periodic derivatives/fractional operators | Library |
| Scalar FEM | Line/Tri/Tet meshes; all seven upstream linear reference-element shapes; Poisson 1D/2D/3D, P1 reaction-diffusion with mixed conditions, transient heat | Library + JSON CLI |
| Structural FEM | Plane stress/strain, component constraints, body force and thermal strain; axial bar modes; coupled piezoelectric bar receive/transmit; elastic/plastic/Neo-Hookean/contact laws | Library + Rust example |
| LBM | Single-precision D2Q9 BGK pull, Guo forcing, moving/stationary/Bouzidi walls, velocity inlet and pressure outlet, Taylor-Green and diagnostics | Library; basic initialization/walls through CLI |
| FVM geometry/operators | General face-based PolyMesh; Cartesian/sheared hexa builders; Gauss/least-squares gradients, divergence, orthogonal/corrected Laplacians, linear/upwind/bounded/TVD interpolation, reusable workspace and CFL bound | Library; generated hexa cases through CLI |
| FVM scalar/flow | Implicit scalar transport, SIMPLE/PISO/PIMPLE, Rhie–Chow pressure flux, non-orthogonal correctors, pressure gauge, Euler/BDF2/CN; structured diffusion and ideal-gas Euler 1D MUSCL/WENO5-Rusanov SSPRK3 | Library + JSON CLI |
| FVM heat/radiation | Solid conduction, explicit fluid heat, material/source callbacks, Boussinesq force, optically thin source, gray exchange and radiosity | Library; constant properties and prescribed velocity through CLI |
| FVM turbulence | SA, k–ε, SST with DES shielding, six-component Reynolds-stress/ε transport; mixing-length/Smagorinsky/WALE closures | Library; prescribed-velocity RANS through CLI |
| Flow diagnostics | Vorticity, Q, lambda2, enstrophy, y+, pressure force, probes and residual callbacks | Library; derived cell fields and step diagnostics through CLI |
| FDTD | Yee Maxwell 1D, PEC/Mur/PML/CPML, conductive media, Debye/Drude/Lorentz ADE, TFSF and lumped parallel RLC | Library; Gaussian reference case through CLI |
| RF | Two-port Z/S/ABCD algebra, cascades, TEM/microstrip/coax/stripline/CPW/TE10 approximations and Touchstone S2P | Library; TEM sweeps through CLI |
| Optics | Ray reflection/refraction/TIR; plane/sphere/conic/asphere sequential traces; paraxial analysis, Sellmeier materials, Gaussian beams, Fresnel/Jones/Stokes and normal thin films | Library; sequential traces through CLI |
| Chemistry | Mass action, Arrhenius/reversibility/orders/third bodies/falloff, ideal phases, binary precipitation and adaptive implicit isothermal batch integration | Library; batch integration through CLI |
| Acoustics | 2D k-space FFT solver, heterogeneous medium, PML, nonlinearity, power-law absorption, source arrays, sensors and time reversal | Library; Gaussian/sensors through CLI |
| Output | Validated ASCII VTU, PVD, CSV, JSON and RI Touchstone S2P | Library + CLI |

The Rust CLI exposes 45 case types with 52 embedded examples. It does not claim compatibility with all 90 fixed reference cases in the upstream dispatcher. Library features are broader than the JSON forms. See [docs/CASE_FORMAT.md](docs/CASE_FORMAT.md), [docs/FVM.md](docs/FVM.md) and [docs/COMPATIBILITY.md](docs/COMPATIBILITY.md).

## Post-processing

Open `heat.pvd`, `flow.pvd`, `solid_heat.pvd`, `thermal.pvd`, `rans.pvd` or `transport.pvd` in ParaView for their time series. `steady_iterations.pvd` uses SIMPLE iteration number, with physical time unchanged. New FVM cases export real hexahedral connectivity and **CellData**. Flow exports velocity, kinematic pressure, vorticity, Q, lambda2 and enstrophy; every step has a diagnostic JSON record. Maxwell E and H fields use separate meshes to retain the Yee half-cell offset. LBM/acoustic/Euler finite-volume fields are **CellData**; FEM and structured nodal diffusion fields are **PointData**. RF sweeps export `.s2p` and complex components in CSV. Optical rays and batch concentrations have domain-specific JSON outputs.

The built-in exporters do not embed ParaView, its filters, remote rendering or its Python server. HDF5/XDMF and OpenFOAM interchange have not been ported. `scripts/validate_exports.py` independently checks XML arrays, associations and PVD references. With `--vtk-reader`, it also compares geometry and fields read through VTK: all 61 CLI VTU files passed with VTK 9.3.1. The ParaView application/filter runtime has not been exercised.

One-way thermal-to-structural transfer is demonstrated with real native solver calls:

```sh
cargo run --release --example thermoelastic_plate -- results/plate
```

## Verification and remaining work

The verification suite runs over 400 Rust tests across 74 integration-test files (408 with default features), including value-by-value comparison with 33 separately compiled C++ oracles (682 reference scenarios, 69,809 numeric outputs plus byte-compared text outputs), and analytical, conservation, patch and serial/Rayon agreement checks. Tests cover the listed small reference problems; they do not establish commercial-solver parity, arbitrary-geometry robustness or hardware/backend equivalence. See [docs/VERIFICATION.md](docs/VERIFICATION.md) and both `parity/goldens*-metadata.json` files.

The 33 partial units are the parts that cannot be built or verified in this offline, safe-Rust crate: SYCL/GPU device residency (LBM, FDTD, FVM, DEM and sparse kernels), MPI transports and runtimes, HDF5 output, `dlopen` of OSDI modules, 64-byte aligned allocation and thread affinity, the exported C ABI of `python/c_api.h`, the Gmsh/VTU readers of `io/mesh_io.hpp`, and C++-ABI-level qualification of some core numerics. Each one has a CPU or transport-agnostic reference path and its limitation recorded in `parity/port-map.json`. Coupled application-level validation (turbulent multiphase, full multiphysics workflows) remains future work.

[docs/COVERAGE.md](docs/COVERAGE.md) lists all 285 public units. [docs/MIGRATION.md](docs/MIGRATION.md) defines the remaining dependency order and acceptance gates. Pending units have no fake solver implementations or success stubs.

`parity/upstream-units.json` retains type/configuration/selector and lexical symbol discovery data. The 4,217 symbol candidates include helpers and private/overloaded call candidates; they are not an AST or a count of independent methods/solvers. All 90 upstream command presets are separately tracked in `parity/upstream-cli-cases.json` without claiming exact Rust CLI compatibility.

## License

MIT, retaining the upstream contributor notice in [LICENSE](LICENSE). Adapted algorithms are identified by the pinned public-unit/source manifest. Vendored dependencies retain their original license files. This project is separate from the original C++ repository and the SolidCraft workspace design artifacts.
