# Verification record

The executable algorithms, not labels or generated source stubs, are tested.

The full `--all-features` suite contains **93 integration tests** across 14 files. It exercises 23 actual CLI case examples, meaningful numerical invariants and 78 C++ reference scenarios containing 14,185 floating-point outputs. The baseline oracle has 32 scenarios / 6,106 outputs, and the new FVM/core/RANS oracle has 46 / 8,079. Fixture metadata records the pinned source revision, compiler, generator hash and fixture hash. `parity/upstream-files.json` retains hashes for the retrieved upstream files.

| Tests | Evidence |
| --- | --- |
| `tests/core.rs` | Sparse duplicate assembly/transposition, malformed CSR rejection, actual residuals for symmetric/nonsymmetric solves, ILU factor action, nonlinear roots, FFT modes and quaternion rotation |
| `tests/advanced_core.rs` | Forward/reverse analytical derivatives, Jv/J-transpose dot identity, actual matrix-free Newton/adjoint, implicit DAE decay/events and generalized eigenvalue/mass orthogonality |
| `tests/fem.rs` | Shape partition and gradient sums, quadrature volumes, exact 1D load, affine triangle/tetrahedron patch tests, assembled/matrix-free agreement, mixed Robin/Neumann field, heat invariance |
| `tests/physics.rs` | LBM uniform mass/velocity and Taylor-Green decay, wall forcing, Bouzidi branches, absorbing-layer energy removal, zero-dispersion limits, disallowed material overlap, Euler conservation/physical states and diffusion convergence |
| `tests/domains.rs` | RF network identities and Touchstone ordering/phase; Snell/TIR/coating/Gaussian focus; kinetic falloff, reactor decay/rollback, phase potentials, precipitation; acoustic spectral attenuation and sensors |
| `tests/structural.rs` | Plane stress/strain patch tests, free thermal expansion, fixed-free modal frequencies, piezoelectric receive/transmit and plastic/hyperelastic constitutive checks |
| `tests/upstream_parity.rs` | Separately compiled C++ values for Poisson/diffusion/heat, 4 LBM scenarios, 9 Maxwell scenarios, 4 Euler scenarios, RF/optics/chemistry/acoustics/time reversal and structural/material cases |
| `tests/cli.rs` | All 23 embedded examples execute and produce listed artifacts; unknown variants/fields fail; existing completed results are retained; association errors are rejected; Rayon and serial LBM agree bit for bit on the tested grid |
| `tests/core_extended.rs` | Block and complex sparse solves, malformed block rejection, independent stencil/CSR agreement, discrete-mode multigrid and large-residual overflow/NaN regression |
| `tests/poly_fvm.rs` | Mesh closure/shear, affine gradients, conservative/bounded face schemes, Rhie–Chow sensitivity, scalar time histories, channel continuity, sheared cavity and external preconditioner callback |
| `tests/thermal_models.rs` | Energy balance, invalid-material rollback, radiative enclosure balance, closure limits, Q/lambda2/vorticity/force/probe/y+ and report callbacks |
| `tests/rans.rs` | C++ SA/k–ε/SST/RSM paths in three time schemes, three DES shields, analytical discrete homogeneous decay, realizability and failed multi-equation rollback |
| `tests/poly_cli.rs` | Real hexahedron centres/volumes match control volumes, malformed new cases and failed linear solves leave no success manifest, SIMPLE completion differs from steady convergence |
| `tests/upstream_fvm_parity.rs` | C++ block/complex/multigrid, Cartesian/sheared geometry/operators, five face schemes, corrected Rhie–Chow, three scalar temporal paths, five flow paths, variable-property heat, solid conduction, closures and radiosity |

Comparisons use method-specific absolute and relative tolerances in the tests. For f32 LBM the C++ comparison tolerances are on the order of 1e-6. Structural/scalar/Maxwell/chemistry/acoustic tests use their own tolerances appropriate to their scale and iterative solves. Small reference comparisons do not constitute mesh/time-step convergence studies across every method, physical experiment validation, commercial parity, GPU qualification or an HPC performance benchmark.

`scripts/validate_exports.py` runs actual CLI examples and independently parses exported XML with Python's standard library. It checks point/cell counts, field arity, finite values, connectivity, cell types and PVD references. Its optional `--vtk-reader` path also checks point coordinates, cell types, field associations, component counts and every field value through `vtkXMLUnstructuredGridReader`. All 61 CLI VTU files passed with VTK 9.3.1. The additional thermoelastic output is checked separately. This does **not** claim a ParaView application/filter test. VTK is optional for export verification and is not a Rust build/runtime dependency.

## Reproduce

```sh
cargo test --locked --all-features
cargo test --locked --release --all-features
python scripts/update_coverage.py --check
python scripts/validate_exports.py --binary target/release/cfd-solve
# Optional, if Python VTK is installed:
python scripts/validate_exports.py --binary target/release/cfd-solve --vtk-reader
cargo run --release --example thermoelastic_plate -- results/plate
```

To regenerate or independently check the C++ fixtures, obtain a checkout of the exact pinned upstream revision. The script verifies header and selected source hashes before compilation:

```sh
python scripts/regenerate_goldens.py --upstream /path/to/cfd_solvers --check
python scripts/regenerate_goldens.py --upstream /path/to/cfd_solvers --group fvm --check
python scripts/regenerate_goldens.py --upstream /path/to/cfd_solvers
```

Only fixture regeneration needs Python and a C++20 compiler. Normal Cargo builds/tests use the stored numeric fixtures. The CI workflow includes Linux/Windows Rust jobs and a separate pinned-upstream oracle job. Writing the workflow does not establish that remote CI has run.

The pinned C++ structured-stencil header references weighted parallel helpers absent from its pinned parallel header. That unit is compared with an independently assembled CSR operator instead of being listed as C++ oracle-qualified. Native multigrid retains the source's cross-only prolongation convention. C++ residual/iteration bookkeeping is not blindly treated as a physical convergence certificate.

The supplied machine-readable `parity/verification.json` records which local checks actually completed. Native dependency source code is vendored for repeatable offline crate resolution; toolchains and platform linkers are not bundled.

## 0.4.0 fan and porous checks

`tests/poly_cli.rs` adds five tests:

- Fan-curve lookups never extrapolate.
- Ambiguous fan boundaries are refused.
- The operating point is checked against an independent fixed-flow run of the same system, for two different curves, with an added baffle and with the total-pressure convention.
- Porous zones are refused when they overlap, carry invalid coefficients or use reserved names.
- A slip-walled porous plug is compared with the analytic Darcy, Forchheimer and combined pressure drops, superficial speeds and dissipated power at 4 and 8 plug cells.

The interface error is first order and decreases under refinement. These are numerical-consistency checks, not physical validation of any fan or filter product.

## 0.5.0 conjugate heat and the offline suite

`tests/poly_cli.rs` adds five tests:

- The energy operator reaches second order (±0.1) on a manufactured solution at 8, 16 and 32 cells.
- Layered conduction with harmonic interfaces matches the series-resistance flow and the
  piecewise-linear profile to 1e-9.
- An adiabatic heated channel advects all power, and its bulk rise matches P/(ρ·cp·Q).
- The temperature rise is exactly linear in power. Pass/warn/fail are classified correctly,
  wall convection takes its share, and invalid thermal inputs are refused.
- `cfd-solve verify` passes and refuses to overwrite an existing folder.

The full suite has 105 passing tests.
