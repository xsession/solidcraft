# Compatibility boundaries

This is an idiomatic Rust reimplementation of selected reference methods. Names, ownership, error handling and storage are not binary/source compatible with C++. Result validity is established for the tests listed in `VERIFICATION.md`, not every possible input.

| Contract | Current state |
| --- | --- |
| Source baseline | Git revision `7cceaaaddafaf2515e01a86406e13d6037fdbf9a`, tree `7eac84f5d0babcc9ee89b94368859ba5e4f0365d` |
| Public-unit inventory | All 285 public headers indexed, including enums/configuration declarations and build evidence |
| Native reference implementations | 53 units; scope stated individually; not full API qualification |
| Partial native coverage | 11 units; gaps stated individually |
| No native implementation | 221 units |
| Upstream CLI | Original 90 fixed reference cases are not reproduced wholesale; new typed Rust runner exposes 18 case types and 23 examples |
| Runtime | Native Rust, no C++ solver runtime or FFI |
| Floating point | f64 for most methods, f32 for D2Q9; toleranced comparisons, not a promise of bit-identical cross-compiler execution |
| Parallelism | Rayon D2Q9 cell updates only; reductions and other domains are serial |
| GPU/SYCL/MPI | Unported; no CPU backend relabelled as a device backend |
| Linux | Built and tested in the recorded environment |
| Windows | Build wrapper and CI job supplied; Windows execution has not been performed in this session |
| SolidCraft integration | Portable library/JSON worker seam only; no SolidCraft source modification or worker adapter is claimed |
| Post-processing | VTU/PVD/CSV/JSON/S2P exports; VTK 9.3.1 XML reader checked; no embedded ParaView filters, HDF5/XDMF or OpenFOAM interchange |
| Persistence | Resolved input/output manifest; checkpoint/restart and distributed controller persistence unported |

Key algorithmic conventions retained include upstream D2Q9 direction ordering and float precision; Fourier sign/normalization and zero Nyquist derivative; Poisson1D midpoint and Poisson2D centroid loading; Tet4 order-2 quadrature; implicit consistent-mass heat; upstream plane-strain thermal factor; bar consistent mass; the mixed-unit piezoelectric absolute pivot convention; Roe-characteristic WENO5 and physical-state fallback; ADE dispersion/TFSF half-step timing; k-space correction and sensor recording order; and step-doubling reactor state commitment.

Rust strengthens some malformed-input checks, replaces exception paths with `Result` in the public numerical APIs, uses owned vectors/closures instead of C++ allocator/SIMD/backend contracts, and changes some assembled Dirichlet systems to symmetric elimination on the full matrix. Sparse vector/storage layout, backend residency and iteration counts are not promised to match C++.

FVM conventions and intentional differences are documented in [FVM.md](FVM.md). In particular, higher-order interpolation exists as operators, while collocated momentum remains the source's upwind operator. Rust rejects unsupported collocated selector values. An external AMG callback is a supplied preconditioner interface, not an implemented AMG solver. RANS equations are native transport solvers but velocity/flux are supplied; automatic eddy-viscosity feedback into momentum is pending.

The upstream source-only files for GSTC, photonic modes, nonlinear polarization, RF integral-equation extensions/surface MoM and FDTD-specific HDF5 remain distinct in the inventory. Their presence in the source is not treated as an upstream compiled-runtime qualification or as a completed Rust port.

Broader commercial workflow parity—mesh import/repair, arbitrary project authoring, adaptive refinement, restart, monitors/far-field/modal reports, optimization, distributed scheduling and coupled solution control—requires separate implementations and validation. None is inferred from a small reference case.
