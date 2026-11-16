# Completion plan and acceptance gates

The requested full Rust clone is **unfinished**. `parity/status.json` is the authoritative per-unit progress record. `native_reference` means the stated numerical reference behaviour exists in Rust, not that every overload, backend, performance contract or application workflow is complete. `native_partial` identifies an intentionally narrower subset. Every unit's `fully_qualified_parity` remains false in this migration release.

## Dependency order

The 0.2.0 checkpoint supplies block/complex sparse methods, serial stencil and Poisson multigrid references; general face-based geometry/operators/schemes; scalar, solid, thermal, collocated and RANS transport; radiation closures/radiosity and flow diagnostics. Remaining work in the first two stages includes distributed infrastructure, built-in AMG, variable-viscosity momentum coupling, moving/adaptive meshes, species/multiphase/reacting systems, import/restart and broader application validation. Detailed scope and intentional source differences are in `FVM.md` and the per-unit ledger.

1. Complete shared infrastructure: block/complex/distributed sparse matrices, multigrid/AMG, mixed precision, grid/SoA/stencil/partition contracts, task graph/profiling and backend selection. Preserve numerical tolerances and explicit breakdown states.
2. Port general polyhedral mesh/operators/FV matrix and runtime schemes. Then pressure-velocity/collocated/projection systems, scalar/species/thermal transport, turbulence, moving/adaptive meshes, VOF, Euler-Euler, reacting/spray/radiation and particle interactions. Keep genuine coupled residuals and conservation budgets.
3. Extend LBM descriptors/pull/in-place/compressed/precision variants, passive scalars, advanced collision/free-surface/particle/voxel workflows. Validate distinct algorithms separately; a BGK kernel does not qualify another collision model.
4. Complete remaining FEM/mixed/axisymmetric/nonlinear transport, magnetostatics, eddy current, elasticity, dynamics, adaptivity and constitutive models. Retain distinct vector/scalar/mixed representations and correct units.
5. Port 3D/cylindrical FDTD, ports/monitors/transmission lines/EMC sources, EM frequency-domain and photonic/metasurface extensions; then N-port/MoM/PEEC, antenna and asymptotic methods. Preserve complex phase, modal conventions and far-field normalization in outputs.
6. Implement the full circuit/SPICE engine and analysis/import/device seams. Keep device models, charge history, nonlinear stabilization, sparse MNA, controlled sources, noise, sensitivity, pole/zero, periodic analysis and external-adapter contracts distinct.
7. Complete remaining optics analysis/nonsequential/polarization/diffraction/interchange/optimization; reaction databases, advanced/aqueous equilibria and reactor networks; detector/radiation/particle/plasma transport; battery/DFN/pack/thermal, TCAD and power networks.
8. Port remaining DG/spectral high-order, kinetic, atomistic, porous/THMC, multibody/structural/Cosserat, and all multiphysics couplings. Transfer consistent geometry/units and expose real convergence/energy balances.
9. Add checkpoint/restart, import formats, workflow/distributed controller, Python bindings and SolidCraft worker adapters. Translate each of the upstream fixed reference cases only after its actual method is implemented.
10. Implement hardware-specific backends after CPU references are trustworthy: communicator/partition/halo/reduction contracts for MPI, device residency/kernel equivalence for SYCL or an explicitly new Rust GPU backend. Replacing SYCL with a different API requires an explicit compatibility decision, not a renamed CPU fallback.

## Gate for a unit marked complete in a future release

- Enumerate all public types, selectors, callbacks, dimensionality, boundaries, outputs and source/build conditions from the pinned header and implementation.
- Port the actual algorithm with explicit unsupported/error paths and no runtime C++ dependency.
- Compile the source equivalent as an isolated C++ oracle where available; compare numerical arrays, conservation, residuals and failure states. Source-only upstream extensions require their own build evidence.
- Add an analytical/convergence/conservation or physically meaningful regression beyond an implementation-mirroring test.
- Verify import/export round trips and representational semantics when that unit owns data interchange.
- Exercise relevant serial/parallel/device/distributed modes independently; record unavailable hardware explicitly.
- Expose the method in the intended library/case/workspace flow and verify the real user path. API implementation and GUI availability are separate claims.
- Update the ledger, limitations and evidence; qualify whole-family parity only after every required unit/method/backend is accounted for.

The complete 285-row inventory and selector/configuration index are shipped to make follow-up work concrete. Empty modules and success placeholders are not used to improve the coverage count.
