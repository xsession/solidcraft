# Native Multiphysics companion host

This application-layer library runs the optional `cfd-solve` companion beside the SolidCraft window thread. It links no numerical solver into the CAD engine. The browser build returns no local process host.

Set `SOLIDCRAFT_CFD_SOLVE` to the native Rust port's executable path, or install `cfd-solve` on `PATH`. Host discovery does not start a process. An invalid explicit path produces an actionable job failure when a run starts.

The adapter accepts native Rust companions **0.2.0**, **0.3.0**, **0.4.0**, **0.5.0**, **0.6.0** and **0.7.0**, upstream baseline revision `7cceaaaddafaf2515e01a86406e13d6037fdbf9a`, and backend `native_rust_serial_cpu`. Cases with nonempty grid obstacles require **0.3.0** or **0.4.0** and the explicit `grid_obstacles_v1` manifest capability. A `fan` requires **0.4.0** or later with `fan_curve_v1`. Nonempty `porous_blocks` require **0.4.0** or later with `porous_blocks_v1`. A `thermal` block requires **0.5.0** or later with `conjugate_heat_v1`. Surface `bodies`, `goals` and thermal `buoyancy` require **0.6.0** or later with `surface_bodies_v1`, `goals_v1` and `natural_convection_v1`. The domain solvers added in 0.7.0 (SPICE, TCAD, battery, antenna, 3-D FDTD, Nédélec cavity, PDN, power network, corrosion, multiphase, spray, elasticity, plasticity, fracture, rods, DEM, PIC, transport, bioheat, Joule heating, Cahn–Hilliard, DG, Vlasov, MD, photonics and THMC) require **0.7.0** with `domain_cases_v1`. Version 0.2.0 remains admissible for legacy cases and retained runs. Live execution and retained evidence share the domain producer validator. A later companion needs a deliberate compatibility update. The C++ named-reference dispatcher is not this protocol.

A run freezes the authored case, creates a private temporary task directory, and invokes `cfd-solve run CASE.json --out DIRECTORY` with separate arguments. It uses no shell. The worker retains stderr up to 16 KiB, cancels and reaps the child on request, and applies a 300-second runtime limit. File count and byte growth are checked during execution. Terminal admission limits are 128 files, 8 MiB per artifact and 16 MiB per run. The raw manifest may use up to the artifact limit, while its parsed compact form must satisfy the domain manifest limit. Authored case input is bounded to 2 MiB.

This companion emits its manifest after execution; it does not supply a live residual or progress stream. Job status is lifecycle status, not a convergence monitor.

Before returning a completed result, the host checks the manifest version, relative artifact paths, regular-file bounds, directory ownership and symbolic-link containment. The domain validator separately checks frozen-case identity, the authored subset of the runner's normalized case, finite metrics, file closure and the retained manifest bytes. The original input case and artifact text travel with the evidence; temporary files are removed after admission. Process success does not establish convergence or physical validation.

Cancellation remains pending until the worker acknowledges termination. Dropping a job requests cancellation without blocking the UI. The native application's shutdown path should cancel and briefly drain pending jobs before exiting, allowing the worker to reap its child.

Run the host tests with:

```sh
cargo test --profile ci -p solidcraft-multiphysics-host
```

Optional integration tests run the real companion. `real_reference_runner_all_templates` covers every shipped template and requires all their capabilities. `real_electronics_enclosure_retains_fluid_only_fields` checks the real obstacle-aware companion's final VTU cell count, solid exclusion, cell-valued velocity, patch report and relative boundary flux balance (below `1e-4` for this bounded demonstration). It verifies the recorded prescribed-flow isothermal case, not fan performance or heat transfer. `real_engine_run_keep_reopen_and_export` creates, executes and retains a case through engine commands, saves and reopens a `.sc` design, and checks the exported solution against `u=x(1-x)`:

```sh
SOLIDCRAFT_TEST_CFD_SOLVE=/path/to/cfd-solve cargo test --profile ci -p solidcraft-multiphysics-host -- --ignored
```

Set `SOLIDCRAFT_TEST_OUTPUT_DIR` to a fresh directory to retain the real-case outputs and the saved design for independent VTK/ParaView checks. This test option is not an application output setting.
