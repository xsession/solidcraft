# Optional numerical companion

`cfd-solvers-rs` is a separate MIT-licensed Cargo project, excluded from the
SolidCraft workspace. It is launched as a child process, not linked into CAD.
Its vendored dependencies and licence files are included for reproducible builds.

From the repository root, run `sh scripts/run-multiphysics.sh` on Linux, or
`scripts\run-multiphysics.bat` on Windows. These build the companion, set
`SOLIDCRAFT_CFD_SOLVE` to its explicit executable path, and build/run SolidCraft.
The normal Rust toolchain and operating-system linker/GUI development requirements
of both projects still apply. SolidCraft requires Rust 1.95 or newer; the companion
retains its own toolchain declaration. These scripts do not install system tools.

In the application select Multiphysics and the `flow-electronics-enclosure`
template. Run, Keep result, then select a VTU artifact in Results and Inspect.
That example uses discretized grid-aligned obstacles and prescribed isothermal
airflow. The `flow-electronics-enclosure-fan` template replaces the prescribed
inlet with a fan pressure-flow curve and adds an exhaust filter as a porous zone.
The solver finds the operating point from the measured system resistance.
General CAD fluid extraction and coupled thermal predictions are not implemented.
`flow-electronics-enclosure-cht` adds conjugate heat transfer: powered, conducting
component blocks (orthotropic PCB) cooled by the computed airflow, with an energy
balance and a pass/warn/fail component table. Run `cfd-solve verify --out DIR` for the
offline self-check suite. `flow-electronics-enclosure-natural` is a sealed box cooled by natural convection. The
0.6.0 companion declares `grid_obstacles_v1`, `fan_curve_v1`, `porous_blocks_v1`,
`conjugate_heat_v1`, `surface_bodies_v1`, `goals_v1` and `natural_convection_v1`.
Nothing in it uses the network. Older companions cannot silently admit a case that needs a
capability they lack.
