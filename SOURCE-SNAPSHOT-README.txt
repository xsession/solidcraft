SolidCraft Multiphysics complete source snapshot — 2026-11-13

The solidcraft/ folder contains the complete repository source tree, including
assets, test fixtures, the Multiphysics integration and the separate MIT-licensed
companion/cfd-solvers-rs 0.6.0 project with vendored dependencies. Git history and
build caches are excluded. This is a full source archive, not a patch archive.

From solidcraft/ run:
  Linux:   sh scripts/run-multiphysics.sh
  Windows: scripts\run-multiphysics.bat
These build the companion, set SOLIDCRAFT_CFD_SOLVE, then build/run SolidCraft.
Install each project's required Rust toolchain and OS native linker/GUI build
prerequisites first. SolidCraft requires Rust1.95 or newer. The companion retains
its own toolchain declaration and is excluded from the outer Cargo workspace.
Windows execution of the launcher has not been tested here; shell syntax passed.

In the app select Multiphysics and flow-electronics-enclosure, Run, Keep result,
then select a VTU artifact in Results and Inspect. Fixed-axis HEX8 slices show
original scalar cell values or vector magnitudes. The example uses prescribed
isothermal flow and coarse grid-aligned PCB/component blocks. General CAD fluid
extraction, fan curves and conjugate heat transfer remain planned.

Fresh verification and actual VTU/PVD/mesh-report output are in verification/.

Fresh tests: 95 solver, 30 domain passed. UI/host type check: passed.
Headless GUI rendering and pointer/keyboard controls now pass at1440x900 and1280x800. Native OS-window capture is not claimed. See docs/multiphysics-gui-verification.md and verification/gui/report.json.
Ready-to-open study: examples/electronics-enclosure-airflow.sc (choose Multiphysics -> Results -> VTU -> Inspect).

Update (fan curve, porous zones, setup records):
- Companion 0.4.0 adds fan_curve_v1 (the operating point follows the measured system
  resistance; curves are never extrapolated) and porous_blocks_v1 (implicit
  Darcy-Forchheimer zones).
- New template: flow-electronics-enclosure-fan. It has a shrouded inlet fan and an
  exhaust filter. The curve and coefficients are illustrative; the flow is isothermal laminar.
- crates/multiphysics adds airflow admission with per-capability producer gating, and
  setup records: named CAD selections, typed boundaries, contacts, orthotropic conduction
  and a coverage audit. The audit blocks unresolved references at Apply.
- Fresh results: 100 solver tests and 34 contract tests passed. The contract tests ran in
  an offline harness; see verification/fan-porous/README.md.
- The prescribed enclosure field is byte-identical to the earlier verification.
- The full SolidCraft workspace and its UI were NOT rebuilt in this increment, because
  crates.io was unavailable. Run cargo xtask ci before merging.

Update (conjugate heat transfer and offline verification, all offline):
- Companion 0.5.0 adds conjugate_heat_v1. A steady energy solve covers fluid and solid
  cells: orthotropic solids with power, convection walls, an energy balance, and a
  pass/warn/fail component table.
- cfd-solve verify --out DIR runs six embedded self-checks with no network.
- New template: flow-electronics-enclosure-cht.
- The contract admits thermal blocks; producer 0.5.0 is required for them.
- Fresh results: 105 solver tests, 6/6 verify checks and 35 contract tests passed. See
  verification/cht-verify/README.md. The full workspace and UI were not rebuilt here.

Update (natural convection, goals, surface bodies, enclosure wizard, all offline):
- Companion 0.6.0 adds:
  - natural_convection_v1: two-way Boussinesq coupling, checked against de Vahl Davis;
  - goals_v1: monitored goals with early stop;
  - surface_bodies_v1: closed triangle surfaces become solids.
- crates/multiphysics adds wizard::enclosure_case, which builds a complete case from
  mm/°C/W inputs and reports every automatic choice.
- New template: flow-electronics-enclosure-natural.
- Fresh results: 109 solver tests, 7/7 verify checks and 37 contract tests passed, plus the
  wizard to companion to admission end-to-end runs. See
  verification/natural-goals-wizard/README.md. The full workspace and UI were not rebuilt here.
