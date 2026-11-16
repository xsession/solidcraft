# Multiphysics GUI verification

Verified on 2026-11-13 with Rust 1.100.0 on Linux.

The Results workspace uses compact slice controls and a contextual run summary.
Raw artifact text is available in a closed disclosure. Results scroll independently
and retain original cell values. Viewing retained output does not require an
installed solver; the status bar distinguishes viewing from starting a new run.

## Evidence and scope

The supplied archive contains verification/gui/ui images at 1440×900 and 1280×800,
logs, a machine-readable report, and examples/electronics-enclosure-airflow.sc.
The images are headless egui raster captures of the actual SolidCraft widgets,
using the real retained VTU data. They are not native OS-window captures or a
recording. The test drives pointer/keyboard events through the existing harness.

The verified sequence selects XZ, edits the layer index, changes velocity to
kinematic pressure and back, visits YZ, and restores XY. It checks that solver
inputs are unchanged, the selected frame coordinate matches retained PVD data,
and the plot, data area and legend fit the central viewport at both resolutions.
The data area is at least 350×200 pixels in this fixture. Source arrays are not
interpolated. The selected field is velocity_m_per_s, not an arbitrary vector.

Host checks cover all 24 templates, the engine run/keep/save/reopen/export
workflow, and a separate enclosure save/reopen/export with exact retained
field text. The fresh enclosure field matches the field used by the UI captures
byte for byte (SHA-256 6fb29f921f03f00eea8f231ad1e76915cc25d22fb76a27173fffbd053f426d8c).

The example remains a grid-aligned primitive, prescribed-flow, isothermal
laminar case. It does not establish mesh independence or physical validation.
General CAD fluid extraction, fan curves and conjugate heat transfer remain planned.

## Reproduce the UI check

From the solidcraft directory in the complete source archive, on Linux:

```sh
export SOLIDCRAFT_ENCLOSURE_RESULT="$(pwd)/../verification/enclosure"
export SOLIDCRAFT_MULTIPHYSICS_EVIDENCE="$(pwd)/../verification/gui/ui"
cargo test -p solidcraft-ui-egui -p solidcraft-multiphysics-host --lib --profile ci \
  --config 'profile.ci.opt-level=0' --config 'profile.ci.debug=0' \
  --config 'profile.ci.package."*".opt-level=0' \
  --config 'profile.ci.package.solidcraft-ui-egui.codegen-units=256' \
  real_enclosure_result_renders_spatial_slices -- --ignored --nocapture --test-threads=1
```

The low-memory compile settings do not change the case model. Build and configure
the companion before invoking ignored host tests, as described in the host README.
Open the saved .sc project, choose Multiphysics, select its study and Results,
then choose a VTU artifact and Inspect. The project contains a numerical study;
it is not a CAD enclosure assembly or an automatic CAD-to-mesh model.


## Orbitable 3D slice acceptance

The latest evidence is in `verification/field3d/report.json`, its logs and `ui/` images.
The prior 2D-only and numerical host evidence remains under `verification/gui/`.
The latest pass used the same retained final enclosure field, with the SHA-256 above.
No numerical solver rerun was needed for this presentation change.

The existing real-data test now selects **3D slice**, orbits using Fusion-style
middle-button plus Shift input, pans with the middle button and zooms with wheel input.
It enables vector glyphs and checks that 1–128 arrows are actually rendered.
**Fit result** restores the domain centre and preserves the chosen orientation.
The test runs with all three Design angle locks enabled. It asserts that the Design
camera, document identity, solver case, angle locks and navigation preferences remain unchanged.
Separate regressions exercise preservation of a pending Design fit and a Render entry-fit token.

Fresh results: 32 portable core tests, 12 workspace tests, the real enclosure input
check, the Render-fit regression and 13 host unit tests passed. One core fixture,
one workspace external-data fixture in the normal run, and four host external-run
checks were skipped; the enclosure UI fixture was then executed explicitly and passed.
Changed Rust files passed the formatting check. The full repository CI was not run.

The 3D captures were inspected at 1440×900 and 1280×800. The 2D capture at 1280×800
was also inspected. Canvas and legend containment passed at both sizes.
These remain headless raster captures of the actual widgets, not native OS-window captures.
The displayed object is a single original cell layer in a domain wireframe.
Vector arrows are samples, not streamlines. This does not add volume rendering,
CAD body overlays or thermal coupling. Use **3D slice**, **Vector glyphs** and
**Fit result** in Results after inspecting the saved example's VTU artifact.
