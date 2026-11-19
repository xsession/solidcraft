# Multiphysics reference menus — verification

Date: 18 November 2026. Base: `50145b4482f531e2b4e99f10deb0c449ac8ef2dc`, SolidCraft 0.44.0.
Environment: Linux, pinned Rust 1.99.1, native egui 0.37.0. The implementation is local source work against the retrieved repository revision.

## Delivered behavior

- Seven grouped menus reuse the existing dispatcher, shared native controls and SolidCraft theme. Unavailable device/geometry/post-processing operations explain their missing adapter.
- The seven-page wizard preserves a transient numerical candidate across Back/Next. Finish creates one template-compatible study through one undoable engine transaction.
- Typed, read-only feature preparation supports an external singleton fan, bounded primitive-grid porous regions, thermal records for existing numerical bodies and outlet volume-flow goals. Drafts require review and explicit Apply.
- Replacing existing prescribed conditions on the fan patch requires explicit acceptance; other patch conditions survive.
- Scalar array members now have editable native controls. Structured fields expose stable RFC 6901 scenario targets.
- Project clone copies the committed setup through its existing command. Tree right-click menus now use the shared menu container.
- One existing VTK tuple conversion was updated to the pinned compiler's safe array-chunk API so strict core clippy passes.

## Completed checks

| Check | Result |
|---|---|
| Core multiphysics library tests | 59 passed; 1 real-VTU fixture test ignored |
| UI-text library tests | 15 passed |
| Native UI and host test build | Passed with low-memory CI profile overrides; existing deprecation/unused-import warnings remain |
| Multiphysics workspace GUI suite | 16 passed; 1 real-enclosure fixture test ignored |
| Native host suite | 14 passed; 4 companion-dependent tests ignored |
| Engine multiphysics job/command suite | Focused transaction check: 1 passed; 0 failed; 0 ignored; full unit-test binary blocked by SIGKILL during compilation |
| Strict core clippy | Passed with `-D warnings` |
| Layers / mesh-decisions / assets / UI-containers | Passed; UI-container gate reports unrelated existing raw-tooltip warnings |
| Formatting and diff whitespace | Passed |
| `xtask lists` | Fails in 13 untouched root locale catalogs; no new pending-key collision |

The root catalogs reported by `lists` are ar, de, en, es, fr, he, ja, ko, pl, uk, vi, zh-TW and zh. They have no working-tree changes. This delivery does not claim a passing full `xtask ci`, wasm build, full application test suite or complete commercial-simulator parity.

## Mouse-driven acceptance observed

`reference_menus_insert_apply_and_clone_without_running`: actual pointer/keyboard events open Insert, exercise an inert unsupported TEC item, create a porous draft, edit its numeric max-bound array, Apply once, and clone the committed study. Preparation preserves document identity/history; Apply adds one undo step, clone adds one, no job starts, and the design camera is retained.

`reference_seven_step_wizard_preserves_edits_and_creates_atomically`: actual pointer/keyboard events open Project → wizard, set a name, choose Flow / channel, edit the method to SIMPLE, navigate Back/Next, and Finish. The final study retains the edited method and name; one undo removes the new study. Until Finish, the authoritative document/history are unchanged.

The four new reference-workflow captures were inspected at 1440×900 and 1280×800. Menus and wizard footer controls are visible. At 1280×800, the tall feature sheet floats over some ribbon controls to keep its footer on screen. Captures come from real headless egui rendering, not an OS window or reconstructed mock-up. They show setup with no installed numerical worker in that test session, not solved airflow.

## Limits and remaining work

Feature values are original illustrative starters, not measured vendor properties. Bounded authoring validation is separate from the companion's full schema and physical checks. The wizard uses explicit SI conventions without implicit unit conversion, and edits only fields exposed by its template.

General exact CAD boundary picking and fluid extraction/lids, local mesh refinement, internal fans, coupled TEC/heat-pipe/PCB/two-resistor/LED/battery devices, contacts, perforations, rotating-region adapters, and several advanced result visualisations remain unavailable in this authoring path. Existing independent battery/optical templates do not establish coupled enclosure device support. Radiation/electrical tree entries similarly do not introduce those models to the primitive flow solver.

See `multiphysics-reference-menu-map.md` for the functional mapping and `multiphysics-reference-inputs.json` for hashed input provenance. Uploaded commercial images, icons, material tables and private model identifiers are excluded from the source ZIP.

## Reproduction

The core and locale tests use `cargo test -p solidcraft-multiphysics --lib` and `cargo test -p solidcraft-ui-text --lib`.

Native acceptance was built with `CARGO_BUILD_JOBS=1`, `--profile ci`, `profile.ci.opt-level=0`, `profile.ci.debug=0`, `profile.ci.package."*".opt-level=0`, and `profile.ci.package.solidcraft-ui-egui.codegen-units=256`. The resulting UI test executable ran with filter `multiphysics_workspace::tests::` and `--test-threads=1`; the host executable ran its full library suite. Set `SOLIDCRAFT_MULTIPHYSICS_EVIDENCE` to an output directory to save headless captures. The engine suite uses filter `multiphysics_jobs::tests::`.

The separate engine unit-test binary was killed with SIGKILL during compilation twice, including a retry with 256 code-generation units. The focused transaction test was then compiled separately, unchanged from its existing repository test file, and linked to the real engine/IO libraries already compiled for the UI (bundled-render-library and profile features). That focused check passed and covers read-only preparation, invalid-case atomicity/ID preservation, undo/redo and design-container save/reopen. The full engine unit suite did not execute. The focused harness source and log are included in the verification bundle.
