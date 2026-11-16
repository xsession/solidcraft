# Multiphysics material library acceptance

Date: 2026-11-13. Status: material import/edit workflow verified; final document-switch regression rebuild pending.

The task exposes the existing internal physical-reference catalog from Multiphysics.
It creates design-owned copies through existing authoritative commands and reuses
the existing material editor. It does not change reference data, the document format,
CAD body assignments, a study's numerical case or any solver implementation.

Changed implementation paths:

- `crates/ui-egui/src/multiphysics_material_library.rs`
- `crates/ui-egui/src/multiphysics_workspace.rs`
- `crates/ui-egui/src/lib.rs`
- `assets/locales/pending/multiphysics-material-library.json`
- `assets/locales/pending/multiphysics-material-library-navigation.json`

The reference panel uses `render.library.list`, `render.library.describe` and
`render.library.physical.prepare`. Preparation is read-only. Its returned command
must be `material.create`, which is then executed explicitly. A visible source-review
checkbox is required before creation. Original values, conditions, sources, omitted
properties and native card receipts remain inspectable.

Design materials reuse `material_editor::rows` and normal `material.create/edit`
commands. Selection follows a material name. A different session UID discards
transient names, review acknowledgement, selection and unfinished material edits.

Independent read-only review found no concrete command, state-isolation or atomicity
issues. The new source passed formatting. The first native mouse/keyboard workflow
passed with 1 test and 0 failures, including import, editing, Undo/Redo, save/reopen,
duplicate refusal and camera/case/body-assignment isolation. Its four headless widget
captures were produced at 1440×900 and 1280×800; the reference view was inspected at
both sizes and the design view at 1440×900. Controls fit, and long details scroll.
The document-switch regression had its registry lookup corrected to use a point
handle instead of a rectangle. Its final rebuild and acceptance run remain pending.

Acceptance tests in the existing workspace test module:

- `material_library_mouse_import_edit_undo_and_reopen_preserve_study`
- `material_library_replacing_document_clears_pending_editor`

The first test searches C11000 by keyboard, selects the reference by mouse, checks
that review is required, imports one copy, checks original card data and one Undo
step, performs Undo/Redo, renames it in the existing editor, saves/reopens the design,
and refuses a duplicate atomically. It checks unchanged study identity, camera,
CAD selection and body assignments. It also captures both tabs at two desktop sizes.
The second test replaces the session while an unsaved material edit exists.

Build uses the existing low-memory native test settings: profile `ci`, optimization
and debug information disabled, dependency optimization disabled and 256 GUI
code-generation units. Test evidence is intended for `verification/materials/` in
the full source archive. Full repository CI, native OS-window capture and new solver
physics are outside this acceptance run.

For the design targets and official-source survey, see
[multiphysics setup architecture](research/multiphysics-setup-architecture.md).
