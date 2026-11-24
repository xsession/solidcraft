# Architecture

Multiphysics is split so that the parts that decide something can be tested without a
window, and the solver can fail without taking SolidCraft with it.

| Piece | Where | What it does |
|---|---|---|
| **Portable model** | `crates/multiphysics` | Studies, cases, the part registry, the component table, Find, grouping, batch reports, result queries. Plain Rust with serialization and hashing only: no geometry, no files, no threads. It builds for the browser too |
| **Commands** | `crates/engine/src/cmd/multiphysics*.rs` | Every action as a command with typed input and output schemas: `multiphysics.model.*`, `multiphysics.parts.*`, `multiphysics.run.*`, `multiphysics.result.*`. Measures parts from exact geometry and hands the portable model what it needs |
| **Jobs** | `crates/engine/src/multiphysics_jobs` | Runs in flight: at most two at once, each in a fresh private folder with bounded output and a five-minute limit, eight kept in the session's history |
| **Process host** | `apps/multiphysics-host` | Starts the companion as a separate process on native hosts |
| **Files** | `crates/io/src/multiphysics.rs` | Case export, result hand-offs for ParaView, batch reports: always into new files or folders, written atomically |
| **Workspace** | `crates/multiphysics-ui` | The tree, panels and plots. Its own crate: editing it rebuilds only it |
| **Viewport** | `crates/ui-egui` | The 3-D view: part colours, ghosts and outlines |
| **Solver companion** | `companion/cfd-solvers-rs` | `cfd-solve`, a separate program and project |

## Everything is a command

A button runs a command, a script runs the same command, and an AI agent over MCP runs it too.
A command called by a program never opens a dialog. Each command has a JSON schema for its
input and its output, so a caller can check its request before sending it.

Commands that change the design (a study's settings, tags, selections) are one undo step
each. View state (colour by, ghosts, the selection) is not saved and is not an undo step.

## Revisions

A study has a revision that goes up with every change. Edits name the revision they were
based on (`expected_revision`), and an edit based on an older one is refused. A person and a
script working on the same study therefore never overwrite each other silently.

## What is kept in the design

The design file holds each study's settings and its **kept runs**: the run's input case and
its SHA-256, the study revision it was made from, the producer's manifest and its output
files (VTU fields, reports, goal histories). A project is bounded to 64 MiB, and a run's output
to 16 MiB: this is a tool for studies of a design, not an archive of large simulations.

## Boundaries

- The CAD kernel is unchanged by Multiphysics; studies read its exact geometry and its display
  meshes.
- Decisions about the model (which part a click means, how big a part is, whether two parts
  are alike) use exact geometry; meshes are for drawing and for the solver's case.
- The companion is trusted only as far as its output is checked: a run whose output does not
  match its frozen input, or arrives for a study that has changed, cannot become a current
  result.
