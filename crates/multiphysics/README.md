# Multiphysics authoring contracts

This independent, safe Rust crate stores studies and portable evidence for the
Multiphysics workspace. It does not invoke a solver, compute geometry, or change
the existing Simulation workspace.

The 23 case templates in `data/templates/` come from the MIT-licensed
`cfd-solvers-rs` 0.2.0 source package. Their original license is retained beside
them. They select 18 bounded numerical case types, not all upstream solver
methods. `data/solver-catalog.json` is the original MIT workspace design v1.4
inventory of 34 broad profile families; `data/DESIGN-LICENSE` retains its notice.
Catalogue presence does not establish an installed or qualified execution path.

Study names, questions and case objects are authored document state. Incomplete
cases remain savable; `preflight` reports bounded structural checks separately
from the companion's detailed schema, numerical and physical validation.
Updates use an expected revision and commit atomically. No-op Apply preserves
the revision. Removal never reuses a study identifier.

Case limits are 2 MiB, depth 32, 200,000 JSON nodes, 32,768 generated mesh
cells/nodes, 20,000 requested steps/iterations and 50 million mesh-step work
units. These bounds protect the authoring interface; the execution host must
also limit runtime, output, memory and detailed method-specific work.

Runs retain the raw submitted case, its hash, study revision, native manifest and
bounded text artifacts. Each run has at most 128 flat ordinary filenames and
16 MiB of artifact text; no individual artifact exceeds 8 MiB. A project holds
at most 64 studies and 12 runs per study. The complete compact serialized project,
including original cases, normalized manifests and escaped artifact text, is
limited to 64 MiB. Every mutation validates that final budget before committing.
These are deliberate portable example limits, not large-dataset post-processing
support. The host validates fresh directory ownership and symlinks before admission.

Admission understands exactly `cfd-solvers-rs/run/1`, producer version 0.2.0,
the pinned upstream revision and `native_rust_serial_cpu`. The normalized
manifest case must preserve every explicitly authored key/value; additional
solver defaults are permitted. Declared D2Q9 single-precision scalars use the
exact f32-to-f64 encoding of the pinned runner. This does not add a numerical
tolerance to other settings. The host independently binds the original
submitted case bytes to the job. This is structural consistency, not a
cryptographic signature or independent proof that a solver executed.

Editing keeps earlier immutable runs as historical evidence. Attaching a new
run requires the current case hash and revision. Execution completion, linear
convergence, steady convergence and physical validation remain distinct. A
manifest cannot establish physical validity, convergence quality, field units,
or a CAD-to-volume-mesh connection.
