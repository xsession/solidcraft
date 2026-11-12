# Multiphysics workspace

Multiphysics adds saved numerical studies, actual companion execution and portable results under a separate workspace and command namespace.
The planned Simulation entry is unchanged. Existing Design and Render defaults are unchanged.

## Author and run

Choose **Multiphysics** in the workspace menu. Choose a template and create a study.
Write its engineering question, edit its case fields or JSON, and choose Apply.
Each study retains its own unapplied draft, including invalid JSON, across study and workspace switches.
Cancel discards the selected draft. Apply, Remove and Keep result use document Undo and Redo.
Save the `.sc` design to retain applied studies and retained results.
Unapplied drafts remain only in the open document session; replacing or closing that document discards them.

Prepare, Physics, Mesh, Couple, Run, Results and Verify are contextual views of the same study.
Mesh exposes actual case discretization settings. Couple explains the unavailable general coupling adapter.
The structural check covers bounded inputs; the numerical companion checks its full case schema and prerequisites at execution.
An empty structural error list is not numerical or physical qualification.

The catalogue contains 34 design profiles. The 23 templates cover 18 bounded case types in the Rust port.
Catalogue presence is separate from executable support. Cases use generated or explicit numerical representations.
The current CAD model is not automatically converted into these cases.

## Native numerical companion

Use `cfd-solve` built from **cfd-solvers-rs 0.2.0**, pinned to upstream revision
`7cceaaaddafaf2515e01a86406e13d6037fdbf9a`.
The companion is a separate process, not a linked dependency or binary shipped by this change.
Another producer version is refused.

```sh
export SOLIDCRAFT_CFD_SOLVE=/absolute/path/to/cfd-solve
cargo run --profile ci -p solidcraft -- --sample
```

```powershell
$env:SOLIDCRAFT_CFD_SOLVE = 'C:\tools\cfd-solve.exe'
cargo run --profile ci -p solidcraft -- --sample
```

An executable named `cfd-solve` on PATH is also discovered. Without it, authoring and case exports remain usable.
Browser builds prepare cases but cannot launch a local process or export a native directory.
Native CLI/MCP sessions inject the configured process host too.

Run freezes committed inputs, their SHA-256 and the study revision. It never silently applies a draft.
At most two workers run concurrently. Each uses a fresh private directory, bounded output and a five-minute timeout.
Cancel remains pending until the worker stops. Closing or replacing a document retires its work.
The native app waits briefly for cancellation acknowledgement at exit.
The transient history retains at most eight jobs; starting a later job can retire the oldest finished entry.
Keep important results in the design before they leave that session history.

## Post-processing and exports

Results displays actual metrics, manifests, artifacts, CSV plots and tables.
Keep result saves admitted current evidence in the design. Later study edits make that evidence historical.
Late or mismatched completions cannot become current results.
Every retained run keeps its original submitted case and hash independently of later edits.
Execution completion, linear convergence, steady convergence and physical validation remain distinct.

Export a case to a new JSON filename or a retained run to a new directory.
The handoff preserves original VTU/PVD/CSV/JSON/Touchstone artifacts, the producer manifest,
`case.json`, a metadata file and a short ParaView guide. Current-study metadata is labelled separately from the captured case.
Open `.pvd` collections or `.vtu` files in ParaView and save a state there.
The exporter does not fabricate `.pvsm` or a spatial field from a CSV plot.

Producer point/cell associations and values remain unchanged. Flow pressure is kinematic pressure in m²/s².
SIMPLE series use iteration, not physical seconds. Maxwell fields may use separate staggered grids.
Generic potential fields have no inferred physical units. Inspect conventions before comparisons.

## Command interface

Every command uses `Session::execute` and typed input/output schemas. Programmatic calls never open dialogs.

| Command | Input | Result |
|---|---|---|
| `multiphysics.catalog` | `{}` | Catalogue, templates, host availability |
| `multiphysics.study.list` | `{}` | Project study/run summaries and case hashes, without case payloads |
| `multiphysics.study.create` | `template`, `name` | Created study authoring state and run summaries |
| `multiphysics.study.update` | `id`, `expected_revision`, `name`, `question`, `case` | Saved study authoring state and run summaries |
| `multiphysics.study.remove` | `id` | Removed identity |
| `multiphysics.preflight` | `id` | Structural issues, separate host status |
| `multiphysics.run.start` | `id` | Asynchronous job |
| `multiphysics.run.status` | `{}` | Recent statuses, without large dataset bytes |
| `multiphysics.run.cancel` | `job_id` | Pending or acknowledged cancellation |
| `multiphysics.run.keep` | `job_id` | Study authoring state and retained-run summaries |
| `multiphysics.case.export` | `id`, `path` | New case file |
| `multiphysics.result.export` | `id`, `run_id`, `directory` | Original result handoff |

Exports respect existing agent path scopes. Scratch computations cannot launch external numerical work.
Projects are bounded to 64 MiB of serialized authoring and evidence; this is deliberately a small-case adapter.
Study command replies include the selected case for editing, but only run identity, stale state,
artifact counts and convergence metadata. Project listings omit numerical cases as well.
Original manifests and artifact text remain in the saved design and native result view;
use case export or result export when their complete bytes are needed.

## Review boundaries

The portable `crates/multiphysics` contract crate sits at layer 1 and uses external serialization and hashing only.
The document owns `multiphysics` plus one pending format note. Engine commands and transient jobs live in their own modules.
IO writes handoffs through the safe writer. The UI module owns workspace presentation and drafts.
`apps/multiphysics-host` owns native processes and stays outside the portable library layers.
No geometry-kernel or numerical solver-repository source is changed.

This is the first native implementation of the v1.4 design. Remaining work includes arbitrary CAD volume meshing,
persistent CAD boundary assignments, coupling-graph execution, scientific spatial field rendering,
the full browser prototype filter pipeline and unported solver/GPU/MPI methods.
Specialist RF, EMC, LED, detector and fan workflows retain the limits of actual connected case adapters.
