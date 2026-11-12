# Native Multiphysics companion host

This application-layer library runs the optional `cfd-solve` companion beside the SolidCraft window thread. It links no numerical solver into the CAD engine. The browser build returns no local process host.

Set `SOLIDCRAFT_CFD_SOLVE` to the native Rust port's executable path, or install `cfd-solve` on `PATH`. Host discovery does not start a process. An invalid explicit path produces an actionable job failure when a run starts.

The adapter accepts only the native Rust port **0.2.0**, upstream revision `7cceaaaddafaf2515e01a86406e13d6037fdbf9a`, and backend `native_rust_serial_cpu`. A later companion needs a deliberate compatibility update. The C++ named-reference dispatcher is not this protocol.

A run freezes the authored case, creates a private temporary task directory, and invokes `cfd-solve run CASE.json --out DIRECTORY` with separate arguments. It uses no shell. The worker retains stderr up to 16 KiB, cancels and reaps the child on request, and applies a 300-second runtime limit. File count and byte growth are checked during execution. Terminal admission limits are 128 files, a 1 MiB manifest, 8 MiB per artifact and 16 MiB per run.

Before returning a completed result, the host checks the manifest version, relative artifact paths, regular-file bounds, directory ownership and symbolic-link containment. The domain validator separately checks frozen-case identity, the authored subset of the runner's normalized case, finite metrics, file closure and the retained manifest bytes. The original input case and artifact text travel with the evidence; temporary files are removed after admission. Process success does not establish convergence or physical validation.

Cancellation remains pending until the worker acknowledges termination. Dropping a job requests cancellation without blocking the UI. The native application's shutdown path should cancel and briefly drain pending jobs before exiting, allowing the worker to reap its child.

Run the host tests with:

```sh
cargo test --profile ci -p solidcraft-multiphysics-host
```

Optional integration tests run the real companion. `real_reference_runner_all_templates` covers every shipped template. `real_engine_run_keep_reopen_and_export` creates, executes and retains a case through engine commands, saves and reopens a `.sc` design, and checks the exported solution against `u=x(1-x)`:

```sh
SOLIDCRAFT_TEST_CFD_SOLVE=/path/to/cfd-solve cargo test --profile ci -p solidcraft-multiphysics-host -- --ignored
```

Set `SOLIDCRAFT_TEST_OUTPUT_DIR` to a fresh directory to retain the real-case outputs and the saved design for independent VTK/ParaView checks. This test option is not an application output setting.
