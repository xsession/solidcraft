# 6. Scripting a study

Every button in these tutorials runs a command, and a command can be sent without the
buttons: from a script over SolidCraft's control channel, from the command line, or from an
AI agent over MCP. A command called this way never opens a dialog.

## The control channel

Start SolidCraft with a control port:

```sh
SOLIDCRAFT_CFD_SOLVE=/path/to/cfd-solve solidcraft --control 7979
```

It listens on `127.0.0.1:7979` for JSON, one request per line. The first request carries the
launch token from `~/.local/share/solidcraft/control/7979.token`:

```json
{"id": 1, "method": "engine.execute", "token": "…",
 "params": {"command": "multiphysics.parts.group", "params": {}}}
```

The reply has `ok` and `result` (or `error`). Never paste the token anywhere public.

## The example, as a script

The example design of these tutorials is built and run by
`docs/multiphysics-site/tools/make_study.py`. Its core, in order:

```python
run("solid.box", {"length": 160, "width": 60, "height": 40, "corner": [0, 0, 0],
                  "operation": "new", "body_name": "Duct"})
run("solid.box", {"length": 170, "width": 52, "height": 32, "center": [80, 30, 20],
                  "operation": "cut"})
# … the board and its parts …
run("multiphysics.parts.group", {})
study = run("multiphysics.model.create",
            {"name": "Board in duct", "analysis": "internal", "cell_mm": 4.0})["study"]
run("multiphysics.model.condition", {"id": sid, "expected_revision": rev(), "name": "Inlet",
    "face": {"body": "Duct", "point": [0, 30, 2]},
    "role": {"kind": "inlet", "speed_m_per_s": 1.0}})
run("multiphysics.model.heat", {"id": sid, "expected_revision": rev(), "enabled": True,
    "inflow_temperature_k": 298.15})
run("multiphysics.model.material", {"id": sid, "expected_revision": rev(), "bodies": ["U1"],
    "conductivity_w_per_m_k": 150.0, "power_w": 2.0, "max_temperature_k": 358.15})
run("multiphysics.model.leave_out_small", {"id": sid, "expected_revision": rev(), "below_mm": 4.0})
job = run("multiphysics.run.start", {"id": sid})
# poll multiphysics.run.status until the job is "succeeded", then:
run("multiphysics.run.keep", {"job_id": job["id"]})
run("file.save_as", {"path": "board-in-duct.sc"})
```

!!! note "Revisions"
    Every edit of a study names the revision it was based on (`expected_revision`). An edit
    based on an older revision is refused, so a script and a person editing the same study
    never overwrite each other's changes silently. Read the revision again
    (`multiphysics.study.list`) after each edit, as `rev()` does.

## From an AI agent

`solidcraft-cli mcp` serves the same commands to an AI agent over the Model Context Protocol,
with the same rules. The agent can build the study, run it, read the component table and
write the batch report; the commands are listed in [Commands](../reference/commands.md).

Next: [climate, humidity and dust](climate-and-dust.md).
