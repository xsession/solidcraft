# Commands

Every Multiphysics action is one of these commands. Each has a JSON schema for its input and
its output; `solidcraft-cli commands` lists them with their parameters, and the
[scripting tutorial](../tutorials/scripting.md) shows how to send one.

**Undo**: *one step* changes the design and is undone in one step; *view* changes only what
is shown (not saved); *none* reads or writes files only.

## Studies of the model

| Command | Does | Undo |
|---|---|---|
| `multiphysics.model.create` | A new study of the model: name, internal or external, bodies or a selection, margin, cell size | one step |
| `multiphysics.model.set` | Change its analysis, bodies or selection, margin, cell size, mesh level, cut cells | one step |
| `multiphysics.model.condition` | Inlet, outlet, wall or symmetry on a face of the model | one step |
| `multiphysics.model.fan` | A fan curve on a face at the domain's boundary | one step |
| `multiphysics.model.heat` | Heat transfer on or off, inflow temperature, fluid properties | one step |
| `multiphysics.model.material` | Conductivity, power and limit of bodies or a selection | one step |
| `multiphysics.model.heat_source` | Power on a face | one step |
| `multiphysics.model.thermal_face` | Temperature, heat flux, convection or adiabatic on a face | one step |
| `multiphysics.model.goal` | A goal on a face, a point, a body or the whole fluid | one step |
| `multiphysics.model.refine` | Finer cells around bodies or faces | one step |
| `multiphysics.model.adapt` | Finer cells where a run's field changes most | one step |
| `multiphysics.model.leave_out_small` | Keep only the study's parts of at least a size | one step |
| `multiphysics.model.build` | Build the case again from the design | one step |
| `multiphysics.model.status` | Is the case built from the design as it is now; the domain | none |
| `multiphysics.model.mesh_study` | Two coarser copies of a study, for a mesh study | one step |
| `multiphysics.model.leak` | The shortest path through the fluid between two faces, if any | none |

## Parts

| Command | Does | Undo |
|---|---|---|
| `multiphysics.parts.list` | Every part with its tags and attributes, every selection with its parts now | none |
| `multiphysics.parts.tag` | Add or remove tags, set attributes | one step |
| `multiphysics.parts.group` | Tag parts by the kind their reference designators name | one step |
| `multiphysics.parts.selection` | Add, replace or remove a saved or query selection | one step |
| `multiphysics.parts.select` | Set the design's selection to bodies or a selection | view |
| `multiphysics.parts.similar` | Parts alike in exact volume and box; optionally select them | view |
| `multiphysics.parts.find` | Parts meeting rules on setup and results; optionally select them | view |
| `multiphysics.parts.ghost` | Draw parts see-through, or solid again | view |
| `multiphysics.view.save` | Keep a result view by name: camera, colour-by, ghosted parts, cut plot, colour range | one step |
| `multiphysics.view.apply` | Put a saved view back: camera, ghosting and colours; the plot settings for the host | view |
| `multiphysics.view.list` | A study's saved result views | none |
| `multiphysics.view.remove` | Drop a saved view | one step |
| `multiphysics.parts.colour` | Colour parts by material, power, tag or verdict, with a legend | view |
| `multiphysics.parts.check` | The checklist of a study | none |
| `multiphysics.parts.table` | The component table, optionally as CSV | none |

## Runs and results

| Command | Does | Undo |
|---|---|---|
| `multiphysics.preflight` | The structural check of a study, and whether a solver is available | none |
| `multiphysics.run.start` | Run a study's committed case | none |
| `multiphysics.run.status` | Runs in flight and finished in this session | none |
| `multiphysics.run.cancel` | Cancel a run | none |
| `multiphysics.run.keep` | Keep a finished run in the design | one step |
| `multiphysics.result.query` | Probe, min/max, line, integral, difference, goals, components, trajectories, isosurface, surface | none |
| `multiphysics.result.inspect` | What a kept VTU or PVD file holds | none |
| `multiphysics.result.convergence` | Grid convergence index of each part over three grids | none |
| `multiphysics.result.compare` | Two runs, of one study or two, part by part: temperatures, change, verdicts | none |
| `multiphysics.result.report` | Batch report over runs: CSV, HTML, pictures | none |
| `multiphysics.result.export` | A run's original files and case, for ParaView, into a new folder | none |
| `multiphysics.case.export` | A study's case as a new JSON file | none |

## Studies in general

| Command | Does | Undo |
|---|---|---|
| `multiphysics.catalog` | Profiles, templates, and whether a solver is available | none |
| `multiphysics.study.list` | Every study with its revision and runs | none |
| `multiphysics.study.create` | A study from a template | one step |
| `multiphysics.study.update` | Apply a study's edited case | one step |
| `multiphysics.study.clone` | Copy a study's committed case and setup (not its runs) | one step |
| `multiphysics.study.remove` | Remove a study | one step |
| `multiphysics.study.compare` | Each study's latest run's metrics side by side | none |
| `multiphysics.study.sweep` | One checked study per value (1 to 16) of a case setting; nothing runs | one step |
| `multiphysics.wizard.enclosure` | A study of an electronics enclosure from its size, cooling and components; nothing runs | one step |
| `multiphysics.cad.bodies` | The design's bodies, for bringing into a case | none |
| `multiphysics.case.component.import` | Design bodies made into parts of a case, returned as a draft to apply | none |
| `multiphysics.case.feature.prepare` | A case feature (fan, porous block, material, goal, device) added to a draft case | none |
| `multiphysics.surface.inspect` | Connectivity report of a triangle surface | none |
