# Theory of operation

These pages explain what happens between the design and the numbers, so that you can judge
what a result means and what it does not.

```text
 the design's bodies        part registry (tags, selections)
          │                           │
          ▼                           ▼
   study of the model ──build──▶ case (JSON: grid, solids, conditions, heat, goals)
          │                           │  frozen with its SHA-256 at Run
          │                           ▼
          │                  solver companion (cfd-solve, its own process)
          │                           │  VTU fields, thermal report, goals, manifest
          ▼                           ▼
   run evidence kept in the design ◀──┘
          │
          ├─▶ component table (setup joined to results per part, margins, verdicts)
          ├─▶ colour by, Find, batch report
          └─▶ cut and surface plots on the model
```

| Page | Question it answers |
|---|---|
| [Architecture](architecture.md) | Which part of SolidCraft does what, and where the boundaries are |
| [From design to case](design-to-case.md) | How bodies become a grid, faces become conditions, and when the case is built again |
| [The solver](solver.md) | What equations are solved, and with which assumptions |
| [Parts and selections](parts.md) | How parts are named, tagged, measured and selected |
| [Results and verdicts](results.md) | How results are kept, matched back to parts and judged |
| [Drawing large assemblies](drawing.md) | How colours, ghosts and outlines stay fast with thousands of parts |
| [What it does not do](limits.md) | The limits to keep in mind |
