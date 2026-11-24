# Glossary

**Analysis (internal, external)**
: Where the fluid is: inside the study's bodies (a duct, an enclosure) or around them, in their
  box grown by a margin.

**Case**
: The JSON the solver reads: grid, solids, surface bodies, conditions, heat, goals, output
  settings. Built from a study of the model; frozen with its SHA-256 when a run starts.

**Cell**
: One box of the grid. Values in field files are mostly per cell.

**Companion**
: `cfd-solve`, the separate solver program SolidCraft starts for each run.

**Component table**
: One row per part of a study: its setup joined to a run's results, with margin and verdict.

**Cut cells**
: An option that keeps the fluid part of cells a surface crosses, instead of making each cell
  wholly solid or fluid.

**Designator (reference designator)**
: The name a board gives a part, such as U12 or R3; its letters say what kind of part it is.

**Domain**
: The box the grid fills: the study's bodies' bounding box, plus the margin for an external
  study.

**Evidence (kept run)**
: A run stored in the design with its input case and hash, the study revision, the manifest and
  the output files.

**Ghost**
: A part drawn see-through, with faint edges and an outline, that clicks pass through.

**Historical run**
: A kept run made from an earlier revision of its study.

**Kinematic pressure**
: Pressure divided by density, in m²/s², as the solver reports it.

**Margin**
: A part's temperature limit less its highest temperature in a run, in K.

**Part**
: A body of the design, by its name.

**Part selection**
: A named set of parts, saved or worked out from rules (a query).

**Revision**
: A study's change counter; edits name the revision they were based on.

**Staircase**
: A body voxelised at cell centres: its faces become steps of cells.

**Study of the model**
: A study whose geometry is the open design's bodies, linked, not copied.

**Verdict**
: pass, warn, fail or no limit: a part's highest temperature judged against its limit.
