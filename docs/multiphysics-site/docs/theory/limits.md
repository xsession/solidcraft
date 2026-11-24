# What it does not do

Know these before trusting a number.

## Physics

- **Laminar flow only.** No turbulence model is applied. Forced convection at the speeds of
  small fans in small ducts is often laminar or transitional; fast flows over large boards
  are not.
- **One-way heat.** The temperature is solved on the settled flow and does not change it:
  no buoyancy, so no natural convection in a study of the model.
- **No radiation.** Hot parts in still air lose a fair share of their heat by radiation; a
  study without it reads high there.
- **Constant properties.** Air and solids keep one density, specific heat and conductivity.

## Geometry

- **Staircase bodies.** Bodies are voxelised at cell centres. Thin walls, slanted faces and
  small gaps are only as good as the cells are small. Check the build's note on how much of
  the bodies' volume the grid holds.
- **Lost parts.** A part thinner than a cell may not be in the grid at all. The checklist says
  which; leave them out or refine around them.
- **Small cases.** At most 32 768 cells, 16 MiB of output per run and 64 MiB per project. The
  example's 4 mm cells are coarse: its temperatures are an illustration of the workflow, not
  a converged answer for that board.
- **Cut cells** keep thin walls and curved faces, but not yet with heat transfer, porous
  media, fans, rotating regions or transferred conditions.

## Process

- **Completed is not converged, and converged is not validated.** Execution, linear
  convergence, steady convergence and physical validation are separate questions. The
  manifest and the goal histories say what happened; whether it describes your product is
  for you to establish, for example by refining the grid until the numbers you care about
  stop changing: a [mesh study](../tutorials/comparing-runs.md#is-the-grid-fine-enough) does
  that on three grids and says how far each part's temperature still moves.
- **Results come from the companion** that ran them. SolidCraft checks their identity and
  bounds, not their physics.

## The workflow

- One fan per study, on a face of the model at the domain's boundary.
- Heat with inlets or outlets on faces inside the domain is not available yet (a wall on such a
  face, to measure a goal there, is).
