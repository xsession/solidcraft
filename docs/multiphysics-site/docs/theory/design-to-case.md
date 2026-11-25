# From design to case

A study of the model keeps a *link* to the design, not a copy of it. The **case**, the JSON the
solver reads, is built from the link and the design whenever the study's settings change,
and again before a run when the design has changed since.

## The link

| Setting | Meaning |
|---|---|
| **Analysis** | *Internal*: the fluid is inside the bodies (a duct, an enclosure). *External*: the fluid is around them |
| **Bodies** or **selection** | The bodies in the study: a list, every body, or a [part selection](parts.md) worked out again at each build |
| **Margin** | External only: space around the bodies on every side |
| **Cell size** or **mesh level** | A typed cell size, or a level from 1 (coarsest) to 7; level 3 is the automatic size, and each level up makes the cells about √2 smaller |
| **Openings** | Conditions on faces: inlets, outlets, walls, symmetry; a fan |
| **Refinements** | Finer cells around picked bodies or faces, or where a run's velocity changes most |
| **Heat** | Heat transfer on or off, the inflow temperature, materials (conductivity, power, limit), heat on faces, thermal conditions |
| **Goals** | Quantities the solver reports and may stop on: flows through faces, pressures, temperatures, global maxima and means |

## The domain and the grid

The **domain** is the study's bodies' exact bounding box: for an internal study as it is; for
an external one grown by the margin. The grid is a tensor-product grid of box cells over the
domain: uniform, or with finer rows where refinements ask. Refining around a part splits every
cell row crossing the part's box, so a refinement runs through the whole domain on the other
two axes. The solver accepts at most 32 768 cells.

Without a typed cell size, the cells are sized from what the bodies ask of the grid. On the
bodies' triangle surfaces, from each face that looks along an axis, one ray goes out into the
fluid and one into the part, along that axis, as the grid sees widths:

- the **narrowest gap**, the shortest way out to another surface, should be **2 cells** across;
- the **thinnest wall**, the shortest way in to the part's own far side, should be **2 cells**.

Both are robust minima: the smallest width under which at least 1 % of the sampled face area
lies, so a chamfer or a thread's end does not set the size of the whole study. A face resting
on another part (a chip on a board) has no gap there. The margin of an external study gets at
least one cell. When that is more than the 32 768 cells allow, the cells grow until it fits,
and the build says which gap or wall is left with fewer cells than it asked for, and how many
it got: refine around it, or [leave out](../tutorials/imported-assemblies.md#leave-out-small-parts)
the parts that cause it. The mesh level then scales the size up or down by √2 a level.

Local refinements come on top of that grid. Refining across a whole section splits every row
crossing it, so on an automatic grid the cells could pass the limit: the base cells then grow,
15 % at a time, until the refinements fit, and the build says by how much. With a typed cell
size the build refuses instead: that size is yours.

In the example the cell size is typed (4 mm): the grid is 40 × 15 × 10 = 6 000 cells.

## Bodies in the grid

Each body is taken as its closed triangle surface, in metres, and voxelised at the cell
centres: a cell whose centre is inside a body is that body's (a *staircase*). Consequences:

- A part thinner than a cell may hold no cell centre and has no cell of its own. With heat it
  is *lumped* (below); without, it is not in the flow. The checklist lists such parts. Refine
  around it, or [leave it out](../tutorials/imported-assemblies.md#leave-out-small-parts).
- Where two parts claim the same cell, the earlier part keeps it, with its heat. The
  checklist lists overlaps.
- A slanted or curved face becomes steps of cells.

### Partial cells

With heat on, the solver measures every cell on a lattice up to 8 times finer per axis and
takes it for what it holds, while the flow still sees the staircase:

- **Power and heat capacity go by solid volume.** A part's power is spread over the sub-cells
  it fills, so its total is exact and its power density does not depend on the cells' sizes
  (on a locally refined grid a staircase would give small cells more power per volume).
- **A cell of two materials conducts like both.** Along each axis its conductivity is the
  mean of two classic bounds over its sub-cells: columns along the axis in series, and layers
  across it in series. A layer of board in a cell conducts across the layer like the board
  and the air in series, along it like them side by side.
- **A part with no cell of its own is lumped.** Its power goes into the cells its volume lies
  in (or the cell at its centre, when it is smaller than a sub-cell), and its temperature in
  the component table is theirs.

The conductance between two cells is in series across the face between them, as before:
the harmonic mean, so heat crossing from a good conductor into a poor one is limited by the
poor one.

### How well the grid holds each part

The checklist measures each part of the study against the grid. Under **How the grid holds
each part** it lists:

| Column | What it is |
|---|---|
| Cells | The cells the part takes |
| Cells thick | Its fewest layers of cells along x, y or z: 0 is lost, 1 a single step |
| Volume held | Its cells' volume as a share of its exact volume, measured on the B-rep |

A part under **2 cells thick**, or whose cells hold more than **25 %** more or less than its
own volume, makes the *capture* line warn: its temperatures are those of a coarse block, so
refine around it. A part with power that takes **no cell** is an error (*lost power*): the run
has nowhere to put that power. And a fluid region the flow **enters but cannot leave**, an
inlet or fan with no outlet beyond it, is an error (*flow path*): parts close the way.

With **Cut cells** on (a flow study without heat), cells the surface crosses keep the fluid
part of their volume instead: each connected piece of fluid in a cut cell is a cell of its
own, with the open part of its faces and a wall face, measured on a lattice up to eight times
finer per axis. A wall thinner than a cell then still separates the fluid on its sides.

## Conditions on faces

Conditions are picked on faces of the model. A **flat face on the domain's boundary** opens
that side of the domain, as with the duct's ends. Any other face becomes a condition on its own
body: the body's grid faces nearest the face's triangles. A heat source on a face puts its
power into the body's cells behind those grid faces.

## When the case is built again

The case is built from the link and the design each time a study setting changes. The build
records a digest of everything it read. Before a run, and in the checklist (*Built from the
design as it is now*), that digest is compared with the design as it is now; if the design
changed (a moved part, a retag that changes a selection), the case is built again first.

## What a build reports

A build reports how much of the bodies' volume the grid holds, the cell counts, notes (a
body without a material conducts like the starter material, 1 W/(m·K); parts thinner than a
cell), and refuses cases the solver would refuse, with the reason.
