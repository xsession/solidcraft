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

The automatic cell size is at most half the thinnest body. In the example (160 × 60 × 40 mm,
4 mm cells) the grid is 40 × 15 × 10 = 6 000 cells.

## Bodies in the grid

Each body is taken as its closed triangle surface, in metres, and voxelised at the cell
centres: a cell whose centre is inside a body is that body's (a *staircase*). Consequences:

- A part thinner than a cell may hold no cell centre and is lost from the grid. The
  checklist lists such parts; the solver refuses a surface body that covers no cell centre.
  Refine around it, or [leave it out](../tutorials/imported-assemblies.md#leave-out-small-parts).
- Where two parts claim the same cell, the earlier part keeps it, with its heat. The
  checklist lists overlaps.
- A slanted or curved face becomes steps of cells.

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
