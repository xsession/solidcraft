# Drawing large assemblies

An electronics assembly can have thousands of parts. Colouring them by verdict, ghosting a
housing or picking a part must not mean sending every triangle to the GPU again.

## A colour table on the GPU

Each body has stable slots in a **palette**, a table of colours on the GPU: one for the body,
one per face, one for its edges and one for its placement. A body's triangles name their
slots, not their colours. Each body also has a slot in a **motion table**, so a body can be
drawn moved without its triangles changing.

Colour by verdict, by power or by tag therefore writes palette rows only: a test renders the
scene on the GPU and checks that colouring parts uploads **no geometry**. Ghosting a part
changes how its triangles are drawn (see-through), which uploads at most that part's sorted
see-through triangles.

The parts that are coloured or ghosted are kept in ordered maps by name, so drawing a scene
looks each part up instead of searching a list per part. With 3 000 parts (a debug build),
a scene rebuild took:

| | Before | After |
|---|---|---|
| No colours | 8.4 ms | 8.2 ms |
| Coloured by setup | 20 ms | 9.8 ms |
| Half the parts ghosted | 27 ms | 11.8 ms |

## Ghosting

A ghosted part is drawn:

- **see-through**, at 30 % opacity, in the plain body colour whatever its appearance. Where the
  device supports it, see-through surfaces are blended without sorting (order-independent
  transparency), so overlapping ghosts do not flicker as the view turns;
- with its **edges at half strength**, so the solid parts seen through it stand out;
- with an **outline** seen from the view (its silhouette), in the view's text colour at 55 %,
  so a curved housing keeps its shape. Outlines are traced from the display mesh for the
  ghosted parts only, a bounded amount of work per frame, so a large scene fills them in over
  a few frames instead of stalling one.

Clicks pass through a ghosted part to what is behind it. That decision, which part a click
means, uses the parts' exact geometry, not the mesh drawn.

## The overview inset

While parts are ghosted, the view's corner shows the whole model, orthographic, from the view's
direction and framed to the model's box. The inset is a second picture of the scenes the view
already holds on the GPU: it borrows the view's uploads for its frame, so nothing is built or
uploaded for it, and it is drawn again only when the view turns or the scenes change. Panning
and zooming move only its frame, which is where the rays through the view's corners meet the
plane through the model's centre square to the view.

The inset is drawn before the view in each frame, so the view's own settings reach the GPU
after the inset's and the view is drawn exactly as without it (a test compares the pixels).
It takes no clicks, and the view does not pick or navigate through it.

## Pictures without a window

The same pipelines draw pictures without a window: the batch report's pictures, the pictures
on these pages, and SolidCraft's own tests, which render scenes offscreen and check them
pixel by pixel.
