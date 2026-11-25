# 5. Imported assemblies

A board brought in from an ECAD or STEP export has hundreds of parts named by their reference
designators. The **Group** section of the Parts panel turns those names into something you can
select, colour and set up by kind.

## Tag by designator

Press **Tag by designator**. Every part whose name starts with a reference designator, one to
three letters and one to five digits, then the end of the name or a character that is
neither (`U12`, `c3 (2)`, `R10_1`), gets the kind of part its letters name as a tag, and the
designator itself as the attribute `designator`:

| Letters | Tag | | Letters | Tag |
|---|---|---|---|---|
| U, IC | `ic` | | J, P | `connector` |
| R | `resistor` | | Y | `crystal` |
| C | `capacitor` | | F | `fuse` |
| L | `inductor` | | T | `transformer` |
| D | `diode` | | K | `relay` |
| Q | `transistor` | | SW | `switch` |
| LED | `led` | | TP | `test_point` |
| FB | `ferrite_bead` | | BT | `battery` |

Names such as `Board` are left as they are, and so are letters with no usual meaning (X is a
crystal on one board and a socket on another). The example gives:
U1, U2 → `ic`; Q1 → `transistor`; J1 → `connector`; C1, C2 → `capacitor`; R1 → `resistor`.
It is one undo step.

??? info "The command"
    `multiphysics.parts.group {}` (every body), or with `bodies` or `selection`.

## Part selections

A **part selection** is a named set of parts that a study, a material or a goal can point at.
It is either *saved* (these bodies) or a *query*, worked out again every time it is used:

```json
{"command": "multiphysics.parts.selection",
 "params": {"name": "Heat sources",
            "members": {"kind": "query", "rules": [{"rule": "tag", "tag": "ic"}]}}}
```

Query rules: *tagged*, *not tagged*, *name matches* (`*` and `?`), an *attribute* compared with
a value (your own, or `volume_mm3` and `size_mm`, measured from the part's exact geometry),
and *in another selection*. A study made of a query follows the design: tag a new part `ic`
and the next build includes it.

## Select similar parts

Pick one part and press **Select similar**: every part whose exact volume and bounding-box
sides each differ by at most 1 % is selected, so every 0603 capacitor, or every copy of a heat
sink, is picked in one go. A part turned a quarter turn about an axis still counts, because
the sides are compared in size order. Mesh bodies have no exact geometry and are never
similar.

??? info "The command"
    `multiphysics.parts.similar {"body": "C1", "tolerance_pct": 1, "select": true}`

## Leave out small parts

*Leave out parts under (mm)* keeps only the study's parts whose largest side is at least that
size. It is how a board with hundreds of passives becomes a study the grid can hold: the
small parts are left out, by a live query of the study's own (its earlier parts, and
`size_mm ≥` the size). A second size replaces the first. Their materials are kept for when
they come back.

## Modules: one setup for every copy

A rack holds the same board four times; a board holds the same regulator twice. Make the
repeated part a **component** and place it as often as it is used; then its Multiphysics setup
needs giving only once.

1. Set up one copy: its chips' materials, power and limits, their tags, a goal on the hottest
   chip, a heat source on a face.
2. Under **Modules** in the Parts panel, pick the component and press **Keep as module**. The
   module records that setup by the component's *own* body names (`U1`, or `Sub/R1` for a body
   of a component inside it), not by the names one copy has in the design (`U1 (Board:2)`).
3. Press **Apply to every copy**. Each placement gets the same materials, tags, attributes,
   goals and face heat, on its own bodies, and the study is built again, in one undo step.
   What a copy already had for those bodies is replaced, so applying again changes nothing
   more; a body a copy does not have is skipped, and the panel lists it.

A module is part of the design. To use it in another design, type a file name and press
**Export**; in the other design, type the file's name and press **Import**, then apply it to a
study there. The other design needs a component of the same name.

??? info "The commands"
    `multiphysics.module.capture {"id": "…", "component": "Board"}` keeps the setup of the
    component's first placement (or `"placement": "Board:2"`); `multiphysics.module.apply
    {"id": "…", "expected_revision": 7, "module": "Board setup"}` gives it to every placement;
    `multiphysics.module.export` / `.import {"path": "…"}` write and read it as
    `solidcraft-multiphysics-module/1` JSON.

Next: [scripting a study](scripting.md).
