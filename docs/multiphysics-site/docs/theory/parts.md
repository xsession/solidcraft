# Parts and selections

A **part** is a body of the design, known by its body name. The design's **part registry**
records what is known about parts beyond their geometry, and names sets of them.

## Tags and attributes

Each part can carry **tags** (single words, at most 32 per part) and **attributes**, a number
or a short text each (at most 32 per part). Two attributes are built in and measured, never
typed:

| Attribute | Meaning |
|---|---|
| `volume_mm3` | The part's volume, from its exact geometry |
| `size_mm` | The largest side of its exact bounding box |

They are measured only when a rule needs them. A mesh body has no exact geometry, so it has
neither.

## Selections

A **part selection** is a named set of parts. It is either:

- **saved**: these bodies, by name, as they were picked; or
- a **query**: the bodies that meet every one of its rules, worked out again each time it is
  used.

| Rule | A part meets it when |
|---|---|
| `tag` / `not_tag` | It carries (or does not carry) the tag |
| `name` | Its name matches a pattern, `*` for any run and `?` for one character, case ignored |
| `attribute` | An attribute (its own or `volume_mm3`, `size_mm`) compares with a value: exists, =, ≠, <, ≤, >, ≥. Numbers compare as numbers and texts as texts |
| `in_selection` | It is in another selection |

A study, a material or a command can name a selection instead of bodies. Because a query is
worked out again at each build, a study of *every part tagged `ic`* follows the design as
parts are added or retagged. A selection a study uses cannot be removed until the study uses
something else.

Bounds keep a hostile design from stalling SolidCraft: 4 096 parts with records, 256
selections, 16 rules in a query, and a name pattern is matched in at most pattern × name steps.

## Reference designators

**Tag by designator** reads a designator from the start of a part's name: one to three letters
and one to five digits, followed by the end of the name or by a character that is neither.
`U12`, `c3 (2)` and `R10_1` carry `U12`, `C3` and `R10`; `Board` carries none. The letters
name the kind of part by the usual class letters (U and IC: `ic`, R: `resistor`, C: `capacitor`,
…); letters with no usual meaning are left alone. The tag and the attribute `designator` are
written in one undo step.

## Similar parts

Two parts are **similar** when their exact volumes and the three sides of their exact bounding
boxes each differ by at most a tolerance (1 % by default, at most 50 %) of the larger value.
The sides are compared in size order, so a part turned a quarter turn about an axis matches.
The search measures bounding boxes first and computes volumes only for parts whose sides
already match.

## Leaving out small parts

*Leave out parts under* writes a query of the study's own, `<study> without small parts`: the
study's parts before (its earlier selection, or its list of bodies saved as `<study> parts`) and
`size_mm ≥` the size. The study then uses it. Calling again replaces the size instead of
stacking another rule. A material on a part the study no longer holds is kept, with a note in
the build, so the part keeps its setup when it comes back.

## The shared selection

The Multiphysics tree, the component table and the 3-D view show one selection, the design's
own. Picking a part in any of them picks it in the others; Find and Select similar put what
they found in it.
