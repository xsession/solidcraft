# Results and verdicts

## A run, from start to kept

1. **Run** freezes the study's committed case, its SHA-256 and the study's revision. A draft
   that has not been applied is never run.
2. The companion runs in a fresh private folder, with bounded output (16 MiB) and a
   five-minute limit. Cancelling waits until it has stopped.
3. When it finishes, its output is checked against what was frozen: the case it reports
   must be the case it was given. A late or mismatched result cannot become current.
4. **Keep result** stores the run in the design as **evidence**: its input case and hash, the
   study revision it was made from, the producer's manifest and its files. Until then it is
   only in the session's history, which keeps the last eight runs.

A kept run becomes **historical** when the study changes afterwards. It keeps its own input
case, so its numbers still mean what they meant; panels and reports mark it as made from an
older setup.

## Matching results to parts

The solver's `thermal_report.json` lists one row per solid: power, lowest, mean and highest
temperature, limit, margin and verdict. Solids are named in the case by their parts' names,
made safe for the solver (letters, digits, `_` and `-`; other characters become `_`, so
`R 2` is `R_2`). The **component table** joins each row back to its part by that name, and
adds the setup: tags, attributes, conductivity, power in the body and on its faces, the
limit. Its total power counts the parts in the study.

The verdict rules are on [The solver](solver.md#verdicts) page. The table, Find, colour by
verdict and the batch report all read the same table, so they always agree.

## Find

A Find is a list of 1 to 16 rules; a part is found when it meets every rule. A rule is either
a [part rule](parts.md#selections) on the setup or a result rule on the table:

- a quantity (power, limit, highest or mean temperature, margin) compared with a number, in
  °C or K; a part without that quantity does not meet it;
- a verdict (pass, warn, fail, no limit).

Result rules need a run with a thermal report. Non-finite numbers are refused.

## Plots

Field files are ASCII VTU grids of box cells. They are parsed off the window's thread, with
bounds (at most 32 768 cells for a slice). Values are shown as the producer wrote them: point
data interpolated, cell data constant per cell, nothing extrapolated outside the mesh.

### Colour scales

All field plots (the 2-D cut, the 3-D cut, the plot on the model) and their legend share one
colour scale:

| Setting | Options |
|---|---|
| Map | Viridis, Turbo, cool–warm, thermal, grayscale; reversed |
| Bands | 0 (smooth) or 2–32 bands of one colour each |
| Range | Data; Locked (typed); 2–98 % (percentiles); Around zero (symmetric) |
| Out of range | Its own colours below and above, chosen per map so they never look like a value on it |

### Where a cut plot goes

A cut is a plane across one axis of the grid, through one layer of cells. The layer is either
typed, or **anchored**: to the cell where the plotted field is largest, or to the centre of
the picked part's exact bounding box, placed in the case's coordinates (metres from the
domain's corner). An anchor is worked out again for every run, field, plane and pick.

## The batch report

A batch report reads the component table of every kept run with a thermal report (or the runs
named), and writes `report.csv`, `report.html` and a picture per run. The HTML has no scripts
and fetches nothing; every text in it is escaped. Pictures are drawn by SolidCraft's own
renderer, with each study part in its verdict's colour and the other bodies see-through.
