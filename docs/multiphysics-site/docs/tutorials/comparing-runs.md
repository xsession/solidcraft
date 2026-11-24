# 4. Comparing runs

One run answers one question. To see what the inlet speed does, run the study again at
0.5 m/s and compare.

## A second run

Pick **Boundary conditions**, pick the inlet face and set its speed to **0.5 m/s**. The study's
revision goes up, and the first run is now *historical*: it was made from the study as it
was. It keeps its own input case, so its results still mean what they meant; the panels mark
it as made from an older setup.

Run again and keep the result. The study now holds both runs.

## The batch report

Under the component table in the Parts panel, **Batch report** takes the component table of
every run that kept a thermal report, oldest first. Type a new folder's name and press
**Write report**. The folder gets:

| File | What it holds |
|---|---|
| `report.html` | A page that needs nothing else: each run's total power, hottest part, least margin and pass / warn / fail counts, then each part's highest temperature per run, coloured by its verdict, and a picture of each run's parts by verdict |
| `report.csv` | Every part's numbers per run: power, limit, highest and mean temperature, margin, verdict |
| `run-1.png`, `run-2.png`, … | Each run's parts in their verdict colours, iso view (when SolidCraft can draw pictures, which it always can in the app) |

The example's report: [report.html](../examples/report/report.html),
[report.csv](../examples/report/report.csv).

| Run | Inlet | Hottest part | Max (°C) | Least margin | Pass | Warn | Fail |
|---|---|---|---|---|---|---|---|
| mp-1-1 *(older setup)* | 1.0 m/s | U2 | 77.1 | Q1, 2.9 K | 4 | 0 | 0 |
| mp-1-2 | 0.5 m/s | U2 | 86.7 | Q1, −7.0 K | 1 | 1 | 2 |

Halving the inlet speed raises every part by 8 to 10 K, and three of the four parts with
limits stop passing. Q1 has the least margin in both runs, though U2 is the hottest: Q1's
limit is 10 K lower.

!!! note "The report never edits the design"
    It reads kept runs only, and writes into a new folder (an existing one is refused, so
    nothing is overwritten). Writing fails as a whole: a folder half written is removed.

??? info "The command"
    `multiphysics.result.report {"id": "…", "directory": "/path/to/new-folder"}`; or
    `"inline": true` to have the CSV and HTML returned instead of written; `run_ids` to pick
    the runs and their order (at most 64).

## Two runs part by part

The report puts every run side by side; **Compare runs**, under it, takes two and lines them up
part by part. Pick the **First run** and the **Second run** and press **Compare**. The runs
can be of this study (two operating points, as here) or of two studies of the model: a study
and its clone with a part moved, or two points of a sweep. Parts are matched by name.

![The two runs compared part by part](../images/compare-runs.png)

Each part's highest temperature in both runs, the change, and its verdict, the largest rise
first. A rise is amber and a fall green; a verdict that got worse is red (**Pass → Fail**),
one that got better green. The line above the table counts the parts warmer, cooler, worse
and better, and gives the change in total power: none here, as only the airflow changed.
Changes under 0.05 K count as none.

| Part | First (°C) | Second (°C) | Change (K) | Verdict |
|---|---|---|---|---|
| Q1 | 72.1 | 82.0 | +9.9 | Pass → Fail |
| U2 | 77.1 | 86.7 | +9.6 | Pass → Fail |
| J1 | 53.3 | 62.5 | +9.3 | Pass |
| U1 | 74.6 | 83.0 | +8.4 | Pass → Warning |

??? info "The command"
    `multiphysics.result.compare {"first": {"id": "…", "run_id": "mp-1-1"}, "second":
    {"id": "…"}}`: a side without `run_id` takes its study's latest run with a thermal report.
    `"csv": true` adds the comparison as CSV, one line per part.

To see the change *where* it happens rather than per part, lock the colour range of the cut
plot ([cut and surface plots](plots.md)) and show the same plane in both runs: one scale, so
the colours mean the same temperature in both.

Next: [imported assemblies](imported-assemblies.md).
