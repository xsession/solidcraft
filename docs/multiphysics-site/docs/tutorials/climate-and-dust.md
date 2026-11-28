# 7. Climate, humidity and dust

You will bring a real board and its printed case in from KiCad and STEP, put them on a shelf in
a dusty shed, and work out what a year of weather does to them: how warm the parts run, where
the air gets humid enough to wet dust, and where the dust settles.

The example is the [Universal Weather Node](https://github.com/franzenhof/universal-weather-node)
v2.0, an open-hardware weather station (hardware under CERN-OHL-P-2.0): a 100 × 100 mm board
with a radio module, RJ12 jacks and screw terminals, in a 3D-printed PETG case with cut-outs for
its connectors and cables.

![The board in its case, the case drawn see-through](../images/climate-design.png)

!!! tip "The script"
    [`tools/climate_example.py`](https://github.com/tangency-labs/solidcraft/blob/main/docs/multiphysics-site/tools/climate_example.py)
    does every step of this page with commands, from the downloaded project to a study ready
    to run, and [`tools/climate_year.py`](https://github.com/tangency-labs/solidcraft/blob/main/docs/multiphysics-site/tools/climate_year.py)
    runs the year and draws the maps. Read them alongside.

## Bring the board and the case in

1. **Open** the board's `.kicad_pcb`. Leave copper and silkscreen out: no flow grid holds a
   35 µm trace. Parts whose 3D models the project ships come in as meshes; the rest need KiCad's
   own model library.
2. **Insert** the case's STEP assembly. This one also holds the project's own blocks for the
   radio module, the jacks and the terminals, which stand in for the parts the board lacks.
3. Move the board onto its standoffs: the assembly puts the board's underside at z = 8.4 mm.
   Name the bodies by their designators (U1 for the module, J3 … J18 for the connectors) so
   [tagging](imported-assemblies.md) and the year's table read well.
4. Draw a plain 100 × 100 × 1.6 mm block where the board is, named **Board**, for the study.
   The KiCad board has some 600 drilled holes and pads: as a surface for the grid that is far
   more than a 3.6 mm cell can use, and more than a study takes.

??? info "The commands"
    `doc.open {"path": "…/universal-weather-node.kicad_pcb", "copper": "none", "silkscreen": false}`,
    `file.insert_step {"path": "…-assembly.step"}`, `occurrence.move`, `body.rename`,
    `solid.box`.

## The study

Make a **Model study**, **External** (the air around the case as well as inside it), with an
8 mm margin, of the case, lid, button, Board, U1 and the connectors. Give it a fluid run of 30 s
that stops once its goals settle.

Turn **Heat** on with **natural convection**, and give the materials: PETG case, lid and button
(0.2 W/(m·K)), the Board (3 W/(m·K) for FR-4 with copper planes, 0.08 W for the regulator and the
rest), U1 (5 W/(m·K), 0.42 W: an ESP32-S3 with Wi-Fi on and LoRa listening) and the connectors
(0.25 W/(m·K)). Add goals on U1's highest and the Board's mean temperature: the run stops when
they settle.

### Close thin walls

The grid holds at most 32 768 cells, which makes them 3.6 mm here, and the case's walls are 1.4
to 2.6 mm. A wall thinner than a cell holds a cell centre only where one happens to fall inside
it: in this case the side walls did, the end walls and nearly all of the lid did not, and the box
was open at both ends and on top.

Tick **Close thin walls** (Domain panel, under the grid). Every grid face that a body's surface
crosses between two air cells then becomes a wall of that body, so the case keeps only its real
openings. A small pocket of air the closed walls seal off (here 6 cells inside the left wall) is
filled with the wall's solid; the run's mesh report lists each one.

??? info "The command"
    `multiphysics.model.set {"id": "…", "expected_revision": 9, "close_thin_walls": true}`

!!! tip "Is it sealed where it should be?"
    Pick a face inside the case and one outside it and press **Find a leak between the two selected faces** (Domain panel): the
    path it draws goes through an opening, or there is none.

## Surroundings

In the **Environment** section, set a **breeze** of 0.1 m/s and keep **Floor under the model**.
The breeze comes in at −X and leaves at +X, the sides along it slip, the top is open: a shelf in
an open shed. The climate below changes the breeze through the day.

??? info "The command"
    `multiphysics.model.surroundings {"id": "…", "expected_revision": 10, "breeze_m_per_s": 0.1, "breeze_from": "-x", "floor": true}`

## Dust

Set **Dust** to 0.5 mg/m³. The sizes default to coarse outdoor dust by mass: 5 % at 1 µm,
10 % at 2.5, 20 % at 5, 35 % at 10 and 30 % at 20 µm, 2 000 kg/m³.

What matters is the dust that gets *into* the case, so count it there: pick the case body, type
its wall, and press **Count what enters the picked body**. The dust is then released where air
crosses into that box. Put the box's sides inside the walls, not in the air on either side,
or air already inside counts as coming in. This case's posts and tabs stand out past its walls,
so its box is given by hand: x −50.5 … 50.5, y 21 … 152, z 1 … 38 mm.

The solver follows each particle through the settled flow: drag on it (with the slip of very
small particles), gravity, the push from warm air towards cold (thermophoresis) and Brownian
motion, until it lands on a part, leaves through an opening or runs out of time.

??? info "The command"
    ```json
    {"command": "multiphysics.model.dust", "params": {"id": "…", "expected_revision": 11,
     "concentration_mg_per_m3": 0.5,
     "entry": {"kind": "box", "low_mm": [-50.5, 21, 1], "high_mm": [50.5, 152, 38]}}}
    ```

## The year

Under **Climate**, type the four seasons: each season's low and high temperature, its dew point
(the air's water, which changes little through a day) and its dust by day. **Set the year**
makes 16 conditions: each season's night, morning, afternoon and evening, 6 hours a day for a
quarter of the year, the warmest at 15:00. The breeze and the dust follow the day: calm and clean
at night (0.05 m/s, half the dust), windy and busy in the afternoon (0.25 m/s, one and a half
times the dust).

| Season | Low °C | High °C | Dew point °C | Dust mg/m³ |
|---|---|---|---|---|
| Winter | −4 | 3 | −5 | 0.1 |
| Spring | 4 | 16 | 3 | 0.4 |
| Summer | 14 | 30 | 13 | 0.8 |
| Autumn | 6 | 16 | 7 | 0.5 |

Then **Make one study per condition**: 16 studies, "Weather node climate · summer afternoon" and
so on, each with the air at its temperature, its breeze and its dust.

??? info "The commands"
    `multiphysics.model.climate {"id": "…", "expected_revision": 12, "seasons": [{"name": "winter", "low_c": -4, "high_c": 3, "dew_point_c": -5, "dust_mg_per_m3": 0.1}, …]}`
    (or `conditions`: your own list, each with `name`, `hours_per_year`, `ambient_c`,
    `dew_point_c`, and optionally `breeze_m_per_s` and `dust_mg_per_m3`), then
    `multiphysics.model.climate_studies {"id": "…"}`.

!!! note "Why a quarter of a day can stand still"
    Each condition is solved as if it lasted: the flow and the temperatures settle to it. The
    case is about 150 g of PETG with some 0.06 m² of outside, so it follows the air with a time
    constant of about 8 minutes, and the air changes by at most 2 K an hour: the walls lag
    the air by about 0.3 K. A day stands for its four quarters well enough.

## Run them

In the app, run the 16 studies (**Batch run**) and keep what you need. A year of a full grid is
more than a design holds (16 MiB a run, 64 MiB in all), so for this size the runs are better
made outside it:

```sh
python3 tools/climate_year.py prepare weather-node.sc "Weather node climate" year/
python3 tools/climate_year.py run year/ --jobs 6 --dt 0.02 --end-time 20
python3 tools/climate_year.py summary year/
```

`prepare` makes the condition studies and writes their cases; `run` runs them with the solver,
six at a time; `summary` adds the year up. In the app a run is kept short, with a cautious time
step; `--dt` and `--end-time` set them outside it, and a run still stops once its goals settle.

## The year's summary

**Year summary** (or `summary`) adds every condition up by its hours.

RESULTS-TABLE

RESULTS-MAPS

### Humid dust

Dust holds salts, and salts take up water well below saturation: from about 60 % relative
humidity at their surface, dust turns into a thin conducting film. The relative humidity at a
surface is the air's water over what the air could hold at the *surface's* temperature, so a
part warmer than the air is drier than the air around it.

The summary gives, for each part, the hours a year its dust spends at or above that limit
(*Humid dust*, weighed by the mass of dust), and the hours its coldest point is below the dew
point, when water condenses on it (*Below dew point*).

RESULTS-HUMID

## What it leaves out

- **Laminar flow**, and a grid of 3.6 mm cells: the openings are two or three cells across. A
  [mesh study](comparing-runs.md#is-the-grid-fine-enough) on a cut-down model (the board's
  end of the case) says how far the numbers still move.
- **Every cut-out open.** The case's connector cut-outs are holes in the CAD; in use a plug fills
  most of each. This is the case with nothing plugged in, the most dust can get in.
- **Dust sticks where it lands, and stays.** No bounce, no wind taking it off again, no
  cleaning. Coarse dust that settles on top of the case outside is not counted (the box counts
  what enters).
- **No sun, no rain, no radiation.** On a shelf in a shed. Sunshine on the lid would warm the
  case by several kelvin by day.
- **Breathing is left out**: the case's air (about 0.4 L) swells and shrinks with the day's
  16 K by 0.02 L a day, against thousands of litres the breeze moves through the cut-outs.
- **A representative year**, not a record: put your site's seasons in.

Previous: [scripting a study](scripting.md).
