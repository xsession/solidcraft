# The solver

Studies of the model run on the companion's `collocated_flow3d` solver, from *cfd-solvers-rs*
0.16.0, a Rust port based on upstream revision `7cceaaad`. It runs as a separate program and
never opens a network connection.

## Flow

The incompressible Navier–Stokes equations are solved by finite volumes on the box grid of
the case, with velocity and pressure stored at cell centres (a *collocated* arrangement):

- Momentum and continuity are coupled by pressure correction (SIMPLE and PISO), with
  Rhie–Chow interpolation of the face fluxes so that the collocated pressure does not
  decouple into a checkerboard.
- No turbulence model is applied: the flow is laminar.
- **Pressure is kinematic**, pressure over density, in m²/s². Multiply by the air's density
  (about 1.2 kg/m³) for pascals. Outlets are given in the same unit.
- Steady cases march in steps to a settled state. Outputs written along the way are numbered
  by step, not by seconds.
- Goals are reported as histories; with goal control, the solver can stop once they have
  settled within a tolerance.

### Fluid and run length

A study of the model has its own **fluid** (air or water at 20 °C, or any kinematic viscosity
and density) and **end time**, under *Fluid and run* in its Mesh panel
(`multiphysics.model.set {"flow": …}`). From them the study works out its run:

- the **time step** is the larger one that keeps both limits of the explicit parts stable on
  the smallest cells: the flow crosses at most 0.4 of a cell per step at twice the fastest
  inflow (Courant number), and the viscous diffusion number stays at or below 0.4;
- the **steps** are as many as reach the end time, up to 20 000 and up to 50 million cell
  steps (cells × steps); past either the build says how much of the end time a run reaches;
- the pressure and the momentum are solved to tolerances **relative to their own size**
  (companion 0.18.0): with real fluids in small cells the solvers' fixed floor would otherwise
  leave them unsolved, and the flow would stay where it started;
- with **Stop once the goals settle** (on by default, at 0.1 %), the run ends early once every
  goal has changed by at most that share between evaluation windows. A study without goals of
  its own watches the fluid's mean speed for this.

A new study starts with air for two seconds. Studies made before these settings existed keep
running with their template's fluid (0.02 m²/s, a stand-in far more viscous than air) and
steps, so their results stay as they were; the panel says so, and picking a fluid runs them as
air or water from then on.

## Heat

With heat on, every solid of the study conducts heat (*conjugate heat transfer*):

- Each part has a conductivity (per axis where needed, for example for a board's copper
  layers), and its power is released evenly through its volume, or on a face.
- The fluid carries heat with the flow it settled into, and conducts it.
- The energy equation is solved **steady and one-way**, on the settled flow: the
  temperature does not change the flow.
- Properties are constant, and there is no buoyancy and no radiation in a study of the model.
- Unlisted outer faces are adiabatic; a face condition can make one a fixed temperature, a
  heat flux or a convection to the outside.

### Verdicts

For each part with a temperature limit, the solver reports its highest and mean temperature
and a verdict:

| Verdict | When |
|---|---|
| **fail** | The highest temperature is above the limit |
| **warn** | The margin (limit − highest) is below 5 % of the rise the limit allows above the inflow temperature |
| **pass** | Otherwise |
| **no limit** | The part has no limit |

In the example at 0.5 m/s, U1's limit of 85 °C allows 60 K above the 25 °C inflow; 5 % is
3 K, and U1's margin is 2.0 K, so it warns.

## What the solver checks

The solver checks the full case against its schema before it starts, and refuses what it
cannot do, with the reason: a surface body that covers no cell centre, an option it does not
support together with another (heat with cut cells), a speed its time step cannot follow.
SolidCraft passes the reason on unchanged.

!!! warning "Checked is not validated"
    A run that completes has solved the equations above on the grid it was given. Whether the
    grid is fine enough, and whether laminar, one-way heat transfer describes your product,
    are for you to judge. See [What it does not do](limits.md).
