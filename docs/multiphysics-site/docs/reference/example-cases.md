# Example cases

The repository's `example_cases` folder holds cases that show what Multiphysics and its solver
can do, each one checked by a test:

- **Studies of the model** (`example_cases/model`): eleven studies built from SolidCraft geometry
  by SolidCraft commands, one for each capability of the 3-D flow solver a study of the model
  uses, with an exact or textbook answer where there is one. A test builds every example and
  checks that its case asks the solver for what the example says it shows; another runs them
  through the solver and judges their results.
- **Case templates** (`example_cases/templates`): the 56 templates, one small complete case for
  each of the solver's families (flow, heat, structures, electromagnetics, plasma, chemistry,
  batteries, circuits and more), run as they come.

| Example | Shows | Last run |
|---|---|---|
| Square duct | air in a study, scaled momentum, goals | −1.8 % from the exact laminar pressure drop |
| Round pipe, cut cells | cut cells | −6.9 % from Hagen–Poiseuille |
| Local refinement | cells split around a body | all the inflow leaves |
| Board, conjugate heat | conjugate heat, partial cells | heat out equals heat in (imbalance 2e-9) |
| Cold-plate stack | heat on a face, a cold plate | energy balance; the stack 15 K above the plate |
| Fan curve | a fan on a face | −10 % from the developed laminar operating point |
| Porous brick | a porous body | −2.4 % from Darcy |
| Jet from a face | a condition on a face inside the domain | the jet's flow leaves, exactly |
| Natural convection | buoyancy | the heat sets the air moving (0.25 m/s) |
| Fan inside a duct | a fan inside the domain | +2.1 % from the developed laminar operating point; all its flow leaves |
| Turning cylinder | a rotating region (the Magnus effect) | 1.1 mPa higher on the side moving against the stream |

Each example is a JSON file: its title, what it shows, the capabilities, the commands and the
checks. The folder's README says how to run them and how to add one; the [scripting
tutorial](../tutorials/scripting.md) shows the same commands sent to a running SolidCraft.
