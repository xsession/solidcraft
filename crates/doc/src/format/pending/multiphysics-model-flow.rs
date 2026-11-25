//! Saves a Multiphysics study of the model's fluid and run length (`flow` on a study's `model`:
//! kinematic viscosity, density, end time and an optional convergence tolerance). New studies
//! get air; a study saved before has none and keeps its template's fluid and steps, so its
//! results are unchanged. The last release warns about the newer format, opens the design and
//! keeps the section it cannot read when saving.
