//! Saves a Multiphysics study of the model's environment: `close_thin_walls`, the air around an
//! external study (`surroundings`: breeze, its side, a floor), the dust in its air (`dust`:
//! concentration, sizes, density, particles, where it is counted from), its year of conditions
//! (`climate`) and, on a study made for one condition, that `condition`. Studies saved before
//! have none of them and build and run as before. The last release warns about the newer
//! format, opens the design and keeps the section it cannot read when saving.
