//! Saves the part registry of a design's Multiphysics project (`parts` on `multiphysics`): tags
//! and attributes on bodies, by body name, and named part selections, each a saved list of
//! bodies or a live query of rules on tags, names and attributes. A model study may take its
//! bodies from a selection (`selection` on a study's `model`). Both are omitted while empty, so
//! older designs are unchanged and open with none. The last release warns about the newer
//! format, opens the design and keeps the section it cannot read when saving.
