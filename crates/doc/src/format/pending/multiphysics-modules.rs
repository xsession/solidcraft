//! Saves a design's Multiphysics modules (`modules` on `multiphysics`): a component's setup by
//! the component's own body names (materials, power, limits, tags, attributes, temperature goals
//! and face heat sources), to give every copy of the component and other designs. Omitted while
//! a design has none, so older designs are unchanged and open with none. The last release warns
//! about the newer format, opens the design and keeps the section it cannot read when saving.
