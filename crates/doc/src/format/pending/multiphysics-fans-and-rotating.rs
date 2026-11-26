//! Saves the fans inside the domain and the rotating regions of a Multiphysics study of the
//! model (`internal_fans` and `rotating` on a study's `model`: the bodies, a fan's direction and
//! curve, a region's axis, speed and margin). Omitted while there are none, so older designs are
//! unchanged. The last release warns about the newer format, opens the design and keeps the
//! section it cannot read when saving.
