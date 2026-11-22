//! Saves the link of a Multiphysics study to the loaded model (`model` on a study): whether the
//! fluid is inside or around the bodies, which bodies of the design it uses, its margin and cell
//! size, the conditions, local refinements and goals on model faces and bodies, and heat transfer
//! (body materials and power, heat sources and thermal conditions on faces; faces by persistent
//! face name), and the digest of the model input
//! its case was last built from. It is omitted for studies of a plain case, so older designs are
//! unchanged and open with none. The last release warns about the newer format, opens the
//! design and keeps the section it cannot read when saving.
