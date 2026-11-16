//! Saves independent Multiphysics studies, bounded case inputs and retained text result artifacts.
//! A study may also carry an optional `setup` (named CAD selections, boundary assignments,
//! contacts and conduction); it is omitted while empty, so designs without one are unchanged.
//! Older designs get an empty collection. The last release warns about the newer format and
//! preserves sections it cannot read when saving. Active jobs, processes and local paths are absent.
