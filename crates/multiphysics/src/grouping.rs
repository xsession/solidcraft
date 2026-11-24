//! Import grouping: making sense of an imported assembly's many parts.
//!
//! - Reference designators: a board's parts arrive named by their designators (U12, C3, R10),
//!   whose letters say what kind of part each is. [`designator`] reads one from a body name and
//!   [`class`] names its kind, so a whole board can be tagged in one go.
//! - Similar parts: parts of one shape (every 0603 capacitor, every copy of a heat sink) are
//!   found from a part's exact volume and bounding box ([`similar`]).
//!
//! Portable: the engine measures each body's [`Shape`] from its exact geometry.
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::{Error, Result};

/// What a designator's letters mean, by the usual class letters for electronic parts. Each
/// kind is a tag the part registry accepts.
const CLASSES: &[(&str, &str)] = &[
    ("BT", "battery"),
    ("C", "capacitor"),
    ("D", "diode"),
    ("F", "fuse"),
    ("FB", "ferrite_bead"),
    ("IC", "ic"),
    ("J", "connector"),
    ("K", "relay"),
    ("L", "inductor"),
    ("LED", "led"),
    ("P", "connector"),
    ("Q", "transistor"),
    ("R", "resistor"),
    ("SW", "switch"),
    ("T", "transformer"),
    ("TP", "test_point"),
    ("U", "ic"),
    ("Y", "crystal"),
];

/// The attribute a grouped part's designator is kept in.
pub const DESIGNATOR: &str = "designator";

/// The reference designator a body name starts with, in capitals: one to three letters and one
/// to five digits, then the end of the name or a character that is neither (`U12`, `c3 (2)`,
/// `R10_1`). A name such as `Board` has none; `Cap1` reads as `CAP1`, a kind [`class`] does
/// not know.
pub fn designator(name: &str) -> Option<String> {
    let letters = name.chars().take_while(char::is_ascii_alphabetic).count();
    if !(1..=3).contains(&letters) {
        return None;
    }
    let rest = name.get(letters..)?;
    let digits = rest.chars().take_while(char::is_ascii_digit).count();
    if !(1..=5).contains(&digits) {
        return None;
    }
    if rest.chars().nth(digits).is_some_and(|c| c.is_ascii_alphanumeric()) {
        return None;
    }
    Some(name.get(..letters + digits)?.to_ascii_uppercase())
}

/// The kind of part a designator's letters name, as a tag (`U12` → `ic`), or none for letters
/// with no usual meaning.
pub fn class(designator: &str) -> Option<&'static str> {
    let letters: String = designator.chars().take_while(char::is_ascii_alphabetic).collect::<String>().to_ascii_uppercase();
    CLASSES.iter().find(|(l, _)| *l == letters).map(|(_, kind)| *kind)
}

/// One part grouped by its designator.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct Grouped {
    pub body: String,
    pub designator: String,
    /// The kind of part, the tag it is given.
    pub class: String,
}

/// The parts among `bodies` whose names carry a designator of a known kind, and the others.
pub fn group(bodies: &[String]) -> (Vec<Grouped>, Vec<String>) {
    let mut grouped = Vec::new();
    let mut other = Vec::new();
    for body in bodies {
        match designator(body).and_then(|d| class(&d).map(|c| (d, c))) {
            Some((designator, class)) => grouped.push(Grouped { body: body.clone(), designator, class: class.to_string() }),
            None => other.push(body.clone()),
        }
    }
    (grouped, other)
}

/// A part's shape as similar parts are compared: its volume and its bounding box's sides,
/// smallest first (so a part turned a quarter turn about an axis is still alike).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Shape {
    pub volume_mm3: f64,
    pub sides_mm: [f64; 3],
}

impl Shape {
    /// A shape from a volume and a box's three sides in any order; refused unless every value
    /// is finite and not negative.
    pub fn new(volume_mm3: f64, sides_mm: [f64; 3]) -> Result<Self> {
        if !volume_mm3.is_finite() || volume_mm3 < 0.0 || sides_mm.iter().any(|s| !s.is_finite() || *s < 0.0) {
            return Err(Error::Invalid("a part's volume and sides must be finite and not negative".into()));
        }
        let mut sides = sides_mm;
        sides.sort_by(f64::total_cmp);
        Ok(Self { volume_mm3, sides_mm: sides })
    }
}

/// Default and largest relative difference of similar parts.
pub const DEFAULT_TOLERANCE: f64 = 0.01;
pub const MAX_TOLERANCE: f64 = 0.5;

/// The parts whose volume and sides each differ from `like`'s by at most `tolerance` of the
/// larger (0.01: one per cent), in the order given.
pub fn similar(like: &Shape, parts: &[(String, Shape)], tolerance: f64) -> Result<Vec<String>> {
    if !tolerance.is_finite() || !(0.0..=MAX_TOLERANCE).contains(&tolerance) {
        return Err(Error::Invalid(format!("the tolerance must be 0 to {MAX_TOLERANCE}")));
    }
    // An absolute floor, so parts measured as equal compare equal however small they are.
    let near = |a: f64, b: f64| (a - b).abs() <= tolerance * a.abs().max(b.abs()) + 1e-9;
    let alike = |s: &Shape| near(s.volume_mm3, like.volume_mm3) && s.sides_mm.iter().zip(like.sides_mm).all(|(a, b)| near(*a, b));
    Ok(parts.iter().filter(|(_, s)| alike(s)).map(|(n, _)| n.clone()).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn designators_and_their_kinds() {
        for (name, want) in [
            ("U12", Some("U12")),
            ("c3 (2)", Some("C3")),
            ("R10_1", Some("R10")),
            ("LED4", Some("LED4")),
            ("Board", None),
            ("Cap1", Some("CAP1")),
            ("U", None),
            ("12", None),
            ("U123456", None),
            ("U1A", None),
            ("", None),
            ("Ü1", None),
        ] {
            assert_eq!(designator(name).as_deref(), want, "{name}");
        }
        assert_eq!(class("U12"), Some("ic"));
        assert_eq!(class("led4"), Some("led"));
        assert_eq!(class("X1"), None, "X is a crystal or a socket: no usual meaning");
        assert_eq!(class("CAP1"), None);
        let names: Vec<String> = ["U1", "C2", "Board", "X1", "TP7"].iter().map(|s| s.to_string()).collect();
        let (grouped, other) = group(&names);
        let kinds: Vec<(&str, &str)> = grouped.iter().map(|g| (g.body.as_str(), g.class.as_str())).collect();
        assert_eq!(kinds, [("U1", "ic"), ("C2", "capacitor"), ("TP7", "test_point")]);
        assert_eq!(other, ["Board", "X1"]);
        for (_, kind) in CLASSES {
            assert!(kind.chars().all(|c| c.is_ascii_lowercase() || c == '_'), "{kind} is a plain tag");
        }
    }

    #[test]
    fn similar_parts_match_by_volume_and_sides_in_any_order() {
        let shape = |v: f64, s: [f64; 3]| Shape::new(v, s).unwrap();
        let like = shape(1.6 * 0.8 * 0.45, [1.6, 0.8, 0.45]);
        let parts = vec![
            ("C1".to_string(), like),
            ("C2".to_string(), shape(1.6 * 0.8 * 0.45, [0.8, 0.45, 1.6])),
            ("C3".to_string(), shape(1.6 * 0.8 * 0.45 * 1.005, [1.6, 0.8, 0.45])),
            ("R1".to_string(), shape(1.0 * 0.5 * 0.35, [1.0, 0.5, 0.35])),
            ("Hollow".to_string(), shape(0.1, [1.6, 0.8, 0.45])),
        ];
        assert_eq!(similar(&like, &parts, DEFAULT_TOLERANCE).unwrap(), ["C1", "C2", "C3"]);
        assert_eq!(similar(&like, &parts, 0.0).unwrap(), ["C1", "C2"], "exactly alike");
    }

    #[test]
    fn hostile_shapes_and_tolerances_are_refused() {
        assert!(Shape::new(f64::NAN, [1.0; 3]).is_err());
        assert!(Shape::new(1.0, [1.0, f64::INFINITY, 1.0]).is_err());
        assert!(Shape::new(-1.0, [1.0; 3]).is_err());
        let like = Shape::new(1.0, [1.0; 3]).unwrap();
        for t in [f64::NAN, -0.1, 0.51, f64::INFINITY] {
            assert!(similar(&like, &[], t).is_err(), "{t}");
        }
    }
}
