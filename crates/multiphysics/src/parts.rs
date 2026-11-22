//! The part registry of a design's Multiphysics project: tags and attributes on the design's
//! bodies, and named selections of them, either a saved list of bodies or a live query that is
//! worked out again against the bodies as they are each time it is used.
//!
//! A part is a body of the design, known by the body name a model study already uses
//! (`model::ModelLink::bodies`). Model studies point at a selection by name
//! (`model::ModelLink::selection`), so a study of "every tagged power part" follows the design
//! as parts are added or retagged instead of holding a list that goes out of date.
//!
//! This module is portable: it never reads geometry. The engine gives it a [`PartFacts`] per
//! body, measured from exact geometry where a query asks for a size or a volume.
use std::collections::{BTreeMap, BTreeSet};

use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::{Error, Result};

/// Most parts with tags or attributes, selections, tags or attributes on one part, and rules in one query.
pub const MAX_PARTS: usize = 4096;
pub const MAX_SELECTIONS: usize = 256;
pub const MAX_TAGS: usize = 32;
pub const MAX_ATTRIBUTES: usize = 32;
pub const MAX_RULES: usize = 16;
/// Built-in attributes every part has, measured from its exact geometry: its volume (mm³) and
/// the largest side of its bounding box (mm). A part's own attributes may not use these names.
pub const VOLUME: &str = "volume_mm3";
pub const SIZE: &str = "size_mm";

/// An attribute's value: a number or a short text.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(untagged)]
pub enum AttrValue {
    Number(f64),
    Text(String),
}

impl AttrValue {
    fn validate(&self, what: &str) -> Result<()> {
        match self {
            Self::Number(n) if n.is_finite() => Ok(()),
            Self::Number(_) => Err(Error::Invalid(format!("{what} must be a finite number"))),
            Self::Text(t) => text(t, what, true),
        }
    }
}

/// What the registry records about one part.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PartInfo {
    #[serde(default, skip_serializing_if = "BTreeSet::is_empty")]
    pub tags: BTreeSet<String>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub attributes: BTreeMap<String, AttrValue>,
}

impl PartInfo {
    pub fn is_empty(&self) -> bool {
        self.tags.is_empty() && self.attributes.is_empty()
    }
}

/// How a rule compares an attribute.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum Compare {
    /// The part has the attribute, whatever its value.
    Exists,
    Eq,
    Ne,
    Lt,
    Le,
    Gt,
    Ge,
}

/// One condition of a query; a part is in the query's selection when it meets every rule.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "rule", rename_all = "snake_case", deny_unknown_fields)]
pub enum Rule {
    /// The part carries the tag.
    Tag { tag: String },
    /// The part does not carry the tag.
    NotTag { tag: String },
    /// The body name matches a pattern where `*` stands for any run of characters and `?` for
    /// one; case is ignored.
    Name { pattern: String },
    /// An attribute (its own, or the built-in `volume_mm3` and `size_mm`) compared with a value.
    /// Numbers compare as numbers and texts as texts; a number never equals a text.
    Attribute {
        key: String,
        compare: Compare,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        value: Option<AttrValue>,
    },
    /// The part is in another selection.
    InSelection { name: String },
}

/// What a selection holds.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Members {
    /// These bodies, by name, as they were picked.
    Saved { bodies: Vec<String> },
    /// The bodies that meet every rule, worked out each time the selection is used.
    Query { rules: Vec<Rule> },
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PartSelection {
    pub name: String,
    pub members: Members,
}

/// Tags, attributes and selections of a project's parts.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
pub struct Registry {
    /// What is recorded about each part, by body name; parts with nothing recorded are absent.
    #[serde(skip_serializing_if = "BTreeMap::is_empty")]
    pub parts: BTreeMap<String, PartInfo>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub selections: Vec<PartSelection>,
}

/// One body as a query sees it.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct PartFacts {
    pub name: String,
    /// The built-in attributes, where the caller measured them (see [`Registry::needs_measure`]).
    pub volume_mm3: Option<f64>,
    pub size_mm: Option<f64>,
}

/// A change to the tags and attributes of some parts.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Retag {
    pub add_tags: Vec<String>,
    pub remove_tags: Vec<String>,
    /// Attributes to set; `None` removes one.
    pub attributes: BTreeMap<String, Option<AttrValue>>,
}

fn text(value: &str, what: &str, allow_spaces: bool) -> Result<()> {
    let bad_char = |c: char| c.is_control() || (!allow_spaces && c.is_whitespace());
    if value.trim().is_empty() || value.len() > 200 || value.chars().any(bad_char) {
        let spaces = if allow_spaces { "" } else { " or spaces" };
        return Err(Error::Invalid(format!("{what} must be 1 to 200 characters without control characters{spaces}")));
    }
    Ok(())
}

fn attribute_key(key: &str) -> Result<()> {
    text(key, "an attribute name", false)?;
    if key == VOLUME || key == SIZE {
        return Err(Error::Invalid(format!("{key} is measured from the part's geometry and cannot be set")));
    }
    Ok(())
}

/// Whether `name` matches `pattern` (`*` any run, `?` one character), ignoring case. Iterative,
/// so a hostile pattern costs at most `pattern × name` steps.
pub fn matches(pattern: &str, name: &str) -> bool {
    let p: Vec<char> = pattern.to_lowercase().chars().collect();
    let n: Vec<char> = name.to_lowercase().chars().collect();
    let (mut i, mut j) = (0usize, 0usize);
    let mut star: Option<(usize, usize)> = None;
    while j < n.len() {
        match p.get(i) {
            Some('*') => {
                star = Some((i, j));
                i += 1;
            }
            Some(&c) if c == '?' || Some(&c) == n.get(j) => {
                i += 1;
                j += 1;
            }
            _ => match star {
                Some((si, sj)) => {
                    i = si + 1;
                    j = sj + 1;
                    star = Some((si, sj + 1));
                }
                None => return false,
            },
        }
    }
    p.get(i..).is_some_and(|rest| rest.iter().all(|&c| c == '*'))
}

fn compare(have: &AttrValue, how: Compare, want: &AttrValue) -> bool {
    use std::cmp::Ordering;
    let order = match (have, want) {
        (AttrValue::Number(a), AttrValue::Number(b)) => a.partial_cmp(b),
        (AttrValue::Text(a), AttrValue::Text(b)) => Some(a.cmp(b)),
        _ => None,
    };
    match (how, order) {
        (Compare::Exists, _) => true,
        (Compare::Ne, None) => true,
        (_, None) => false,
        (Compare::Eq, Some(o)) => o == Ordering::Equal,
        (Compare::Ne, Some(o)) => o != Ordering::Equal,
        (Compare::Lt, Some(o)) => o == Ordering::Less,
        (Compare::Le, Some(o)) => o != Ordering::Greater,
        (Compare::Gt, Some(o)) => o == Ordering::Greater,
        (Compare::Ge, Some(o)) => o != Ordering::Less,
    }
}

impl Registry {
    pub fn is_empty(&self) -> bool {
        self.parts.is_empty() && self.selections.is_empty()
    }

    pub fn selection(&self, name: &str) -> Option<&PartSelection> {
        self.selections.iter().find(|s| s.name == name)
    }

    /// Bounded structure: names, limits, and selections that refer only to earlier selections,
    /// so a query never loops.
    pub fn validate(&self) -> Result<()> {
        if self.parts.len() > MAX_PARTS {
            return Err(Error::Limit(format!("the registry records at most {MAX_PARTS} parts")));
        }
        for (body, info) in &self.parts {
            text(body, "a body name", true)?;
            if info.is_empty() {
                return Err(Error::Invalid(format!("{body} is recorded with no tags or attributes")));
            }
            if info.tags.len() > MAX_TAGS || info.attributes.len() > MAX_ATTRIBUTES {
                return Err(Error::Limit(format!("a part has at most {MAX_TAGS} tags and {MAX_ATTRIBUTES} attributes")));
            }
            for tag in &info.tags {
                text(tag, "a tag", false)?;
            }
            for (key, value) in &info.attributes {
                attribute_key(key)?;
                value.validate("an attribute value")?;
            }
        }
        if self.selections.len() > MAX_SELECTIONS {
            return Err(Error::Limit(format!("a project has at most {MAX_SELECTIONS} part selections")));
        }
        let mut earlier = BTreeSet::new();
        for selection in &self.selections {
            text(&selection.name, "a selection name", true)?;
            validate_members(&selection.members, &earlier)?;
            if !earlier.insert(selection.name.as_str()) {
                return Err(Error::Invalid(format!("two selections are named {}", selection.name)));
            }
        }
        Ok(())
    }

    /// Change the tags and attributes of `bodies`; a part left with nothing recorded is dropped.
    pub fn retag(&mut self, bodies: &[String], change: &Retag) -> Result<()> {
        if bodies.is_empty() || bodies.len() > MAX_PARTS {
            return Err(Error::Invalid(format!("name 1 to {MAX_PARTS} bodies")));
        }
        for tag in change.add_tags.iter().chain(&change.remove_tags) {
            text(tag, "a tag", false)?;
        }
        for (key, value) in &change.attributes {
            attribute_key(key)?;
            if let Some(value) = value {
                value.validate("an attribute value")?;
            }
        }
        let mut next = self.clone();
        for body in bodies {
            text(body, "a body name", true)?;
            let info = next.parts.entry(body.clone()).or_default();
            for tag in &change.remove_tags {
                info.tags.remove(tag);
            }
            info.tags.extend(change.add_tags.iter().cloned());
            for (key, value) in &change.attributes {
                match value {
                    Some(value) => info.attributes.insert(key.clone(), value.clone()),
                    None => info.attributes.remove(key),
                };
            }
            if info.is_empty() {
                next.parts.remove(body);
            }
        }
        next.validate()?;
        *self = next;
        Ok(())
    }

    /// Add or replace the selection called `name`, keeping its place; `None` removes it. A
    /// selection a later one refers to cannot be removed or made to refer to a later one.
    pub fn set_selection(&mut self, name: &str, members: Option<Members>) -> Result<()> {
        let mut next = self.clone();
        let at = next.selections.iter().position(|s| s.name == name);
        match (at, members) {
            (Some(at), Some(members)) => {
                if let Some(s) = next.selections.get_mut(at) {
                    s.members = members;
                }
            }
            (None, Some(members)) => next.selections.push(PartSelection { name: name.to_string(), members }),
            (Some(at), None) => {
                next.selections.remove(at);
            }
            (None, None) => return Err(Error::NotFound(format!("there is no part selection called {name}"))),
        }
        next.validate()?;
        *self = next;
        Ok(())
    }

    /// Whether resolving the selection called `name` reads the built-in measured attributes, so
    /// the caller knows to measure each body first.
    pub fn needs_measure(&self, name: &str) -> bool {
        let mut seen = BTreeSet::new();
        self.needs_measure_inner(name, &mut seen)
    }

    fn needs_measure_inner<'a>(&'a self, name: &'a str, seen: &mut BTreeSet<&'a str>) -> bool {
        if !seen.insert(name) {
            return false;
        }
        let Some(PartSelection { members: Members::Query { rules }, .. }) = self.selection(name) else { return false };
        rules.iter().any(|rule| match rule {
            Rule::Attribute { key, .. } => key == VOLUME || key == SIZE,
            Rule::InSelection { name } => self.needs_measure_inner(name, seen),
            Rule::Tag { .. } | Rule::NotTag { .. } | Rule::Name { .. } => false,
        })
    }

    /// The bodies of `bodies` that the selection called `name` holds, in `bodies`' order. A
    /// saved selection's bodies that the design no longer has are listed in the second value.
    pub fn resolve(&self, name: &str, bodies: &[PartFacts]) -> Result<(Vec<String>, Vec<String>)> {
        let selection = self.selection(name).ok_or_else(|| Error::NotFound(format!("there is no part selection called {name}")))?;
        let mut depth = 0;
        let held = self.holds(selection, bodies, &mut depth)?;
        let chosen = bodies.iter().zip(&held).filter(|(_, h)| **h).map(|(b, _)| b.name.clone()).collect();
        let missing = match &selection.members {
            Members::Saved { bodies: saved } => saved.iter().filter(|s| !bodies.iter().any(|b| &&b.name == s)).cloned().collect(),
            Members::Query { .. } => vec![],
        };
        Ok((chosen, missing))
    }

    fn holds(&self, selection: &PartSelection, bodies: &[PartFacts], depth: &mut usize) -> Result<Vec<bool>> {
        *depth += 1;
        if *depth > MAX_SELECTIONS {
            return Err(Error::Limit("part selections refer to each other too deeply".into()));
        }
        let held = match &selection.members {
            Members::Saved { bodies: saved } => bodies.iter().map(|b| saved.contains(&b.name)).collect(),
            Members::Query { rules } => {
                let mut held = vec![true; bodies.len()];
                for rule in rules {
                    let meets: Vec<bool> = match rule {
                        Rule::InSelection { name } => {
                            let other = self.selection(name).ok_or_else(|| Error::NotFound(format!("there is no part selection called {name}")))?;
                            self.holds(other, bodies, depth)?
                        }
                        _ => bodies.iter().map(|b| self.meets(rule, b)).collect(),
                    };
                    for (h, m) in held.iter_mut().zip(meets) {
                        *h = *h && m;
                    }
                }
                held
            }
        };
        *depth -= 1;
        Ok(held)
    }

    fn meets(&self, rule: &Rule, body: &PartFacts) -> bool {
        let info = self.parts.get(&body.name);
        match rule {
            Rule::Tag { tag } => info.is_some_and(|i| i.tags.contains(tag)),
            Rule::NotTag { tag } => !info.is_some_and(|i| i.tags.contains(tag)),
            Rule::Name { pattern } => matches(pattern, &body.name),
            Rule::Attribute { key, compare: how, value } => {
                let have = match key.as_str() {
                    VOLUME => body.volume_mm3.map(AttrValue::Number),
                    SIZE => body.size_mm.map(AttrValue::Number),
                    _ => info.and_then(|i| i.attributes.get(key)).cloned(),
                };
                match (have, value) {
                    (None, _) => false,
                    (Some(_), _) if *how == Compare::Exists => true,
                    (Some(have), Some(want)) => compare(&have, *how, want),
                    (Some(_), None) => false,
                }
            }
            Rule::InSelection { .. } => false,
        }
    }
}

fn validate_members(members: &Members, earlier: &BTreeSet<&str>) -> Result<()> {
    match members {
        Members::Saved { bodies } => {
            if bodies.len() > MAX_PARTS {
                return Err(Error::Limit(format!("a saved selection holds at most {MAX_PARTS} bodies")));
            }
            for body in bodies {
                text(body, "a body name", true)?;
            }
        }
        Members::Query { rules } => {
            if rules.is_empty() || rules.len() > MAX_RULES {
                return Err(Error::Invalid(format!("a query has 1 to {MAX_RULES} rules")));
            }
            for rule in rules {
                match rule {
                    Rule::Tag { tag } | Rule::NotTag { tag } => text(tag, "a tag", false)?,
                    Rule::Name { pattern } => text(pattern, "a name pattern", true)?,
                    Rule::Attribute { key, compare: how, value } => {
                        text(key, "an attribute name", false)?;
                        match (how, value) {
                            (Compare::Exists, _) => {}
                            (_, Some(value)) => value.validate("a compared value")?,
                            (_, None) => return Err(Error::Invalid(format!("the rule on {key} needs a value to compare with"))),
                        }
                    }
                    Rule::InSelection { name } => {
                        if !earlier.contains(name.as_str()) {
                            return Err(Error::Invalid(format!("a selection may only use selections listed before it, and none called {name} is")));
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn facts(names: &[(&str, f64)]) -> Vec<PartFacts> {
        names.iter().map(|(n, v)| PartFacts { name: (*n).into(), volume_mm3: Some(*v), size_mm: Some(v.cbrt()) }).collect()
    }

    fn retag(tags: &[&str]) -> Retag {
        Retag { add_tags: tags.iter().map(|t| (*t).into()).collect(), ..Retag::default() }
    }

    #[test]
    fn patterns_match_like_a_file_name() {
        assert!(matches("U*", "u12"));
        assert!(matches("R?_*", "R1_top"));
        assert!(!matches("R?_*", "R12_top"));
        assert!(matches("*cap*", "decoupling CAP 3"));
        assert!(matches("**", ""));
        assert!(!matches("a", ""));
        // A hostile pattern stays cheap.
        assert!(!matches(&"*a".repeat(100), &"a".repeat(99)));
    }

    #[test]
    fn a_query_follows_tags_attributes_and_size() {
        let mut r = Registry::default();
        r.retag(&["U1".into(), "U2".into()], &retag(&["power"])).ok();
        let mut watts = Retag::default();
        watts.attributes.insert("power_w".into(), Some(AttrValue::Number(4.5)));
        r.retag(&["U1".into()], &watts).ok();
        let rules = vec![
            Rule::Tag { tag: "power".into() },
            Rule::Attribute { key: "power_w".into(), compare: Compare::Gt, value: Some(AttrValue::Number(1.0)) },
        ];
        assert_eq!(r.set_selection("hot", Some(Members::Query { rules })), Ok(()));
        let small = vec![Rule::Attribute { key: SIZE.into(), compare: Compare::Lt, value: Some(AttrValue::Number(3.0)) }];
        assert_eq!(r.set_selection("small", Some(Members::Query { rules: small })), Ok(()));
        let bodies = facts(&[("U1", 1000.0), ("U2", 1000.0), ("C1", 8.0)]);
        assert_eq!(r.resolve("hot", &bodies), Ok((vec!["U1".to_string()], vec![])));
        assert_eq!(r.resolve("small", &bodies).map(|r| r.0), Ok(vec!["C1".to_string()]));
        assert!(r.needs_measure("small") && !r.needs_measure("hot"));
        // Retagging changes what the live query holds.
        r.retag(&["U1".into()], &Retag { remove_tags: vec!["power".into()], ..Retag::default() }).ok();
        assert_eq!(r.resolve("hot", &bodies).map(|r| r.0), Ok(vec![]));
    }

    #[test]
    fn saved_selections_report_bodies_the_design_lost() {
        let mut r = Registry::default();
        let members = Members::Saved { bodies: vec!["A".into(), "Gone".into()] };
        assert_eq!(r.set_selection("pick", Some(members)), Ok(()));
        assert_eq!(r.resolve("pick", &facts(&[("A", 1.0), ("B", 1.0)])), Ok((vec!["A".to_string()], vec!["Gone".to_string()])));
    }

    #[test]
    fn selections_only_refer_to_earlier_ones() {
        let mut r = Registry::default();
        let uses = |name: &str| Members::Query { rules: vec![Rule::InSelection { name: name.into() }] };
        assert!(r.set_selection("a", Some(uses("a"))).is_err());
        assert_eq!(r.set_selection("base", Some(Members::Saved { bodies: vec!["X".into()] })), Ok(()));
        assert_eq!(r.set_selection("derived", Some(uses("base"))), Ok(()));
        assert_eq!(r.resolve("derived", &facts(&[("X", 1.0), ("Y", 1.0)])).map(|r| r.0), Ok(vec!["X".to_string()]));
        // Removing what a later selection uses is refused, and nothing changes.
        let before = r.clone();
        assert!(r.set_selection("base", None).is_err());
        assert_eq!(r, before);
    }

    #[test]
    fn hostile_registry_input_is_refused() {
        let mut r = Registry::default();
        assert!(r.retag(&["A".into()], &retag(&["two words"])).is_err());
        assert!(r.retag(&["A".into()], &retag(&[""])).is_err());
        let mut nan = Retag::default();
        nan.attributes.insert("k".into(), Some(AttrValue::Number(f64::NAN)));
        assert!(r.retag(&["A".into()], &nan).is_err());
        let mut builtin = Retag::default();
        builtin.attributes.insert(VOLUME.into(), Some(AttrValue::Number(1.0)));
        assert!(r.retag(&["A".into()], &builtin).is_err());
        assert!(r.retag(&[], &retag(&["t"])).is_err());
        assert!(r.set_selection("q", Some(Members::Query { rules: vec![] })).is_err());
        let no_value = vec![Rule::Attribute { key: "k".into(), compare: Compare::Lt, value: None }];
        assert!(r.set_selection("q", Some(Members::Query { rules: no_value })).is_err());
        assert!(r.set_selection("missing", None).is_err());
        assert!(r.is_empty());
        // A tag removed from every part leaves nothing recorded.
        r.retag(&["A".into()], &retag(&["t"])).ok();
        r.retag(&["A".into()], &Retag { remove_tags: vec!["t".into()], ..Retag::default() }).ok();
        assert!(r.is_empty());
    }
}
