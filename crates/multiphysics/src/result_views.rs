//! Saved result views of a study: how its results were being looked at, kept by name so the same
//! picture can be brought back after another run, by someone else, or by a script. A view holds
//! the camera, what the parts are coloured by, which parts are ghosted, the cut plot (field,
//! direction and what its plane goes through) and the colour scale and range.
//!
//! A view names parts by body name and fields by name, never by index, so it still means the
//! same after the design or the run changes; a part the design no longer has is skipped when the
//! view is applied, and said so. Views are not setup: saving one does not make runs stale.
//!
//! Portable: data and checks only; the engine applies a view to the session, the host to its plots.
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};

use crate::results::{ColorScale, MAX_BANDS, RangeMode};
use crate::setup_check::ColourBy;
use crate::{Error, Result, Study};

/// Most views one study keeps.
pub const MAX_VIEWS: usize = 64;
/// Most characters in a view's name.
pub const MAX_NAME: usize = 128;
/// Most ghosted parts one view lists.
pub const MAX_GHOSTED: usize = 100_000;

/// An orbit camera, as the design's named views keep one: the point it turns about (mm), its
/// heading and elevation (radians), how far it stands off, its field of view (0: parallel) and roll.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ViewCamera {
    pub target: [f64; 3],
    pub yaw: f64,
    pub pitch: f64,
    pub distance: f64,
    #[serde(default)]
    pub fov: f64,
    #[serde(default)]
    pub roll: f64,
}

/// What a cut plot's plane goes through.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum PlotAnchor {
    /// The layer of cells given.
    #[default]
    Layer,
    /// The centre of the part named.
    Part,
    /// The cell where the field is largest.
    Maximum,
}

/// A cut plot on the model.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct PlotView {
    /// The field's name in the run's field file, such as `temperature`.
    pub field: String,
    /// The plane's normal: 0 for x, 1 for y, 2 for z.
    pub axis: u8,
    #[serde(default)]
    pub anchor: PlotAnchor,
    /// The layer of cells, for `layer`.
    #[serde(default)]
    pub layer: usize,
    /// The part, for `part`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub part: Option<String>,
    /// Velocity arrows on the plot.
    #[serde(default)]
    pub vectors: bool,
}

/// One saved view.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ResultView {
    pub name: String,
    /// None keeps the camera where it is.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub camera: Option<ViewCamera>,
    #[serde(default)]
    pub colour_by: ColourBy,
    /// Parts drawn see-through, by body name.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub ghosted: Vec<String>,
    /// The cut plot; none for no plot.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub plot: Option<PlotView>,
    /// The model's wetted faces coloured by the run (surface plots).
    #[serde(default)]
    pub surfaces: bool,
    #[serde(default)]
    pub scale: ColorScale,
    #[serde(default)]
    pub range: RangeMode,
    /// The low and high values of a `locked` range.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub locked_range: Option<[f64; 2]>,
}

fn invalid(message: impl Into<String>) -> Error {
    Error::Invalid(message.into())
}

/// Refuses a view that could not be drawn as saved: a blank or long name, a camera that is not
/// finite or has no distance, a plot with no field or a bad axis, or a locked range that is not
/// a finite low below a high.
pub fn check(view: &ResultView) -> Result<()> {
    let name = view.name.trim();
    if name.is_empty() || name.chars().count() > MAX_NAME {
        return Err(invalid(format!("a view's name is 1 to {MAX_NAME} characters")));
    }
    if let Some(c) = &view.camera {
        let numbers = [c.target[0], c.target[1], c.target[2], c.yaw, c.pitch, c.distance, c.fov, c.roll];
        if !numbers.iter().all(|x| x.is_finite()) || c.distance <= 0.0 || c.fov < 0.0 {
            return Err(invalid("the camera must be finite, with a positive distance and a field of view of at least 0"));
        }
        if c.pitch.abs() > std::f64::consts::FRAC_PI_2 + 1e-9 {
            return Err(invalid("the camera's pitch must be within ±90°"));
        }
    }
    if view.ghosted.len() > MAX_GHOSTED {
        return Err(invalid(format!("a view ghosts at most {MAX_GHOSTED} parts")));
    }
    if let Some(plot) = &view.plot {
        if plot.field.trim().is_empty() {
            return Err(invalid("a cut plot names its field"));
        }
        if plot.axis > 2 {
            return Err(invalid("a cut plot's axis is 0 (x), 1 (y) or 2 (z)"));
        }
        if plot.anchor == PlotAnchor::Part && plot.part.as_deref().is_none_or(|p| p.trim().is_empty()) {
            return Err(invalid("a cut plot through a part names the part"));
        }
    }
    if view.scale.bands > MAX_BANDS || view.scale.bands == 1 {
        return Err(invalid(format!("a colour scale has 0 (continuous) or 2 to {MAX_BANDS} bands")));
    }
    match (view.range, view.locked_range) {
        (RangeMode::Locked, None) => return Err(invalid("a locked range gives its low and high values")),
        (_, Some([low, high])) if !(low.is_finite() && high.is_finite() && low < high) => {
            return Err(invalid("a locked range is a finite low below a finite high"));
        }
        _ => {}
    }
    Ok(())
}

/// Keeps `view` in `study`, in place of one of the same name (names compare trimmed); refused
/// past [`MAX_VIEWS`]. Returns whether it replaced one.
pub fn save(study: &mut Study, mut view: ResultView) -> Result<bool> {
    view.name = view.name.trim().to_string();
    check(&view)?;
    let mut seen = std::collections::BTreeSet::new();
    view.ghosted.retain(|g| seen.insert(g.clone()));
    if let Some(slot) = study.views.iter_mut().find(|v| v.name == view.name) {
        *slot = view;
        return Ok(true);
    }
    if study.views.len() >= MAX_VIEWS {
        return Err(invalid(format!("a study keeps at most {MAX_VIEWS} views")));
    }
    study.views.push(view);
    Ok(false)
}

/// Drops the view named `name`; an error if there is none.
pub fn remove(study: &mut Study, name: &str) -> Result<ResultView> {
    let at =
        study.views.iter().position(|v| v.name == name.trim()).ok_or_else(|| invalid(format!("study {} has no view named {name}", study.name)))?;
    Ok(study.views.remove(at))
}

/// The view named `name`.
pub fn find<'a>(study: &'a Study, name: &str) -> Result<&'a ResultView> {
    study.views.iter().find(|v| v.name == name.trim()).ok_or_else(|| invalid(format!("study {} has no view named {name}", study.name)))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn study() -> Study {
        crate::Project::default().create("flow-channel", "Board").unwrap()
    }

    fn view(name: &str) -> ResultView {
        ResultView {
            name: name.into(),
            camera: Some(ViewCamera { target: [10.0, 0.0, 5.0], yaw: 0.6, pitch: 0.4, distance: 300.0, fov: 0.5, roll: 0.0 }),
            colour_by: ColourBy::Verdict,
            ghosted: vec!["Duct".into(), "Duct".into()],
            plot: Some(PlotView { field: "temperature".into(), axis: 2, anchor: PlotAnchor::Maximum, layer: 0, part: None, vectors: false }),
            surfaces: false,
            scale: ColorScale::default(),
            range: RangeMode::Locked,
            locked_range: Some([25.0, 90.0]),
        }
    }

    #[test]
    fn a_view_is_kept_by_name_and_replaced_by_the_same_name() {
        let mut s = study();
        assert!(!save(&mut s, view(" Hot spots ")).unwrap());
        assert_eq!(s.views[0].name, "Hot spots");
        assert_eq!(s.views[0].ghosted, ["Duct"], "a part ghosted twice is listed once");
        let mut other = view("Hot spots");
        other.colour_by = ColourBy::Power;
        assert!(save(&mut s, other).unwrap());
        assert_eq!(s.views.len(), 1);
        assert_eq!(find(&s, "Hot spots").unwrap().colour_by, ColourBy::Power);
        assert_eq!(remove(&mut s, "Hot spots").unwrap().name, "Hot spots");
        assert!(remove(&mut s, "Hot spots").is_err());
        // Saved and read back as the design keeps it.
        save(&mut s, view("Side")).unwrap();
        let text = serde_json::to_string(&s).unwrap();
        let back: Study = serde_json::from_str(&text).unwrap();
        assert_eq!(back.views, s.views);
    }

    #[test]
    fn views_that_cannot_be_drawn_are_refused() {
        let mut s = study();
        let bad = |f: &dyn Fn(&mut ResultView)| {
            let mut v = view("v");
            f(&mut v);
            check(&v).is_err()
        };
        assert!(bad(&|v| v.name = "  ".into()));
        assert!(bad(&|v| v.name = "x".repeat(MAX_NAME + 1)));
        assert!(bad(&|v| v.camera.as_mut().unwrap().distance = 0.0));
        assert!(bad(&|v| v.camera.as_mut().unwrap().yaw = f64::NAN));
        assert!(bad(&|v| v.camera.as_mut().unwrap().pitch = 2.0));
        assert!(bad(&|v| v.plot.as_mut().unwrap().axis = 3));
        assert!(bad(&|v| v.plot.as_mut().unwrap().field = String::new()));
        assert!(bad(&|v| v.plot.as_mut().unwrap().anchor = PlotAnchor::Part));
        assert!(bad(&|v| v.locked_range = None));
        assert!(bad(&|v| v.locked_range = Some([5.0, 5.0])));
        assert!(bad(&|v| v.scale.bands = 1));
        assert!(!bad(&|_| {}));
        for i in 0..MAX_VIEWS {
            save(&mut s, view(&format!("v{i}"))).unwrap();
        }
        assert!(save(&mut s, view("one more")).is_err());
        assert!(save(&mut s, view("v0")).is_ok(), "replacing one still works when full");
    }

    /// A study saved before views existed reads with none, and one without views saves as before.
    #[test]
    fn older_studies_read_with_no_views_and_save_unchanged() {
        let s = study();
        let text = serde_json::to_string(&s).unwrap();
        assert!(!text.contains("\"views\""));
        let back: Study = serde_json::from_str(&text).unwrap();
        assert!(back.views.is_empty());
    }
}
