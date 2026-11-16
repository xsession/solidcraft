//! Bounded three-dimensional presentation of one genuine cell slice.
//!
//! This is a planar scientific layer and a wireframe domain, not a volume renderer,
//! CAD body display or particle-path integrator. Display normalization leaves the
//! retained coordinates and original field arrays untouched.
use std::fmt;
use std::sync::Arc;

use egui::{Align2, FontId, Pos2, Sense, Shape, Stroke};
use solidcore_record::Vec3;
use solidcraft_engine::multiphysics::vtk::{BoxGrid, SliceCell3d};
use solidcraft_theme::Tokens;
use solidcraft_viewport::vec::Bounds;
use solidcraft_viewport::{Camera, StandardView};

use crate::SolveApp;
use crate::ui_kit::BackingTooltipExt as _;
use crate::viewport::ViewportState;

const MAX_GLYPHS: usize = 128;
const DISPLAY_SPAN: f64 = 100.0;

#[derive(Clone, Debug)]
struct SliceCache {
    key: (usize, usize, u64),
    cells: Arc<Vec<SliceCell3d>>,
    range: (f64, f64),
}

/// Presentation state only. The parent recreates it for each inspected dataset.
/// Cloning retains the view but never duplicates an in-progress pointer gesture.
pub(crate) struct State3d {
    pub(crate) camera: Camera,
    pub(crate) glyphs: bool,
    /// Set by the parent each frame so 2D and 3D share one colour scale.
    pub(crate) color_map: solidcraft_engine::multiphysics::results::ColorMap,
    pub(crate) color_range: Option<(f64, f64)>,
    initialized: bool,
    viewport: Box<ViewportState>,
    cache: Option<SliceCache>,
}

impl Default for State3d {
    fn default() -> Self {
        Self {
            camera: Camera::default(),
            glyphs: false,
            color_map: Default::default(),
            color_range: None,
            initialized: false,
            viewport: Box::default(),
            cache: None,
        }
    }
}

impl Clone for State3d {
    fn clone(&self) -> Self {
        Self {
            camera: self.camera,
            glyphs: self.glyphs,
            color_map: self.color_map,
            color_range: self.color_range,
            initialized: self.initialized,
            viewport: Box::default(),
            cache: self.cache.clone(),
        }
    }
}

impl fmt::Debug for State3d {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("State3d")
            .field("camera", &self.camera)
            .field("glyphs", &self.glyphs)
            .field("initialized", &self.initialized)
            .finish_non_exhaustive()
    }
}

#[derive(Clone, Copy)]
struct DisplayTransform {
    center: [f64; 3],
    span: f64,
    bounds: [[f64; 2]; 3],
}

impl DisplayTransform {
    fn new(grid: &BoxGrid) -> Option<Self> {
        let mut center = [0.0; 3];
        let mut span = 0.0_f64;
        for (axis, range) in grid.bounds.iter().enumerate() {
            let width = range[1] - range[0];
            if !range[0].is_finite() || !width.is_finite() || width <= 0.0 {
                return None;
            }
            *center.get_mut(axis)? = range[0] + width * 0.5;
            span = span.max(width);
        }
        (span > 0.0 && span.is_finite()).then_some(Self { center, span, bounds: grid.bounds })
    }

    fn point(self, point: [f64; 3]) -> Option<Vec3> {
        // Divide before multiplying: a subnormal producer length must not create
        // an infinite scale factor. This transform is display-only.
        let mut out = [0.0; 3];
        for (axis, value) in point.into_iter().enumerate() {
            let v = (value - *self.center.get(axis)?) / self.span * DISPLAY_SPAN;
            if !v.is_finite() {
                return None;
            }
            *out.get_mut(axis)? = v;
        }
        Some(Vec3::new(*out.first()?, *out.get(1)?, *out.get(2)?))
    }

    fn corners(self) -> Option<[Vec3; 8]> {
        let mut out = [Vec3::ZERO; 8];
        for (index, point) in out.iter_mut().enumerate() {
            let x = self.bounds.first()?.get(usize::from(index & 1 != 0))?;
            let y = self.bounds.get(1)?.get(usize::from(index & 2 != 0))?;
            let z = self.bounds.get(2)?.get(usize::from(index & 4 != 0))?;
            *point = self.point([*x, *y, *z])?;
        }
        Some(out)
    }

    fn camera_bounds(self) -> Option<Bounds> {
        let mut bounds = Bounds::EMPTY;
        for point in self.corners()? {
            bounds.add(point);
        }
        Some(bounds)
    }
}

fn invalid() -> String {
    solidcraft_ui_text::tr("multiphysics.field_3d_invalid").into()
}

fn fit(state: &mut State3d, bounds: &Bounds, aspect: f64) {
    state.camera.fov = 0.0;
    if !state.initialized {
        state.camera.set_view(StandardView::Iso);
    }
    if state.camera.fit_aspect(bounds, aspect) {
        state.initialized = true;
    }
}

/// Keep the existing navigation engine and all its button/modifier bindings,
/// while isolating the data view's camera, gestures and caches from CAD.
fn navigate(app: &mut SolveApp, ui: &mut egui::Ui, state: &mut State3d, rect: egui::Rect, response: &egui::Response) {
    let original_camera = std::mem::replace(&mut app.cam, state.camera);
    let original_animation = app.cam_anim.take();
    let original_pending_fit = app.pending_fit.take();
    let original_locks = (app.ui.view_angles_lock_yaw, app.ui.view_angles_lock_pitch, app.ui.view_angles_lock_roll);
    app.ui.view_angles_lock_yaw = false;
    app.ui.view_angles_lock_pitch = false;
    app.ui.view_angles_lock_roll = false;
    let original_perspective = std::mem::replace(&mut app.ui.perspective, false);
    let original_auto_ortho = std::mem::replace(&mut app.ui.auto_ortho_faces, false);
    let original_pivot = std::mem::replace(&mut app.preferences.orbit_around, crate::prefs::OrbitAround::ViewCentre);
    let original_zoom_cursor = std::mem::replace(&mut app.preferences.zoom_to_cursor, false);
    let original_zoom_surface = std::mem::replace(&mut app.preferences.zoom_stop_at_model, false);
    let original_readout = std::mem::replace(&mut app.preferences.show_view_angles, false);
    // CAD Home/Fit resolves CAD bounds. This data view has an explicit Fit action
    // instead; MiddleDoubleClick::of reads this setting directly for every preset.
    let original_middle = std::mem::replace(&mut app.preferences.middle_double_click, "off".into());
    std::mem::swap(&mut app.viewport, state.viewport.as_mut());
    app.viewport.rect = Some(rect);
    // This view draws its raw-coordinate domain axes, not the Design compass or
    // navigation bar. An explicit empty overlay layout prevents invisible
    // navigation exclusions in the upper-right and lower canvas corners.
    app.viewport.layout = Some(crate::viewport_layout::Layout {
        viewport: rect,
        drawable: rect,
        middle: rect,
        content: rect,
        axis_view: egui::Rect::NOTHING,
        navigation: egui::Rect::NOTHING,
        confirm: egui::Rect::NOTHING,
        usable: None,
    });
    crate::render_workspace::preserve_entry_fit_during(app, |app| {
        crate::navigation::handle(app, ui, rect, response, crate::navigation::NavCaps { orbit: true, pick: false });
    });
    state.camera = app.cam;
    std::mem::swap(&mut app.viewport, state.viewport.as_mut());
    app.cam = original_camera;
    app.cam_anim = original_animation;
    app.pending_fit = original_pending_fit;
    (app.ui.view_angles_lock_yaw, app.ui.view_angles_lock_pitch, app.ui.view_angles_lock_roll) = original_locks;
    app.ui.perspective = original_perspective;
    app.ui.auto_ortho_faces = original_auto_ortho;
    app.preferences.orbit_around = original_pivot;
    app.preferences.zoom_to_cursor = original_zoom_cursor;
    app.preferences.zoom_stop_at_model = original_zoom_surface;
    app.preferences.show_view_angles = original_readout;
    app.preferences.middle_double_click = original_middle;
}

/// Draw a single selected plane and return the number of visible source cells.
/// The parent supplies total remaining height before its common field legend.
pub(crate) fn show(
    app: &mut SolveApp,
    ui: &mut egui::Ui,
    state: &mut State3d,
    grid: &BoxGrid,
    field: usize,
    normal: usize,
    coordinate: f64,
    height: f32,
) -> Result<usize, String> {
    if !height.is_finite() || height <= 0.0 || normal > 2 || !coordinate.is_finite() {
        return Err(invalid());
    }
    let transform = DisplayTransform::new(grid).ok_or_else(invalid)?;
    let bounds = transform.camera_bounds().ok_or_else(invalid)?;
    let corners = transform.corners().ok_or_else(invalid)?;
    let key = (field, normal, coordinate.to_bits());
    if state.cache.as_ref().is_none_or(|cached| cached.key != key) {
        let cells = grid.slice3d(field, normal, coordinate).map_err(|error| error.to_string())?;
        let range = (0..grid.cells.len())
            .filter_map(|cell| grid.value(field, cell))
            .fold(None, |range, value| {
                Some(match range {
                    None => (value, value),
                    Some((low, high)) => (value.min(low), value.max(high)),
                })
            })
            .ok_or_else(invalid)?;
        state.cache = Some(SliceCache { key, cells: Arc::new(cells), range });
    }
    let cached = state.cache.as_ref().ok_or_else(invalid)?.clone();
    let t = Tokens::get();
    let top = ui.next_widget_position().y;
    let mut requested_fit = false;
    ui.horizontal_wrapped(|ui| {
        let response = crate::ui_kit::button(ui, solidcraft_ui_text::tr("multiphysics.field_3d_fit"), crate::ui_kit::Button::Normal);
        requested_fit = response.clicked();
        crate::scenario::publish_target("multiphysics:3d_fit", response.rect);
        response.backing_on_hover_text(solidcraft_ui_text::tr("multiphysics.field_3d_hint"));
        let vector = grid.fields.get(field).is_some_and(|field| field.components == 3);
        let response = ui
            .add_enabled_ui(vector, |ui| crate::ui_kit::checkbox(ui, solidcraft_ui_text::tr("multiphysics.field_3d_glyphs"), &mut state.glyphs))
            .inner;
        crate::scenario::publish_target("multiphysics:3d_glyphs", response.rect);
        response.backing_on_hover_text(solidcraft_ui_text::tr("multiphysics.field_3d_scale"));
    });
    let used = (ui.next_widget_position().y - top).max(0.0);
    let visible_width = (ui.clip_rect().right() - ui.next_widget_position().x).max(1.0);
    let size = egui::vec2(ui.available_width().min(visible_width).max(1.0), (height - used).clamp(64.0, 2048.0));
    let (rect, response) = ui.allocate_exact_size(size, Sense::click_and_drag());
    crate::scenario::publish_rect("multiphysics:field_slice", rect);
    crate::scenario::publish_target("multiphysics:3d_canvas", rect);
    if !state.camera.is_valid() {
        state.camera = Camera::default();
        state.initialized = false;
    }
    if !state.initialized || requested_fit {
        fit(state, &bounds, f64::from(rect.width() / rect.height().max(1.0)));
    }
    navigate(app, ui, state, rect, &response);
    let painter = ui.painter().with_clip_rect(ui.clip_rect().intersect(rect));
    painter.rect_filled(rect, t.control_radius, t.control_fill);
    let projection = crate::viewport::Proj {
        cam: state.camera,
        vp: state.camera.view_proj_bounds(f64::from(rect.width() / rect.height().max(1.0)), &bounds),
        rect,
    };
    let project = |point: Vec3| projection.to_screen(point).filter(|point| point.is_finite());
    let plane_point = transform
        .point({
            let mut p = transform.center;
            if let Some(value) = p.get_mut(normal) {
                *value = coordinate;
            }
            p
        })
        .ok_or_else(invalid)?;
    let component = |point: Vec3| match normal {
        0 => point.x,
        1 => point.y,
        _ => point.z,
    };
    let plane = component(plane_point);
    let side = component(state.camera.eye()) - plane;
    let mut behind = Vec::new();
    let mut front = Vec::new();
    for (index, a) in corners.iter().enumerate() {
        for axis in 0..3 {
            let other = index ^ (1 << axis);
            if other <= index {
                continue;
            }
            let Some(b) = corners.get(other) else { continue };
            let da = component(*a) - plane;
            let db = component(*b) - plane;
            if da * db < 0.0 {
                let ratio = da / (da - db);
                let cross = *a + (*b - *a) * ratio;
                for (x, y) in [(*a, cross), (cross, *b)] {
                    let far = (component((x + y) * 0.5) - plane) * side < 0.0;
                    if far {
                        behind.push((x, y));
                    } else {
                        front.push((x, y));
                    }
                }
            } else if (da + db) * side < 0.0 {
                behind.push((*a, *b));
            } else {
                front.push((*a, *b));
            }
        }
    }
    let draw_line = |a: Vec3, b: Vec3| {
        if let (Some(a), Some(b)) = (project(a), project(b)) {
            painter.line_segment([a, b], Stroke::new(1.0, t.text_dim));
        }
    };
    for (a, b) in behind {
        draw_line(a, b);
    }
    let (low, high) = state.color_range.unwrap_or(cached.range);
    let scale = low.abs().max(high.abs()).max(1.0);
    let mut visible = 0;
    for cell in cached.cells.iter() {
        let points: Option<Vec<Pos2>> = cell.corners.iter().map(|&p| transform.point(p).and_then(project)).collect();
        let Some(points) = points else { continue };
        if !egui::Rect::from_points(&points).intersects(rect) {
            continue;
        }
        let fraction = if high > low { ((cell.value / scale - low / scale) / (high / scale - low / scale)).clamp(0.0, 1.0) as f32 } else { 0.5 };
        let [r, g, b] = state.color_map.color(f64::from(fraction));
        let color = egui::Color32::from_rgb(r, g, b);
        let mut mesh = egui::Mesh::default();
        for point in points {
            mesh.colored_vertex(point, color);
        }
        mesh.add_triangle(0, 1, 2);
        mesh.add_triangle(0, 2, 3);
        painter.add(Shape::mesh(mesh));
        visible += 1;
    }
    for (a, b) in front {
        draw_line(a, b);
    }
    let mut glyphs = 0;
    // Glyph length follows the data, not a user-locked colour range.
    let high = cached.range.1;
    if state.glyphs && high > 0.0 && grid.fields.get(field).is_some_and(|field| field.components == 3) {
        let stride = cached.cells.len().div_ceil(MAX_GLYPHS).max(1);
        for cell in cached.cells.iter().step_by(stride) {
            let Some(vector) = cell.vector else { continue };
            let Some(start) = transform.point(cell.center) else { continue };
            let end = start + Vec3::new(vector[0] / high, vector[1] / high, vector[2] / high) * 6.0;
            let (Some(a), Some(b)) = (project(start), project(end)) else { continue };
            if !rect.contains(a) {
                continue;
            }
            let delta = b - a;
            let length = delta.length();
            if !length.is_finite() || length < 1.0 {
                continue;
            }
            let direction = delta / length;
            let perpendicular = egui::vec2(-direction.y, direction.x);
            let head = length.min(t.space_4);
            let stroke = Stroke::new(1.0, t.text);
            painter.line_segment([a, b], stroke);
            painter.line_segment([b, b - direction * head + perpendicular * head * 0.5], stroke);
            painter.line_segment([b, b - direction * head - perpendicular * head * 0.5], stroke);
            glyphs += 1;
        }
    }
    crate::scenario::publish_count("multiphysics:3d_glyph_count", glyphs as f64);
    if let Some(origin) = corners.first().copied() {
        for (index, label) in [(1, "X"), (2, "Y"), (4, "Z")] {
            if let Some(point) = corners.get(index).copied().and_then(project) {
                painter.text(point, Align2::LEFT_BOTTOM, label, FontId::proportional(t.font_small), t.text_dim);
            }
        }
        if let Some(point) = project(origin) {
            crate::scenario::publish_rect("multiphysics:3d_origin", egui::Rect::from_center_size(point, egui::Vec2::ZERO));
        }
    }
    Ok(visible)
}
