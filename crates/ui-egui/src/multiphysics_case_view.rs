//! The Multiphysics graphics window for a 3-D airflow case: the computational domain, the
//! parts in it and the role of each domain face, drawn from the case being edited.
//!
//! It works the way a flow-analysis add-in does in the CAD window:
//! - Computational Domain shows an arrow on each of the six domain faces. Dragging one moves
//!   that face outward or inward by whole cells; the panel's buttons do the same one cell at
//!   a time.
//! - Solid materials and the thermal devices pick parts (obstacles, porous blocks, surface
//!   bodies). Click picks one, Ctrl+click adds or removes, clicking empty space clears.
//! - Boundary conditions and fans pick domain faces. Clicking again at the same spot steps to
//!   the next thing under the pointer (the face behind the model, for instance).
//!
//! Every change goes through `case_scene` and lands in the project's draft, never straight
//! into the document: Apply in the panel commits it, as with any other edit.
use solidcraft_engine::geom::Vec3Geom as _;
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use egui::{Align2, Color32, FontId, Pos2, Sense, Shape, Stroke};
use serde_json::Value;
use solidcore_record::Vec3;
use solidcraft_engine::multiphysics::case_scene::{self as cs, PATCHES, PartKind, PatchRole, Scene, Target};
use solidcraft_theme::Tokens;
use solidcraft_viewport::vec::Bounds;

use crate::SolveApp;
use crate::multiphysics_field_view::State3d;
use crate::ui_kit::BackingTooltipExt as _;

/// The display frame spans this many units, like the results view.
const DISPLAY_SPAN: f64 = 100.0;
/// Arrow length as a fraction of the display span.
const ARROW: f64 = 0.16;
/// How close (in points) the pointer must be to an arrow to grab it.
const GRAB: f32 = 9.0;

/// What clicking in the view picks.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Mode {
    /// Domain arrows; clicks pick nothing.
    Domain,
    /// Obstacles, porous blocks and surface bodies.
    Parts,
    /// The six domain faces.
    Faces,
}

/// One picked thing.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) enum Pick {
    Part(String),
    /// A domain face, as a [`PATCHES`] index.
    Patch(usize),
}

impl Pick {
    pub(crate) fn label(&self) -> String {
        match self {
            Self::Part(name) => name.clone(),
            Self::Patch(index) => PATCHES.get(*index).copied().unwrap_or("?").to_string(),
        }
    }
}

/// The display frame, fixed when the view is fitted so edits do not make the picture jump.
#[derive(Clone, Copy, Debug)]
struct Frame {
    center: [f64; 3],
    span: f64,
}

#[derive(Clone, Copy, Debug)]
struct Drag {
    patch: usize,
    start: Pos2,
    /// Screen motion for one metre along the face's outward normal.
    per_metre: egui::Vec2,
    cells: i64,
}

#[derive(Default)]
pub(crate) struct CaseView {
    view: State3d,
    frame: Option<Frame>,
    /// Added to case coordinates before display. A minimum face moved outward shifts the
    /// case frame; this keeps the parts still on screen.
    shift: [f64; 3],
    scene: Option<(u64, Result<Arc<Scene>, String>)>,
    drag: Option<Drag>,
    /// The last click and which of the hits under it was taken (for stepping through them).
    cycle: Option<(Pos2, usize)>,
    pub(crate) selection: Vec<Pick>,
    /// Show the CAD model instead of the case.
    pub(crate) cad: bool,
    /// Draw the basic mesh: the grid's cell lines on the domain faces.
    pub(crate) mesh: bool,
    /// Draw the mesh regardless (set each frame while the Mesh node is open).
    pub(crate) force_mesh: bool,
    /// Sealed fluid pockets to outline (bounds in metres), and whether to draw them (set each
    /// frame while the Computational Domain node is open).
    pub(crate) pockets: Vec<([f64; 3], [f64; 3])>,
    pub(crate) show_pockets: bool,
    /// The project the selection and frame belong to.
    owner: Option<String>,
}

impl CaseView {
    /// Keep the selection only while the same project is shown.
    pub(crate) fn own(&mut self, project: &str) {
        if self.owner.as_deref() != Some(project) {
            self.reset();
            self.owner = Some(project.to_string());
        }
    }

    /// Forget the selection and the frame (another project was picked).
    pub(crate) fn reset(&mut self) {
        self.selection.clear();
        self.frame = None;
        self.shift = [0.0; 3];
        self.drag = None;
        self.cycle = None;
    }

    pub(crate) fn parts(&self) -> Vec<String> {
        self.selection.iter().filter_map(|pick| if let Pick::Part(name) = pick { Some(name.clone()) } else { None }).collect()
    }

    pub(crate) fn patches(&self) -> Vec<usize> {
        self.selection.iter().filter_map(|pick| if let Pick::Patch(index) = pick { Some(*index) } else { None }).collect()
    }

    /// A minimum face moved by `cells`: keep the parts where they were on screen.
    pub(crate) fn moved(&mut self, patch: usize, cells: i64, cell: f64) {
        if patch.is_multiple_of(2)
            && let Some(slot) = self.shift.get_mut(patch / 2)
        {
            *slot += cells as f64 * cell;
        }
    }
}

/// What the view asks the workspace to do with the draft.
#[derive(Default)]
pub(crate) struct Outcome {
    /// A domain drag finished: (patch, cells).
    pub(crate) move_face: Option<(usize, i64)>,
    /// A click changed the selection this frame.
    pub(crate) picked: bool,
}

/// The scene of `case_text`, cached by its text.
fn scene(view: &mut CaseView, case_text: &str) -> Result<Arc<Scene>, String> {
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    case_text.hash(&mut hasher);
    let key = hasher.finish();
    if let Some((cached, scene)) = &view.scene
        && *cached == key
    {
        return scene.clone();
    }
    let built = serde_json::from_str::<Value>(case_text)
        .map_err(|error| error.to_string())
        .and_then(|case| cs::scene(&case).map(Arc::new).map_err(|error| error.to_string()));
    view.scene = Some((key, built.clone()));
    built
}

/// Whether `case_text` is a case this view can draw.
pub(crate) fn drawable(view: &mut CaseView, case_text: &str) -> bool {
    scene(view, case_text).is_ok()
}

fn tr(key: &str) -> &str {
    solidcraft_ui_text::tr(key)
}

fn role_key(role: PatchRole) -> &'static str {
    match role {
        PatchRole::Inlet => "multiphysics.case_role_inlet",
        PatchRole::Outlet => "multiphysics.case_role_outlet",
        PatchRole::Wall => "multiphysics.case_role_wall",
        PatchRole::Symmetry => "multiphysics.case_role_symmetry",
        PatchRole::Fan => "multiphysics.case_role_fan",
        PatchRole::Default => "multiphysics.case_role_default",
    }
}

fn role_color(t: &Tokens, role: PatchRole) -> Color32 {
    match role {
        PatchRole::Inlet => t.surface_accent,
        PatchRole::Outlet => t.form_success,
        PatchRole::Wall => t.text_dim,
        PatchRole::Symmetry => t.sketch_accent,
        PatchRole::Fan => t.render_accent,
        PatchRole::Default => t.border,
    }
}

/// The four corners of box face `face` (a [`PATCHES`] index), wound around the face.
fn box_face(min: [f64; 3], max: [f64; 3], face: usize) -> [[f64; 3]; 4] {
    let axis = face / 2;
    let fixed = if face.is_multiple_of(2) { min[axis] } else { max[axis] };
    let [u, v] = cs::face_axes(axis);
    let corner = |a: f64, b: f64| {
        let mut p = [0.0; 3];
        p[axis] = fixed;
        p[u] = a;
        p[v] = b;
        p
    };
    [corner(min[u], min[v]), corner(max[u], min[v]), corner(max[u], max[v]), corner(min[u], max[v])]
}

fn face_center(min: [f64; 3], max: [f64; 3], face: usize) -> [f64; 3] {
    let axis = face / 2;
    let mut c = [0, 1, 2].map(|a| (min[a] + max[a]) * 0.5);
    c[axis] = if face.is_multiple_of(2) { min[axis] } else { max[axis] };
    c
}

fn outward(face: usize) -> [f64; 3] {
    let mut n = [0.0; 3];
    n[face / 2] = if face.is_multiple_of(2) { -1.0 } else { 1.0 };
    n
}

/// Distance from `p` to segment `a`–`b`.
fn to_segment(p: Pos2, a: Pos2, b: Pos2) -> f32 {
    let ab = b - a;
    let len = ab.length_sq();
    let s = if len > 0.0 { ((p - a).dot(ab) / len).clamp(0.0, 1.0) } else { 0.0 };
    (p - (a + ab * s)).length()
}

struct Polygon {
    depth: f64,
    points: Vec<Pos2>,
    fill: Color32,
    stroke: Stroke,
}

/// Draw the case and handle picking and domain drags. `None` when the case cannot be drawn
/// (the caller shows the CAD window instead).
pub(crate) fn show(app: &mut SolveApp, ui: &mut egui::Ui, view: &mut CaseView, case_text: &str, mode: Mode) -> Option<Outcome> {
    let t = Tokens::get();
    let scene = scene(view, case_text).ok()?;
    let mut outcome = Outcome::default();
    let mut refit = false;
    ui.horizontal(|ui| {
        let response = crate::ui_kit::button(ui, tr("multiphysics.case_fit"), crate::ui_kit::Button::Normal);
        crate::scenario::publish_target("multiphysics:case_fit", response.rect);
        refit = response.clicked();
        let response = crate::ui_kit::checkbox(ui, tr("multiphysics.case_show_cad"), &mut view.cad);
        crate::scenario::publish_target("multiphysics:case_show_cad", response.rect);
        let response = crate::ui_kit::checkbox(ui, tr("multiphysics.tree_show_basic_mesh"), &mut view.mesh);
        crate::scenario::publish_target("multiphysics:case_show_mesh", response.rect);
        if scene.truncated {
            ui.small(tr("multiphysics.case_truncated"));
        }
    });
    let size = ui.available_size().max(egui::vec2(64.0, 64.0));
    let (rect, response) = ui.allocate_exact_size(size, Sense::click_and_drag());
    crate::scenario::publish_target("multiphysics:case_canvas", rect);
    if view.frame.is_none() || refit {
        view.frame = Some(Frame { center: scene.lengths.map(|l| l * 0.5), span: scene.lengths.iter().copied().fold(0.0, f64::max).max(1e-9) });
        view.shift = [0.0; 3];
        refit = true;
    }
    let frame = view.frame?;
    let shift = view.shift;
    let display = move |p: [f64; 3]| {
        Vec3::new(
            (p[0] - shift[0] - frame.center[0]) / frame.span * DISPLAY_SPAN,
            (p[1] - shift[1] - frame.center[1]) / frame.span * DISPLAY_SPAN,
            (p[2] - shift[2] - frame.center[2]) / frame.span * DISPLAY_SPAN,
        )
    };
    let to_case = move |p: Vec3| {
        [
            p.x / DISPLAY_SPAN * frame.span + frame.center[0] + shift[0],
            p.y / DISPLAY_SPAN * frame.span + frame.center[1] + shift[1],
            p.z / DISPLAY_SPAN * frame.span + frame.center[2] + shift[2],
        ]
    };
    let mut bounds = Bounds::EMPTY;
    for corner in 0..8 {
        let p = [0, 1, 2].map(|axis| if corner & (1 << axis) != 0 { scene.lengths[axis] } else { 0.0 });
        bounds.add(display(p));
    }
    // Margin for the arrows.
    let pad = Vec3::new(DISPLAY_SPAN * ARROW * 1.5, DISPLAY_SPAN * ARROW * 1.5, DISPLAY_SPAN * ARROW * 1.5);
    bounds.add(bounds.min - pad);
    bounds.add(bounds.max + pad);
    // Arrows from last frame's camera decide whether the pointer is on a handle before the
    // navigation engine sees the press.
    let pointer = response.hover_pos();
    let arrows = |proj: &crate::viewport::Proj| -> Vec<(usize, Pos2, Pos2)> {
        (0..6)
            .filter_map(|face| {
                let c = face_center([0.0; 3], scene.lengths, face);
                let n = outward(face);
                let tip = [0, 1, 2].map(|a| c[a] + n[a] * frame.span * ARROW);
                Some((face, proj.to_screen(display(c))?, proj.to_screen(display(tip))?))
            })
            .collect()
    };
    let before = crate::multiphysics_field_view::projection(&view.view, rect, &bounds);
    let grab_at = |p: Pos2| arrows(&before).into_iter().find(|(_, a, b)| to_segment(p, *a, *b) < GRAB).map(|(face, ..)| face);
    let domain = mode == Mode::Domain;
    let grabbed = pointer.filter(|_| domain).and_then(grab_at);
    // The arrow under the press, so a drag that leaves the arrow quickly still holds it.
    let press = ui.input(|i| i.pointer.press_origin());
    let pressed_on = press.filter(|_| domain).and_then(grab_at);
    let pressing = ui.input(|i| i.pointer.primary_down());
    let navigate_now = view.drag.is_none() && !(pressed_on.is_some() && pressing);
    let proj = crate::multiphysics_field_view::camera(app, ui, &mut view.view, rect, &response, &bounds, refit, navigate_now);
    let eye = view.view.camera.eye();
    let project = |p: [f64; 3]| proj.to_screen(display(p)).filter(|p| p.is_finite());
    let depth = |points: &[[f64; 3]]| {
        let n = points.len().max(1) as f64;
        let mut c = [0.0; 3];
        for p in points {
            for a in 0..3 {
                c[a] += p[a] / n;
            }
        }
        let d = display(c) - eye;
        d.dot(d)
    };

    // ---- picking ----
    let hits_at = |p: Pos2| {
        let (o, d) = proj.ray(p);
        let origin = to_case(o);
        let direction = [d.x, d.y, d.z];
        cs::pick(&scene, origin, direction)
            .into_iter()
            .filter(|hit| match (&hit.target, mode) {
                (Target::Part { .. }, Mode::Parts) | (Target::Patch { .. }, Mode::Faces) => true,
                _ => false,
            })
            .collect::<Vec<_>>()
    };
    let as_pick = |target: &Target| match target {
        Target::Part { name, .. } => Pick::Part(name.clone()),
        Target::Patch { index } => Pick::Patch(*index),
    };
    let hover = match (pointer, view.drag) {
        (Some(p), None) if mode != Mode::Domain => hits_at(p).first().map(|hit| as_pick(&hit.target)),
        _ => None,
    };
    if response.clicked()
        && mode != Mode::Domain
        && let Some(p) = response.interact_pointer_pos()
    {
        let hits = hits_at(p);
        let index = match view.cycle {
            Some((last, index)) if last.distance(p) < 4.0 && !hits.is_empty() => (index + 1) % hits.len(),
            _ => 0,
        };
        view.cycle = Some((p, index));
        let before_click = view.selection.clone();
        let add = ui.input(|i| i.modifiers.command);
        match hits.get(index).map(|hit| as_pick(&hit.target)) {
            Some(pick) if add => {
                if let Some(at) = view.selection.iter().position(|old| *old == pick) {
                    view.selection.remove(at);
                } else {
                    view.selection.push(pick);
                }
            }
            Some(pick) => view.selection = vec![pick],
            None if !add => view.selection.clear(),
            None => {}
        }
        outcome.picked = view.selection != before_click;
    }

    // ---- domain drag ----
    if mode == Mode::Domain {
        if response.drag_started()
            && let (Some(face), Some(start)) = (pressed_on, press)
        {
            let c = face_center([0.0; 3], scene.lengths, face);
            let n = outward(face);
            let probe = 0.01 * frame.span;
            let ahead = [0, 1, 2].map(|a| c[a] + n[a] * probe);
            if let (Some(a), Some(b)) = (project(c), project(ahead)) {
                let per_metre = (b - a) / probe as f32;
                if per_metre.length_sq() > 1e-6 {
                    view.drag = Some(Drag { patch: face, start, per_metre, cells: 0 });
                }
            }
        }
        if let Some(drag) = view.drag.as_mut()
            && let Some(p) = response.interact_pointer_pos().or(pointer)
        {
            let metres = f64::from((p - drag.start).dot(drag.per_metre) / drag.per_metre.length_sq());
            let cell = scene.edge_cell(drag.patch);
            let count = scene.dimensions.get(drag.patch / 2).copied().unwrap_or(1) as i64;
            if cell > 0.0 && metres.is_finite() {
                drag.cells = ((metres / cell).round() as i64).clamp(1 - count, cs::MAX_AXIS_CELLS as i64 - count);
            }
        }
        if response.drag_stopped()
            && let Some(drag) = view.drag.take()
            && drag.cells != 0
        {
            outcome.move_face = Some((drag.patch, drag.cells));
        }
    } else {
        view.drag = None;
    }

    // ---- drawing ----
    let painter = ui.painter_at(rect);
    painter.rect_filled(rect, t.control_radius, t.control_fill);
    let selected = |pick: &Pick| view.selection.contains(pick);
    let mut polygons = Vec::new();
    let mut push = |points: &[[f64; 3]], fill: Color32, stroke: Stroke| {
        let projected: Option<Vec<Pos2>> = points.iter().map(|p| project(*p)).collect();
        if let Some(projected) = projected {
            polygons.push(Polygon { depth: depth(points), points: projected, fill, stroke });
        }
    };
    for (face, role) in scene.patches.iter().enumerate() {
        let pick = Pick::Patch(face);
        let (fill, stroke) = if selected(&pick) {
            (t.sel_face.gamma_multiply(0.55), Stroke::new(2.0, t.sel_edge))
        } else if hover.as_ref() == Some(&pick) {
            (t.hover.gamma_multiply(0.45), Stroke::new(1.5, t.hover_edge))
        } else {
            (role_color(&t, *role).gamma_multiply(if mode == Mode::Faces { 0.22 } else { 0.08 }), Stroke::new(1.0, t.text_dim))
        };
        push(box_face([0.0; 3], scene.lengths, face).as_slice(), fill, stroke);
    }
    for part in &scene.parts {
        let pick = Pick::Part(part.name.clone());
        let base = match part.kind {
            PartKind::Porous => t.surface_soft,
            _ if part.solid => t.warning,
            PartKind::Body => t.form_fill,
            PartKind::Obstacle => t.text_dim,
        };
        let (fill, stroke) = if selected(&pick) {
            (t.sel_face, Stroke::new(1.5, t.sel_edge))
        } else if hover.as_ref() == Some(&pick) {
            (t.hover, Stroke::new(1.5, t.hover_edge))
        } else {
            (base.gamma_multiply(if mode == Mode::Faces { 0.35 } else { 0.8 }), Stroke::new(1.0, t.text.gamma_multiply(0.6)))
        };
        if part.kind == PartKind::Body {
            for triangle in &part.triangles {
                push(triangle.as_slice(), fill, Stroke::NONE);
            }
        } else {
            for face in 0..6 {
                push(box_face(part.min, part.max, face).as_slice(), fill, stroke);
            }
        }
    }
    polygons.sort_by(|a, b| b.depth.total_cmp(&a.depth));
    for polygon in polygons {
        if polygon.points.iter().any(|p| rect.expand(rect.width()).contains(*p)) {
            for shape in crate::viewport::projected_polygon(&polygon.points, polygon.fill, polygon.stroke) {
                painter.add(shape);
            }
        }
    }
    // Basic mesh: every cell plane where it meets the domain faces (thinned to 64 per axis).
    if view.mesh || view.force_mesh {
        let stroke = Stroke::new(0.5, t.text_dim.gamma_multiply(0.6));
        let mut lines = 0usize;
        for face in 0..6usize {
            let axis = face / 2;
            let fixed = if face.is_multiple_of(2) { 0.0 } else { scene.lengths[axis] };
            let [u, v] = cs::face_axes(axis);
            for (along, across) in [(u, v), (v, u)] {
                let nodes = &scene.nodes[along];
                let stride = nodes.len().div_ceil(64).max(1);
                for &at in nodes.iter().step_by(stride) {
                    let mut a = [0.0; 3];
                    a[axis] = fixed;
                    a[along] = at;
                    let mut b = a;
                    b[across] = scene.lengths[across];
                    if let (Some(a), Some(b)) = (project(a), project(b)) {
                        painter.line_segment([a, b], stroke);
                        lines += 1;
                    }
                }
            }
        }
        crate::scenario::publish_count("multiphysics:basic_mesh_lines", lines as f64);
    }
    // Face labels with their role, so a wall and an outlet can be told apart without colour.
    for (face, role) in scene.patches.iter().enumerate() {
        if let Some(at) = project(face_center([0.0; 3], scene.lengths, face)) {
            let name = PATCHES.get(face).copied().unwrap_or("");
            let text = format!("{name} · {}", tr(role_key(*role)));
            painter.text(at, Align2::CENTER_CENTER, text, FontId::proportional(t.font_small), t.text);
            crate::scenario::publish_target(&format!("multiphysics:case_face:{name}"), egui::Rect::from_center_size(at, egui::vec2(GRAB, GRAB)));
        }
    }
    // Sealed fluid pockets found by the fluid volume check: warning-coloured wire boxes.
    if view.show_pockets {
        for (n, (min, max)) in view.pockets.iter().enumerate() {
            let stroke = Stroke::new(2.0, t.warning);
            let corner = |i: usize| [0, 1, 2].map(|a| if i >> a & 1 == 0 { min[a] } else { max[a] });
            for i in 0..8usize {
                for a in 0..3 {
                    if i >> a & 1 == 0
                        && let (Some(p), Some(q)) = (project(corner(i)), project(corner(i | 1 << a)))
                    {
                        painter.line_segment([p, q], stroke);
                    }
                }
            }
            if let Some(at) = project([0, 1, 2].map(|a| 0.5 * (min[a] + max[a]))) {
                crate::scenario::publish_target(
                    &format!("multiphysics:fluid_pocket:{}", n + 1),
                    egui::Rect::from_center_size(at, egui::vec2(GRAB, GRAB)),
                );
            }
        }
    }
    // Internal fans: a wire box and an arrow through it along the blowing direction.
    for fan in &scene.fans {
        let stroke = Stroke::new(1.5, t.accent);
        let corner = |i: usize| [0, 1, 2].map(|a| if i >> a & 1 == 0 { fan.min[a] } else { fan.max[a] });
        for i in 0..8usize {
            for a in 0..3 {
                if i >> a & 1 == 0
                    && let (Some(p), Some(q)) = (project(corner(i)), project(corner(i | 1 << a)))
                {
                    painter.line_segment([p, q], stroke);
                }
            }
        }
        let centre = [0, 1, 2].map(|a| 0.5 * (fan.min[a] + fan.max[a]));
        let half = 0.5 * (0..3).map(|a| (fan.max[a] - fan.min[a]) * fan.direction[a].abs()).sum::<f64>() + 0.5 * scene.cell(0).min(scene.cell(1));
        let end = |sign: f64| [0, 1, 2].map(|a| centre[a] + sign * half * fan.direction[a]);
        if let (Some(a), Some(b)) = (project(end(-1.0)), project(end(1.0))) {
            let delta = b - a;
            let length = delta.length().max(1.0);
            let dir = delta / length;
            let side = egui::vec2(-dir.y, dir.x);
            let head = (length * 0.3).clamp(6.0, 14.0);
            painter.line_segment([a, b], Stroke::new(2.0, t.accent));
            painter.add(Shape::convex_polygon(
                vec![b, b - dir * head + side * head * 0.45, b - dir * head - side * head * 0.45],
                t.accent,
                Stroke::NONE,
            ));
            painter.text(b, Align2::LEFT_BOTTOM, &fan.name, FontId::proportional(t.font_small), t.text);
            crate::scenario::publish_target(
                &format!("multiphysics:internal_fan:{}", fan.name),
                egui::Rect::from_center_size(b, egui::vec2(GRAB, GRAB)),
            );
        }
    }
    // Rotating regions: a wire box, the axis through it and the speed.
    for region in &scene.rotating {
        let stroke = Stroke::new(1.5, t.accent);
        let corner = |i: usize| [0, 1, 2].map(|a| if i >> a & 1 == 0 { region.min[a] } else { region.max[a] });
        for i in 0..8usize {
            for a in 0..3 {
                if i >> a & 1 == 0
                    && let (Some(p), Some(q)) = (project(corner(i)), project(corner(i | 1 << a)))
                {
                    painter.line_segment([p, q], stroke);
                }
            }
        }
        let norm = region.axis.iter().map(|c| c * c).sum::<f64>().sqrt();
        if norm > 0.0 {
            let axis = region.axis.map(|c| c / norm);
            // Half the box's extent along the axis, either side of the axis point.
            let half = 0.5 * (0..3).map(|a| (region.max[a] - region.min[a]) * axis[a].abs()).sum::<f64>();
            let end = |sign: f64| [0, 1, 2].map(|a| region.axis_point[a] + sign * half * axis[a]);
            if let (Some(p), Some(q)) = (project(end(-1.0)), project(end(1.0))) {
                painter.line_segment([p, q], Stroke::new(2.0, t.accent));
                let text = format!("{} · {} rpm", region.name, region.rpm);
                painter.text(q, Align2::LEFT_BOTTOM, text, FontId::proportional(t.font_small), t.text);
                crate::scenario::publish_target(
                    &format!("multiphysics:rotating_region:{}", region.name),
                    egui::Rect::from_center_size(q, egui::vec2(GRAB, GRAB)),
                );
            }
        }
    }
    for part in &scene.parts {
        if let Some(at) = project([0, 1, 2].map(|a| (part.min[a] + part.max[a]) * 0.5)) {
            crate::scenario::publish_target(
                &format!("multiphysics:case_part:{}", part.name),
                egui::Rect::from_center_size(at, egui::vec2(GRAB, GRAB)),
            );
            if mode == Mode::Parts && (selected(&Pick::Part(part.name.clone())) || hover == Some(Pick::Part(part.name.clone()))) {
                painter.text(at, Align2::CENTER_BOTTOM, &part.name, FontId::proportional(t.font_small), t.text);
            }
        }
    }
    if mode == Mode::Domain {
        let live = arrows(&proj);
        for (face, a, b) in &live {
            let active = view.drag.is_some_and(|d| d.patch == *face) || (view.drag.is_none() && grabbed == Some(*face));
            let color = if active { t.accent } else { t.text };
            let delta = *b - *a;
            let length = delta.length().max(1.0);
            let dir = delta / length;
            let side = egui::vec2(-dir.y, dir.x);
            let head = (length * 0.3).clamp(6.0, 14.0);
            painter.line_segment([*a, *b], Stroke::new(if active { 3.0 } else { 2.0 }, color));
            painter.add(Shape::convex_polygon(
                vec![*b, *b - dir * head + side * head * 0.45, *b - dir * head - side * head * 0.45],
                color,
                Stroke::NONE,
            ));
            let name = PATCHES.get(*face).copied().unwrap_or("");
            crate::scenario::publish_target(&format!("multiphysics:domain_arrow:{name}"), egui::Rect::from_two_pos(*a, *b).expand(GRAB));
            crate::scenario::publish_target(
                &format!("multiphysics:domain_arrow_base:{name}"),
                egui::Rect::from_center_size(*a, egui::vec2(2.0, 2.0)),
            );
            crate::scenario::publish_target(
                &format!("multiphysics:domain_arrow_tip:{name}"),
                egui::Rect::from_center_size(*b - dir * head * 0.5, egui::vec2(2.0, 2.0)),
            );
        }
        if let Some(drag) = view.drag {
            let offset = drag.cells as f64 * scene.edge_cell(drag.patch);
            let n = outward(drag.patch);
            let ghost: Vec<[f64; 3]> =
                box_face([0.0; 3], scene.lengths, drag.patch).iter().map(|p| [0, 1, 2].map(|a| p[a] + n[a] * offset)).collect();
            let points: Option<Vec<Pos2>> = ghost.iter().map(|p| project(*p)).collect();
            if let Some(points) = points {
                for shape in crate::viewport::projected_polygon(&points, t.accent_soft.gamma_multiply(0.5), Stroke::new(1.5, t.accent)) {
                    painter.add(shape);
                }
            }
            if let Some(p) = pointer {
                let text = solidcraft_ui_text::format(
                    "multiphysics.case_drag_cells",
                    &[("cells", &format!("{:+}", drag.cells)), ("metres", &format!("{:.4}", offset))],
                );
                painter.text(p + egui::vec2(12.0, -12.0), Align2::LEFT_BOTTOM, text, FontId::proportional(t.font_small), t.text);
            }
            crate::scenario::publish_count("multiphysics:domain_drag_cells", drag.cells as f64);
        }
    }
    let hint = tr(match mode {
        Mode::Domain => "multiphysics.case_hint_domain",
        Mode::Parts => "multiphysics.case_hint_parts",
        Mode::Faces => "multiphysics.case_hint_faces",
    });
    let galley = painter.layout(hint.to_owned(), FontId::proportional(t.font_small), t.text_dim, rect.width() - 2.0 * t.space_12);
    painter.galley(egui::pos2(rect.left() + t.space_12, rect.bottom() - t.space_8 - galley.size().y), galley, t.text_dim);
    crate::scenario::publish_count("multiphysics:case_selection", view.selection.len() as f64);
    if let Some(pick) = &hover {
        let _ = response.clone().backing_on_hover_text(pick.label());
    }
    Some(outcome)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn box_faces_wind_on_their_plane_and_arrows_point_out() {
        let (min, max) = ([0.0, 0.0, 0.0], [1.0, 2.0, 3.0]);
        for face in 0..6 {
            let corners = box_face(min, max, face);
            let axis = face / 2;
            let plane = if face.is_multiple_of(2) { min[axis] } else { max[axis] };
            assert!(corners.iter().all(|c| c[axis] == plane));
            let c = face_center(min, max, face);
            assert_eq!(c[axis], plane);
            assert_eq!(outward(face)[axis], if face.is_multiple_of(2) { -1.0 } else { 1.0 });
        }
        assert!((to_segment(egui::pos2(5.0, 3.0), egui::pos2(0.0, 0.0), egui::pos2(10.0, 0.0)) - 3.0).abs() < 1e-6);
        assert!((to_segment(egui::pos2(-4.0, 3.0), egui::pos2(0.0, 0.0), egui::pos2(10.0, 0.0)) - 5.0).abs() < 1e-6);
    }
}
