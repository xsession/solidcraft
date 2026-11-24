//! Results of a study of the loaded model drawn on the model in the design's own 3-D view: the
//! cut plot on its plane, the surface plot (each wetted face coloured by the fluid cell next to
//! it), velocity arrows on the cut plot, a flow tool's trajectories or isosurface, and a probe
//! where the user clicks. The model and its camera are the design's; only the
//! result is painted over them. The display mesh is used for display and for where a click
//! lands, never to decide the model.

use serde_json::json;
use solidcore_record::Vec3;
use solidcraft_engine::multiphysics::vtk::BoxGrid;
use solidcraft_theme::Tokens;

use super::{ViewState, post, scale_color, tr};
use crate::SolveApp;

/// More bins than this and the lookup grid gets coarser.
const MAX_BINS: f64 = 2_097_152.0;
/// More display triangles than this and the surface plot is not drawn.
const MAX_TRIANGLES: usize = 400_000;
/// Velocity arrows drawn on a cut plot at most (cells are skipped evenly beyond this).
const MAX_ARROWS: usize = 600;
/// Trajectory segments and isosurface triangles drawn per frame at most.
const MAX_SEGMENTS: usize = 80_000;
const MAX_OVERLAY_TRIANGLES: usize = 60_000;

/// A display triangle of the model in the design (mm) and its unit outward normal.
type Triangle = ([[f64; 3]; 3], [f64; 3]);

/// The model's display triangles for a result: the wetted ones with the fluid cell next to each.
#[derive(Debug)]
pub(super) struct Surface {
    key: (usize, usize, usize),
    wet: Vec<(Triangle, usize)>,
    dry: Vec<Triangle>,
    too_many: bool,
}

/// Where the user probed: the point in the design (mm), in the case (m) and its fluid cell.
#[derive(Clone, Debug, PartialEq)]
pub(super) struct Probe {
    pub at_mm: [f64; 3],
    pub case_m: [f64; 3],
    pub cell: Option<usize>,
}

/// What is drawn on the model.
pub(super) enum Plot {
    /// The cut plot: plane normal axis, its in-plane axes and its coordinate (case m).
    Slice { axis: usize, along: (usize, usize), coordinate: f64 },
    /// The surface plot on the wetted faces.
    Surfaces,
    /// Only the flow tool's trajectories or isosurface.
    Flow,
}

/// The result's velocity field: a three-component field, preferring one named for velocity.
pub(super) fn vector_field(grid: &BoxGrid) -> Option<usize> {
    let vectors = || grid.fields.iter().enumerate().filter(|(_, f)| f.components == 3);
    vectors().find(|(_, f)| f.name.contains("velocity")).or_else(|| vectors().next()).map(|(index, _)| index)
}

fn vector(grid: &BoxGrid, field: usize, cell: usize) -> Option<[f64; 3]> {
    let start = cell.checked_mul(3)?;
    let tuple = grid.fields.get(field)?.values.get(start..start.checked_add(3)?)?;
    let v = [tuple[0], tuple[1], tuple[2]];
    v.iter().all(|c| c.is_finite()).then_some(v)
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn add_scaled(a: [f64; 3], b: [f64; 3], k: f64) -> [f64; 3] {
    [a[0] + b[0] * k, a[1] + b[1] * k, a[2] + b[2] * k]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn arr(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn vec(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}

#[cfg(test)]
impl Surface {
    pub(super) fn wet_count(&self) -> usize {
        self.wet.len()
    }
    pub(super) fn dry_count(&self) -> usize {
        self.dry.len()
    }
}

/// Finds the result cell holding a point: the cells binned on a uniform grid over the result.
pub(super) struct Bins<'a> {
    grid: &'a BoxGrid,
    low: [f64; 3],
    size: [f64; 3],
    count: [usize; 3],
    cells: Vec<Vec<u32>>,
}

impl<'a> Bins<'a> {
    pub(super) fn new(grid: &'a BoxGrid) -> Bins<'a> {
        let smallest = smallest_cell(grid).max(f64::MIN_POSITIVE);
        let low = [grid.bounds[0][0], grid.bounds[1][0], grid.bounds[2][0]];
        let extent = [0, 1, 2].map(|a| (grid.bounds[a][1] - grid.bounds[a][0]).max(smallest));
        let mut count = extent.map(|e| ((e / smallest).ceil()).clamp(1.0, MAX_BINS));
        let total = count[0] * count[1] * count[2];
        if total > MAX_BINS {
            let shrink = (total / MAX_BINS).cbrt();
            count = count.map(|c| (c / shrink).floor().max(1.0));
        }
        let count = count.map(|c| c as usize);
        let size = [0, 1, 2].map(|a| extent[a] / count[a] as f64);
        let mut cells = vec![Vec::new(); count[0] * count[1] * count[2]];
        for (index, cell) in grid.cells.iter().enumerate() {
            let range = [0, 1, 2].map(|a| {
                let first = ((cell[a][0] - low[a]) / size[a]).floor().clamp(0.0, (count[a] - 1) as f64) as usize;
                let last = ((cell[a][1] - low[a]) / size[a]).floor().clamp(0.0, (count[a] - 1) as f64) as usize;
                (first, last.max(first))
            });
            for k in range[2].0..=range[2].1 {
                for j in range[1].0..=range[1].1 {
                    for i in range[0].0..=range[0].1 {
                        if let Some(bin) = cells.get_mut((k * count[1] + j) * count[0] + i) {
                            bin.push(index as u32);
                        }
                    }
                }
            }
        }
        Bins { grid, low, size, count, cells }
    }

    /// The cell holding `p` (case m), if any.
    pub(super) fn find(&self, p: [f64; 3]) -> Option<usize> {
        let mut at = [0usize; 3];
        for a in 0..3 {
            let i = ((p[a] - self.low[a]) / self.size[a]).floor();
            if !i.is_finite() || i < 0.0 || i >= self.count[a] as f64 {
                return None;
            }
            at[a] = i as usize;
        }
        let bin = self.cells.get((at[2] * self.count[1] + at[1]) * self.count[0] + at[0])?;
        bin.iter().map(|&c| c as usize).find(|&c| {
            self.grid.cells.get(c).is_some_and(|cell| (0..3).all(|a| {
                let slack = (cell[a][1] - cell[a][0]).abs() * 1e-9;
                p[a] >= cell[a][0] - slack && p[a] <= cell[a][1] + slack
            }))
        })
    }
}

fn smallest_cell(grid: &BoxGrid) -> f64 {
    grid.cells.iter().flat_map(|cell| cell.iter().map(|r| r[1] - r[0])).filter(|e| *e > 0.0).fold(f64::INFINITY, f64::min)
}

/// The fluid cell next to a wetted triangle: tried at half, one and one and a half of the
/// smallest cell off the face along its outward normal.
fn cell_beside(bins: &Bins, origin: [f64; 3], centre: [f64; 3], normal: [f64; 3], step_mm: f64) -> Option<usize> {
    [0.5, 1.0, 1.5].into_iter().find_map(|k| bins.find(to_case(add_scaled(centre, normal, step_mm * k), origin)))
}

fn to_case(p_mm: [f64; 3], origin: [f64; 3]) -> [f64; 3] {
    [0, 1, 2].map(|a| (p_mm[a] - origin[a]) / 1000.0)
}

/// The model's display triangles split into wetted (with their cell) and dry.
fn surface(app: &SolveApp, id: &str, grid: &BoxGrid, origin: [f64; 3], key: (usize, usize, usize)) -> Surface {
    let bodies = app.session.doc.multiphysics.study(id).and_then(|s| s.model.as_ref()).map(|m| m.bodies.clone()).unwrap_or_default();
    let world = crate::open_meshing::shown_world(app);
    let shown: Vec<_> = world.bodies.iter().filter(|b| (bodies.is_empty() || bodies.contains(&b.name)) && app.body_visible(&b.name)).collect();
    let meshes: Vec<_> = shown.iter().map(|b| b.mesh()).collect();
    let mut surface = Surface { key, wet: Vec::new(), dry: Vec::new(), too_many: false };
    if meshes.iter().map(|m| m.triangles.len()).sum::<usize>() > MAX_TRIANGLES {
        surface.too_many = true;
        return surface;
    }
    let bins = Bins::new(grid);
    let step_mm = smallest_cell(grid) * 1000.0;
    for mesh in &meshes {
        for t in &mesh.triangles {
            let Some(points) = mesh.tri(t) else { continue };
            let p = points.map(arr);
            let n = cross(sub(p[1], p[0]), sub(p[2], p[0]));
            let length = dot(n, n).sqrt();
            if !(length > 0.0 && length.is_finite()) {
                continue;
            }
            let normal = n.map(|c| c / length);
            let centre = [0, 1, 2].map(|a| (p[0][a] + p[1][a] + p[2][a]) / 3.0);
            match step_mm.is_finite().then(|| cell_beside(&bins, origin, centre, normal, step_mm)).flatten() {
                Some(cell) => surface.wet.push(((p, normal), cell)),
                None => surface.dry.push((p, normal)),
            }
        }
    }
    surface
}

/// Where the study's domain starts in the design (mm), read once per document and study.
pub(super) fn origin(app: &mut SolveApp, view: &mut ViewState, id: &str) -> Option<[f64; 3]> {
    let key = std::sync::Arc::as_ptr(&app.session.doc) as usize;
    if view.model_origin.as_ref().is_none_or(|(k, study, _)| *k != key || study != id) {
        let origin = app
            .run("multiphysics.model.status", json!({"id":id}))
            .ok()
            .and_then(|reply| serde_json::from_value::<[[f64; 3]; 2]>(reply["domain_mm"].clone()).ok())
            .map(|[low, _]| low);
        view.model_origin = Some((key, id.to_string(), origin));
    }
    view.model_origin.as_ref().and_then(|(_, _, origin)| *origin)
}

/// The nearest wetted triangle a ray from the camera meets from its fluid side.
fn ray_hit(wet: &[(Triangle, usize)], o: [f64; 3], d: [f64; 3]) -> Option<(f64, usize)> {
    let mut best: Option<(f64, usize)> = None;
    for (index, ((p, normal), _)) in wet.iter().enumerate() {
        if dot(*normal, d) >= 0.0 {
            continue;
        }
        let (e1, e2) = (sub(p[1], p[0]), sub(p[2], p[0]));
        let h = cross(d, e2);
        let det = dot(e1, h);
        if det.abs() < 1e-300 {
            continue;
        }
        let s = sub(o, p[0]);
        let u = dot(s, h) / det;
        let q = cross(s, e1);
        let v = dot(d, q) / det;
        let t = dot(e2, q) / det;
        if u >= 0.0 && v >= 0.0 && u + v <= 1.0 && t > 0.0 && best.is_none_or(|(near, _)| t < near) {
            best = Some((t, index));
        }
    }
    best
}

/// Draws `plot` over the model in `rect`, takes a click as a probe, and returns how many cells or
/// faces it coloured.
pub(super) fn show(app: &mut SolveApp, ui: &mut egui::Ui, view: &mut ViewState, grid: &BoxGrid, rect: egui::Rect, plot: Plot, (low, high): (f64, f64)) -> usize {
    let Some(id) = view.selected.clone() else { return 0 };
    let Some(origin) = origin(app, view, &id) else { return 0 };
    // Only the part of the view on screen: the model's picture is fitted to what is on screen,
    // so a view scrolled half out of the panel would otherwise not line up with the result.
    let full = rect;
    let rect = rect.intersect(ui.clip_rect());
    if rect.width() < 8.0 || rect.height() < 8.0 {
        return 0;
    }
    let mut child = ui.new_child(egui::UiBuilder::new().max_rect(rect));
    crate::viewport::show(app, &mut child);
    // The design's camera is framed for the design view: frame the model once in this one.
    // Framed once the view is (nearly) all on screen, so the fit uses its real size.
    if view.model_framed.as_deref() != Some(id.as_str()) && rect.height() >= full.height() * 0.9 {
        view.model_framed = Some(id.clone());
        let _ = app.run("view.fit", json!({}));
        ui.ctx().request_repaint();
    }
    let proj = crate::viewport::projection(app, rect);
    let painter = ui.painter_at(rect);
    let (field, scale) = (view.slice_field, view.color_scale);
    let colour = |cell: usize| grid.value(field, cell).map(|value| scale_color(scale, value, low, high));
    let to_screen = |p: [f64; 3]| proj.to_screen(vec(p));
    let (eye, ahead) = proj.ray(rect.center());
    let (eye, ahead) = (arr(eye), arr(ahead));
    let clicked = ui.input(|i| if i.pointer.primary_clicked() { i.pointer.interact_pos() } else { None }).filter(|p| rect.contains(*p));
    let mut visible = 0;
    match plot {
        Plot::Slice { axis, along: (horizontal, vertical), coordinate } => {
            for (cell_index, cell) in grid.cells.iter().enumerate() {
                if !grid.intersects(cell_index, axis, coordinate) {
                    continue;
                }
                let Some(color) = colour(cell_index) else { continue };
                let (Some(x), Some(y)) = (cell.get(horizontal), cell.get(vertical)) else { continue };
                let corner = |a: f64, b: f64| {
                    let mut p = [0.0; 3];
                    p[axis] = coordinate;
                    p[horizontal] = a;
                    p[vertical] = b;
                    to_screen([0, 1, 2].map(|k| p[k] * 1000.0 + origin[k]))
                };
                let points: Option<Vec<egui::Pos2>> = [corner(x[0], y[0]), corner(x[1], y[0]), corner(x[1], y[1]), corner(x[0], y[1])].into_iter().collect();
                let Some(points) = points else { continue };
                let middle = points.iter().fold(egui::Vec2::ZERO, |sum, p| sum + p.to_vec2()) / points.len() as f32;
                crate::scenario::publish_target("multiphysics:slice_probe_target", egui::Rect::from_center_size(middle.to_pos2(), egui::vec2(2.0, 2.0)));
                for shape in crate::viewport::projected_polygon(&points, color.gamma_multiply(0.9), egui::Stroke::NONE) {
                    painter.add(shape);
                }
                visible += 1;
            }
            crate::scenario::publish_count("multiphysics:slice_on_model_cells", visible as f64);
            if view.model_vectors {
                let arrows = arrows(&painter, grid, axis, coordinate, origin, &to_screen);
                crate::scenario::publish_count("multiphysics:model_vectors_drawn", arrows as f64);
            }
            if let Some(pos) = clicked {
                let (o, d) = proj.ray(pos);
                let (o, d) = (arr(o), arr(d));
                let plane = coordinate * 1000.0 + origin[axis];
                if d[axis].abs() > 1e-12 {
                    let t = (plane - o[axis]) / d[axis];
                    let at_mm = add_scaled(o, d, t);
                    let case_m = to_case(at_mm, origin);
                    let cell = Bins::new(grid).find(case_m);
                    view.probe_point = case_m;
                    view.model_probe = Some(Probe { at_mm, case_m, cell });
                }
            }
        }
        Plot::Surfaces => {
            let key = (std::sync::Arc::as_ptr(&app.session.doc) as usize, std::ptr::from_ref(grid) as usize, grid.cells.len());
            if view.model_surface.as_ref().is_none_or(|s| s.key != key) {
                view.model_surface = Some(std::sync::Arc::new(surface(app, &id, grid, origin, key)));
            }
            let Some(surface) = view.model_surface.clone() else { return 0 };
            if surface.too_many {
                ui.colored_label(Tokens::get().warning, tr("multiphysics.surface_on_model_too_big"));
                return 0;
            }
            let shell = Tokens::get().control_fill.gamma_multiply(0.55);
            // Far to near, so nearer faces paint over farther ones; the dry faces are see-through.
            let mut order: Vec<(f64, [[f64; 3]; 3], egui::Color32, bool)> = Vec::new();
            for ((p, normal), cell) in &surface.wet {
                if dot(*normal, ahead) < 0.0
                    && let Some(color) = colour(*cell)
                {
                    let centre = [0, 1, 2].map(|a| (p[0][a] + p[1][a] + p[2][a]) / 3.0);
                    order.push((dot(sub(centre, eye), ahead), *p, color, true));
                }
            }
            for (p, normal) in &surface.dry {
                if dot(*normal, ahead) < 0.0 {
                    let centre = [0, 1, 2].map(|a| (p[0][a] + p[1][a] + p[2][a]) / 3.0);
                    order.push((dot(sub(centre, eye), ahead), *p, shell, false));
                }
            }
            order.sort_by(|a, b| b.0.total_cmp(&a.0));
            let mut mesh = egui::Mesh::default();
            let mut nearest: Option<egui::Pos2> = None;
            for (_, p, color, wet) in &order {
                let Some(points) = p.iter().map(|q| to_screen(*q)).collect::<Option<Vec<egui::Pos2>>>() else { continue };
                let base = mesh.vertices.len() as u32;
                for point in &points {
                    mesh.colored_vertex(*point, *color);
                }
                mesh.add_triangle(base, base + 1, base + 2);
                if *wet {
                    visible += 1;
                    nearest = Some(egui::pos2((points[0].x + points[1].x + points[2].x) / 3.0, (points[0].y + points[1].y + points[2].y) / 3.0));
                }
            }
            painter.add(egui::Shape::mesh(mesh));
            if let Some(point) = nearest {
                crate::scenario::publish_target("multiphysics:model_probe_target", egui::Rect::from_center_size(point, egui::vec2(2.0, 2.0)));
            }
            crate::scenario::publish_count("multiphysics:surface_on_model_faces", visible as f64);
            if let Some(pos) = clicked {
                let (o, d) = proj.ray(pos);
                let (o, d) = (arr(o), arr(d));
                let step_mm = smallest_cell(grid) * 1000.0;
                view.model_probe = ray_hit(&surface.wet, o, d).and_then(|(t, index)| {
                    let ((_, normal), cell) = surface.wet.get(index)?;
                    let at_mm = add_scaled(o, d, t);
                    Some(Probe { at_mm, case_m: to_case(add_scaled(at_mm, *normal, step_mm * 0.5), origin), cell: Some(*cell) })
                });
                if let Some(probe) = &view.model_probe {
                    view.probe_point = probe.case_m;
                }
            }
        }
        Plot::Flow => {}
    }
    if !matches!(plot, Plot::Surfaces)
        && let Some(overlay) = view.field_3d.overlay.clone()
    {
        // The flow tool's own scale, unless the plot is the flow itself (then the shown range).
        let (low, high) = if matches!(plot, Plot::Flow) { (low, high) } else { overlay.range };
        let drawn = flow_overlay(&painter, &overlay, origin, (eye, ahead), (low, high), scale, &to_screen);
        crate::scenario::publish_count("multiphysics:model_flow_drawn", drawn as f64);
        if matches!(plot, Plot::Flow) {
            visible = drawn;
        }
    }
    if !matches!(plot, Plot::Flow) {
        probe_marker(ui, &painter, view, grid, &|p| to_screen(p));
    }
    visible
}

/// The flow tool's lines (or isosurface) on the graphics window's model, while the tool is
/// open in the property panel beside it.
pub(super) fn flow_on_graphics(app: &mut SolveApp, ui: &egui::Ui, view: &mut ViewState, rect: egui::Rect) {
    let Some(overlay) = view.field_3d.overlay.clone() else { return };
    let Some(id) = view.selected.clone() else { return };
    let Some(origin) = origin(app, view, &id) else { return };
    // The design's camera may not show the model: frame it once, as the view in the panel does.
    if view.model_framed.as_deref() != Some(id.as_str()) {
        view.model_framed = Some(id.clone());
        let _ = app.run("view.fit", json!({}));
        ui.ctx().request_repaint();
    }
    let proj = crate::viewport::projection(app, rect);
    let (eye, ahead) = proj.ray(rect.center());
    let painter = ui.painter_at(rect);
    let to_screen = |p: [f64; 3]| proj.to_screen(vec(p));
    let drawn = flow_overlay(&painter, &overlay, origin, (arr(eye), arr(ahead)), view.model_flow_range, view.color_scale, &to_screen);
    crate::scenario::publish_count("multiphysics:model_flow_drawn", drawn as f64);
}

/// Arrows of the in-plane velocity at the cut plot's cells, scaled to the largest; returns how
/// many were drawn.
fn arrows(
    painter: &egui::Painter,
    grid: &BoxGrid,
    axis: usize,
    coordinate: f64,
    origin: [f64; 3],
    to_screen: &dyn Fn([f64; 3]) -> Option<egui::Pos2>,
) -> usize {
    let Some(field) = vector_field(grid) else { return 0 };
    let mut cells: Vec<(usize, [f64; 3])> = Vec::new();
    for cell in 0..grid.cells.len() {
        if grid.intersects(cell, axis, coordinate)
            && let Some(mut v) = vector(grid, field, cell)
        {
            v[axis] = 0.0;
            cells.push((cell, v));
        }
    }
    let largest = cells.iter().map(|(_, v)| dot(*v, *v).sqrt()).fold(0.0, f64::max);
    if largest <= 0.0 || !largest.is_finite() {
        return 0;
    }
    let stride = cells.len().div_ceil(MAX_ARROWS).max(1);
    // Skipped cells give the drawn arrows more room.
    let room = (stride as f64).sqrt().min(3.0) * 0.9;
    let color = Tokens::get().text;
    let mut drawn = 0;
    for (cell, v) in cells.iter().step_by(stride) {
        let Some(bounds) = grid.cells.get(*cell) else { continue };
        let size = (0..3).filter(|a| *a != axis).map(|a| bounds[a][1] - bounds[a][0]).fold(f64::INFINITY, f64::min);
        let mut centre = [0, 1, 2].map(|a| (bounds[a][0] + bounds[a][1]) * 0.5);
        centre[axis] = coordinate;
        let length = dot(*v, *v).sqrt();
        if length <= largest * 1e-3 {
            continue;
        }
        let tip = add_scaled(centre, *v, size * room / largest);
        let (Some(a), Some(b)) = (to_screen(add_scaled(origin, centre, 1000.0)), to_screen(add_scaled(origin, tip, 1000.0))) else { continue };
        let along = b - a;
        let pixels = along.length();
        if pixels < 1.0 {
            continue;
        }
        let stroke = egui::Stroke::new(1.5, color);
        painter.line_segment([a, b], stroke);
        let head = (pixels * 0.35).min(7.0);
        let back = -along / pixels * head;
        for turn in [0.45_f32, -0.45] {
            let (sin, cos) = turn.sin_cos();
            let side = egui::vec2(back.x * cos - back.y * sin, back.x * sin + back.y * cos);
            painter.line_segment([b, b + side], stroke);
        }
        drawn += 1;
    }
    drawn
}

/// The flow tool's trajectories (coloured by speed) and isosurface triangles over the model;
/// returns how many segments and triangles were drawn.
fn flow_overlay(
    painter: &egui::Painter,
    overlay: &crate::multiphysics_field_view::Overlay,
    origin: [f64; 3],
    (eye, ahead): ([f64; 3], [f64; 3]),
    (low, high): (f64, f64),
    scale: post::ColorScale,
    to_screen: &dyn Fn([f64; 3]) -> Option<egui::Pos2>,
) -> usize {
    let colour = |value: f64| scale_color(scale, value, low, high);
    let at = |p: [f64; 3]| add_scaled(origin, p, 1000.0);
    let mut drawn = 0;
    let stride = overlay.triangles.len().div_ceil(MAX_OVERLAY_TRIANGLES).max(1);
    let mut order: Vec<(f64, &([[f64; 3]; 3], [f64; 3]))> = overlay
        .triangles
        .iter()
        .step_by(stride)
        .map(|t| {
            let centre = [0, 1, 2].map(|a| (t.0[0][a] + t.0[1][a] + t.0[2][a]) / 3.0);
            (dot(sub(at(centre), eye), ahead), t)
        })
        .collect();
    order.sort_by(|a, b| b.0.total_cmp(&a.0));
    let mut mesh = egui::Mesh::default();
    for (_, (corners, values)) in order {
        let Some(points) = corners.iter().map(|q| to_screen(at(*q))).collect::<Option<Vec<egui::Pos2>>>() else { continue };
        let base = mesh.vertices.len() as u32;
        for (point, value) in points.iter().zip(values) {
            mesh.colored_vertex(*point, colour(*value).gamma_multiply(0.85));
        }
        mesh.add_triangle(base, base + 1, base + 2);
        drawn += 1;
    }
    painter.add(egui::Shape::mesh(mesh));
    let segments: usize = overlay.lines.iter().map(|line| line.len().saturating_sub(1)).sum();
    let stride = segments.div_ceil(MAX_SEGMENTS).max(1);
    for line in &overlay.lines {
        let points: Vec<(Option<egui::Pos2>, f64)> = line.iter().step_by(stride).map(|&(p, v)| (to_screen(at(p)), v)).collect();
        for pair in points.windows(2) {
            if let [(Some(a), va), (Some(b), vb)] = pair {
                painter.line_segment([*a, *b], egui::Stroke::new(2.0, colour((va + vb) * 0.5)));
                drawn += 1;
            }
        }
    }
    drawn
}

/// The probe's marker on the model and its value below the view.
fn probe_marker(ui: &mut egui::Ui, painter: &egui::Painter, view: &ViewState, grid: &BoxGrid, to_screen: &dyn Fn([f64; 3]) -> Option<egui::Pos2>) {
    let Some(probe) = &view.model_probe else { return };
    let t = Tokens::get();
    let field = grid.fields.get(view.slice_field).map(|f| f.name.as_str()).unwrap_or_default();
    let unit = post::unit_hint(field).map(|u| format!(" {u}")).unwrap_or_default();
    let value = probe.cell.and_then(|cell| grid.value(view.slice_field, cell));
    let text = match value {
        Some(value) => format!("{field} = {value:.4e}{unit}"),
        None => tr("multiphysics.probe_no_fluid").to_string(),
    };
    if let Some(at) = to_screen(probe.at_mm) {
        painter.circle(at, 5.0, t.accent_soft, egui::Stroke::new(2.0, t.accent));
        painter.text(at + egui::vec2(9.0, -9.0), egui::Align2::LEFT_BOTTOM, &text, egui::FontId::proportional(t.font_small), t.text);
    }
    let response = ui.monospace(format!(
        "{} ({:.2}, {:.2}, {:.2}) mm: {text}",
        tr("multiphysics.probe_on_model"),
        probe.at_mm[0],
        probe.at_mm[1],
        probe.at_mm[2]
    ));
    crate::scenario::publish_rect("multiphysics:model_probe", response.rect);
    if let Some(value) = value {
        crate::scenario::publish_count("multiphysics:model_probe_value", value);
    }
}
