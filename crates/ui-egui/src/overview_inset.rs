//! The overview inset: while parts are drawn see-through (ghosted), a small picture in the view's
//! lower right corner shows the whole model from the view's own direction, with a frame round
//! what the view itself shows. Zoomed in on a few parts among ghosts, where they sit in the
//! assembly stays in sight.
//!
//! The inset is a second picture of the scenes the view already holds on the GPU
//! ([`solidcraft_viewport_egui::InsetCallback`]), framed to the model's box, drawn with the view's
//! own uploads: nothing is built for it. It changes only when the view turns or the scenes change, so panning and zooming the view
//! only move the frame. It takes no input of its own; the view does not pick or navigate
//! through it either ([`crate::navigation`]).

use egui::{Pos2, Rect, Sense, Shape, Stroke, StrokeKind, pos2, vec2};
use solidcore_record::Vec3;
use solidcraft_engine::geom::Aabb3;
use solidcraft_engine::geom::prelude::*;
use solidcraft_theme::Tokens;
use solidcraft_viewport::Camera;

use crate::SolveApp;
use crate::ui_kit::BackingTooltipExt as _;
use crate::viewport::{FrameQuality, Proj};

/// The GPU slot of the inset's picture: far above the slots the views of the main and the native
/// windows count up to ([`crate::viewport::ViewportState::scene_slot`]).
const SLOT: u64 = u64::MAX - 1;
/// The inset's size, in points.
const SIZE: egui::Vec2 = vec2(240.0, 160.0);
/// A view smaller than this each way has no room for the inset.
const ROOM: egui::Vec2 = vec2(560.0, 400.0);
/// The model's box fills this much of the inset, so its edges are not cut.
const FILL: f64 = 0.86;

/// Whether the view shows the inset: parts are ghosted and the view is the window's only one.
fn wanted(app: &SolveApp) -> bool {
    !app.session.ghosted_bodies.is_empty() && app.viewport.gpu.is_some() && app.viewport.tile.is_none() && !app.viewport.card
}

/// Where the inset goes in `view`: its lower right corner, above the panels at the bottom.
pub(crate) fn place(app: &SolveApp, view: Rect) -> Option<Rect> {
    if view.width() < ROOM.x || view.height() < ROOM.y {
        return None;
    }
    let free = crate::shell::bottom_readout_rect(app, view);
    let gap = Tokens::get().space_8 * 1.5;
    let rect = Rect::from_min_size(pos2(free.right() - gap - SIZE.x, free.bottom() - gap - SIZE.y), SIZE);
    view.contains_rect(rect).then_some(rect)
}

/// The inset's camera: the view's direction, orthographic, framing `bounds` in a picture of
/// `aspect`.
pub(crate) fn camera(view: &Camera, bounds: &Aabb3, aspect: f64) -> Camera {
    let mut cam = *view;
    cam.fov = 0.0;
    cam.target = bounds.center();
    let (right, up, _) = cam.basis();
    let (mut wide, mut high) = (0.0f64, 0.0f64);
    for i in 0..8 {
        let corner = Vec3::new(
            if i & 1 == 0 { bounds.min.x } else { bounds.max.x },
            if i & 2 == 0 { bounds.min.y } else { bounds.max.y },
            if i & 4 == 0 { bounds.min.z } else { bounds.max.z },
        ) - cam.target;
        wide = wide.max(corner.dot(right).abs());
        high = high.max(corner.dot(up).abs());
    }
    let half_height = high.max(wide / aspect.max(1e-3)).max(1e-3) / FILL;
    // An orthographic camera shows `0.4 × distance` above and below its target (Camera::half_height).
    cam.distance = half_height / 0.4;
    cam
}

/// What `view` shows, outlined in `inset`: its corners' rays met with the plane through the
/// model's centre square to the view, in the inset's picture.
pub(crate) fn footprint(view: &Proj, inset: &Proj, centre: Vec3) -> Option<[Pos2; 4]> {
    let normal = view.cam.back();
    let r = view.rect;
    let corner = |p: Pos2| {
        let (origin, dir) = view.ray(p);
        let along = dir.dot(normal);
        if along.abs() < 1e-9 {
            return None;
        }
        let t = (centre - origin).dot(normal) / along;
        let at = origin + dir * t;
        inset.to_screen(at).filter(|p| p.is_finite())
    };
    Some([corner(r.left_top())?, corner(r.right_top())?, corner(r.right_bottom())?, corner(r.left_bottom())?])
}

/// The inset's frame this frame, prepared by the view ([`prepare`]) for the workspace to draw
/// ([`paint`]).
pub(crate) struct Pending {
    rect: Rect,
    inset: Proj,
    frame: std::sync::Arc<solidcraft_viewport_gpu::view::Frame>,
}

/// Called by the window's view before it adds its own picture to `painter`: when the inset is
/// wanted over `view`, its picture is prepared first, with the view's scenes, so the view's own
/// writes to the GPU come after it ([`solidcraft_viewport_egui::InsetCallback`]).
pub(crate) fn prepare(app: &mut SolveApp, ctx: &egui::Context, painter: &egui::Painter, view: Rect) {
    app.viewport.inset_pending = None;
    let bounds = crate::viewport::model_bounds(app);
    let Some(rect) = place(app, view).filter(|_| wanted(app) && !bounds.is_empty() && bounds.diagonal().is_finite()) else { return };
    let aspect = f64::from(rect.width() / rect.height());
    let cam = camera(&app.cam, &bounds, aspect);
    let vp = cam.view_proj_bounds(aspect, &solidcraft_engine::render::view_bounds(&bounds));
    let ppp = ctx.pixels_per_point();
    let quality = FrameQuality { ssaa: 2.0, ao_samples: solidcraft_viewport_gpu::view::AO_SAMPLES, edge_px: solidcraft_viewport_gpu::view::EDGE_PX };
    let frame = std::sync::Arc::new(crate::viewport::model_frame(app, ctx, &cam, &vp, [rect.width() * ppp, rect.height() * ppp], quality, None));
    painter.add(egui_wgpu::Callback::new_paint_callback(
        rect,
        solidcraft_viewport_egui::InsetCallback { slot: SLOT, frame: frame.clone(), draws: false },
    ));
    app.viewport.inset_pending = Some(Pending { rect, inset: Proj { cam, vp, rect, bounds: Some(bounds) }, frame });
}

/// Draws the inset the view prepared over `view`, after the view and what is drawn over it.
pub fn paint(app: &mut SolveApp, ui: &mut egui::Ui, view: Rect) {
    let Some(Pending { rect, inset, frame }) = app.viewport.inset_pending.take() else { return };
    app.viewport.inset = Some(rect);
    let t = Tokens::get();
    let painter = ui.painter_at(rect);
    solidcraft_viewport_egui::paint_background(&painter, rect, t.viewport_top, t.viewport_bottom);
    painter.add(egui_wgpu::Callback::new_paint_callback(rect, solidcraft_viewport_egui::InsetCallback { slot: SLOT, frame, draws: true }));
    let main = crate::viewport::projection(app, view);
    if let Some(corners) = inset.bounds.and_then(|bounds| footprint(&main, &inset, bounds.center())) {
        painter.add(Shape::closed_line(corners.to_vec(), Stroke::new(1.5, t.accent)));
    }
    // The card's edge: drawn by an outer painter, so the clip does not halve it.
    ui.painter().rect_stroke(rect, egui::CornerRadius::same(t.radius_md), Stroke::new(1.0, t.border), StrokeKind::Outside);
    ui.interact(rect, ui.scope_id().with("overview_inset"), Sense::hover())
        .backing_on_hover_text(solidcraft_ui_text::tr("viewport.overview_inset_tip"));
    crate::scenario::publish_rect("viewport:overview_inset", rect);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cube(size: f64) -> Aabb3 {
        Aabb3 { min: Vec3::new(0.0, 0.0, 0.0), max: Vec3::new(size, size * 0.5, size * 0.25) }
    }

    fn proj(cam: Camera, rect: Rect, bounds: &Aabb3) -> Proj {
        let aspect = f64::from(rect.width() / rect.height());
        Proj { cam, vp: cam.view_proj_bounds(aspect, &solidcraft_engine::render::view_bounds(bounds)), rect, bounds: Some(*bounds) }
    }

    #[test]
    fn the_inset_frames_the_whole_model_from_the_views_direction() {
        let bounds = cube(400.0);
        let rect = Rect::from_min_size(pos2(10.0, 20.0), SIZE);
        let mut view = Camera { distance: 30.0, target: Vec3::new(5.0, 5.0, 5.0), fov: 0.6, ..Camera::default() };
        for (yaw, pitch) in [(0.3, 0.4), (2.0, -0.8), (-1.2, 1.3)] {
            view.yaw = yaw;
            view.pitch = pitch;
            let cam = camera(&view, &bounds, f64::from(SIZE.x / SIZE.y));
            assert_eq!((cam.yaw, cam.pitch, cam.fov), (yaw, pitch, 0.0));
            let inset = proj(cam, rect, &bounds);
            let mut fill = 0.0f32;
            for i in 0..8 {
                let corner =
                    Vec3::new(if i & 1 == 0 { 0.0 } else { 400.0 }, if i & 2 == 0 { 0.0 } else { 200.0 }, if i & 4 == 0 { 0.0 } else { 100.0 });
                let p = inset.to_screen(corner).unwrap();
                assert!(rect.contains(p), "{p:?} outside {rect:?}");
                let off = p - rect.center();
                fill = fill.max((off.x / (rect.width() / 2.0)).abs()).max((off.y / (rect.height() / 2.0)).abs());
            }
            // The model fills the inset (FILL of it one way), not a speck of it.
            assert!((fill - FILL as f32).abs() < 0.01, "{fill}");
        }
    }

    #[test]
    fn the_frame_shows_what_the_view_shows() {
        let bounds = cube(400.0);
        let screen = Rect::from_min_size(Pos2::ZERO, vec2(800.0, 600.0));
        let inset_rect = Rect::from_min_size(pos2(500.0, 400.0), SIZE);
        // Zoomed in on the model's centre: a small frame round the inset's middle.
        let view = Camera { target: bounds.center(), distance: 40.0, ..Camera::default() };
        let main = proj(view, screen, &bounds);
        let inset = proj(camera(&view, &bounds, 1.5), inset_rect, &bounds);
        let corners = footprint(&main, &inset, bounds.center()).unwrap();
        let framed = Rect::from_points(&corners);
        assert!(inset_rect.contains_rect(framed));
        assert!(framed.width() < inset_rect.width() * 0.3);
        assert!((framed.center() - inset_rect.center()).length() < 2.0, "{:?}", framed.center());
        // Panned to a side: the frame moves with it the same way.
        let panned = Camera { target: bounds.center() + view.basis().0 * 100.0, ..view };
        let moved = Rect::from_points(&footprint(&proj(panned, screen, &bounds), &inset, bounds.center()).unwrap());
        assert!(moved.center().x > framed.center().x + 10.0);
        assert!((moved.width() - framed.width()).abs() < 0.5);
    }

    #[test]
    fn a_small_view_has_no_inset() {
        let app = SolveApp::new(solidcraft_engine::Session::default(), crate::Services::default());
        assert!(place(&app, Rect::from_min_size(Pos2::ZERO, vec2(500.0, 900.0))).is_none());
        let view = Rect::from_min_size(pos2(300.0, 100.0), vec2(900.0, 700.0));
        let rect = place(&app, view).unwrap();
        assert!(view.contains_rect(rect));
        assert!(rect.right() > view.center().x && rect.bottom() > view.center().y);
    }

    /// The inset is a second picture of the view's own scenes, prepared before the view's: it
    /// shows the model, and the view beside it is drawn exactly as without it.
    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn the_inset_shows_the_model_and_leaves_the_view_as_it_was() {
        use serde_json::json;
        let mut app = SolveApp::new(solidcraft_engine::Session::default(), crate::Services::default());
        let size = [900u32, 640u32];
        let Some(mut gpu) = crate::gpu::EguiFixture::new(&mut app, size) else { return };
        app.ui.ground_shadow = false;
        app.set_origin_shown(false);
        app.run("solid.box", json!({"length": 40, "width": 30, "height": 20, "corner": [-80, 0, 0], "body_name": "A"})).unwrap();
        app.run("solid.box", json!({"length": 40, "width": 30, "height": 20, "corner": [40, 0, 0], "body_name": "B"})).unwrap();
        app.fit_view();
        let ctx = crate::fonts::test_context();
        let view = std::cell::Cell::new(Rect::NOTHING);
        let mut time = 1.0;
        // The view turns a little every frame, as when orbiting: it is drawn anew each time, not
        // shown from its cached picture.
        let start = app.cam;
        let mut draw = |app: &mut SolveApp| {
            let mut pixels = Vec::new();
            for turn in 0..8 {
                time += 0.5;
                app.cam.yaw = start.yaw + 0.01 * f64::from(turn);
                pixels = gpu.frame_at(&ctx, time, vec![], |ui| {
                    let rect = ui.max_rect();
                    view.set(rect);
                    crate::viewport::show(app, ui);
                    paint(app, ui, rect);
                });
            }
            pixels
        };
        let pixel = |pixels: &[u8], p: Pos2| {
            let at = (p.y as usize * size[0] as usize + p.x as usize) * 4;
            [pixels[at], pixels[at + 1], pixels[at + 2]]
        };
        let plain = draw(&mut app);
        assert!(app.viewport.inset.is_none(), "no inset while nothing is ghosted");
        app.session.ghosted_bodies.insert("A".into());
        let ghosted = draw(&mut app);
        let rect = app.viewport.inset.expect("an inset while a part is ghosted");
        // The solid box B in the view: the same pixels with the inset as without.
        let main = crate::viewport::projection(&app, view.get());
        let b = main.to_screen(Vec3::new(60.0, 15.0, 10.0)).unwrap();
        assert!(!rect.expand(4.0).contains(b));
        for dx in [-6.0, 0.0, 6.0] {
            for dy in [-6.0, 0.0, 6.0] {
                let p = b + vec2(dx, dy);
                assert_eq!(pixel(&plain, p), pixel(&ghosted, p), "the view's pixel at {p:?}");
            }
        }
        // B in the inset: the model, not the inset's background.
        let bounds = crate::viewport::model_bounds(&app);
        let aspect = f64::from(rect.width() / rect.height());
        let inset = proj(camera(&app.cam, &bounds, aspect), rect, &bounds);
        let in_inset = inset.to_screen(Vec3::new(60.0, 15.0, 10.0)).unwrap();
        let background = pixel(&ghosted, rect.left_top() + vec2(4.0, rect.height() / 2.0));
        let drawn = pixel(&ghosted, in_inset);
        let differs = drawn.iter().zip(background).map(|(a, b)| i32::from(*a).abs_diff(i32::from(b))).max().unwrap_or(0);
        assert!(differs > 12, "the inset shows box B: {drawn:?} against the background {background:?}");
    }
}
