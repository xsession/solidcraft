//! Ruled tables and tree rows shared by the Multiphysics database, wizard and project tree: a
//! tinted column header, group bands, striped Parameter | Value rows with a column rule, and
//! clickable multi-column item rows. Only presentation; callers own the values and handles.
use egui::{Rect, Sense, pos2, vec2};
use solidcraft_theme::Tokens;

/// Share of the width given to the parameter column of a two-column grid.
const LABEL_SHARE: f32 = 0.46;

fn text_in(painter: &egui::Painter, rect: Rect, text: &str, font: egui::FontId, color: egui::Color32) {
    let galley = painter.layout_no_wrap(text.to_owned(), font, color);
    let y = rect.center().y - galley.size().y / 2.0;
    painter.with_clip_rect(rect.intersect(painter.clip_rect())).galley(pos2(rect.left(), y), galley, color);
}

fn rule(painter: &egui::Painter, rect: Rect) {
    let t = Tokens::get();
    painter.hline(rect.x_range(), rect.bottom() - 0.5, egui::Stroke::new(1.0, t.dialog_rule));
}

/// Column edges for `shares` (fractions summing to about one) across `rect`.
fn columns(rect: Rect, shares: &[f32]) -> Vec<Rect> {
    let mut x = rect.left();
    let total: f32 = shares.iter().sum::<f32>().max(f32::EPSILON);
    shares
        .iter()
        .map(|share| {
            let width = rect.width() * share / total;
            let cell = Rect::from_min_max(pos2(x, rect.top()), pos2(x + width, rect.bottom()));
            x += width;
            cell
        })
        .collect()
}

/// The tinted header of a table: one caption per column.
pub(crate) fn header(ui: &mut egui::Ui, captions: &[&str], shares: &[f32]) {
    let t = Tokens::get();
    let (rect, _) = ui.allocate_exact_size(vec2(ui.available_width(), crate::ui_kit::row_height()), Sense::hover());
    let painter = ui.painter();
    painter.rect_filled(rect, 0.0, t.panel_header);
    let cells = columns(rect, shares);
    for (index, (caption, cell)) in captions.iter().zip(&cells).enumerate() {
        if index > 0 {
            painter.vline(cell.left(), cell.y_range(), egui::Stroke::new(1.0, t.dialog_rule));
        }
        text_in(painter, cell.shrink2(vec2(t.space_8, 0.0)), caption, crate::ui_kit::type_face::semibold(t.font_small), t.text_dim);
    }
    rule(painter, rect);
}

/// A group band spanning the table, such as "Thermodynamic parameters".
pub(crate) fn band(ui: &mut egui::Ui, text: &str) {
    let t = Tokens::get();
    let (rect, _) = ui.allocate_exact_size(vec2(ui.available_width(), crate::ui_kit::row_height()), Sense::hover());
    let painter = ui.painter();
    painter.rect_filled(rect, 0.0, t.control_fill);
    text_in(painter, rect.shrink2(vec2(t.space_8, 0.0)), text, crate::ui_kit::type_face::semibold(t.font_body), t.text);
    rule(painter, rect);
}

/// A sub-heading row inside a group (one fan, one boundary), indented by `depth`.
pub(crate) fn heading(ui: &mut egui::Ui, depth: usize, text: &str) {
    let t = Tokens::get();
    let (rect, _) = ui.allocate_exact_size(vec2(ui.available_width(), crate::ui_kit::row_height()), Sense::hover());
    let painter = ui.painter();
    let left = rect.left() + t.space_8 + depth as f32 * t.space_14;
    text_in(painter, Rect::from_min_max(pos2(left, rect.top()), rect.max), text, crate::ui_kit::type_face::semibold(t.font_body), t.text);
    rule(painter, rect);
}

/// One Parameter | Value row: the label in the left column, `value` drawn in the right one.
/// Rows alternate their fill when `stripe` is set; the response covers the whole row.
pub(crate) fn row<R>(ui: &mut egui::Ui, depth: usize, label: &str, stripe: bool, value: impl FnOnce(&mut egui::Ui) -> R) -> (egui::Response, R) {
    let t = Tokens::get();
    let height = crate::ui_kit::row_height().max(t.control_height + t.space_4);
    let (rect, response) = ui.allocate_exact_size(vec2(ui.available_width(), height), Sense::hover());
    if stripe {
        ui.painter().rect_filled(rect, 0.0, t.panel);
    }
    let split = rect.left() + rect.width() * LABEL_SHARE;
    let label_rect = Rect::from_min_max(pos2(rect.left() + t.space_8 + depth as f32 * t.space_14, rect.top()), pos2(split - t.space_4, rect.bottom()));
    text_in(ui.painter(), label_rect, label, egui::FontId::proportional(t.font_body), t.text);
    ui.painter().vline(split, rect.y_range(), egui::Stroke::new(1.0, t.dialog_rule));
    rule(ui.painter(), rect);
    let cell = Rect::from_min_max(pos2(split + t.space_6, rect.top() + t.space_2), pos2(rect.right() - t.space_6, rect.bottom() - t.space_2));
    let inner = ui.scope_builder(egui::UiBuilder::new().max_rect(cell).layout(egui::Layout::left_to_right(egui::Align::Center)), |ui| {
        ui.set_clip_rect(cell.intersect(ui.clip_rect()));
        value(ui)
    });
    let full = ui.painter().layout_no_wrap(label.to_owned(), egui::FontId::proportional(t.font_body), t.text).size().x;
    let response = if label_rect.width() < full { response.on_hover_text(label) } else { response };
    (response, inner.inner)
}

/// One clickable item row of a multi-column table; `selected` fills it with the accent.
pub(crate) fn item(ui: &mut egui::Ui, cells: &[&str], shares: &[f32], selected: bool, stripe: bool) -> egui::Response {
    let t = Tokens::get();
    let (rect, response) = ui.allocate_exact_size(vec2(ui.available_width(), crate::ui_kit::row_height()), Sense::click());
    let painter = ui.painter();
    if selected {
        painter.rect_filled(rect, 0.0, t.accent_soft);
    } else if response.hovered() {
        painter.rect_filled(rect, 0.0, t.control_fill);
    } else if stripe {
        painter.rect_filled(rect, 0.0, t.panel);
    }
    let color = if selected { t.selected_ink } else { t.text };
    for (index, (text, cell)) in cells.iter().zip(columns(rect, shares)).enumerate() {
        if index > 0 {
            painter.vline(cell.left(), cell.y_range(), egui::Stroke::new(1.0, t.dialog_rule));
        }
        text_in(painter, cell.shrink2(vec2(t.space_8, 0.0)), text, egui::FontId::proportional(t.font_body), color);
    }
    rule(painter, rect);
    response.widget_info(|| egui::WidgetInfo::labeled(egui::Role::Button, true, cells.first().copied().unwrap_or("")));
    response
}

/// A framed pane with a tinted title strip, as the database tree and the wizard's step list use.
pub(crate) fn pane<R>(ui: &mut egui::Ui, rect: Rect, title: &str, add: impl FnOnce(&mut egui::Ui) -> R) -> R {
    let t = Tokens::get();
    ui.painter().rect_filled(rect, t.control_radius, t.panel);
    ui.painter().rect_stroke(rect, t.control_radius, egui::Stroke::new(1.0, t.border), egui::StrokeKind::Inside);
    let strip = Rect::from_min_size(rect.min, vec2(rect.width(), crate::ui_kit::row_height()));
    ui.painter().rect_filled(strip.shrink(1.0), t.control_radius, t.panel_header);
    text_in(ui.painter(), strip.shrink2(vec2(t.space_8, 0.0)), title, crate::ui_kit::type_face::semibold(t.font_small), t.text_dim);
    let body = Rect::from_min_max(pos2(rect.left() + t.space_2, strip.bottom() + t.space_2), rect.max - vec2(t.space_2, t.space_2));
    ui.scope_builder(egui::UiBuilder::new().max_rect(body).layout(egui::Layout::top_down(egui::Align::Min)), |ui| {
        ui.set_clip_rect(body.intersect(ui.clip_rect()));
        add(ui)
    })
    .inner
}

pub(crate) fn disclosure_left(depth: usize) -> f32 {
    Tokens::get().space_4 + depth as f32 * Tokens::get().space_16
}

/// One tree row: indent, an optional disclosure triangle, the icon and the label. `bold` marks
/// the active project, as flow-analysis trees do.
#[allow(clippy::too_many_arguments)]
pub(crate) fn tree_row(ui: &mut egui::Ui, depth: usize, open: Option<bool>, icon: &str, text: &str, selected: bool, enabled: bool, bold: bool) -> egui::Response {
    let t = Tokens::get();
    let height = crate::ui_kit::row_height();
    let (rect, response) = ui.allocate_exact_size(vec2(ui.available_width(), height), if enabled { Sense::click() } else { Sense::hover() });
    let painter = ui.painter();
    if selected {
        painter.rect_filled(rect, t.control_radius, t.accent_soft);
    } else if response.hovered() && enabled {
        painter.rect_filled(rect, t.control_radius, t.control_fill);
    }
    let mut x = rect.left() + disclosure_left(depth);
    let glyph = t.space_12;
    if let Some(open) = open {
        let disclosure = Rect::from_center_size(pos2(x + glyph / 2.0, rect.center().y), vec2(glyph, glyph));
        crate::icons::paint_ink(painter, disclosure, if open { "disclosure_open" } else { "disclosure_closed" }, t.text_dim);
    }
    x += glyph + t.space_4;
    let icon_size = t.space_16;
    let icon_rect = Rect::from_center_size(pos2(x + icon_size / 2.0, rect.center().y), vec2(icon_size, icon_size));
    crate::icons::paint_icon(painter, icon_rect, icon, if enabled { 1.0 } else { 0.4 });
    x += icon_size + t.space_6;
    let font = if bold { crate::ui_kit::type_face::semibold(t.font_body) } else { egui::FontId::proportional(t.font_body) };
    let color = if !enabled {
        t.text_dim
    } else if selected {
        t.selected_ink
    } else {
        t.text
    };
    let room = (rect.right() - x - t.space_4).max(0.0);
    let galley = painter.layout(text.to_owned(), font, color, room);
    let galley_height = galley.size().y;
    painter.with_clip_rect(rect).galley(pos2(x, rect.center().y - galley_height / 2.0), galley, color);
    response.widget_info(|| egui::WidgetInfo::labeled(egui::Role::Button, enabled, text));
    response
}
