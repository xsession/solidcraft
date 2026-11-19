//! Rotating regions (`rotating_regions_v1`, MRF): a closed box lying wholly in a region
//! settles into exact solid-body rotation with p = ½Ω²r²; a rotor in a partial region stirs
//! the fluid in its direction of turn; bad regions are refused.
use cfd_solvers::fvm::{PressureBoundaryType, VelocityBoundaryType};
use cfd_solvers::{case::rotating::RotatingRegion, case::*, core::Vec3};
use std::sync::atomic::{AtomicUsize, Ordering};
static NEXT: AtomicUsize = AtomicUsize::new(0);

struct Scratch(std::path::PathBuf);
impl Scratch {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!("cfd-rotating-{}-{}", std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
        std::fs::create_dir_all(&p).unwrap();
        Self(p)
    }
}
impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// One radian per second.
const RPM_ONE_RAD: f64 = 60.0 / std::f64::consts::TAU;

/// A closed n × n × 1 box of side `side`, no-slip walls, slip front and back.
fn closed_box(n: usize, side: f64, steps: usize, dt: f64) -> CollocatedFlowCase {
    let mut case = CollocatedFlowCase {
        mesh: PolyMeshSpec { dimensions: [n, n, 1], lengths: Vec3::new(side, side, side / n as f64), xy_shear: 0.0, coordinates: None },
        method: FlowMethod::Piso,
        output: PolyOutput { steps, output_interval: steps },
        ..Default::default()
    };
    case.config.dt = dt;
    case.config.density = 1.0;
    case.config.kinematic_viscosity = 1e-3;
    case.config.pressure_tolerance = 1e-10;
    case.config.pressure_iterations = 4000;
    case.config.piso_neighbour_update = true;
    case.config.balanced_body_force = true;
    for (name, kind) in [
        ("left", VelocityBoundaryType::FixedValue),
        ("right", VelocityBoundaryType::FixedValue),
        ("bottom", VelocityBoundaryType::FixedValue),
        ("top", VelocityBoundaryType::FixedValue),
        ("front", VelocityBoundaryType::Slip),
        ("back", VelocityBoundaryType::Slip),
    ] {
        case.velocity_boundaries.push(VelocityPatch { patch: name.into(), kind, value: Vec3::ZERO });
    }
    let _ = PressureBoundaryType::ZeroGradient;
    case
}

fn region(name: &str, min: [usize; 3], max: [usize; 3], centre: Vec3, rpm: f64) -> RotatingRegion {
    RotatingRegion { name: name.into(), min, max, axis_point: centre, axis: Vec3::new(0.0, 0.0, 1.0), rpm }
}

/// Cell velocities and pressures of the last VTU frame, in cell order.
fn fields(dir: &std::path::Path) -> (Vec<Vec3>, Vec<f64>) {
    let mut frames: Vec<_> = std::fs::read_dir(dir)
        .unwrap()
        .filter_map(|e| e.ok().map(|e| e.file_name().to_string_lossy().to_string()))
        .filter(|n| n.ends_with(".vtu"))
        .collect();
    frames.sort();
    let text = std::fs::read_to_string(dir.join(frames.last().unwrap())).unwrap();
    let array = |name: &str| -> Vec<f64> {
        let start = text.find(&format!("Name=\"{name}\"")).unwrap();
        let body = &text[start..];
        let body = &body[body.find('>').unwrap() + 1..body.find("</DataArray>").unwrap()];
        body.split_whitespace().map(|v| v.parse().unwrap()).collect()
    };
    let u = array("velocity_m_per_s").chunks(3).map(|c| Vec3::new(c[0], c[1], c[2])).collect();
    (u, array("kinematic_pressure_m2_per_s2"))
}

#[test]
fn a_box_turning_with_its_frame_settles_into_solid_body_rotation() {
    let scratch = Scratch::new();
    let (n, side) = (16, 0.1);
    let centre = Vec3::new(0.05, 0.05, 0.0);
    // Ω = 1 rad/s, dt 0.02 s (diffusion number 0.5), 1500 steps: 30 s, three viscous
    // times L²/ν. With balanced body forces solid-body rotation is reproduced to round-off;
    // with the cell Gauss gradient it was off by 0.49 % of Ω R on this grid.
    let mut case = closed_box(n, side, 1500, 0.02);
    case.rotating_regions.push(region("box", [0, 0, 0], [n, n, 1], centre, RPM_ONE_RAD));
    let manifest = Case::CollocatedFlow3d(case).run(&scratch.0).unwrap();
    assert!(manifest.capabilities.iter().any(|c| c == "rotating_regions_v1"));
    let (u, p) = fields(&scratch.0);
    let h = side / n as f64;
    let centre_of = |c: usize| Vec3::new((c % n) as f64 * h + h / 2.0, (c / n) as f64 * h + h / 2.0, 0.0) - centre;
    let mut worst: f64 = 0.0;
    for (c, v) in u.iter().enumerate() {
        let r = centre_of(c);
        let exact = Vec3::new(-r.y, r.x, 0.0);
        worst = worst.max((*v - exact).norm());
    }
    let tip = 0.05; // Ω R
    eprintln!("solid-body rotation: max velocity error {:.3e} of Ω R", worst / tip);
    assert!(worst < 1e-6 * tip, "{worst}");
    // p = a + ½ Ω² r²: least-squares fit over the cells off the wall, whose
    // zero-gradient pressure boundary is only first order.
    let (mut sx, mut sy, mut sxx, mut sxy, mut k) = (0.0, 0.0, 0.0, 0.0, 0.0);
    for c in 0..n * n {
        let (i, j) = (c % n, c / n);
        if i == 0 || j == 0 || i == n - 1 || j == n - 1 {
            continue;
        }
        let r = centre_of(c);
        let x = r.x * r.x + r.y * r.y;
        sx += x;
        sy += p[c];
        sxx += x * x;
        sxy += x * p[c];
        k += 1.0;
    }
    let slope = (k * sxy - sx * sy) / (k * sxx - sx * sx);
    eprintln!("pressure curvature {slope:.5}, exact 0.5");
    assert!((slope - 0.5).abs() < 1e-6, "{slope}");
    let report: serde_json::Value = serde_json::from_slice(&std::fs::read(scratch.0.join("rotating_regions.json")).unwrap()).unwrap();
    let swirl = report["regions"][0]["swirl_ratio"].as_f64().unwrap();
    assert!((swirl - 1.0).abs() < 0.01, "{report}");
}

#[test]
fn a_rotor_in_a_partial_region_stirs_the_fluid_its_way() {
    let scratch = Scratch::new();
    let (n, side) = (16, 0.16);
    let centre = Vec3::new(0.08, 0.08, 0.0);
    let mut swirls = vec![];
    for (name, rpm) in [("forward", 2.0 * RPM_ONE_RAD), ("reverse", -2.0 * RPM_ONE_RAD)] {
        let mut case = closed_box(n, side, 1500, 0.01);
        // A 6 × 2 bar through the centre, inside an 8 × 8 rotating region.
        case.obstacles.push(GridObstacle { name: "rotor".into(), min: [5, 7, 0], max: [11, 9, 1] });
        case.rotating_regions.push(region("hub", [4, 4, 0], [12, 12, 1], centre, rpm));
        let dir = scratch.0.join(name);
        Case::CollocatedFlow3d(case).run(&dir).unwrap();
        let report: serde_json::Value = serde_json::from_slice(&std::fs::read(dir.join("rotating_regions.json")).unwrap()).unwrap();
        assert_eq!(report["regions"][0]["cells"], 64 - 12);
        let mesh: serde_json::Value = serde_json::from_slice(&std::fs::read(dir.join("mesh_report.json")).unwrap()).unwrap();
        assert!(mesh["relative_boundary_flux_imbalance"].as_f64().unwrap() < 1e-6, "{mesh}");
        swirls.push(report["regions"][0]["swirl_ratio"].as_f64().unwrap());
    }
    eprintln!("swirl ratios {swirls:?}");
    // The fluid turns with the rotor, slower than it, whichever way it spins.
    assert!(swirls.iter().all(|s| *s > 0.05 && *s < 1.0), "{swirls:?}");
    assert!((swirls[0] - swirls[1]).abs() < 0.05, "symmetric: {swirls:?}");
}

#[test]
fn bad_rotating_regions_are_refused() {
    let scratch = Scratch::new();
    let centre = Vec3::new(0.05, 0.05, 0.0);
    let run = |r: Vec<RotatingRegion>, dt: f64, name: &str| {
        let mut case = closed_box(10, 0.1, 2, dt);
        case.rotating_regions = r;
        Case::CollocatedFlow3d(case).run(&scratch.0.join(name))
    };
    assert!(run(vec![region("a", [0, 0, 0], [10, 10, 1], centre, RPM_ONE_RAD)], 0.1, "ok").is_ok());
    assert!(run(vec![region("a", [0, 0, 0], [11, 10, 1], centre, 1.0)], 0.1, "outside").is_err());
    assert!(run(vec![region("a", [3, 0, 0], [3, 10, 1], centre, 1.0)], 0.1, "empty").is_err());
    let mut flat = region("a", [0, 0, 0], [5, 5, 1], centre, 1.0);
    flat.axis = Vec3::ZERO;
    assert!(run(vec![flat], 0.1, "axis").is_err());
    assert!(run(vec![region("a", [0, 0, 0], [5, 5, 1], centre, 1.0), region("b", [4, 4, 0], [8, 8, 1], centre, 1.0)], 0.1, "overlap").is_err());
    assert!(run(vec![region("a", [0, 0, 0], [5, 5, 1], centre, 1.0), region("a", [6, 6, 0], [8, 8, 1], centre, 1.0)], 0.1, "names").is_err());
    // 3 rad/s × 0.1 s = 0.3 rad per step: too coarse for the explicit frame term.
    assert!(run(vec![region("a", [0, 0, 0], [5, 5, 1], centre, 3.0 * RPM_ONE_RAD)], 0.1, "fast").is_err());
    // JSON round trip; unknown keys are refused.
    let mut case = closed_box(10, 0.1, 2, 0.1);
    case.rotating_regions.push(region("a", [0, 0, 0], [5, 5, 1], centre, 100.0));
    let json = serde_json::to_value(&case).unwrap();
    assert_eq!(json["rotating_regions"][0]["rpm"], 100.0);
    let mut extra = json.clone();
    extra["rotating_regions"][0]["slip"] = serde_json::json!(true);
    assert!(serde_json::from_value::<CollocatedFlowCase>(extra).is_err());
}
