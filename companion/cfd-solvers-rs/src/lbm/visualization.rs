//! LBM post-processing and lightweight viewers: legacy VTK, live JSON state, vorticity,
//! Q-criterion, slices, streamlines, PGM and a self-contained HTML viewer, ported from upstream
//! `cfd/solvers/lbm/visualization` (MIT, cfd_solvers contributors). Numbers are written with
//! C++ iostream `%.Ng` semantics so the files are byte-identical to upstream output.
// Index loops mirror the upstream lattice stencils direction by direction.
#![allow(clippy::needless_range_loop)]
use crate::{require, Result, SolverError};
use std::{fmt::Write as _, path::Path};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vector2 {
    pub x: f64,
    pub y: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Vector3 {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// `printf("%.{precision}g")`, as used by iostreams with `setprecision` and default float format.
pub fn format_general(value: f64, precision: usize) -> String {
    if value.is_nan() {
        return "nan".into();
    }
    if value.is_infinite() {
        return if value > 0.0 {
            "inf".into()
        } else {
            "-inf".into()
        };
    }
    let p = precision.max(1);
    if value == 0.0 {
        return if value.is_sign_negative() {
            "-0".into()
        } else {
            "0".into()
        };
    }
    // Round to p significant digits first; the exponent of the rounded value decides the style.
    let sci = format!("{:.*e}", p - 1, value);
    let (mantissa, exponent) = sci.split_once('e').unwrap_or((&sci, "0"));
    let exp: i32 = exponent.parse().unwrap_or(0);
    let strip = |s: &str| -> String {
        if s.contains('.') {
            s.trim_end_matches('0').trim_end_matches('.').to_string()
        } else {
            s.to_string()
        }
    };
    if exp < -4 || exp >= p as i32 {
        let sign = if exp < 0 { '-' } else { '+' };
        format!("{}e{}{:02}", strip(mantissa), sign, exp.abs())
    } else {
        let decimals = (p as i32 - 1 - exp).max(0) as usize;
        strip(&format!("{:.*}", decimals, value))
    }
}

fn g17(v: f64) -> String {
    format_general(v, 17)
}

fn write_text(path: &Path, text: &str, what: &str) -> Result<()> {
    std::fs::write(path, text).map_err(|e| SolverError::Io(format!("cannot create {what}: {e}")))
}

fn create_parent(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
        std::fs::create_dir_all(parent)?;
    }
    Ok(())
}

/// Legacy ASCII VTK structured points (2-D, unit spacing) with one scalar and optional velocity.
pub fn write_vtk_structured_2d(
    path: impl AsRef<Path>,
    nx: usize,
    ny: usize,
    scalar: &[f64],
    velocity: &[Vector2],
    scalar_name: &str,
) -> Result<()> {
    require(
        nx * ny == scalar.len() && (velocity.is_empty() || velocity.len() == scalar.len()),
        "LBM VTK field size mismatch",
    )?;
    let mut o = format!(
        "# vtk DataFile Version 3.0\ncfd_solvers LBM\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS {nx} {ny} 1\nORIGIN 0 0 0\nSPACING 1 1 1\nPOINT_DATA {}\nSCALARS {scalar_name} double 1\nLOOKUP_TABLE default\n",
        scalar.len()
    );
    for v in scalar {
        let _ = writeln!(o, "{}", g17(*v));
    }
    if !velocity.is_empty() {
        o.push_str("VECTORS velocity double\n");
        for v in velocity {
            let _ = writeln!(o, "{} {} 0", g17(v.x), g17(v.y));
        }
    }
    write_text(path.as_ref(), &o, "LBM VTK")
}

/// Legacy ASCII VTK structured points (3-D) with optional velocity and solid mask.
#[allow(clippy::too_many_arguments)]
pub fn write_vtk_structured_3d(
    path: impl AsRef<Path>,
    nx: usize,
    ny: usize,
    nz: usize,
    scalar: &[f64],
    velocity: &[Vector3],
    solid: &[u8],
    scalar_name: &str,
) -> Result<()> {
    let cells = nx * ny * nz;
    require(
        nx > 0
            && ny > 0
            && nz > 0
            && scalar.len() == cells
            && (velocity.is_empty() || velocity.len() == cells)
            && (solid.is_empty() || solid.len() == cells),
        "LBM 3-D VTK field size mismatch",
    )?;
    let path = path.as_ref();
    create_parent(path)?;
    let mut o = format!(
        "# vtk DataFile Version 3.0\ncfd_solvers LBM 3-D\nASCII\nDATASET STRUCTURED_POINTS\nDIMENSIONS {nx} {ny} {nz}\nORIGIN 0 0 0\nSPACING 1 1 1\nPOINT_DATA {cells}\nSCALARS {scalar_name} double 1\nLOOKUP_TABLE default\n"
    );
    for v in scalar {
        let _ = writeln!(o, "{}", g17(*v));
    }
    if !velocity.is_empty() {
        o.push_str("VECTORS velocity double\n");
        for v in velocity {
            let _ = writeln!(o, "{} {} {}", g17(v.x), g17(v.y), g17(v.z));
        }
    }
    if !solid.is_empty() {
        o.push_str("SCALARS solid unsigned_char 1\nLOOKUP_TABLE default\n");
        for v in solid {
            let _ = writeln!(o, "{v}");
        }
    }
    write_text(path, &o, "3-D LBM VTK")
}

/// Live-viewer JSON state (values at 9 significant digits), published atomically via a `.tmp` rename.
#[allow(clippy::too_many_arguments)]
pub fn write_json_state_3d(
    path: impl AsRef<Path>,
    nx: usize,
    ny: usize,
    nz: usize,
    step: u64,
    scalar: &[f64],
    solid: &[u8],
    scalar_name: &str,
) -> Result<()> {
    let cells = nx * ny * nz;
    require(
        nx > 0
            && ny > 0
            && nz > 0
            && scalar.len() == cells
            && (solid.is_empty() || solid.len() == cells),
        "LBM 3-D live state size mismatch",
    )?;
    let path = path.as_ref();
    create_parent(path)?;
    let values: Vec<String> = scalar.iter().map(|v| format_general(*v, 9)).collect();
    let mut o = format!("{{\n  \"dimensions\": [{nx}, {ny}, {nz}],\n  \"step\": {step},\n  \"field\": \"{scalar_name}\",\n  \"values\": [{}]", values.join(", "));
    if !solid.is_empty() {
        let s: Vec<String> = solid.iter().map(|v| v.to_string()).collect();
        let _ = write!(o, ",\n  \"solid\": [{}]", s.join(", "));
    }
    o.push_str("\n}\n");
    let mut tmp = path.as_os_str().to_owned();
    tmp.push(".tmp");
    let tmp = std::path::PathBuf::from(tmp);
    write_text(&tmp, &o, "3-D LBM live state")?;
    if std::fs::rename(&tmp, path).is_err() {
        let _ = std::fs::remove_file(path);
        std::fs::rename(&tmp, path)
            .map_err(|e| SolverError::Io(format!("cannot publish 3-D LBM live state: {e}")))?;
    }
    Ok(())
}

/// Interior central-difference vorticity `dv/dx - du/dy`; boundary rows/columns stay zero.
pub fn vorticity_2d(
    nx: usize,
    ny: usize,
    velocity: &[Vector2],
    dx: f64,
    dy: f64,
) -> Result<Vec<f64>> {
    require(
        nx >= 3 && ny >= 3 && velocity.len() == nx * ny && dx > 0.0 && dy > 0.0,
        "invalid LBM vorticity grid",
    )?;
    let id = |x: usize, y: usize| y * nx + x;
    let mut w = vec![0.0; nx * ny];
    for y in 1..ny - 1 {
        for x in 1..nx - 1 {
            let dvdx = (velocity[id(x + 1, y)].y - velocity[id(x - 1, y)].y) / (2.0 * dx);
            let dudy = (velocity[id(x, y + 1)].x - velocity[id(x, y - 1)].x) / (2.0 * dy);
            w[id(x, y)] = dvdx - dudy;
        }
    }
    Ok(w)
}

/// Interior Q-criterion `0.5 (|Omega|^2 - |S|^2)`; boundary cells stay zero.
#[allow(clippy::too_many_arguments)]
pub fn q_criterion_3d(
    nx: usize,
    ny: usize,
    nz: usize,
    velocity: &[Vector3],
    dx: f64,
    dy: f64,
    dz: f64,
) -> Result<Vec<f64>> {
    require(
        nx >= 3
            && ny >= 3
            && nz >= 3
            && velocity.len() == nx * ny * nz
            && dx > 0.0
            && dy > 0.0
            && dz > 0.0,
        "invalid LBM Q-criterion grid",
    )?;
    let at = |x: usize, y: usize, z: usize| velocity[(z * ny + y) * nx + x];
    let mut q = vec![0.0; nx * ny * nz];
    for z in 1..nz - 1 {
        for y in 1..ny - 1 {
            for x in 1..nx - 1 {
                let (xp, xm, yp, ym, zp, zm) = (
                    at(x + 1, y, z),
                    at(x - 1, y, z),
                    at(x, y + 1, z),
                    at(x, y - 1, z),
                    at(x, y, z + 1),
                    at(x, y, z - 1),
                );
                let g = [
                    [
                        (xp.x - xm.x) / (2.0 * dx),
                        (yp.x - ym.x) / (2.0 * dy),
                        (zp.x - zm.x) / (2.0 * dz),
                    ],
                    [
                        (xp.y - xm.y) / (2.0 * dx),
                        (yp.y - ym.y) / (2.0 * dy),
                        (zp.y - zm.y) / (2.0 * dz),
                    ],
                    [
                        (xp.z - xm.z) / (2.0 * dx),
                        (yp.z - ym.z) / (2.0 * dy),
                        (zp.z - zm.z) / (2.0 * dz),
                    ],
                ];
                let (mut s2, mut o2) = (0.0, 0.0);
                for i in 0..3 {
                    for j in 0..3 {
                        let sij = 0.5 * (g[i][j] + g[j][i]);
                        let oij = 0.5 * (g[i][j] - g[j][i]);
                        s2 += sij * sij;
                        o2 += oij * oij;
                    }
                }
                q[(z * ny + y) * nx + x] = 0.5 * (o2 - s2);
            }
        }
    }
    Ok(q)
}

/// Row `y` of a row-major `nx * ny` field.
pub fn horizontal_slice(nx: usize, ny: usize, field: &[f64], y: usize) -> Result<Vec<f64>> {
    require(field.len() == nx * ny && y < ny, "invalid LBM slice")?;
    Ok(field[y * nx..(y + 1) * nx].to_vec())
}

/// Fixed-length unit steps along nearest-lower-cell velocity until leaving the domain, hitting a
/// stagnant cell or reaching `max_points`.
pub fn streamline_2d(
    nx: usize,
    ny: usize,
    velocity: &[Vector2],
    start: Vector2,
    step: f64,
    max_points: usize,
) -> Result<Vec<Vector2>> {
    require(
        velocity.len() == nx * ny && step > 0.0 && max_points > 0,
        "invalid streamline controls",
    )?;
    let sample = |p: Vector2| {
        if p.x < 0.0 || p.y < 0.0 || p.x > (nx - 1) as f64 || p.y > (ny - 1) as f64 {
            return Vector2::default();
        }
        let x = p.x.clamp(0.0, (nx - 1) as f64) as usize;
        let y = p.y.clamp(0.0, (ny - 1) as f64) as usize;
        velocity[y * nx + x]
    };
    let mut out = Vec::new();
    let mut p = start;
    for _ in 0..max_points {
        if p.x < 0.0 || p.y < 0.0 || p.x >= nx as f64 || p.y >= ny as f64 {
            break;
        }
        out.push(p);
        let v = sample(p);
        let m = v.x.hypot(v.y);
        if m < 1e-14 {
            break;
        }
        p.x += step * v.x / m;
        p.y += step * v.y / m;
    }
    Ok(out)
}

/// ASCII PGM (P2) with min/max scaling to 0..255.
pub fn write_pgm_scalar(path: impl AsRef<Path>, nx: usize, ny: usize, field: &[f64]) -> Result<()> {
    require(
        field.len() == nx * ny && !field.is_empty(),
        "invalid PGM field",
    )?;
    let mn = field.iter().copied().fold(f64::INFINITY, f64::min);
    let mx = field.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let span = mx - mn;
    let mut o = format!("P2\n{nx} {ny}\n255\n");
    for v in field {
        let q = if span > 0.0 {
            (255.0 * (v - mn) / span).round() as i64
        } else {
            0
        };
        let _ = write!(o, "{} ", q.clamp(0, 255));
    }
    write_text(path.as_ref(), &o, "PGM")
}

const HTML_HEAD: &str = "<!doctype html><meta charset=\"utf-8\"><title>cfd_solvers LBM viewer</title>\n<style>html,body{margin:0;height:100%;background:#111;color:#ddd;font:14px sans-serif;overflow:hidden}#hud{position:fixed;left:12px;top:10px;background:#000b;padding:8px;border-radius:6px}canvas{width:100%;height:100%;image-rendering:pixelated;cursor:grab}</style>\n<canvas id=\"view\"></canvas><div id=\"hud\">wheel: zoom \u{b7} drag: pan<br><span id=\"value\"></span></div><script>\n";

const HTML_SCRIPT: &str = r#"const c=document.getElementById('view'),g=c.getContext('2d'),v=document.getElementById('value');let zoom=1,ox=0,oy=0,drag=false,lx=0,ly=0;
const lo=Math.min(...field),hi=Math.max(...field),span=hi-lo||1,off=document.createElement('canvas');off.width=nx;off.height=ny;const og=off.getContext('2d'),im=og.createImageData(nx,ny);
for(let i=0;i<field.length;i++){const t=Math.max(0,Math.min(1,(field[i]-lo)/span)),r=Math.round(255*t),b=Math.round(255*(1-t)),q=4*i;im.data[q]=r;im.data[q+1]=Math.round(255*(1-Math.abs(2*t-1)));im.data[q+2]=b;im.data[q+3]=255;}og.putImageData(im,0,0);
function resize(){c.width=innerWidth*devicePixelRatio;c.height=innerHeight*devicePixelRatio;draw()}function draw(){g.setTransform(1,0,0,1,0,0);g.clearRect(0,0,c.width,c.height);const base=Math.min(c.width/nx,c.height/ny)*zoom,w=nx*base,h=ny*base,x=(c.width-w)/2+ox,y=(c.height-h)/2+oy;g.imageSmoothingEnabled=false;g.drawImage(off,x,y,w,h);c._map={x,y,w,h}}
c.addEventListener('wheel',e=>{e.preventDefault();zoom*=Math.exp(-e.deltaY*.001);zoom=Math.max(.1,Math.min(100,zoom));draw()},{passive:false});c.addEventListener('pointerdown',e=>{drag=true;lx=e.clientX*devicePixelRatio;ly=e.clientY*devicePixelRatio;c.setPointerCapture(e.pointerId)});c.addEventListener('pointermove',e=>{const x=e.clientX*devicePixelRatio,y=e.clientY*devicePixelRatio;if(drag){ox+=x-lx;oy+=y-ly;lx=x;ly=y;draw()}const m=c._map,ix=Math.floor((x-m.x)/m.w*nx),iy=Math.floor((y-m.y)/m.h*ny);v.textContent=(ix>=0&&iy>=0&&ix<nx&&iy<ny)?`(${ix}, ${iy}) = ${field[iy*nx+ix]}`:''});c.addEventListener('pointerup',()=>drag=false);addEventListener('resize',resize);resize();
</script>"#;

/// Offline single-file HTML viewer with zoom, pan and value read-out.
pub fn write_html_scalar_viewer(
    path: impl AsRef<Path>,
    nx: usize,
    ny: usize,
    field: &[f64],
) -> Result<()> {
    require(
        field.len() == nx * ny && !field.is_empty(),
        "invalid HTML viewer field",
    )?;
    let values: Vec<String> = field.iter().map(|v| g17(*v)).collect();
    let o = format!(
        "{HTML_HEAD}const nx={nx},ny={ny},field=[{}];\n{HTML_SCRIPT}",
        values.join(",")
    );
    write_text(path.as_ref(), &o, "LBM HTML viewer")
}
