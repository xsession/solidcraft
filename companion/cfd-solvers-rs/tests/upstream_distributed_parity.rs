//! Value-by-value comparison with the C++ upstream distributed/{cartesian, halo_grid,
//! selective_halo, device_assignment}. Goldens come from `parity/oracle_distributed.cpp`.
#![allow(clippy::needless_range_loop)]
use cfd_solvers::{
    distributed::*,
    lbm::{D3Q19Descriptor as D19, D3Q27Descriptor as D27},
};
use serde_json::Value;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-distributed.json")).unwrap()
}

#[track_caller]
fn same(name: &str, rust: &[f64]) {
    let g = goldens();
    let cpp: Vec<f64> = g[name]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect();
    assert_eq!(rust, cpp.as_slice(), "{name}");
}

fn opt(r: Option<usize>) -> f64 {
    r.map_or(-1.0, |v| v as f64)
}

#[test]
fn cartesian_topology() {
    let mut g = vec![];
    for (e, r, d) in [
        (Extent3::new(64, 64, 64), 8, 3),
        (Extent3::new(120, 40, 10), 12, 3),
        (Extent3::new(100, 50, 1), 6, 2),
        (Extent3::new(7, 5, 3), 30, 3),
        (Extent3::new(256, 32, 32), 16, 3),
    ] {
        let p = choose_process_grid(e, r, d).unwrap();
        g.extend([p.x as f64, p.y as f64, p.z as f64]);
    }
    let pg = Extent3::new(3, 2, 2);
    for r in 0..12 {
        let c = coords_from_rank(r, pg).unwrap();
        g.extend([
            c.x as f64,
            c.y as f64,
            c.z as f64,
            rank_from_coords(c, pg).unwrap() as f64,
        ]);
    }
    for c in 0..4 {
        let (b, n) = partition_axis(10, c, 4).unwrap();
        g.extend([b as f64, n as f64]);
    }
    let brick = make_brick(Extent3::new(10, 7, 5), pg, 7).unwrap();
    g.extend(
        [
            brick.begin.x,
            brick.begin.y,
            brick.begin.z,
            brick.extent.x,
            brick.extent.y,
            brick.extent.z,
        ]
        .map(|v| v as f64),
    );
    for o in neighbor_offsets() {
        g.extend([
            neighbor_tag(o) as f64,
            opt(neighbor_rank(&brick, o, [true; 3]).unwrap()),
            opt(neighbor_rank(&brick, o, [false, true, false]).unwrap()),
        ]);
    }
    same("cartesian", &g);
}

#[test]
fn halo_and_selective_packing() {
    let mut h = HaloSoA3D::<f64, 2>::new(Extent3::new(3, 2, 2), 1).unwrap();
    for c in 0..2 {
        for z in 0..4 {
            for y in 0..4 {
                for x in 0..5 {
                    h.set(
                        c,
                        x,
                        y,
                        z,
                        1000.0 * c as f64 + 100.0 * z as f64 + 10.0 * y as f64 + x as f64,
                    )
                    .unwrap();
                }
            }
        }
    }
    let mut hv = vec![];
    for o in neighbor_offsets() {
        let p = h.pack(o);
        hv.push(p.len() as f64);
        hv.extend(p);
    }
    let mut h2 = HaloSoA3D::<f64, 2>::new(Extent3::new(3, 2, 2), 1).unwrap();
    for o in neighbor_offsets() {
        let p: Vec<f64> = h.pack(o.opposite()).iter().map(|x| x + 0.5).collect();
        h2.unpack(o, &p).unwrap();
    }
    for c in 0..2 {
        for z in 0..4 {
            for y in 0..4 {
                for x in 0..5 {
                    hv.push(h2.get(c, x, y, z).unwrap());
                }
            }
        }
    }
    same("halo", &hv);

    let mut s = HaloSoA3D::<f32, 19>::new(Extent3::new(3, 2, 2), 1).unwrap();
    for c in 0..19 {
        for z in 0..4 {
            for y in 0..4 {
                for x in 0..5 {
                    s.set(c, x, y, z, (c * 100 + z * 20 + y * 5 + x) as f32)
                        .unwrap();
                }
            }
        }
    }
    let mut sv = vec![];
    for o in neighbor_offsets() {
        let cp = crossing_populations::<D19>(o);
        sv.push(cp.len() as f64);
        sv.extend(cp.iter().map(|d| *d as f64));
        sv.extend(
            pack_selective::<D19, f32, 19>(&s, o)
                .unwrap()
                .iter()
                .map(|v| *v as f64),
        );
        let e = selective_receive_entry::<D19>(s.interior_extent(), o, 1);
        sv.extend([e.direction, e.x, e.y, e.z].map(|v| v as f64));
    }
    let mut s2 = HaloSoA3D::<f32, 19>::new(Extent3::new(3, 2, 2), 1).unwrap();
    for o in neighbor_offsets() {
        let p = pack_selective::<D19, f32, 19>(&s, o.opposite()).unwrap();
        unpack_selective::<D19, f32, 19>(&mut s2, o, &p).unwrap();
    }
    for c in 0..19 {
        for z in 0..4 {
            for y in 0..4 {
                for x in 0..5 {
                    sv.push(s2.get(c, x, y, z).unwrap() as f64);
                }
            }
        }
    }
    for st in [
        selective_halo_traffic::<D19>(Extent3::new(32, 16, 8), 4),
        full_halo_traffic::<D19>(Extent3::new(32, 16, 8), 4),
        selective_halo_traffic::<D27>(Extent3::new(10, 10, 10), 8),
    ] {
        sv.extend([st.messages as f64, st.values as f64, st.bytes as f64]);
    }
    same("selective", &sv);
}

#[test]
fn device_assignment() {
    let mut dv = vec![];
    for (r, n, v) in [(0, 4, 2), (3, 4, 2), (1, 2, 8), (2, 3, 0)] {
        let a = assign_device(r, n, v).unwrap();
        dv.extend([
            a.device_ordinal as f64,
            if a.oversubscribed { 1.0 } else { 0.0 },
            a.visible_devices as f64,
        ]);
    }
    for (r, n, v, m) in [(0, 3, 8, 0), (2, 3, 8, 0), (1, 3, 8, 2), (3, 4, 2, 0)] {
        let gr = assign_device_group(r, n, v, m).unwrap();
        dv.push(gr.device_ordinals.len() as f64);
        dv.extend(gr.device_ordinals.iter().map(|o| *o as f64));
        dv.push(if gr.oversubscribed { 1.0 } else { 0.0 });
        for sl in split_rank_work(103, &gr) {
            dv.extend([sl.device_ordinal as f64, sl.begin as f64, sl.end as f64]);
        }
    }
    same("devices", &dv);
}
