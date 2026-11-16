//! Value-by-value comparison with the C++ upstream solvers/lbm/distributed_pull. Goldens come from
//! `parity/oracle_lbm_dist.cpp`.
use cfd_solvers::{
    distributed::{choose_process_grid, make_brick, Extent3},
    lbm::*,
};
use serde_json::Value;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-lbm-dist.json")).unwrap()
}

#[track_caller]
fn close(name: &str, rust: &[f64], rel: f64, abs: f64) {
    let g = goldens();
    let cpp: Vec<f64> = g[name]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect();
    assert_eq!(rust.len(), cpp.len(), "{name}: length");
    for (i, (r, c)) in rust.iter().zip(&cpp).enumerate() {
        assert!(
            (r - c).abs() <= rel * c.abs() + abs,
            "{name}[{i}]: rust {r} vs C++ {c}"
        );
    }
}

fn run<D: LatticeDescriptor, const Q: usize>(
    global: Extent3,
    ranks: usize,
    selective: bool,
) -> Vec<f64> {
    let cfg = DistributedLbmConfig {
        global,
        tau: 0.7,
        acceleration_x: 1e-5,
        acceleration_y: -2e-6,
        acceleration_z: 0.0,
    };
    let grid = choose_process_grid(global, ranks, D::DIMENSIONS as u32).unwrap();
    let mut blocks: Vec<DistributedPullBlock<D, Q>> = (0..grid.cells())
        .map(|r| {
            let mut b =
                DistributedPullBlock::new(cfg, make_brick(global, grid, r).unwrap()).unwrap();
            b.initialize_taylor_green(0.04).unwrap();
            b
        })
        .collect();
    let mut out = vec![grid.x as f64, grid.y as f64, grid.z as f64];
    for _ in 0..6 {
        if selective {
            virtual_selective_step(&mut blocks, [true; 3]).unwrap();
        } else {
            virtual_step(&mut blocks, [true; 3]).unwrap();
        }
    }
    out.push(blocks.iter().map(|b| b.local_mass().unwrap()).sum());
    let m = gather_virtual_macroscopic(&blocks).unwrap();
    for i in 0..m.rho.len() {
        out.extend([
            m.rho[i] as f64,
            m.ux[i] as f64,
            m.uy[i] as f64,
            m.uz[i] as f64,
        ]);
    }
    out
}

#[test]
fn distributed_pull_matches_upstream() {
    // f32 kernels: bit-for-bit up to libm sin/cos differences in the initializer.
    close(
        "d3q19_full",
        &run::<D3Q19Descriptor, 19>(Extent3::new(8, 6, 4), 4, false),
        2e-6,
        1e-8,
    );
    close(
        "d3q19_sel",
        &run::<D3Q19Descriptor, 19>(Extent3::new(8, 6, 4), 4, true),
        2e-6,
        1e-8,
    );
    close(
        "d3q27_sel",
        &run::<D3Q27Descriptor, 27>(Extent3::new(6, 5, 4), 3, true),
        2e-6,
        1e-8,
    );
    close(
        "d2q9_full",
        &run::<D2Q9InPlaceDescriptor, 9>(Extent3::new(12, 10, 1), 4, false),
        2e-6,
        1e-8,
    );
}

#[test]
fn checkpoint_layout_matches_upstream() {
    let global = Extent3::new(4, 3, 2);
    let cfg = DistributedLbmConfig {
        global,
        tau: 0.8,
        acceleration_x: 1e-4,
        acceleration_y: 0.0,
        acceleration_z: 0.0,
    };
    let mut b = DistributedPullBlock::<D3Q19Descriptor, 19>::new(
        cfg,
        make_brick(global, Extent3::new(1, 1, 1), 0).unwrap(),
    )
    .unwrap();
    b.initialize_taylor_green(0.02).unwrap();
    b.compute_all().unwrap();
    b.finish_step();
    let path = std::env::temp_dir().join(format!("cfd_lbm_dist_{}.chk", std::process::id()));
    b.save_checkpoint(&path).unwrap();
    let bytes = std::fs::read(&path).unwrap();
    let g = goldens();
    let cpp: Vec<u8> = g["checkpoint"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap() as u8)
        .collect();
    assert_eq!(bytes.len(), cpp.len());
    // Header and step are exact; payload floats agree to f32 rounding of libm.
    assert_eq!(&bytes[..100], &cpp[..100]);
    for (r, c) in bytes[100..].chunks(4).zip(cpp[100..].chunks(4)) {
        let (a, bb) = (
            f32::from_le_bytes(r.try_into().unwrap()),
            f32::from_le_bytes(c.try_into().unwrap()),
        );
        assert!((a - bb).abs() <= 1e-6 * bb.abs());
    }
    let mut c = DistributedPullBlock::<D3Q19Descriptor, 19>::new(
        cfg,
        make_brick(global, Extent3::new(1, 1, 1), 0).unwrap(),
    )
    .unwrap();
    c.load_checkpoint(&path).unwrap();
    assert_eq!(c.time_step(), 1);
    assert_eq!(
        c.compute_macroscopic().unwrap(),
        b.compute_macroscopic().unwrap()
    );
    let other = DistributedLbmConfig { tau: 0.9, ..cfg };
    let mut d = DistributedPullBlock::<D3Q19Descriptor, 19>::new(
        other,
        make_brick(global, Extent3::new(1, 1, 1), 0).unwrap(),
    )
    .unwrap();
    assert!(d.load_checkpoint(&path).is_err());
    std::fs::remove_file(&path).ok();
}

#[test]
fn decomposition_is_invisible() {
    let one = run::<D3Q19Descriptor, 19>(Extent3::new(8, 6, 4), 1, false);
    let four = run::<D3Q19Descriptor, 19>(Extent3::new(8, 6, 4), 4, true);
    for (a, b) in one[3..].iter().zip(&four[3..]) {
        assert!((a - b).abs() <= 1e-6 * a.abs().max(1e-6));
    }
}
