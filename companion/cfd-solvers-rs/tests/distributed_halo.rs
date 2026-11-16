//! In-process halo exchange equals a periodic single-domain halo fill.
use cfd_solvers::distributed::*;

#[test]
fn in_process_exchange_matches_periodic_wrap() {
    let global = Extent3::new(6, 5, 4);
    let grid = choose_process_grid(global, 4, 3).unwrap();
    let bricks: Vec<Brick> = (0..grid.cells())
        .map(|r| make_brick(global, grid, r).unwrap())
        .collect();
    let value = |x: usize, y: usize, z: usize, c: usize| (c * 1000 + z * 100 + y * 10 + x) as f64;
    let mut grids: Vec<HaloSoA3D<f64, 2>> = bricks
        .iter()
        .map(|b| {
            let mut g = HaloSoA3D::new(b.extent, 1).unwrap();
            for c in 0..2 {
                for z in 0..b.extent.z {
                    for y in 0..b.extent.y {
                        for x in 0..b.extent.x {
                            g.set(
                                c,
                                x + 1,
                                y + 1,
                                z + 1,
                                value(b.begin.x + x, b.begin.y + y, b.begin.z + z, c),
                            )
                            .unwrap();
                        }
                    }
                }
            }
            g
        })
        .collect();
    exchange_halos_in_process(&mut grids, &bricks, [true; 3]).unwrap();
    let wrap = |v: i64, n: usize| v.rem_euclid(n as i64) as usize;
    for (g, b) in grids.iter().zip(&bricks) {
        let t = g.total_extent();
        for c in 0..2 {
            for z in 0..t.z {
                for y in 0..t.y {
                    for x in 0..t.x {
                        let gx = wrap(b.begin.x as i64 + x as i64 - 1, global.x);
                        let gy = wrap(b.begin.y as i64 + y as i64 - 1, global.y);
                        let gz = wrap(b.begin.z as i64 + z as i64 - 1, global.z);
                        assert_eq!(g.get(c, x, y, z).unwrap(), value(gx, gy, gz, c));
                    }
                }
            }
        }
    }
    assert!(HaloSoA3D::<f64, 1>::new(Extent3::new(2, 2, 2), 2).is_err());
    assert!(choose_process_grid(Extent3::new(2, 2, 2), 27, 3).is_err());
}
