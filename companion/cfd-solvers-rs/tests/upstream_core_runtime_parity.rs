//! Value-by-value comparison with the C++ upstream core runtime helpers and
//! `distributed/repartition`. Goldens come from `parity/oracle_core_runtime.cpp`.
use cfd_solvers::{core::*, distributed::*};
use serde_json::Value;

fn goldens() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-core-runtime.json")).unwrap()
}

#[track_caller]
fn close(name: &str, rust: &[f64], rel: f64, abs: f64) {
    let g = goldens();
    let cpp: Vec<f64> = g[name]
        .as_array()
        .unwrap_or_else(|| panic!("{name} missing"))
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

fn lap(n: usize) -> CsrMatrix {
    let (mut off, mut col, mut val) = (vec![0], vec![], vec![]);
    for r in 0..n {
        if r > 0 {
            col.push(r - 1);
            val.push(-1.0 - 0.01 * r as f64);
        }
        col.push(r);
        val.push(4.0 + 0.1 * r as f64);
        if r + 1 < n {
            col.push(r + 1);
            val.push(-1.2);
        }
        if r + 3 < n {
            col.push(r + 3);
            val.push(0.3);
        }
        off.push(col.len());
    }
    CsrMatrix::new(n, n, off, col, val).unwrap()
}

#[test]
fn runtime_helpers_match_upstream() {
    let p: Vec<f64> = (0..4)
        .flat_map(|r| {
            let s = partition_1d(23, r, 4, 2).unwrap();
            [s.begin as f64, s.end as f64, s.halo as f64]
        })
        .collect();
    close("partition_1d", &p, 0.0, 0.0);
    let a = lap(10);
    let part = DistributedCsrPartition::new(a.clone(), 3, 7).unwrap();
    let halo: Vec<f64> = part
        .halo_columns()
        .iter()
        .map(|c| 0.1 * *c as f64 + 1.0)
        .collect();
    let mut y = vec![0.0; 4];
    part.multiply(&[1.0, -2.0, 0.5, 3.0], &halo, &mut y)
        .unwrap();
    let mut d: Vec<f64> = part.halo_columns().iter().map(|c| *c as f64).collect();
    d.extend(y);
    close("dcsr", &d, 1e-15, 0.0);
    let mut asm = DistributedCsrAssembler::new(6, 6, vec![0, 2, 4, 6]).unwrap();
    for k in 0..20 {
        asm.add((k * 5) % 6, (k * 7) % 6, 0.1 * k as f64 - 0.7, k % 3)
            .unwrap();
    }
    let mut out = Vec::new();
    for r in 0..3 {
        let b = asm.build_owned(r).unwrap();
        out.extend(b.matrix.values());
        out.extend(b.matrix.columns().iter().map(|c| *c as f64));
    }
    let g = asm.build_global_reference().unwrap();
    out.extend(g.values());
    out.extend(g.row_offsets().iter().map(|o| *o as f64));
    close("assembly", &out, 1e-15, 1e-16);
    let b: Vec<f64> = (0..10).map(|i| (f64::from(i) + 1.0).sin()).collect();
    let mut x = vec![0.0; 10];
    let rr = mixed_precision_iterative_refinement(&a, &b, &mut x, 12, 30, 1e-12).unwrap();
    x.extend([
        f64::from(u8::from(rr.converged)),
        rr.outer_iterations as f64,
    ]);
    let g = goldens();
    let cpp: Vec<f64> = g["refine"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect();
    for i in 0..12 {
        assert!(
            (x[i] - cpp[i]).abs() <= 1e-13 * cpp[i].abs().max(1.0),
            "refine[{i}]"
        );
    }
    assert!(rr.relative_residual < 1e-12 && cpp[12] < 1e-12);

    let els = [
        [0, 1, 4],
        [1, 5, 4],
        [1, 2, 5],
        [2, 6, 5],
        [2, 3, 6],
        [3, 7, 6],
        [4, 5, 8],
        [5, 9, 8],
    ];
    let plan = partition_tri3_mesh(&els, 10, 3, &[1.0, 2.0, 1.0, 1.0, 3.0, 1.0, 1.0, 2.0]).unwrap();
    let mut mp: Vec<f64> = plan.element_boundaries.iter().map(|v| *v as f64).collect();
    for r in &plan.partitions {
        for (tag, list) in [
            (-1.0, &r.owned_nodes),
            (-2.0, &r.halo_nodes),
            (-3.0, &r.shared_owned_nodes),
        ] {
            mp.push(tag);
            mp.extend(list.iter().map(|n| *n as f64));
        }
    }
    close("mesh_partition", &mp, 0.0, 0.0);
    let rp = plan_contiguous_repartition(&[1.0, 5.0, 2.0, 2.0, 0.5, 7.0, 1.0, 1.0, 3.0, 2.0], 4, 2)
        .unwrap();
    let mut rv: Vec<f64> = rp.boundaries.iter().map(|v| *v as f64).collect();
    rv.extend(&rp.loads);
    rv.extend([rp.average_load, rp.imbalance]);
    let z = plan_contiguous_repartition(&[0.0; 7], 3, 1).unwrap();
    rv.extend(z.boundaries.iter().map(|v| *v as f64));
    rv.push(z.imbalance);
    for s in build_migration_segments(&[0, 3, 7, 10], &[0, 4, 5, 10]).unwrap() {
        rv.extend(
            [
                s.source_part,
                s.destination_part,
                s.global_begin,
                s.global_end,
            ]
            .map(|v| v as f64),
        );
    }
    for (c, p) in [(1.4, 1.1), (1.1, 1.0), (1.3, 1.28)] {
        rv.push(f64::from(u8::from(
            repartition_recommended(c, p, 1.15, 0.05).unwrap(),
        )));
    }
    close("repartition", &rv, 1e-15, 0.0);

    let mut op = DistributedReverseOperator::new(5).unwrap();
    op.add_partition(ReverseAdPartition {
        global_indices: vec![0, 1, 2],
        objective: Box::new(|_, v| {
            Ok(v[0].clone() * v[1].clone() + v[2].clone() * v[2].clone() * 0.5)
        }),
    })
    .unwrap();
    op.add_partition(ReverseAdPartition {
        global_indices: vec![2, 3, 4, 0],
        objective: Box::new(|_, v| {
            let d = v[0].clone() - v[1].clone();
            Ok((d.clone() * d) / (v[2].clone() + 2.0) + v[3].clone() * 3.0)
        }),
    })
    .unwrap();
    let mut gg = vec![0.0; 5];
    let ro = op.evaluate(&[0.3, -1.2, 2.0, 0.7, 1.5], &mut gg).unwrap();
    gg.extend([ro.objective, ro.partitions as f64]);
    close("reverse_ad", &gg, 1e-14, 0.0);
    close(
        "cpu_list",
        &parse_cpu_list(" 0-3, 8 ,6,10-11,2")
            .unwrap()
            .iter()
            .map(|c| f64::from(*c))
            .collect::<Vec<_>>(),
        0.0,
        0.0,
    );

    let mut db = KernelTuningDatabase::default();
    let dev = DeviceTuningKey {
        vendor: "Acme\tGPU".into(),
        name: "X1".into(),
        driver: "1.2".into(),
    };
    let cands = [
        KernelTuningCandidate {
            work_group_size: 64,
            vector_width: 1,
            fused_stages: 1,
        },
        KernelTuningCandidate {
            work_group_size: 128,
            vector_width: 2,
            fused_stages: 1,
        },
        KernelTuningCandidate {
            work_group_size: 256,
            vector_width: 4,
            fused_stages: 2,
        },
    ];
    db.autotune(&dev, "spmv", &cands, &mut |c| {
        1e-3 / c.work_group_size as f64 + 1e-6 * c.fused_stages as f64
    })
    .unwrap();
    db.upsert(KernelTuningRecord {
        device: DeviceTuningKey {
            vendor: "B".into(),
            name: "Y".into(),
            driver: "0".into(),
        },
        kernel_group: "stencil".into(),
        candidate: KernelTuningCandidate {
            work_group_size: 32,
            vector_width: 1,
            fused_stages: 3,
        },
        score_seconds: 0.00123456789,
    })
    .unwrap();
    assert_eq!(db.to_tsv(), goldens()["tuning_tsv"].as_str().unwrap());
    assert_eq!(
        device_tuning_key_string(&dev),
        goldens()["tuning_key"].as_str().unwrap()
    );
}
