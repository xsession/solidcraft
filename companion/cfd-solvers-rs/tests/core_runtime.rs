//! Behavioural checks for the core runtime helpers.
use cfd_solvers::{core::*, distributed::*};
use std::sync::{
    atomic::{AtomicUsize, Ordering},
    Arc, Mutex,
};

#[test]
fn backends_parallel_helpers_and_profiler() {
    assert!(backend_available(RuntimeBackend::Serial) && !backend_available(RuntimeBackend::Mpi));
    assert_eq!(select_backend("serial").unwrap(), RuntimeBackend::Serial);
    assert!(matches!(
        select_backend("mpi"),
        Err(cfd_solvers::SolverError::Unsupported(_))
    ));
    assert!(parse_backend("quantum").is_err());
    assert!(compiled_backends().len() == 4);
    assert_eq!(parallel_sum(100, |i| i as f64), 4950.0);
    assert_eq!(parallel_max(10, |i| -(i as f64)), 0.0);
    let hits = AtomicUsize::new(0);
    parallel_for(37, |_| {
        hits.fetch_add(1, Ordering::Relaxed);
    });
    assert_eq!(hits.load(Ordering::Relaxed), 37);
    let p = Profiler::default();
    {
        let _a = ScopedProfile::new("assemble", &p);
    }
    p.record("assemble", 2.0);
    p.record("solve", 1.0);
    let s = p.snapshot();
    assert_eq!(
        (s[0].name.as_str(), s[0].calls, s[0].maximum_seconds),
        ("assemble", 2, 2.0)
    );
    p.set_enabled(false);
    p.record("ignored", 1.0);
    assert_eq!(p.snapshot().len(), 2);
    let mut stats = DeviceTransferStats::default();
    stats.record_host_to_device(10);
    stats.record_device_to_device(5);
    assert_eq!(
        (stats.explicit_transfer_bytes(), stats.host_transfer_bytes()),
        (15, 10)
    );
}

#[test]
fn arenas_packs_and_soa() {
    let mut arena = ScratchArena::<f64>::new(8);
    let a = arena.allocate(3).unwrap();
    let b = arena.allocate(5).unwrap();
    arena.slice_mut(b.clone())[4] = 2.0;
    assert_eq!(
        (arena.used(), arena.slice(a)[0], arena.slice(b)[4]),
        (8, 0.0, 2.0)
    );
    assert!(arena.allocate(1).is_err() && arena.reserve(16).is_err());
    arena.reset();
    arena.reserve(16).unwrap();
    let p = SimdPack::<f64, 4> {
        lane: [1.0, 2.0, 3.0, 4.0],
    };
    assert_eq!(((p + p) * 0.5 - p).sum(), 0.0);
    assert_eq!((p * p).sum(), 30.0);
    let mut soa = StaticSoA::<f64, 3>::new(4);
    soa[(2, 1)] = 5.0;
    assert_eq!(soa.component(2).unwrap()[1], 5.0);
    assert!(soa.component(3).is_err());
}

#[test]
fn task_graph_respects_dependencies() {
    for asynchronous in [false, true] {
        let log = Arc::new(Mutex::new(Vec::new()));
        let mut g = TaskGraph::default();
        let push = |name: &'static str| {
            let log = Arc::clone(&log);
            Box::new(move || {
                log.lock().unwrap().push(name);
                Ok(())
            }) as Task
        };
        let a = g.add("a", push("a"), vec![]).unwrap();
        let b = g.add("b", push("b"), vec![a]).unwrap();
        let c = g.add("c", push("c"), vec![a]).unwrap();
        g.add("d", push("d"), vec![b, c]).unwrap();
        g.run(asynchronous).unwrap();
        let log = log.lock().unwrap();
        assert_eq!((log[0], log[3], log.len()), ("a", "d", 4));
        assert!(g.add("bad", push("x"), vec![9]).is_err());
    }
    let mut g = TaskGraph::default();
    g.add(
        "fail",
        Box::new(|| Err(cfd_solvers::SolverError::NumericalFailure("boom".into()))),
        vec![],
    )
    .unwrap();
    assert!(g.run(true).is_err());
}

#[test]
fn tuning_database_round_trip_and_validation() {
    let mut db = KernelTuningDatabase::default();
    let dev = DeviceTuningKey {
        vendor: "v".into(),
        name: "n".into(),
        driver: "d".into(),
    };
    let rec = |s| KernelTuningRecord {
        device: dev.clone(),
        kernel_group: "k".into(),
        candidate: KernelTuningCandidate::default(),
        score_seconds: s,
    };
    db.upsert(rec(2.0)).unwrap();
    db.upsert(rec(3.0)).unwrap();
    assert_eq!(db.lookup(&dev, "k").unwrap().score_seconds, 2.0);
    let dir = std::env::temp_dir().join(format!("cfd-tune-{}", std::process::id()));
    let path = dir.join("sub/tune.tsv");
    db.save(&path).unwrap();
    assert_eq!(KernelTuningDatabase::load(&path).unwrap(), db);
    assert!(KernelTuningDatabase::load(dir.join("missing.tsv"))
        .unwrap()
        .records()
        .is_empty());
    assert!(KernelTuningDatabase::from_tsv("a\tb\n").is_err());
    assert!(KernelTuningDatabase::from_tsv("a\tb\tc\tk\t0\t1\t1\t1.0\n").is_err());
    assert!(db.upsert(rec(f64::NAN)).is_err());
    let _ = std::fs::remove_dir_all(dir);
}

#[test]
fn numa_planning_on_synthetic_topology() {
    let nodes = vec![
        NumaNodeInfo {
            id: 0,
            cpus: vec![0, 1],
        },
        NumaNodeInfo {
            id: 1,
            cpus: vec![4, 5, 6],
        },
    ];
    let spread =
        plan_numa_thread_placement_on(&nodes, 4, NumaPlacementPolicy::Spread, None).unwrap();
    assert_eq!(
        (spread.thread_cpus, spread.thread_nodes),
        (vec![0, 4, 1, 5], vec![0, 1, 0, 1])
    );
    let compact =
        plan_numa_thread_placement_on(&nodes, 0, NumaPlacementPolicy::Compact, None).unwrap();
    assert_eq!(compact.thread_cpus, vec![0, 1, 4, 5, 6]);
    let explicit =
        plan_numa_thread_placement_on(&nodes, 4, NumaPlacementPolicy::ExplicitNode, Some(1))
            .unwrap();
    assert_eq!(explicit.thread_cpus, vec![4, 5, 6, 4]);
    assert!(
        plan_numa_thread_placement_on(&nodes, 4, NumaPlacementPolicy::ExplicitNode, Some(7))
            .is_err()
    );
    assert!(
        !plan_numa_thread_placement(2, NumaPlacementPolicy::Inherit, None)
            .unwrap()
            .binds_threads()
    );
    assert!(!discover_numa_nodes().is_empty());
    assert!(parse_cpu_list("3-1").is_err() && parse_cpu_list("1,,2").is_err());
    assert_eq!(parse_cpu_list("0-2,").unwrap(), vec![0, 1, 2]);
    let mut v = vec![0.0; 5];
    numa_first_touch(&mut v, &1.5);
    assert!(v.iter().all(|x| *x == 1.5));
}

#[test]
fn distributed_linear_algebra_reproduces_serial() {
    let n: usize = 12;
    let (mut off, mut col, mut val) = (vec![0], vec![], vec![]);
    for r in 0..n {
        for c in [r.wrapping_sub(1), r, r + 1] {
            if c < n {
                col.push(c);
                val.push(if c == r { 4.0 } else { -1.0 });
            }
        }
        off.push(col.len());
    }
    let a = CsrMatrix::new(n, n, off, col, val).unwrap();
    let x: Vec<f64> = (0..n).map(|i| (i as f64).cos()).collect();
    let mut serial = vec![0.0; n];
    a.multiply(&x, &mut serial).unwrap();
    let plan = plan_contiguous_repartition(&[1.0; 12], 3, 1).unwrap();
    for p in 0..3 {
        let (b, e) = (plan.boundaries[p], plan.boundaries[p + 1]);
        let part = DistributedCsrPartition::new(a.clone(), b, e).unwrap();
        let halo: Vec<f64> = part.halo_columns().iter().map(|c| x[*c]).collect();
        let mut y = vec![0.0; e - b];
        part.multiply(&x[b..e], &halo, &mut y).unwrap();
        assert_eq!(y, serial[b..e]);
    }
    let gpu = SyclCsrLinearAlgebra::new(&a);
    let mut sol = vec![0.0; n];
    assert!(
        gpu.conjugate_gradient(&serial, &mut sol, 200, 1e-12)
            .unwrap()
            .converged
    );
    assert!(sol.iter().zip(&x).all(|(a, b)| (a - b).abs() < 1e-9));
    let mut sol = vec![0.0; n];
    assert!(
        gpu.bicgstab(&serial, &mut sol, 200, 1e-6)
            .unwrap()
            .converged
    );
    let mut x0 = vec![0.0; n];
    assert!(
        mixed_precision_iterative_refinement(&a, &serial, &mut x0, 12, 50, 1e-12)
            .unwrap()
            .converged
    );
    let amg = AmgPreconditionerAdapter::new(Box::new(|r, z| {
        z.copy_from_slice(r);
        Ok(())
    }));
    let mut z = vec![0.0; 2];
    amg.apply(&[1.0, 2.0], &mut z).unwrap();
    assert_eq!(z, vec![1.0, 2.0]);
    assert!(amg.apply(&[1.0], &mut z).is_err());
    assert_eq!(build_migration_segments(&[0, 5], &[0, 5]).unwrap(), vec![]);
}
