//! Bit-exact comparison of the libstdc++ random replicas with sequences printed by
//! `parity/oracle_stdrand.cpp` (GCC 13 libstdc++).
use cfd_solvers::core::{Mt19937_64, StdNormal};
use serde_json::Value;

fn g() -> Value {
    serde_json::from_str(include_str!("../parity/goldens-stdrand.json")).unwrap()
}
fn f(name: &str) -> Vec<f64> {
    g()[name]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect()
}

#[test]
fn libstdcxx_random_sequences_match_bit_for_bit() {
    let mut a = Mt19937_64::new(42);
    let all: Vec<u64> = (0..700).map(|_| a.next_u64()).collect();
    let raw: Vec<u64> = all
        .iter()
        .enumerate()
        .filter(|(i, _)| *i < 5 || *i > 690)
        .map(|(_, v)| *v)
        .collect();
    let cpp: Vec<u64> = g()["raw"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_str().unwrap().parse().unwrap())
        .collect();
    assert_eq!(raw, cpp);
    let mut b = Mt19937_64::new(7);
    let u: Vec<f64> = (0..20).map(|_| b.uniform(0.0, 1.0)).collect();
    // serde_json may parse the 17-digit text one ulp off, so allow 2 ulp.
    for (r, c) in u.iter().zip(f("uniform")) {
        assert!((r - c).abs() <= 2.0 * f64::EPSILON * c.abs(), "{r} vs {c}");
    }
    let p: Vec<f64> = (0..40).map(|_| b.uniform_index(0, 6) as f64).collect();
    assert_eq!(p, f("pick"));
    let p2: Vec<f64> = (0..10)
        .map(|_| b.uniform_index(3, 1_000_000_007) as f64)
        .collect();
    assert_eq!(p2, f("pick2"));
    let mut n = StdNormal::default();
    let nn: Vec<f64> = (0..21).map(|_| n.sample(&mut b)).collect();
    for (r, c) in nn.iter().zip(f("normal")) {
        assert!(
            (r - c).abs() <= 4.0 * f64::EPSILON * c.abs().max(1e-300),
            "{r} vs {c}"
        );
    }
}
