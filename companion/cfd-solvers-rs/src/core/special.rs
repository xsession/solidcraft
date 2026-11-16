//! Error functions to near double precision without a libm dependency: a Maclaurin series for
//! small arguments and a modified-Lentz continued fraction for the complementary tail.

const FRAC_2_SQRT_PI: f64 = std::f64::consts::FRAC_2_SQRT_PI;

fn erf_series(x: f64) -> f64 {
    // erf(x) = 2/sqrt(pi) * sum_n (-1)^n x^(2n+1) / (n! (2n+1))
    let x2 = x * x;
    let mut term = x;
    let mut sum = x;
    for n in 1..200 {
        term *= -x2 / n as f64;
        let add = term / (2 * n + 1) as f64;
        sum += add;
        if add.abs() <= 1e-17 * sum.abs() {
            break;
        }
    }
    FRAC_2_SQRT_PI * sum
}

/// `erfc(x)` for `x >= 0.5` via the continued fraction
/// `sqrt(pi) e^{x^2} erfc(x) = 1/(x + (1/2)/(x + 1/(x + (3/2)/(x + ...))))`.
fn erfc_tail(x: f64) -> f64 {
    let tiny = 1e-300;
    let mut f = x;
    let mut c = x;
    let mut d = 0.0;
    for n in 1..20000 {
        let a = n as f64 * 0.5;
        d = x + a * d;
        if d.abs() < tiny {
            d = tiny;
        }
        c = x + a / c;
        if c.abs() < tiny {
            c = tiny;
        }
        d = 1.0 / d;
        let delta = c * d;
        f *= delta;
        if (delta - 1.0).abs() < 1e-16 {
            break;
        }
    }
    (-x * x).exp() / (f * std::f64::consts::PI.sqrt())
}

/// Complementary error function.
pub fn erfc(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x < 0.0 {
        return 2.0 - erfc(-x);
    }
    if x < 0.5 {
        1.0 - erf_series(x)
    } else if x > 27.3 {
        0.0
    } else {
        erfc_tail(x)
    }
}

/// Error function.
pub fn erf(x: f64) -> f64 {
    if x.is_nan() {
        return f64::NAN;
    }
    if x.abs() < 2.0 {
        erf_series(x)
    } else {
        x.signum() * (1.0 - erfc(x.abs()))
    }
}
