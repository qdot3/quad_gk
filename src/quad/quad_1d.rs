use std::ops::Range;

use crate::{Integral, QuadConfig, QuadGKCoef};

/// 1d numerical integration by Gauss-Kronrod quadrature rule.
pub fn quad_gk_1d<const NX: usize>(
    integrand: impl Fn(f64) -> f64,
    start: f64,
    end: f64,
    coef: impl QuadGKCoef<{ NX }>,
) -> Integral {
    //quad_gk_1d_inner(integrand, start, end, &coef, 10_000, 1e-10)
    todo!()
}

pub fn quadgk_1d<const N: usize>(
    integrand: impl Fn(f64) -> f64,
    range: Range<f64>,
    config: QuadConfig<N>,
) -> Integral {
    let QuadConfig {
        rel_tol,
        max_interval_count,
        coef,
    } = config;
    let Range {start, end} = range;

    let mut interval_count = 1;
    loop {
        let interval_width = (end - start) / interval_count as f64;
        let integral: Integral = (0..interval_count)
            .map(|i| {
                let r = 0.5 * interval_width;
                let c = start + i as f64 * interval_width + r;

                let mut eval_g = 0.0;
                let mut eval_gk = 0.0;

                for (node, w_gk, w_g) in coef {
                    let eval = integrand(r * node + c);

                    eval_gk += eval * w_gk;
                    if let Some(w_g) = w_g {
                        eval_g += eval * w_g
                    }
                }

                Integral {
                    value: eval_gk * r,
                    ref_value: eval_g * r,
                }
            })
            .sum();

        if integral.is_ok_rel(rel_tol) {
            return integral;
        }

        if interval_count < max_interval_count {
            if interval_count < max_interval_count >> 1 {
                interval_count <<= 1
            } else {
                interval_count = max_interval_count
            }
        } else {
            return integral;
        }
    }
}

#[macro_export]
macro_rules! quadgk_1d {
    ($integrand:expr, $range:expr $(, $set:ident = $val:expr)* $(,)?) => {{
        let config = QuadConfig::<0>::new()$(.$set($val))*;
        quadgk_1d($integrand, $range, config)
    }}
}

#[cfg(test)]
mod tests {
    use float_cmp::assert_approx_eq;

    use crate::*;

    fn base<const N: usize>(
        integral: f64,
        integrand: impl Fn(f64) -> f64,
        start: f64,
        end: f64,
        coef: impl QuadGKCoef<{ N }>,
    ) {
        let eval = quad_gk_1d(integrand, start, end, coef);
        println!("{:?}", eval);

        assert_approx_eq!(f64, eval.value, integral, epsilon = 1e-6);
        assert_approx_eq!(f64, eval.ref_value, integral, epsilon = 1e-6);
        assert!(eval.is_ok_rel(1e-10))
    }

    #[test]
    fn identity() {
        base(3.0 * 4.0, |_| 3.0, -2.0, 2.0, GK11)
    }

    #[test]
    fn damping_oscillation() {
        base(0.500031396154354, |x| (-x).exp() * x.sin(), 0.0, 10.0, GK21)
    }

    #[test]
    fn test_macro() {
        let integral = quadgk_1d!(
            |x| x.sin(),
            0.0..100.0,
            rel_tol = 1e-6,
            max_interval_count = 1_000_000,
            coef = GKCoef::<11>::coef(),
        );
        println!("{:?}", integral);
        assert!(integral.is_ok_rel(1e-6))
    }
}
