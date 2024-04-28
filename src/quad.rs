//! Gauss-Kronrod quadrature rule from 1d to 4d.
use std::{ops::Range, sync::Arc};

use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::{Integral, GKConfig};

pub fn quad_gk<const N: usize>(
    integrand: Arc<impl Fn(f64) -> f64 + Send + Sync>,
    range: Range<f64>,
    config: GKConfig<N>,
) -> Integral {
    let GKConfig {
        rel_tol,
        max_interval_count,
        coef,
    } = config;
    let Range { start, end } = range;

    let mut interval_count = 1;
    loop {
        let interval_width = (end - start) / interval_count as f64;
        let integral: Integral = (0..interval_count)
            .into_par_iter()
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

        if integral.is_ok(rel_tol) {
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
macro_rules! quad_gk {
    ($integrand:expr, $range:expr $(, $set:ident = $val:expr)* $(,)?) => {{
        let config = GKConfig::<0>::new()$(.$set($val))*;
        quad_gk($integrand, $range, config)
    }}
}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    use crate::*;

    #[test]
    fn test_macro() {
        let integral = quad_gk!(
            Arc::new(|x: f64| x.sin()),
            0.0..100.0,
            rel_tol = 1e-6,
            max_interval_count = 1_000_000,
            coef = GK11,
        );
        println!("{:?}", integral);
        assert!(integral.is_ok(1e-6))
    }
}
