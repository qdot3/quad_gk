use crate::{Integral, QuadGKCoef};

/// 1d numerical integration by Gauss-Kronrod quadrature rule.
pub fn quad_gk_1d<const N: usize>(
    integrand: impl Fn(f64) -> f64,
    start: f64,
    end: f64,
    coef: impl QuadGKCoef<{ N }>,
) -> Integral {
    let r = 0.5 * (end - start);
    let c = 0.5 * (end + start);

    let mut eval_g = 0.0;
    let mut eval_gk = 0.0;

    for (node, w_gk, w_g) in coef.quad_gk_coef() {
        let eval = integrand(r * node + c);

        eval_gk += eval * w_gk;
        if w_g.is_some() {
            eval_g += eval * w_g.unwrap()
        }
    }

    Integral {
        value: eval_gk * r,
        ref_value: eval_g * r,
    }
}

#[cfg(test)]
mod tests {
    use float_cmp::assert_approx_eq;

    use crate::{GK11, GK21};

    use super::*;

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
        assert!(eval.is_ok_by(1e-6, 1))
    }

    #[test]
    fn identity() {
        base(3.0 * 4.0, |_| 3.0, -2.0, 2.0, GK11)
    }

    #[test]
    fn damping_oscillation() {
        base(0.500031396154354, |x| (-x).exp() * x.sin(), 0.0, 10.0, GK21)
    }
}
