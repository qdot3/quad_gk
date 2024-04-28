use crate::{quadgk_1d, Integral, QuadConfig, QuadGKCoef};

/// 2d numerical integration by Gauss-Kronrod quadrature rule.
pub fn quad_gk_2d<const NX: usize, const NY: usize>(
    integrand: impl Fn(f64, f64) -> f64,
    start: (f64, impl Fn(f64) -> f64),
    end: (f64, impl Fn(f64) -> f64),
    coef: (impl QuadGKCoef<{ NX }>, impl QuadGKCoef<{ NY }>),
) -> Integral {
    let mut eval_gk = 0.0;
    let mut eval_g = 0.0;

    let rx = 0.5 * (end.0 - start.0);
    let cx = 0.5 * (end.0 + start.0);

    for (node_x, w_gk_x, w_g_x) in coef.0.quad_gk_coef() {
        let x = rx * node_x + cx;

        let ry = 0.5 * ((end.1)(x) - (start.1)(x));
        let cy = 0.5 * ((end.1)(x) + (start.1)(x));

        let mut eval_gk_y = 0.0;
        let mut eval_g_y = 0.0;
        for (node_y, w_gk_y, w_g_y) in coef.1.quad_gk_coef() {
            let eval = integrand(x, ry * node_y + cy);

            eval_gk_y += eval * w_gk_y;
            if w_g_x.is_some() && w_g_y.is_some() {
                eval_g_y += eval * w_g_y.unwrap()
            }
        }
        eval_gk += ry * eval_gk_y * w_gk_x;
        if w_g_x.is_some() {
            eval_g += ry * eval_g_y * w_g_x.unwrap();
        }
    }

    Integral {
        value: eval_gk * rx,
        ref_value: eval_g * rx,
    }
}

#[cfg(test)]
mod tests {

    use float_cmp::assert_approx_eq;

    use crate::{GK11, GK21, GK31};

    use super::*;

    fn base<const NX: usize, const NY: usize>(
        integral: f64,
        integrand: impl Fn(f64, f64) -> f64,
        start: (f64, impl Fn(f64) -> f64),
        end: (f64, impl Fn(f64) -> f64),
        coef: (impl QuadGKCoef<{ NX }>, impl QuadGKCoef<{ NY }>),
    ) {
        let eval = quad_gk_2d(integrand, start, end, coef);
        println!("{:?}", eval);

        assert_approx_eq!(f64, eval.value, integral, epsilon = 1e-6);
        assert_approx_eq!(f64, eval.ref_value, integral, epsilon = 1e-6);
        assert!(eval.is_ok_by(1e-6, 1))
    }

    #[test]
    fn idnetity() {
        base(
            2.0 * 3.0 * 5.0,
            |_, _| 5.0,
            (0.0, |_| 0.0),
            (2.0, |_| 3.0),
            (GK21, GK31),
        )
    }

    #[test]
    fn log() {
        base(
            5.58832107266053356560471,
            |x, y| (x.powi(2) + y).ln(),
            (1.0, |x| 2.0 * x),
            (2.0, |_| 6.0),
            (GK11, GK11),
        )
    }
}
