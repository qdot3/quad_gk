use crate::{Integral, QuadGKCoef};

/// 4d numerical integration by Gauss-Kronrod quadrature rule.
pub fn quad_gk_4d<const NX: usize, const NY: usize, const NZ: usize, const NW: usize>(
    integrand: impl Fn(f64, f64, f64, f64) -> f64,
    start: (
        f64,
        impl Fn(f64) -> f64,
        impl Fn(f64, f64) -> f64,
        impl Fn(f64, f64, f64) -> f64,
    ),
    end: (
        f64,
        impl Fn(f64) -> f64,
        impl Fn(f64, f64) -> f64,
        impl Fn(f64, f64, f64) -> f64,
    ),
    coef: (
        impl QuadGKCoef<{ NX }>,
        impl QuadGKCoef<{ NY }>,
        impl QuadGKCoef<{ NZ }>,
        impl QuadGKCoef<{ NW }>,
    ),
) -> Integral {
    let mut eval_gk = 0.0;
    let mut eval_g = 0.0;

    let rx = 0.5 * (end.0 - start.0);
    let cx = 0.5 * (end.0 + start.0);
    for (node_x, w_gk_x, w_g_x) in coef.0.quad_gk_coef() {
        let x = rx * node_x + cx;

        let mut eval_gk_y = 0.0;
        let mut eval_g_y = 0.0;

        let ry = 0.5 * ((end.1)(x) - (start.1)(x));
        let cy = 0.5 * ((end.1)(x) + (start.1)(x));
        for (node_y, w_gk_y, w_g_y) in coef.1.quad_gk_coef() {
            let y = ry * node_y + cy;

            let mut eval_gk_z = 0.0;
            let mut eval_g_z = 0.0;

            let rz = 0.5 * ((end.2)(x, y) - (start.2)(x, y));
            let cz = 0.5 * ((end.2)(x, y) + (start.2)(x, y));
            for (node_z, w_gk_z, w_g_z) in coef.2.quad_gk_coef() {
                let z = rz * node_z + cz;

                let mut eval_gk_w = 0.0;
                let mut eval_g_w = 0.0;

                let rw = 0.5 * ((end.3)(x, y, z) - (start.3)(x, y, z));
                let cw = 0.5 * ((end.3)(x, y, z) + (start.3)(x, y, z));
                for (node_w, w_gk_w, w_g_w) in coef.3.quad_gk_coef() {
                    let eval = integrand(x, y, z, rw * node_w + cw);

                    eval_gk_w += eval * w_gk_w;
                    if w_g_x.is_some() && w_g_y.is_some() && w_g_z.is_some() && w_g_w.is_some() {
                        eval_g_w += eval * w_g_w.unwrap();
                    }
                }
                eval_gk_z += rw * eval_gk_w * w_gk_z;
                if w_g_x.is_some() && w_g_y.is_some() && w_g_z.is_some() {
                    eval_g_z += rw * eval_g_w * w_g_z.unwrap();
                }
            }
            eval_gk_y += rz * eval_gk_z * w_gk_y;
            if w_g_x.is_some() && w_g_y.is_some() {
                eval_g_y += rz * eval_g_z * w_g_y.unwrap();
            }
        }
        eval_gk += ry * eval_gk_y * w_gk_x;
        if w_g_x.is_some() {
            eval_g += ry * eval_g_y * w_g_x.unwrap()
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

    use crate::{GK11, GK71, GK81, GK91};

    use super::*;

    fn base<const NX: usize, const NY: usize, const NZ: usize, const NW: usize>(
        integral: f64,
        integrand: impl Fn(f64, f64, f64, f64) -> f64,
        start: (
            f64,
            impl Fn(f64) -> f64,
            impl Fn(f64, f64) -> f64,
            impl Fn(f64, f64, f64) -> f64,
        ),
        end: (
            f64,
            impl Fn(f64) -> f64,
            impl Fn(f64, f64) -> f64,
            impl Fn(f64, f64, f64) -> f64,
        ),
        coef: (
            impl QuadGKCoef<{ NX }>,
            impl QuadGKCoef<{ NY }>,
            impl QuadGKCoef<{ NZ }>,
            impl QuadGKCoef<{ NW }>,
        ),
    ) {
        let eval = quad_gk_4d(integrand, start, end, coef);
        println!("{:?}", eval);

        assert_approx_eq!(f64, eval.value, integral, epsilon = 1e-6);
        assert_approx_eq!(f64, eval.ref_value, integral, epsilon = 1e-6);
        assert!(eval.is_ok_by(1e-6, 1))
    }

    #[test]
    fn idnetity() {
        base(
            2.0 * 3.0 * 5.0 * 7.0 * 11.0,
            |_, _, _, _| 11.0,
            (0.0, |_| 0.0, |_, _| 0.0, |_, _, _| 0.0),
            (2.0, |_| 3.0, |_, _| 5.0, |_, _, _| 7.0),
            (GK71, GK81, GK91, GK11),
        )
    }
}
