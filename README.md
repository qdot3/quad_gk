# quad_gk

Numerical Integration library based on Gauss Kronrod quadrature rule.

# Example
```
use quad_gk::{quad_gk_2d, Integral, GK21};

let integral: Integral = quad_gk_2d(
    |x, y| (-x.powi(2) - y.powi(4) / 5.0),
    (0.0, |_| 0.0),
    (5.0, |x| x.powf(1.5)),
    (GK21, GK21),
);

// estimated approximation error is smaller than 1e-6 or 1ULP.
assert!(integral.is_ok_by(1e-6, 1))
```