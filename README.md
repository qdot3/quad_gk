# quad_gk

Fast and precise numerical integration library based on Gauss Kronrod quadrature rule.

## Basic usage
```rust
use std::sync::Arc;

use quad_gk::*;

let integral: Integral = quad_gk!(
    Arc::new(|x: f64| x.sin() * (-x).exp()),
    0.0..100.0,
);

// estimated relative error in numerical integration is smaller than 1e-6 by default.
assert!(integral.is_ok(1e-6))
```

## More precise calculation
```rust
use std::sync::Arc;

use quad_gk::*;

let integral: Integral = quad_gk!(
    Arc::new(|x: f64| x.sin() * (-x).exp()),
    0.0..100.0,
    rel_tol=1e-14,
    max_interval_count=1_000_000,
    coef=GK91,
);

// estimated relative error in numerical integration is smaller than 1e-14!
assert!(integral.is_ok(1e-14))
```

