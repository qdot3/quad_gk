use std::ops::{Add, AddAssign};

use float_cmp::{approx_eq, F64Margin};

/// Integral with reference.
///
/// You can estimate approximation error by `|value - ref_value|`.
#[derive(Debug, Clone)]
pub struct Integral {
    /// Integral by Gauss-Kronrod quadrature rule.
    pub value: f64,
    /// Integral by Gauss quadrature rule.
    pub ref_value: f64,
}

impl Integral {
    /// Retuns true, if approximation error is estimated to be small enough.
    pub fn is_ok(&self) -> bool {
        approx_eq!(f64, self.value, self.ref_value)
    }

    /// Returns true, if estimated approximation error meets given condition.
    pub fn is_ok_by(&self, epsilon: f64, ulps: i64) -> bool {
        approx_eq!(f64, self.value, self.ref_value, F64Margin { epsilon, ulps })
    }
}

impl Add for Integral {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            value: self.value + rhs.value,
            ref_value: self.ref_value + rhs.ref_value,
        }
    }
}

impl AddAssign for Integral {
    fn add_assign(&mut self, rhs: Self) {
        *self = Self {
            value: self.value + rhs.value,
            ref_value: self.ref_value + rhs.ref_value,
        }
    }
}
