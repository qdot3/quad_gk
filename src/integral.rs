use std::{
    iter::Sum,
    ops::{Add, AddAssign},
};

/// Integral and reference value with low precision.
#[derive(Debug, Clone)]
pub struct Integral {
    /// Integral by (2N+1)-th order Gauss-Kronrod quadrature rule.
    pub value: f64,
    /// Integral by N-th order Gauss quadrature rule.
    pub ref_value: f64,
}

impl Integral {
    // Retuns true, if calculation meets given condition.
    pub fn is_ok(&self, rel_tol: f64) -> bool {
        (self.value - self.ref_value).abs() <= rel_tol * self.value.abs()
    }

    // Returns estimated relative error in calculation.
    pub fn rel_err(&self) -> f64 {
        (self.value - self.ref_value).abs() / self.value.abs()
    }
}

impl Sum for Integral {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(
            Integral {
                value: 0.0,
                ref_value: 0.0,
            },
            |a, b| a + b,
        )
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
