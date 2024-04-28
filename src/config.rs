//! Configuration in numerical integration.

use crate::GK11;

pub struct GKConfig<const N: usize> {
    // required relative tolerance
    pub rel_tol: f64,
    // maximum number of interval divisions
    pub max_interval_count: usize,
    // integral points and coefficients
    pub coef: [(f64, f64, Option<f64>); N],
}

impl<const N: usize> GKConfig<N> {
    // Default configuration.
    pub fn new() -> GKConfig<11> {
        GKConfig {
            rel_tol: 1e-6,
            max_interval_count: 100,
            coef: GK11,
        }
    }

    pub fn rel_tol(self, rel_tol: f64) -> Self {
        Self { rel_tol, ..self }
    }

    pub fn max_interval_count(self, max_interval_count: usize) -> Self {
        Self {
            max_interval_count,
            ..self
        }
    }

    pub fn coef<const M: usize>(self, coef: [(f64, f64, Option<f64>); M]) -> GKConfig<M> {
        GKConfig {
            rel_tol: self.rel_tol,
            max_interval_count: self.max_interval_count,
            coef,
        }
    }
}
