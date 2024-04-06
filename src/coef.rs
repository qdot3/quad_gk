//! Integral points and weights in Gauss-Krnrod quadrature rules.

mod gk11;
mod gk21;
mod gk31;
mod gk41;
mod gk51;
mod gk61;
mod gk71;
mod gk81;
mod gk91;

pub use gk11::GK11;
pub use gk21::GK21;
pub use gk31::GK31;
pub use gk41::GK41;
pub use gk51::GK51;
pub use gk61::GK61;
pub use gk71::GK71;
pub use gk81::GK81;
pub use gk91::GK91;

/// Trait to provide nodes and weights of Gauss Kronrod quadrature rule.
pub trait QuadGKCoef<const NX: usize> {
    /// Returns an array of `(Node, Weight_GK, Weight_G)`, where `Node` is normalized to [-1, 1].
    fn quad_gk_coef(&self) -> [(f64, f64, Option<f64>); NX];
}
