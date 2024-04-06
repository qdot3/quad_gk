//! Gauss-Kronrod quadrature rule from 1d to 4d.

mod quad_1d;
mod quad_2d;
mod quad_3d;
mod quad_4d;

pub use quad_1d::quad_gk_1d;
pub use quad_2d::quad_gk_2d;
pub use quad_3d::quad_gk_3d;
pub use quad_4d::quad_gk_4d;
