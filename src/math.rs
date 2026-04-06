#[cfg(feature = "std")]
mod std_math {
    pub fn ceil(x: f32) -> f32 {
        x.ceil()
    }

    pub fn f32_ln(x: f32) -> f32 {
        x.ln()
    }

    pub fn f32_powf(x: f32, n: f32) -> f32 {
        x.powf(n)
    }

    pub fn f32_sqrt(x: f32) -> f32 {
        x.sqrt()
    }

    pub fn f64_ln(x: f64) -> f64 {
        x.ln()
    }

    pub fn f64_powf(x: f64, n: f64) -> f64 {
        x.powf(n)
    }

    pub fn f64_sqrt(x: f64) -> f64 {
        x.sqrt()
    }
}

#[cfg(not(feature = "std"))]
mod core_math {
    pub use core::f32::math::{ceil, sqrt as f32_sqrt};
    pub use core::f64::math::sqrt as f64_sqrt;

    pub fn f32_ln(_x: f32) -> f32 {
        todo!()
    }

    pub fn f32_powf(_x: f32, _n: f32) -> f32 {
        todo!()
    }

    pub fn f64_ln(_x: f64) -> f64 {
        todo!()
    }

    pub fn f64_powf(_x: f64, _n: f64) -> f64 {
        todo!()
    }
}

#[cfg(feature = "std")]
pub use std_math::*;

#[cfg(not(feature = "std"))]
pub use core_math::*;
