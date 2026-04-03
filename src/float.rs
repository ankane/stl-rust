use core::iter::Sum;
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};

// TODO seal?
#[doc(hidden)]
pub trait Float:
    Add<Output = Self>
    + AddAssign
    + Copy
    + Div<Output = Self>
    + DivAssign
    + Mul<Output = Self>
    + MulAssign
    + PartialOrd
    + Sub<Output = Self>
    + SubAssign
    + Sum
{
    fn from_f64(x: f64) -> Self;
    fn from_usize(x: usize) -> Self;
    fn one() -> Self;
    fn zero() -> Self;

    fn abs(&self) -> Self;
    fn as_f64(&self) -> f64;
    #[cfg(feature = "std")]
    fn ln(&self) -> Self;
    fn max(&self, x: Self) -> Self;
    #[cfg(feature = "std")]
    fn powf(&self, x: Self) -> Self;
    fn sqrt(&self) -> Self;
}

impl Float for f32 {
    fn from_f64(x: f64) -> Self {
        x as Self
    }

    fn from_usize(x: usize) -> Self {
        x as Self
    }

    fn one() -> Self {
        1.0
    }

    fn zero() -> Self {
        0.0
    }

    fn abs(&self) -> Self {
        f32::abs(*self)
    }

    fn as_f64(&self) -> f64 {
        *self as f64
    }

    #[cfg(feature = "std")]
    fn ln(&self) -> Self {
        f32::ln(*self)
    }

    fn max(&self, x: Self) -> Self {
        f32::max(*self, x)
    }

    #[cfg(feature = "std")]
    fn powf(&self, x: Self) -> Self {
        f32::powf(*self, x)
    }

    #[cfg(feature = "std")]
    fn sqrt(&self) -> Self {
        f32::sqrt(*self)
    }

    #[cfg(not(feature = "std"))]
    fn sqrt(&self) -> Self {
        core::f32::math::sqrt(*self)
    }
}

impl Float for f64 {
    fn from_f64(x: f64) -> Self {
        x
    }

    fn from_usize(x: usize) -> Self {
        x as Self
    }

    fn one() -> Self {
        1.0
    }

    fn zero() -> Self {
        0.0
    }

    fn abs(&self) -> Self {
        f64::abs(*self)
    }

    fn as_f64(&self) -> f64 {
        *self
    }

    #[cfg(feature = "std")]
    fn ln(&self) -> Self {
        f64::ln(*self)
    }

    fn max(&self, x: Self) -> Self {
        f64::max(*self, x)
    }

    #[cfg(feature = "std")]
    fn powf(&self, x: Self) -> Self {
        f64::powf(*self, x)
    }

    #[cfg(feature = "std")]
    fn sqrt(&self) -> Self {
        f64::sqrt(*self)
    }

    #[cfg(not(feature = "std"))]
    fn sqrt(&self) -> Self {
        core::f64::math::sqrt(*self)
    }
}
