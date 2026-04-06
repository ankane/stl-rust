use crate::math;
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
    fn ln(&self) -> Self;
    fn max(&self, x: Self) -> Self;
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

    fn ln(&self) -> Self {
        math::f32_ln(*self)
    }

    fn max(&self, x: Self) -> Self {
        f32::max(*self, x)
    }

    fn powf(&self, x: Self) -> Self {
        math::f32_powf(*self, x)
    }

    fn sqrt(&self) -> Self {
        math::f32_sqrt(*self)
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

    fn ln(&self) -> Self {
        math::f64_ln(*self)
    }

    fn max(&self, x: Self) -> Self {
        f64::max(*self, x)
    }

    fn powf(&self, x: Self) -> Self {
        math::f64_powf(*self, x)
    }

    fn sqrt(&self) -> Self {
        math::f64_sqrt(*self)
    }
}
