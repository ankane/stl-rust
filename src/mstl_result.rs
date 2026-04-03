use super::stl_result::strength;
use super::Float;
use alloc::vec::Vec;

/// A MSTL result.
#[derive(Clone, Debug)]
pub struct MstlResult<T: Float = f32> {
    pub(crate) seasonal: Vec<Vec<T>>,
    pub(crate) trend: Vec<T>,
    pub(crate) remainder: Vec<T>,
}

impl<T: Float> MstlResult<T> {
    /// Returns the seasonal components.
    pub fn seasonal(&self) -> &[Vec<T>] {
        &self.seasonal[..]
    }

    /// Returns the trend component.
    pub fn trend(&self) -> &[T] {
        &self.trend
    }

    /// Returns the remainder.
    pub fn remainder(&self) -> &[T] {
        &self.remainder
    }

    /// Returns the seasonal strength.
    pub fn seasonal_strength(&self) -> Vec<f64> {
        self.seasonal()
            .iter()
            .map(|s| strength(s, self.remainder()))
            .collect()
    }

    /// Returns the trend strength.
    pub fn trend_strength(&self) -> f64 {
        strength(self.trend(), self.remainder())
    }

    /// Consumes the result, returning the seasonal components, trend component, and remainder.
    pub fn into_parts(self) -> (Vec<Vec<T>>, Vec<T>, Vec<T>) {
        (self.seasonal, self.trend, self.remainder)
    }
}
