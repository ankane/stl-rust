use super::Float;

/// A STL result.
#[derive(Clone, Debug)]
pub struct StlResult<T: Float = f32> {
    pub(crate) seasonal: Vec<T>,
    pub(crate) trend: Vec<T>,
    pub(crate) remainder: Vec<T>,
    pub(crate) weights: Vec<T>,
}

fn var<T: Float>(series: &[T]) -> f64 {
    let mean = series.iter().map(|v| (*v).as_f64()).sum::<f64>() / series.len() as f64;
    series
        .iter()
        .map(|v| ((*v).as_f64() - mean).powf(2.0))
        .sum::<f64>()
        / (series.len() as f64 - 1.0)
}

pub(crate) fn strength<T: Float>(component: &[T], remainder: &[T]) -> f64 {
    let sr = component
        .iter()
        .zip(remainder)
        .map(|(a, b)| *a + *b)
        .collect::<Vec<T>>();
    (1.0 - var(remainder) / var(&sr)).max(0.0)
}

impl<T: Float> StlResult<T> {
    /// Returns the seasonal component.
    pub fn seasonal(&self) -> &[T] {
        &self.seasonal
    }

    /// Returns the trend component.
    pub fn trend(&self) -> &[T] {
        &self.trend
    }

    /// Returns the remainder.
    pub fn remainder(&self) -> &[T] {
        &self.remainder
    }

    /// Returns the weights.
    pub fn weights(&self) -> &[T] {
        &self.weights
    }

    /// Returns the seasonal strength.
    pub fn seasonal_strength(&self) -> f64 {
        strength(self.seasonal(), self.remainder())
    }

    /// Returns the trend strength.
    pub fn trend_strength(&self) -> f64 {
        strength(self.trend(), self.remainder())
    }

    /// Consumes the result, returning the seasonal component, trend component, remainder, and weights.
    pub fn into_parts(self) -> (Vec<T>, Vec<T>, Vec<T>, Vec<T>) {
        (self.seasonal, self.trend, self.remainder, self.weights)
    }
}
