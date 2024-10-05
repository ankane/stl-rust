use super::stl_result::strength;

/// A MSTL result.
#[derive(Clone, Debug)]
pub struct MstlResult {
    pub(crate) seasonal: Vec<Vec<f32>>,
    pub(crate) trend: Vec<f32>,
    pub(crate) remainder: Vec<f32>,
}

impl MstlResult {
    /// Returns the seasonal components.
    pub fn seasonal(&self) -> &[Vec<f32>] {
        &self.seasonal[..]
    }

    /// Returns the trend component.
    pub fn trend(&self) -> &[f32] {
        &self.trend
    }

    /// Returns the remainder.
    pub fn remainder(&self) -> &[f32] {
        &self.remainder
    }

    /// Returns the seasonal strength.
    pub fn seasonal_strength(&self) -> Vec<f32> {
        self.seasonal()
            .iter()
            .map(|s| strength(s, self.remainder()))
            .collect()
    }

    /// Returns the trend strength.
    pub fn trend_strength(&self) -> f32 {
        strength(self.trend(), self.remainder())
    }

    /// Consumes the result, returning the seasonal components, trend component, and remainder.
    pub fn into_parts(self) -> (Vec<Vec<f32>>, Vec<f32>, Vec<f32>) {
        (self.seasonal, self.trend, self.remainder)
    }
}
