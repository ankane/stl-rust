use super::stl_result::strength;

#[derive(Clone, Debug)]
pub struct MstlResult {
    pub(crate) seasonal: Vec<Vec<f32>>,
    pub(crate) trend: Vec<f32>,
    pub(crate) remainder: Vec<f32>,
}

impl MstlResult {
    pub fn seasonal(&self) -> &[Vec<f32>] {
        &self.seasonal[..]
    }

    pub fn trend(&self) -> &[f32] {
        &self.trend
    }

    pub fn remainder(&self) -> &[f32] {
        &self.remainder
    }

    pub fn seasonal_strength(&self) -> Vec<f32> {
        self.seasonal().iter().map(|s| strength(s, self.remainder())).collect()
    }

    pub fn trend_strength(&self) -> f32 {
        strength(self.trend(), self.remainder())
    }
}
