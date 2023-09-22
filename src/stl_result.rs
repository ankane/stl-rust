#[derive(Clone, Debug)]
pub struct StlResult {
    pub(crate) seasonal: Vec<f32>,
    pub(crate) trend: Vec<f32>,
    pub(crate) remainder: Vec<f32>,
    pub(crate) weights: Vec<f32>
}

fn var(series: &[f32]) -> f32 {
    let mean = series.iter().sum::<f32>() / series.len() as f32;
    series.iter().map(|v| (v - mean).powf(2.0)).sum::<f32>() / (series.len() as f32 - 1.0)
}

pub(crate) fn strength(component: &[f32], remainder: &[f32]) -> f32 {
    let sr = component.iter().zip(remainder).map(|(a, b)| a + b).collect::<Vec<f32>>();
    (1.0 - var(remainder) / var(&sr)).max(0.0)
}

impl StlResult {
    pub fn seasonal(&self) -> &[f32] {
        &self.seasonal
    }

    pub fn trend(&self) -> &[f32] {
        &self.trend
    }

    pub fn remainder(&self) -> &[f32] {
        &self.remainder
    }

    pub fn weights(&self) -> &[f32] {
        &self.weights
    }

    pub fn seasonal_strength(&self) -> f32 {
        strength(self.seasonal(), self.remainder())
    }

    pub fn trend_strength(&self) -> f32 {
        strength(self.trend(), self.remainder())
    }

    pub fn into_parts(self) -> (Vec<f32>, Vec<f32>, Vec<f32>, Vec<f32>) {
        (self.seasonal, self.trend, self.remainder, self.weights)
    }
}
