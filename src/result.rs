#[derive(Clone, Debug)]
pub struct StlResult {
    pub(crate) seasonal: Vec<f32>,
    pub(crate) trend: Vec<f32>,
    pub(crate) remainder: Vec<f32>,
    pub(crate) weights: Vec<f32>
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
        let sr = self.seasonal().iter().zip(self.remainder()).map(|(a, b)| a + b).collect::<Vec<f32>>();
        (1.0 - self.var(self.remainder()) / self.var(&sr)).max(0.0)
    }

    pub fn trend_strength(&self) -> f32 {
        let tr = self.trend().iter().zip(self.remainder()).map(|(a, b)| a + b).collect::<Vec<f32>>();
        (1.0 - self.var(self.remainder()) / self.var(&tr)).max(0.0)
    }

    fn var(&self, series: &[f32]) -> f32 {
        let mean = series.iter().sum::<f32>() / series.len() as f32;
        series.iter().map(|v| (v - mean).powf(2.0)).sum::<f32>() / (series.len() as f32 - 1.0)
    }
}
