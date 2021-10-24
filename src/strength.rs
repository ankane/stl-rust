// https://otexts.com/fpp2/seasonal-strength.html

use crate::StlResult;

pub fn seasonal_strength(result: &StlResult) -> f32 {
    let sr = result.seasonal().iter().zip(result.remainder()).map(|(a, b)| a + b).collect::<Vec<f32>>();
    (1.0 - var(result.remainder()) / var(&sr)).max(0.0)
}

pub fn trend_strength(result: &StlResult) -> f32 {
    let tr = result.trend().iter().zip(result.remainder()).map(|(a, b)| a + b).collect::<Vec<f32>>();
    (1.0 - var(result.remainder()) / var(&tr)).max(0.0)
}

fn var(series: &[f32]) -> f32 {
    let mean = series.iter().sum::<f32>() / series.len() as f32;
    series.iter().map(|v| (v - mean).powf(2.0)).sum::<f32>() / (series.len() as f32 - 1.0)
}
