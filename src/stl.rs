use super::{Error, StlParams, StlResult};

/// Seasonal-trend decomposition using Loess (STL).
pub struct Stl;

impl Stl {
    /// Decomposes a time series.
    pub fn fit(series: &[f32], period: usize) -> Result<StlResult, Error> {
        StlParams::new().fit(series, period)
    }

    /// Creates a new set of parameters.
    pub fn params() -> StlParams {
        StlParams::new()
    }
}

#[cfg(test)]
mod tests {
    use crate::{Error, Stl};

    fn assert_in_delta(exp: f32, act: f32) {
        assert!((exp - act).abs() < 0.001);
    }

    fn assert_elements_in_delta(exp: &[f32], act: &[f32]) {
        assert_eq!(exp.len(), act.len());
        for i in 0..exp.len() {
            assert_in_delta(exp[i], act[i]);
        }
    }

    fn generate_series() -> Vec<f32> {
        return vec![
            5.0, 9.0, 2.0, 9.0, 0.0, 6.0, 3.0, 8.0, 5.0, 8.0, 7.0, 8.0, 8.0, 0.0, 2.0, 5.0, 0.0,
            5.0, 6.0, 7.0, 3.0, 6.0, 1.0, 4.0, 4.0, 4.0, 3.0, 7.0, 5.0, 8.0,
        ];
    }

    #[test]
    fn test_works() {
        let result = Stl::fit(&generate_series(), 7).unwrap();
        assert_elements_in_delta(
            &[0.36926576, 0.75655484, -1.3324139, 1.9553658, -0.6044802],
            &result.seasonal()[..5],
        );
        assert_elements_in_delta(
            &[4.804099, 4.9097075, 5.015316, 5.16045, 5.305584],
            &result.trend()[..5],
        );
        assert_elements_in_delta(
            &[-0.17336464, 3.3337379, -1.6829021, 1.8841844, -4.7011037],
            &result.remainder()[..5],
        );
        assert_elements_in_delta(&[1.0, 1.0, 1.0, 1.0, 1.0], &result.weights()[..5]);
    }

    #[test]
    fn test_robust() {
        let result = Stl::params()
            .robust(true)
            .fit(&generate_series(), 7)
            .unwrap();
        assert_elements_in_delta(
            &[0.14922355, 0.47939026, -1.833231, 1.7411387, 0.8200711],
            &result.seasonal()[..5],
        );
        assert_elements_in_delta(
            &[5.397365, 5.4745436, 5.5517216, 5.6499176, 5.748114],
            &result.trend()[..5],
        );
        assert_elements_in_delta(
            &[-0.5465884, 3.0460663, -1.7184906, 1.6089439, -6.5681853],
            &result.remainder()[..5],
        );
        assert_elements_in_delta(
            &[0.99374926, 0.8129377, 0.9385952, 0.9458036, 0.29742217],
            &result.weights()[..5],
        );
    }

    #[test]
    fn test_into_parts() {
        let result = Stl::fit(&generate_series(), 7).unwrap();
        let (seasonal, trend, remainder, weights) = result.into_parts();
        assert_elements_in_delta(
            &[0.36926576, 0.75655484, -1.3324139, 1.9553658, -0.6044802],
            &seasonal[..5],
        );
        assert_elements_in_delta(
            &[4.804099, 4.9097075, 5.015316, 5.16045, 5.305584],
            &trend[..5],
        );
        assert_elements_in_delta(
            &[-0.17336464, 3.3337379, -1.6829021, 1.8841844, -4.7011037],
            &remainder[..5],
        );
        assert_elements_in_delta(&[1.0, 1.0, 1.0, 1.0, 1.0], &weights[..5]);
    }

    #[test]
    fn test_too_few_periods() {
        let result = Stl::params().fit(&generate_series(), 16);
        assert_eq!(
            result.unwrap_err(),
            Error::Series("series has less than two periods".to_string())
        );
    }

    #[test]
    fn test_bad_seasonal_degree() {
        let result = Stl::params().seasonal_degree(2).fit(&generate_series(), 7);
        assert_eq!(
            result.unwrap_err(),
            Error::Parameter("seasonal_degree must be 0 or 1".to_string())
        );
    }

    #[test]
    fn test_seasonal_strength() {
        let result = Stl::fit(&generate_series(), 7).unwrap();
        assert_in_delta(0.284111676315015, result.seasonal_strength());
    }

    #[test]
    fn test_seasonal_strength_max() {
        let series = (0..30).map(|v| (v % 7) as f32).collect::<Vec<f32>>();
        let result = Stl::fit(&series, 7).unwrap();
        assert_in_delta(1.0, result.seasonal_strength());
    }

    #[test]
    fn test_trend_strength() {
        let result = Stl::fit(&generate_series(), 7).unwrap();
        assert_in_delta(0.16384245231864702, result.trend_strength());
    }

    #[test]
    fn test_trend_strength_max() {
        let series = (0..30).map(|v| v as f32).collect::<Vec<f32>>();
        let result = Stl::fit(&series, 7).unwrap();
        assert_in_delta(1.0, result.trend_strength());
    }
}
