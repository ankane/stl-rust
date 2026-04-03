use super::{Error, Float, StlParams};

#[cfg(feature = "std")]
use super::StlResult;

/// Seasonal-trend decomposition using Loess (STL).
pub struct Stl;

impl Stl {
    /// Decomposes a time series.
    #[cfg(feature = "std")]
    pub fn fit<T: Float>(series: &[T], period: usize) -> Result<StlResult<T>, Error> {
        StlParams::new().fit(series, period)
    }

    /// Decomposes a time series.
    #[cfg(not(feature = "std"))]
    pub fn fit<T: Float>(
        series: &[T],
        period: usize,
        seasonal: &mut [T],
        trend: &mut [T],
        weights: &mut [T],
        work: &mut [T],
    ) -> Result<(), Error> {
        StlParams::new().fit(series, period, seasonal, trend, weights, work)
    }

    /// Creates a new set of parameters.
    pub fn params() -> StlParams {
        StlParams::new()
    }
}

#[cfg(test)]
#[cfg(feature = "std")]
mod tests {
    use crate::{Error, Float, Stl};

    fn assert_in_delta<T: Float>(exp: T, act: T) {
        assert!((exp - act).abs() < T::from_f64(0.001));
    }

    fn assert_elements_in_delta<T: Float>(exp: &[T], act: &[T]) {
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
    fn test_f32() {
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
    fn test_f64() {
        let series: Vec<f64> = vec![
            5.0, 9.0, 2.0, 9.0, 0.0, 6.0, 3.0, 8.0, 5.0, 8.0, 7.0, 8.0, 8.0, 0.0, 2.0, 5.0, 0.0,
            5.0, 6.0, 7.0, 3.0, 6.0, 1.0, 4.0, 4.0, 4.0, 3.0, 7.0, 5.0, 8.0,
        ];
        let result = Stl::fit(&series, 7).unwrap();
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
    fn test_period_one() {
        let result = Stl::fit(&generate_series(), 1);
        let err = result.unwrap_err();
        assert_eq!(err, Error::Period);
        assert_eq!(err.to_string(), "period must be at least 2");
    }

    #[test]
    fn test_too_few_periods() {
        let result = Stl::params().fit(&generate_series(), 16);
        let err = result.unwrap_err();
        assert_eq!(err, Error::Series);
        assert_eq!(err.to_string(), "series has less than two periods");
    }

    #[test]
    fn test_bad_seasonal_degree() {
        let result = Stl::params().seasonal_degree(2).fit(&generate_series(), 7);
        let err = result.unwrap_err();
        assert_eq!(err, Error::SeasonalDegree);
        assert_eq!(err.to_string(), "seasonal_degree must be 0 or 1");
    }

    #[test]
    fn test_bad_trend_degree() {
        let result = Stl::params().trend_degree(2).fit(&generate_series(), 7);
        let err = result.unwrap_err();
        assert_eq!(err, Error::TrendDegree);
        assert_eq!(err.to_string(), "trend_degree must be 0 or 1");
    }

    #[test]
    fn test_bad_low_pass_degree() {
        let result = Stl::params().low_pass_degree(2).fit(&generate_series(), 7);
        let err = result.unwrap_err();
        assert_eq!(err, Error::LowPassDegree);
        assert_eq!(err.to_string(), "low_pass_degree must be 0 or 1");
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

#[cfg(test)]
#[cfg(not(feature = "std"))]
mod tests {
    use crate::{Float, Stl};

    fn assert_in_delta<T: Float>(exp: T, act: T) {
        assert!((exp - act).abs() < T::from_f64(0.001));
    }

    fn assert_elements_in_delta<T: Float>(exp: &[T], act: &[T]) {
        assert_eq!(exp.len(), act.len());
        for i in 0..exp.len() {
            assert_in_delta(exp[i], act[i]);
        }
    }

    #[test]
    fn test_works() {
        let series = [
            5.0, 9.0, 2.0, 9.0, 0.0, 6.0, 3.0, 8.0, 5.0, 8.0, 7.0, 8.0, 8.0, 0.0, 2.0, 5.0, 0.0,
            5.0, 6.0, 7.0, 3.0, 6.0, 1.0, 4.0, 4.0, 4.0, 3.0, 7.0, 5.0, 8.0,
        ];
        // make non-zero to test correctness
        let mut seasonal = [1.0; 30];
        let mut trend = [2.0; 30];
        let mut weights = [3.0; 30];
        let mut work = [4.0; (30 + 2 * 7) * 5];

        Stl::fit(
            &series,
            7,
            &mut seasonal,
            &mut trend,
            &mut weights,
            &mut work,
        )
        .unwrap();

        assert_elements_in_delta(
            &[0.36926576, 0.75655484, -1.3324139, 1.9553658, -0.6044802],
            &seasonal[..5],
        );
        assert_elements_in_delta(
            &[4.804099, 4.9097075, 5.015316, 5.16045, 5.305584],
            &trend[..5],
        );
        assert_elements_in_delta(&[1.0, 1.0, 1.0, 1.0, 1.0], &weights[..5]);
    }

    #[test]
    fn test_robust() {
        let series = [
            5.0, 9.0, 2.0, 9.0, 0.0, 6.0, 3.0, 8.0, 5.0, 8.0, 7.0, 8.0, 8.0, 0.0, 2.0, 5.0, 0.0,
            5.0, 6.0, 7.0, 3.0, 6.0, 1.0, 4.0, 4.0, 4.0, 3.0, 7.0, 5.0, 8.0,
        ];
        // make non-zero to test correctness
        let mut seasonal = [1.0; 30];
        let mut trend = [2.0; 30];
        let mut weights = [3.0; 30];
        let mut work = [4.0; (30 + 2 * 7) * 5];

        Stl::params()
            .robust(true)
            .fit(
                &series,
                7,
                &mut seasonal,
                &mut trend,
                &mut weights,
                &mut work,
            )
            .unwrap();

        assert_elements_in_delta(
            &[0.14922355, 0.47939026, -1.833231, 1.7411387, 0.8200711],
            &seasonal[..5],
        );
        assert_elements_in_delta(
            &[5.397365, 5.4745436, 5.5517216, 5.6499176, 5.748114],
            &trend[..5],
        );
        assert_elements_in_delta(
            &[0.99374926, 0.8129377, 0.9385952, 0.9458036, 0.29742217],
            &weights[..5],
        );
    }
}
