use super::{Error, MstlParams, MstlResult};

/// Multiple seasonal-trend decomposition using Loess (MSTL).
pub struct Mstl;

impl Mstl {
    /// Decomposes a time series.
    pub fn fit(series: &[f32], periods: &[usize]) -> Result<MstlResult, Error> {
        MstlParams::new().fit(series, periods)
    }

    /// Creates a new set of parameters.
    pub fn params() -> MstlParams {
        MstlParams::new()
    }
}

#[cfg(test)]
mod tests {
    use crate::{Error, Mstl, Stl};

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
        let result = Mstl::fit(&generate_series(), &[6, 10]).unwrap();
        assert_elements_in_delta(
            &[0.28318232, 0.70529824, -1.980384, 2.1643379, -2.3356874],
            &result.seasonal()[0][..5],
        );
        assert_elements_in_delta(
            &[1.4130436, 1.6048906, 0.050958008, -1.8706754, -1.7704514],
            &result.seasonal()[1][..5],
        );
        assert_elements_in_delta(
            &[5.139485, 5.223691, 5.3078976, 5.387292, 5.4666862],
            &result.trend()[..5],
        );
        assert_elements_in_delta(
            &[-1.835711, 1.4661198, -1.3784716, 3.319045, -1.3605475],
            &result.remainder()[..5],
        );
    }

    #[test]
    fn test_into_parts() {
        let result = Mstl::fit(&generate_series(), &[6, 10]).unwrap();
        let (seasonal, trend, remainder) = result.into_parts();
        assert_elements_in_delta(
            &[0.28318232, 0.70529824, -1.980384, 2.1643379, -2.3356874],
            &seasonal[0][..5],
        );
        assert_elements_in_delta(
            &[1.4130436, 1.6048906, 0.050958008, -1.8706754, -1.7704514],
            &seasonal[1][..5],
        );
        assert_elements_in_delta(
            &[5.139485, 5.223691, 5.3078976, 5.387292, 5.4666862],
            &trend[..5],
        );
        assert_elements_in_delta(
            &[-1.835711, 1.4661198, -1.3784716, 3.319045, -1.3605475],
            &remainder[..5],
        );
    }

    #[test]
    fn test_unsorted_periods() {
        let result = Mstl::fit(&generate_series(), &[10, 6]).unwrap();
        assert_elements_in_delta(
            &[1.4130436, 1.6048906, 0.050958008, -1.8706754, -1.7704514],
            &result.seasonal()[0][..5],
        );
        assert_elements_in_delta(
            &[0.28318232, 0.70529824, -1.980384, 2.1643379, -2.3356874],
            &result.seasonal()[1][..5],
        );
        assert_elements_in_delta(
            &[5.139485, 5.223691, 5.3078976, 5.387292, 5.4666862],
            &result.trend()[..5],
        );
        assert_elements_in_delta(
            &[-1.835711, 1.4661198, -1.3784716, 3.319045, -1.3605475],
            &result.remainder()[..5],
        );
    }

    #[test]
    fn test_lambda() {
        let result = Mstl::params()
            .lambda(0.5)
            .fit(&generate_series(), &[6, 10])
            .unwrap();
        assert_elements_in_delta(
            &[0.43371448, 0.10503793, -0.7178911, 1.2356076, -1.8253292],
            &result.seasonal()[0][..5],
        );
        assert_elements_in_delta(
            &[1.0437742, 0.8650516, 0.07303603, -1.428663, -1.1990008],
            &result.seasonal()[1][..5],
        );
        assert_elements_in_delta(
            &[2.0748303, 2.1291165, 2.1834028, 2.2330272, 2.2826517],
            &result.trend()[..5],
        );
        assert_elements_in_delta(
            &[-1.0801829, 0.900794, -0.7101207, 1.9600279, -1.2583216],
            &result.remainder()[..5],
        );
    }

    #[test]
    fn test_lambda_zero() {
        let series: Vec<f32> = generate_series().iter().map(|&v| v + 1.0).collect();
        let result = Mstl::params().lambda(0.0).fit(&series, &[6, 10]).unwrap();
        assert_elements_in_delta(
            &[0.18727916, 0.029921893, -0.2716494, 0.47748315, -0.7320051],
            &result.seasonal()[0][..5],
        );
        assert_elements_in_delta(
            &[
                0.42725056,
                0.32145387,
                -0.019030934,
                -0.56607914,
                -0.46765903,
            ],
            &result.seasonal()[1][..5],
        );
        assert_elements_in_delta(
            &[1.592807, 1.6144379, 1.6360688, 1.6559447, 1.6758206],
            &result.trend()[..5],
        );
        assert_elements_in_delta(
            &[-0.41557717, 0.33677137, -0.24677622, 0.7352363, -0.47615635],
            &result.remainder()[..5],
        );
    }

    #[test]
    fn test_lambda_out_of_range() {
        let result = Mstl::params().lambda(2.0).fit(&generate_series(), &[6, 10]);
        assert_eq!(
            result.unwrap_err(),
            Error::Parameter("lambda must be between 0 and 1".to_string())
        );
    }

    #[test]
    fn test_empty_periods() {
        let periods: Vec<usize> = Vec::new();
        let result = Mstl::fit(&generate_series(), &periods);
        assert_eq!(
            result.unwrap_err(),
            Error::Parameter("periods must not be empty".to_string())
        );
    }

    #[test]
    fn test_period_one() {
        let result = Mstl::fit(&generate_series(), &[1]);
        assert_eq!(
            result.unwrap_err(),
            Error::Parameter("periods must be at least 2".to_string())
        );
    }

    #[test]
    fn test_too_few_periods() {
        let result = Mstl::fit(&generate_series(), &[16]);
        assert_eq!(
            result.unwrap_err(),
            Error::Series("series has less than two periods".to_string())
        );
    }

    #[test]
    fn test_seasonal_strength() {
        let mut stl_params = Stl::params();
        stl_params.seasonal_length(7);
        let result = Mstl::params()
            .stl_params(stl_params)
            .fit(&generate_series(), &[7])
            .unwrap();
        assert_in_delta(0.284111676315015, result.seasonal_strength()[0]);
    }

    #[test]
    fn test_seasonal_strength_max() {
        let series = (0..30).map(|v| (v % 7) as f32).collect::<Vec<f32>>();
        let mut stl_params = Stl::params();
        stl_params.seasonal_length(7);
        let result = Mstl::params()
            .stl_params(stl_params)
            .fit(&series, &[7])
            .unwrap();
        assert_in_delta(1.0, result.seasonal_strength()[0]);
    }

    #[test]
    fn test_trend_strength() {
        let mut stl_params = Stl::params();
        stl_params.seasonal_length(7);
        let result = Mstl::params()
            .stl_params(stl_params)
            .fit(&generate_series(), &[7])
            .unwrap();
        assert_in_delta(0.16384245231864702, result.trend_strength());
    }

    #[test]
    fn test_trend_strength_max() {
        let series = (0..30).map(|v| v as f32).collect::<Vec<f32>>();
        let mut stl_params = Stl::params();
        stl_params.seasonal_length(7);
        let result = Mstl::params()
            .stl_params(stl_params)
            .fit(&series, &[7])
            .unwrap();
        assert_in_delta(1.0, result.trend_strength());
    }
}
