use super::{Error, MstlParams, MstlResult};

pub struct Mstl;

impl Mstl {
    pub fn fit(series: &[f32], periods: &[usize]) -> Result<MstlResult, Error> {
        MstlParams::new().fit(series, periods)
    }

    pub fn params() -> MstlParams {
        MstlParams::new()
    }
}

#[cfg(test)]
mod tests {
    use crate::{Error, Mstl};

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
            5.0, 9.0, 2.0, 9.0, 0.0, 6.0, 3.0, 8.0, 5.0, 8.0,
            7.0, 8.0, 8.0, 0.0, 2.0, 5.0, 0.0, 5.0, 6.0, 7.0,
            3.0, 6.0, 1.0, 4.0, 4.0, 4.0, 3.0, 7.0, 5.0, 8.0
        ];
    }

    #[test]
    fn test_works() {
        let result = Mstl::fit(&generate_series(), &[6, 10]).unwrap();
        assert_elements_in_delta(&[0.28318232, 0.70529824, -1.980384, 2.1643379, -2.3356874], &result.seasonal()[0][..5]);
        assert_elements_in_delta(&[1.4130436, 1.6048906, 0.050958008, -1.8706754, -1.7704514], &result.seasonal()[1][..5]);
        assert_elements_in_delta(&[5.139485, 5.223691, 5.3078976, 5.387292, 5.4666862], &result.trend()[..5]);
        assert_elements_in_delta(&[-1.835711, 1.4661198, -1.3784716, 3.319045, -1.3605475], &result.remainder()[..5]);
    }

    #[test]
    fn test_unsorted_periods() {
        let result = Mstl::fit(&generate_series(), &[10, 6]).unwrap();
        assert_elements_in_delta(&[1.4130436, 1.6048906, 0.050958008, -1.8706754, -1.7704514], &result.seasonal()[0][..5]);
        assert_elements_in_delta(&[0.28318232, 0.70529824, -1.980384, 2.1643379, -2.3356874], &result.seasonal()[1][..5]);
        assert_elements_in_delta(&[5.139485, 5.223691, 5.3078976, 5.387292, 5.4666862], &result.trend()[..5]);
        assert_elements_in_delta(&[-1.835711, 1.4661198, -1.3784716, 3.319045, -1.3605475], &result.remainder()[..5]);
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
}
