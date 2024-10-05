use super::mstl_impl::mstl;
use super::{Error, MstlResult, StlParams};

/// A set of MSTL parameters.
#[derive(Clone, Debug)]
pub struct MstlParams {
    iterate: usize,
    lambda: Option<f32>,
    swin: Option<Vec<usize>>,
    stl_params: StlParams,
}

impl MstlParams {
    /// Creates a new set of parameters.
    pub fn new() -> Self {
        Self {
            iterate: 2,
            lambda: None,
            swin: None,
            stl_params: StlParams::new(),
        }
    }

    /// Sets the number of iterations.
    pub fn iterations(&mut self, iterations: usize) -> &mut Self {
        self.iterate = iterations;
        self
    }

    /// Sets lambda for Box-Cox transformation.
    pub fn lambda(&mut self, lambda: f32) -> &mut Self {
        self.lambda = Some(lambda);
        self
    }

    /// Sets the lengths of the seasonal smoothers.
    pub fn seasonal_lengths(&mut self, lengths: &[usize]) -> &mut Self {
        self.swin = Some(lengths.to_vec());
        self
    }

    /// Sets the STL parameters.
    pub fn stl_params(&mut self, stl_params: StlParams) -> &mut Self {
        self.stl_params = stl_params;
        self
    }

    /// Decomposes a time series.
    pub fn fit(&self, series: &[f32], periods: &[usize]) -> Result<MstlResult, Error> {
        // return error to be consistent with stl
        // and ensure seasonal is always same length as periods
        if periods.iter().any(|&v| v < 2) {
            return Err(Error::Parameter("periods must be at least 2".to_string()));
        }

        // return error to be consistent with stl
        // and ensure seasonal is always same length as periods
        for np in periods {
            if series.len() < np * 2 {
                return Err(Error::Series(
                    "series has less than two periods".to_string(),
                ));
            }
        }

        if let Some(lambda) = self.lambda {
            if !(0.0..=1.0).contains(&lambda) {
                return Err(Error::Parameter(
                    "lambda must be between 0 and 1".to_string(),
                ));
            }
        }

        if let Some(swin) = &self.swin {
            if swin.len() != periods.len() {
                return Err(Error::Parameter(
                    "seasonal_lengths must have the same length as periods".to_string(),
                ));
            }
        }

        let (trend, remainder, seasonal) = mstl(
            series,
            periods,
            self.iterate,
            self.lambda,
            &self.swin,
            &self.stl_params,
        )?;

        Ok(MstlResult {
            seasonal,
            trend,
            remainder,
        })
    }
}

impl Default for MstlParams {
    fn default() -> Self {
        Self::new()
    }
}
