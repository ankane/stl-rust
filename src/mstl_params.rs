use alloc::{vec, vec::Vec};

use super::mstl_impl::mstl;
use super::{Error, Float, MstlResult, StlParams};

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
    pub fn fit<T: Float>(&self, series: &[T], periods: &[usize]) -> Result<MstlResult<T>, Error> {
        if periods.is_empty() {
            // TODO use Friedman's Super Smoother for trend
            return Err(Error::EmptyPeriods);
        }

        if periods.iter().any(|&v| v < 2) {
            return Err(Error::Period);
        }

        for np in periods {
            if series.len() / 2 < *np {
                return Err(Error::Series);
            }
        }

        if let Some(lambda) = self.lambda {
            if !(0.0..=1.0).contains(&lambda) {
                return Err(Error::Lambda);
            }
        }

        if let Some(swin) = &self.swin {
            if swin.len() != periods.len() {
                return Err(Error::SeasonalLengths);
            }
        }

        let n = series.len();
        let mut seasonal = Vec::with_capacity(periods.len());
        for _ in 0..periods.len() {
            seasonal.push(vec![T::zero(); n]);
        }
        let mut trend = vec![T::zero(); n];
        let mut remainder = vec![T::zero(); n];
        let mut weights = vec![T::zero(); n];
        let mut work = vec![T::zero(); (n + 2 * periods.iter().max().unwrap()) * 5];

        mstl(
            series,
            periods,
            self.iterate,
            self.lambda,
            &self.swin,
            &self.stl_params,
            &mut seasonal,
            &mut trend,
            &mut remainder,
            &mut weights,
            &mut work,
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
