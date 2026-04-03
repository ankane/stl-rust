use super::stl_impl::stl;
use super::{Error, Float};

#[cfg(feature = "alloc")]
use super::StlResult;

#[cfg(feature = "alloc")]
use alloc::{vec, vec::Vec};

#[cfg(feature = "std")]
fn ceil(x: f32) -> f32 {
    x.ceil()
}

#[cfg(not(feature = "std"))]
use core::f32::math::ceil;

/// A set of STL parameters.
#[derive(Clone, Debug)]
pub struct StlParams {
    pub(crate) ns: Option<usize>,
    nt: Option<usize>,
    nl: Option<usize>,
    isdeg: i32,
    itdeg: i32,
    ildeg: Option<i32>,
    nsjump: Option<usize>,
    ntjump: Option<usize>,
    nljump: Option<usize>,
    ni: Option<usize>,
    no: Option<usize>,
    robust: bool,
}

impl StlParams {
    /// Creates a new set of parameters.
    pub fn new() -> Self {
        Self {
            ns: None,
            nt: None,
            nl: None,
            isdeg: 0,
            itdeg: 1,
            ildeg: None,
            nsjump: None,
            ntjump: None,
            nljump: None,
            ni: None,
            no: None,
            robust: false,
        }
    }

    /// Sets the length of the seasonal smoother.
    pub fn seasonal_length(&mut self, length: usize) -> &mut Self {
        self.ns = Some(length);
        self
    }

    /// Sets the length of the trend smoother.
    pub fn trend_length(&mut self, length: usize) -> &mut Self {
        self.nt = Some(length);
        self
    }

    /// Sets the length of the low-pass filter.
    pub fn low_pass_length(&mut self, length: usize) -> &mut Self {
        self.nl = Some(length);
        self
    }

    /// Sets the degree of locally-fitted polynomial in seasonal smoothing.
    pub fn seasonal_degree(&mut self, degree: i32) -> &mut Self {
        self.isdeg = degree;
        self
    }

    /// Sets the degree of locally-fitted polynomial in trend smoothing.
    pub fn trend_degree(&mut self, degree: i32) -> &mut Self {
        self.itdeg = degree;
        self
    }

    /// Sets the degree of locally-fitted polynomial in low-pass smoothing.
    pub fn low_pass_degree(&mut self, degree: i32) -> &mut Self {
        self.ildeg = Some(degree);
        self
    }

    /// Sets the skipping value for seasonal smoothing.
    pub fn seasonal_jump(&mut self, jump: usize) -> &mut Self {
        self.nsjump = Some(jump);
        self
    }

    /// Sets the skipping value for trend smoothing.
    pub fn trend_jump(&mut self, jump: usize) -> &mut Self {
        self.ntjump = Some(jump);
        self
    }

    /// Sets the skipping value for low-pass smoothing.
    pub fn low_pass_jump(&mut self, jump: usize) -> &mut Self {
        self.nljump = Some(jump);
        self
    }

    /// Sets the number of loops for updating the seasonal and trend components.
    pub fn inner_loops(&mut self, loops: usize) -> &mut Self {
        self.ni = Some(loops);
        self
    }

    /// Sets the number of iterations of robust fitting.
    pub fn outer_loops(&mut self, loops: usize) -> &mut Self {
        self.no = Some(loops);
        self
    }

    /// Sets whether robustness iterations are to be used.
    pub fn robust(&mut self, robust: bool) -> &mut Self {
        self.robust = robust;
        self
    }

    /// Decomposes a time series.
    #[cfg(feature = "alloc")]
    pub fn fit<T: Float>(&self, series: &[T], period: usize) -> Result<StlResult<T>, Error> {
        let n = series.len();
        let np = period;

        // check before allocating
        if n / 2 < np {
            return Err(Error::Series);
        }
        let np = np.max(2);

        let mut seasonal = vec![T::zero(); n];
        let mut trend = vec![T::zero(); n];
        let mut weights = vec![T::zero(); n];
        let mut work = vec![T::zero(); (n + 2 * np) * 5];

        self.fit_impl(
            series,
            period,
            &mut seasonal,
            &mut trend,
            &mut weights,
            &mut work,
        )?;

        let mut remainder = Vec::with_capacity(n);
        for i in 0..n {
            remainder.push(series[i] - seasonal[i] - trend[i]);
        }

        Ok(StlResult {
            seasonal,
            trend,
            remainder,
            weights,
        })
    }

    /// Decomposes a time series with zero allocations.
    #[cfg(not(feature = "alloc"))]
    pub fn fit_zero<T: Float>(
        &self,
        series: &[T],
        period: usize,
        seasonal: &mut [T],
        trend: &mut [T],
        weights: &mut [T],
        work: &mut [T],
    ) -> Result<(), Error> {
        let n = series.len();
        let np = period;

        if n / 2 < np {
            return Err(Error::Series);
        }
        let np = np.max(2);

        debug_assert!(seasonal.len() >= n);
        debug_assert!(trend.len() >= n);
        debug_assert!(weights.len() >= n);
        debug_assert!(work.len() >= (n + 2 * np) * 5);

        self.fit_impl(series, period, seasonal, trend, weights, work)
    }

    pub(crate) fn fit_impl<T: Float>(
        &self,
        series: &[T],
        period: usize,
        seasonal: &mut [T],
        trend: &mut [T],
        weights: &mut [T],
        work: &mut [T],
    ) -> Result<(), Error> {
        if period < 2 {
            return Err(Error::Period);
        }

        let np = period;
        let ns = self.ns.unwrap_or(np);

        let isdeg = self.isdeg;
        let itdeg = self.itdeg;
        let ildeg = self.ildeg.unwrap_or(itdeg);

        let mut newns = ns.max(3);
        if newns % 2 == 0 {
            newns += 1;
        }

        let newnp = np;
        let mut nt = ceil((1.5 * newnp as f32) / (1.0 - 1.5 / newns as f32)) as usize;
        nt = self.nt.unwrap_or(nt);
        nt = nt.max(3);
        if nt % 2 == 0 {
            nt += 1;
        }

        let mut nl = self.nl.unwrap_or(newnp);
        if nl % 2 == 0 && self.nl.is_none() {
            nl += 1;
        }

        let ni = self.ni.unwrap_or(if self.robust { 1 } else { 2 });
        let no = self.no.unwrap_or(if self.robust { 15 } else { 0 });

        let nsjump = self.nsjump.unwrap_or(ceil((newns as f32) / 10.0) as usize);
        let ntjump = self.ntjump.unwrap_or(ceil((nt as f32) / 10.0) as usize);
        let nljump = self.nljump.unwrap_or(ceil((nl as f32) / 10.0) as usize);

        if isdeg != 0 && isdeg != 1 {
            return Err(Error::SeasonalDegree);
        }
        if itdeg != 0 && itdeg != 1 {
            return Err(Error::TrendDegree);
        }
        if ildeg != 0 && ildeg != 1 {
            return Err(Error::LowPassDegree);
        }

        debug_assert!(newnp >= 2);
        debug_assert!(newns % 2 == 1);
        debug_assert!(newns >= 3);
        debug_assert!(nt % 2 == 1);
        debug_assert!(nt >= 3);
        debug_assert!(nl % 2 == 1);
        debug_assert!(nl >= 3);

        trend.fill(T::zero());

        stl(
            series, newnp, newns, nt, nl, isdeg, itdeg, ildeg, nsjump, ntjump, nljump, ni, no,
            weights, seasonal, trend, work,
        );

        Ok(())
    }
}

impl Default for StlParams {
    fn default() -> Self {
        Self::new()
    }
}
