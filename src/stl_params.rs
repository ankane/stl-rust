use super::stl_impl::stl;
use super::Error;

#[cfg(not(feature = "no_std"))]
use super::StlResult;

#[cfg(feature = "no_std")]
use libm::ceilf;

#[cfg(not(feature = "no_std"))]
fn ceilf(x: f32) -> f32 {
    x.ceil()
}

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
    #[cfg(not(feature = "no_std"))]
    pub fn fit(&self, series: &[f32], period: usize) -> Result<StlResult, Error> {
        let n = series.len();
        let np = period;

        // check before allocating
        if n / 2 < np {
            return Err(Error::Series);
        }
        let np = np.max(2);

        let mut seasonal = vec![0.0; n];
        let mut trend = vec![0.0; n];
        let mut weights = vec![0.0; n];
        let mut work = vec![0.0; (n + 2 * np) * 5];

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

    /// Decomposes a time series.
    #[cfg(feature = "no_std")]
    pub fn fit(
        &self,
        series: &[f32],
        period: usize,
        seasonal: &mut [f32],
        trend: &mut [f32],
        weights: &mut [f32],
        work: &mut [f32],
    ) -> Result<(), Error> {
        if series.len() / 2 < period {
            return Err(Error::Series);
        }
        self.fit_impl(series, period, seasonal, trend, weights, work)
    }

    fn fit_impl(
        &self,
        series: &[f32],
        period: usize,
        seasonal: &mut [f32],
        trend: &mut [f32],
        weights: &mut [f32],
        work: &mut [f32],
    ) -> Result<(), Error> {
        let np = period;
        let ns = self.ns.unwrap_or(np);

        let isdeg = self.isdeg;
        let itdeg = self.itdeg;

        let ildeg = self.ildeg.unwrap_or(itdeg);
        let mut newns = ns.max(3);
        if newns % 2 == 0 {
            newns += 1;
        }

        let newnp = np.max(2);
        let mut nt = ceilf((1.5 * newnp as f32) / (1.0 - 1.5 / newns as f32)) as usize;
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

        let nsjump = self.nsjump.unwrap_or(ceilf((newns as f32) / 10.0) as usize);
        let ntjump = self.ntjump.unwrap_or(ceilf((nt as f32) / 10.0) as usize);
        let nljump = self.nljump.unwrap_or(ceilf((nl as f32) / 10.0) as usize);

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

#[cfg(test)]
#[cfg(feature = "no_std")]
mod tests {
    use crate::StlParams;

    fn assert_in_delta(exp: f32, act: f32) {
        assert!((exp - act).abs() < 0.001);
    }

    fn assert_elements_in_delta(exp: &[f32], act: &[f32]) {
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

        let mut seasonal = [0.0; 30];
        let mut trend = [0.0; 30];
        let mut weights = [0.0; 30];
        let mut work = [0.0; (30 + 2 * 7) * 5];

        StlParams::new()
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
            &[0.36926576, 0.75655484, -1.3324139, 1.9553658, -0.6044802],
            &seasonal[..5],
        );
        assert_elements_in_delta(
            &[4.804099, 4.9097075, 5.015316, 5.16045, 5.305584],
            &trend[..5],
        );
        assert_elements_in_delta(&[1.0, 1.0, 1.0, 1.0, 1.0], &weights[..5]);
    }
}
