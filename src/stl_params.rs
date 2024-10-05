use super::stl_impl::stl;
use super::{Error, StlResult};

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
    pub fn fit(&self, series: &[f32], period: usize) -> Result<StlResult, Error> {
        let y = series;
        let np = period;
        let n = y.len();

        if n < np * 2 {
            return Err(Error::Series(
                "series has less than two periods".to_string(),
            ));
        }

        let ns = self.ns.unwrap_or(np);

        let isdeg = self.isdeg;
        let itdeg = self.itdeg;

        let mut rw = vec![0.0; n];
        let mut season = vec![0.0; n];
        let mut trend = vec![0.0; n];

        let ildeg = self.ildeg.unwrap_or(itdeg);
        let mut newns = ns.max(3);
        if newns % 2 == 0 {
            newns += 1;
        }

        let newnp = np.max(2);
        let mut nt = ((1.5 * newnp as f32) / (1.0 - 1.5 / newns as f32)).ceil() as usize;
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

        let nsjump = self
            .nsjump
            .unwrap_or(((newns as f32) / 10.0).ceil() as usize);
        let ntjump = self.ntjump.unwrap_or(((nt as f32) / 10.0).ceil() as usize);
        let nljump = self.nljump.unwrap_or(((nl as f32) / 10.0).ceil() as usize);

        if newns < 3 {
            return Err(Error::Parameter(
                "seasonal_length must be at least 3".to_string(),
            ));
        }
        if nt < 3 {
            return Err(Error::Parameter(
                "trend_length must be at least 3".to_string(),
            ));
        }
        if nl < 3 {
            return Err(Error::Parameter(
                "low_pass_length must be at least 3".to_string(),
            ));
        }
        if newnp < 2 {
            return Err(Error::Parameter("period must be at least 2".to_string()));
        }

        if isdeg != 0 && isdeg != 1 {
            return Err(Error::Parameter(
                "seasonal_degree must be 0 or 1".to_string(),
            ));
        }
        if itdeg != 0 && itdeg != 1 {
            return Err(Error::Parameter("trend_degree must be 0 or 1".to_string()));
        }
        if ildeg != 0 && ildeg != 1 {
            return Err(Error::Parameter(
                "low_pass_degree must be 0 or 1".to_string(),
            ));
        }

        if newns % 2 != 1 {
            return Err(Error::Parameter("seasonal_length must be odd".to_string()));
        }
        if nt % 2 != 1 {
            return Err(Error::Parameter("trend_length must be odd".to_string()));
        }
        if nl % 2 != 1 {
            return Err(Error::Parameter("low_pass_length must be odd".to_string()));
        }

        stl(
            y,
            n,
            newnp,
            newns,
            nt,
            nl,
            isdeg,
            itdeg,
            ildeg,
            nsjump,
            ntjump,
            nljump,
            ni,
            no,
            &mut rw,
            &mut season,
            &mut trend,
        );

        let mut remainder = Vec::with_capacity(n);
        for i in 0..n {
            remainder.push(y[i] - season[i] - trend[i]);
        }

        Ok(StlResult {
            seasonal: season,
            trend,
            remainder,
            weights: rw,
        })
    }
}

impl Default for StlParams {
    fn default() -> Self {
        Self::new()
    }
}
