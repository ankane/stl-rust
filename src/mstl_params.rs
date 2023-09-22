// Bandara, K., Hyndman, R. J., & Bergmeir, C. (2021).
// MSTL: A Seasonal-Trend Decomposition Algorithm for Time Series with Multiple Seasonal Patterns.
// arXiv:2107.13462 [stat.AP]. https://doi.org/10.48550/arXiv.2107.13462

use super::{Error, MstlResult, StlParams};

pub struct MstlParams {
    iterate: usize,
    lambda: Option<f32>,
    stl_params: StlParams,
}

impl MstlParams {
    pub fn new() -> Self {
        Self {
            iterate: 2,
            lambda: None,
            stl_params: StlParams::new(),
        }
    }

    pub fn iterations(&mut self, iterate: usize) -> &mut Self {
        self.iterate = iterate;
        self
    }

    pub fn lambda(&mut self, lambda: f32) -> &mut Self {
        self.lambda = Some(lambda);
        self
    }

    pub fn stl_params(&mut self, stl_params: StlParams) -> &mut Self {
        self.stl_params = stl_params;
        self
    }

    pub fn fit(&self, series: &[f32], periods: &[usize]) -> Result<MstlResult, Error> {
        let x = series;
        let seas_ids = periods;
        let k = x.len();

        // return error to be consistent with stl
        // and ensure seasonal is always same length as seas_ids
        if seas_ids.iter().any(|&v| v < 2) {
            return Err(Error::Parameter("periods must be at least 2".to_string()));
        }

        // return error to be consistent with stl
        // and ensure seasonal is always same length as seas_ids
        for np in seas_ids {
            if k < np * 2 {
                return Err(Error::Series("series has less than two periods".to_string()));
            }
        }

        // keep track of indices instead of sorting seas_ids
        // so order is preserved with seasonality
        let mut indices: Vec<usize> = (0..seas_ids.len()).collect();
        indices.sort_by_key(|&i| &seas_ids[i]);

        let mut iterate = self.iterate;
        if seas_ids.len() == 1 {
            iterate = 1;
        }

        let mut seasonality = Vec::with_capacity(seas_ids.len());
        let mut trend = Vec::new();

        let mut deseas = if let Some(lambda) = self.lambda {
            box_cox(x, lambda)
        } else {
            x.to_vec()
        };

        if !seas_ids.is_empty() {
            for _ in 0..seas_ids.len() {
                seasonality.push(Vec::new());
            }

            for j in 0..iterate {
                for (i, &idx) in indices.iter().enumerate() {
                    let np = seas_ids[idx];

                    if j > 0 {
                        for (d, s) in deseas.iter_mut().zip(&seasonality[idx]) {
                            *d += s;
                        }
                    }

                    // TODO add seasonal_lengths param
                    let fit = if self.stl_params.ns.is_some() {
                        self.stl_params.fit(&deseas, np)?
                    } else {
                        let seasonal_length = 7 + 4 * (i + 1);
                        self.stl_params.clone().seasonal_length(seasonal_length).fit(&deseas, np)?
                    };

                    (seasonality[idx], trend, _, _) = fit.into_parts();

                    for (d, s) in deseas.iter_mut().zip(&seasonality[idx]) {
                        *d -= s;
                    }
                }
            }
        } else {
            // TODO use Friedman's Super Smoother for trend
            return Err(Error::Parameter("periods must not be empty".to_string()));
        }

        let mut remainder = Vec::with_capacity(k);
        for i in 0..k {
            remainder.push(deseas[i] - trend[i]);
        }

        Ok(MstlResult {
            seasonal: seasonality,
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

fn box_cox(y: &[f32], lambda: f32) -> Vec<f32> {
    if lambda > 0.0 {
        y.iter().map(|yi| (yi.powf(lambda) - 1.0) / lambda).collect()
    } else {
        y.iter().map(|yi| yi.ln()).collect()
    }
}
