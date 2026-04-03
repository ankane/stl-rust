// Bandara, K., Hyndman, R. J., & Bergmeir, C. (2021).
// MSTL: A Seasonal-Trend Decomposition Algorithm for Time Series with Multiple Seasonal Patterns.
// arXiv:2107.13462 [stat.AP]. https://doi.org/10.48550/arXiv.2107.13462

use super::{Error, StlParams};

#[allow(clippy::too_many_arguments)]
pub fn mstl(
    x: &[f32],
    seas_ids: &[usize],
    iterate: usize,
    lambda: Option<f32>,
    swin: &Option<Vec<usize>>,
    stl_params: &StlParams,
    trend: &mut [f32],
    remainder: &mut [f32],
) -> Result<Vec<Vec<f32>>, Error> {
    let k = x.len();

    // keep track of indices instead of sorting seas_ids
    // so order is preserved with seasonality
    let mut indices: Vec<usize> = (0..seas_ids.len()).collect();
    indices.sort_by_key(|&i| &seas_ids[i]);

    let mut iterate = iterate;
    if seas_ids.len() == 1 {
        iterate = 1;
    }

    let mut seasonality = Vec::with_capacity(seas_ids.len());

    let mut deseas = if let Some(lam) = lambda {
        box_cox(x, lam)
    } else {
        x.to_vec()
    };

    if !seas_ids.is_empty() {
        for _ in 0..seas_ids.len() {
            seasonality.push(vec![0.0; k]);
        }

        for j in 0..iterate {
            for (i, &idx) in indices.iter().enumerate() {
                if j > 0 {
                    for (d, s) in deseas.iter_mut().zip(&seasonality[idx]) {
                        *d += s;
                    }
                }

                let mut params = stl_params.clone();
                if let Some(sw) = &swin {
                    params.seasonal_length(sw[idx]);
                } else if stl_params.ns.is_none() {
                    params.seasonal_length(7 + 4 * (i + 1));
                }

                // TODO confirm needed
                seasonality[idx].fill(0.0);
                trend.fill(0.0);
                let mut weights = vec![0.0; k];
                let mut work = vec![0.0; (k + 2 * seas_ids[idx]) * 5];

                params.fit_impl(
                    &deseas,
                    seas_ids[idx],
                    &mut seasonality[idx],
                    trend,
                    &mut weights,
                    &mut work,
                )?;

                for (d, s) in deseas.iter_mut().zip(&seasonality[idx]) {
                    *d -= s;
                }
            }
        }
    } else {
        // TODO use Friedman's Super Smoother for trend
        return Err(Error::EmptyPeriods);
    }

    for i in 0..k {
        remainder[i] = deseas[i] - trend[i];
    }

    Ok(seasonality)
}

fn box_cox(y: &[f32], lambda: f32) -> Vec<f32> {
    if lambda != 0.0 {
        y.iter()
            .map(|yi| (yi.powf(lambda) - 1.0) / lambda)
            .collect()
    } else {
        y.iter().map(|yi| yi.ln()).collect()
    }
}
