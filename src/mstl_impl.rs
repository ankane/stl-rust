// Bandara, K., Hyndman, R. J., & Bergmeir, C. (2021).
// MSTL: A Seasonal-Trend Decomposition Algorithm for Time Series with Multiple Seasonal Patterns.
// arXiv:2107.13462 [stat.AP]. https://doi.org/10.48550/arXiv.2107.13462

use super::{Error, Float, StlParams};

#[cfg(feature = "alloc")]
use alloc::vec::Vec;

#[allow(clippy::too_many_arguments)]
pub fn mstl<T: Float>(
    x: &[T],
    seas_ids: &[usize],
    iterate: usize,
    lambda: Option<f32>,
    swin: &Option<Vec<usize>>,
    stl_params: &StlParams,
    seasonality: &mut [Vec<T>],
    trend: &mut [T],
    remainder: &mut [T],
    weights: &mut [T],
    work: &mut [T],
) -> Result<(), Error> {
    let k = x.len();

    // keep track of indices instead of sorting seas_ids
    // so order is preserved with seasonality
    let mut indices: Vec<usize> = (0..seas_ids.len()).collect();
    indices.sort_by_key(|&i| &seas_ids[i]);

    let mut iterate = iterate;
    if seas_ids.len() == 1 {
        iterate = 1;
    }

    let mut deseas = if let Some(lam) = lambda {
        box_cox(x, T::from_f64(lam as f64))
    } else {
        x.to_vec()
    };

    for j in 0..iterate {
        for (i, &idx) in indices.iter().enumerate() {
            if j > 0 {
                for (d, s) in deseas.iter_mut().zip(&seasonality[idx]) {
                    *d += *s;
                }
            }

            let mut params = stl_params.clone();
            if let Some(sw) = &swin {
                params.seasonal_length(sw[idx]);
            } else if stl_params.ns.is_none() {
                params.seasonal_length(7 + 4 * (i + 1));
            }

            params.fit_impl(
                &deseas,
                seas_ids[idx],
                &mut seasonality[idx],
                trend,
                weights,
                work,
            )?;

            for (d, s) in deseas.iter_mut().zip(&seasonality[idx]) {
                *d -= *s;
            }
        }
    }

    for i in 0..k {
        remainder[i] = deseas[i] - trend[i];
    }

    Ok(())
}

fn box_cox<T: Float>(y: &[T], lambda: T) -> Vec<T> {
    if lambda != T::zero() {
        y.iter()
            .map(|yi| ((*yi).powf(lambda) - T::one()) / lambda)
            .collect()
    } else {
        y.iter().map(|yi| (*yi).ln()).collect()
    }
}
