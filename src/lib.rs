#![doc = include_str!("../README.md")]
#![cfg_attr(feature = "no_std", no_std)]
#![cfg_attr(feature = "no_std", feature(core_float_math))]

mod error;
mod stl;
mod stl_impl;
mod stl_params;

#[cfg(not(feature = "no_std"))]
mod mstl;
#[cfg(not(feature = "no_std"))]
mod mstl_impl;
#[cfg(not(feature = "no_std"))]
mod mstl_params;
#[cfg(not(feature = "no_std"))]
mod mstl_result;
#[cfg(not(feature = "no_std"))]
mod stl_result;

pub use error::Error;
pub use stl::Stl;
pub use stl_params::StlParams;

#[cfg(not(feature = "no_std"))]
pub use {mstl::Mstl, mstl_params::MstlParams, mstl_result::MstlResult, stl_result::StlResult};

/// Creates a new set of STL parameters.
pub fn params() -> StlParams {
    StlParams::new()
}
