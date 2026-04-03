#![doc = include_str!("../README.md")]
#![cfg_attr(not(feature = "std"), no_std)]
#![cfg_attr(not(feature = "std"), feature(core_float_math))]

mod error;
mod stl;
mod stl_impl;
mod stl_params;

#[cfg(feature = "std")]
mod mstl;
#[cfg(feature = "std")]
mod mstl_impl;
#[cfg(feature = "std")]
mod mstl_params;
#[cfg(feature = "std")]
mod mstl_result;
#[cfg(feature = "std")]
mod stl_result;

pub use error::Error;
pub use stl::Stl;
pub use stl_params::StlParams;

#[cfg(feature = "std")]
pub use {mstl::Mstl, mstl_params::MstlParams, mstl_result::MstlResult, stl_result::StlResult};

/// Creates a new set of STL parameters.
pub fn params() -> StlParams {
    StlParams::new()
}
