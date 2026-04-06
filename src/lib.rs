#![doc = include_str!("../README.md")]
#![cfg_attr(not(feature = "std"), no_std)]
#![cfg_attr(not(feature = "std"), feature(core_float_math))]

#[cfg(feature = "alloc")]
extern crate alloc;

mod error;
mod float;
mod math;
mod stl;
mod stl_impl;
mod stl_params;

#[cfg(feature = "alloc")]
mod mstl;
#[cfg(feature = "alloc")]
mod mstl_impl;
#[cfg(feature = "alloc")]
mod mstl_params;
#[cfg(feature = "alloc")]
mod mstl_result;
#[cfg(feature = "alloc")]
mod stl_result;

pub use error::Error;
pub use float::Float;
pub use stl::Stl;
pub use stl_params::StlParams;

#[cfg(feature = "alloc")]
pub use {mstl::Mstl, mstl_params::MstlParams, mstl_result::MstlResult, stl_result::StlResult};

/// Creates a new set of STL parameters.
pub fn params() -> StlParams {
    StlParams::new()
}
