#![doc = include_str!("../README.md")]

mod error;
mod mstl;
mod mstl_impl;
mod mstl_params;
mod mstl_result;
mod stl;
mod stl_impl;
mod stl_params;
mod stl_result;

pub use error::Error;
pub use mstl::Mstl;
pub use mstl_params::MstlParams;
pub use mstl_result::MstlResult;
pub use stl::Stl;
pub use stl_params::StlParams;
pub use stl_result::StlResult;

/// Creates a new set of STL parameters.
pub fn params() -> StlParams {
    StlParams::new()
}
