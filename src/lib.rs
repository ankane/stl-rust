//! Seasonal-trend decomposition for Rust
//!
//! [View the docs](https://github.com/ankane/stl-rust)

mod error;
mod mstl;
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

pub fn params() -> StlParams {
    StlParams::new()
}
