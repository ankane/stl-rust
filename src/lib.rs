//! Seasonal-trend decomposition for Rust
//!
//! [View the docs](https://github.com/ankane/stl-rust)

mod error;
mod params;
mod result;
mod stl;
mod stl_impl;

pub use error::Error;
pub use params::StlParams;
pub use result::StlResult;
pub use stl::Stl;

pub fn params() -> StlParams {
    StlParams::new()
}
