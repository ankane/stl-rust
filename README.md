# STL Rust

Seasonal-trend decomposition for Rust

:tada: Zero dependencies

[![Build Status](https://github.com/ankane/stl-rust/workflows/build/badge.svg?branch=master)](https://github.com/ankane/stl-rust/actions)

## Installation

Add this line to your application’s `Cargo.toml` under `[dependencies]`:

```toml
stlrs = "0.2"
```

## Getting Started

Decompose a time series

```rust
let series = vec![
    5.0, 9.0, 2.0, 9.0, 0.0, 6.0, 3.0, 8.0, 5.0, 8.0,
    7.0, 8.0, 8.0, 0.0, 2.0, 5.0, 0.0, 5.0, 6.0, 7.0,
    3.0, 6.0, 1.0, 4.0, 4.0, 4.0, 3.0, 7.0, 5.0, 8.0
];
let period = 7; // period of the seasonal component

let res = stlrs::params().fit(&series, period).unwrap();
```

Get the components

```rust
res.seasonal();
res.trend();
res.remainder();
```

## Robustness

Use robustness iterations

```rust
let res = stlrs::params().robust(true).fit(&series, period).unwrap();
```

Get robustness weights

```rust
res.weights();
```

## Multiple Seasonality

Specify multiple periods [unreleased]

```rust
use stlrs::Mstl;

let periods = [6, 10];
let res = Mstl::fit(&series, &periods).unwrap();
```

## Parameters

Set STL parameters

```rust
stlrs::params()
    .seasonal_length(7)     // length of the seasonal smoother
    .trend_length(15)       // length of the trend smoother
    .low_pass_length(7)     // length of the low-pass filter
    .seasonal_degree(0)     // degree of locally-fitted polynomial in seasonal smoothing
    .trend_degree(1)        // degree of locally-fitted polynomial in trend smoothing
    .low_pass_degree(1)     // degree of locally-fitted polynomial in low-pass smoothing
    .seasonal_jump(1)       // skipping value for seasonal smoothing
    .trend_jump(2)          // skipping value for trend smoothing
    .low_pass_jump(1)       // skipping value for low-pass smoothing
    .inner_loops(2)         // number of loops for updating the seasonal and trend components
    .outer_loops(0)         // number of iterations of robust fitting
    .robust(false)          // if robustness iterations are to be used
```

Set MSTL parameters [unreleased]

```rust
Mstl::params()
    .iterations(2)                 // number of iterations
    .lambda(0.5)                   // lambda for Box-Cox transformation
    .seasonal_lengths(&[11, 15])   // lengths of the seasonal smoothers
    .stl_params(Stl::params())     // STL params
```

## Strength

Get the seasonal strength

```rust
res.seasonal_strength();
```

Get the trend strength

```rust
res.trend_strength();
```

## Credits

This library was ported from the [Fortran implementation](https://www.netlib.org/a/stl).

## References

- [STL: A Seasonal-Trend Decomposition Procedure Based on Loess](https://www.scb.se/contentassets/ca21efb41fee47d293bbee5bf7be7fb3/stl-a-seasonal-trend-decomposition-procedure-based-on-loess.pdf)
- [MSTL: A Seasonal-Trend Decomposition Algorithm for Time Series with Multiple Seasonal Patterns](https://arxiv.org/pdf/2107.13462.pdf)
- [Measuring strength of trend and seasonality](https://otexts.com/fpp2/seasonal-strength.html)

## History

View the [changelog](https://github.com/ankane/stl-rust/blob/master/CHANGELOG.md)

## Contributing

Everyone is encouraged to help improve this project. Here are a few ways you can help:

- [Report bugs](https://github.com/ankane/stl-rust/issues)
- Fix bugs and [submit pull requests](https://github.com/ankane/stl-rust/pulls)
- Write, clarify, or fix documentation
- Suggest or add new features

To get started with development:

```sh
git clone https://github.com/ankane/stl-rust.git
cd stl-rust
cargo test
```
