use core::error;
use core::fmt;

/// An error.
#[derive(Debug, Eq, PartialEq)]
#[non_exhaustive]
pub enum Error {
    Series,
    SeasonalDegree,
    TrendDegree,
    LowPassDegree,
    Periods,
    Lambda,
    SeasonalLengths,
    EmptyPeriods,
}

impl error::Error for Error {}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Error::Series => write!(f, "series has less than two periods"),
            Error::SeasonalDegree => write!(f, "seasonal_degree must be 0 or 1"),
            Error::TrendDegree => write!(f, "trend_degree must be 0 or 1"),
            Error::LowPassDegree => write!(f, "low_pass_degree must be 0 or 1"),
            Error::Periods => write!(f, "periods must be at least 2"),
            Error::Lambda => write!(f, "lambda must be between 0 and 1"),
            Error::SeasonalLengths => {
                write!(f, "seasonal_lengths must have the same length as periods")
            }
            Error::EmptyPeriods => write!(f, "periods must not be empty"),
        }
    }
}
