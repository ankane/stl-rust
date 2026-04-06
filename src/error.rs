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
    Period,
    Lambda,
    SeasonalLengths,
    EmptyPeriods,
}

impl error::Error for Error {}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Error::Series => f.write_str("series has less than two periods"),
            Error::SeasonalDegree => f.write_str("seasonal_degree must be 0 or 1"),
            Error::TrendDegree => f.write_str("trend_degree must be 0 or 1"),
            Error::LowPassDegree => f.write_str("low_pass_degree must be 0 or 1"),
            Error::Period => f.write_str("period must be at least 2"),
            Error::Lambda => f.write_str("lambda must be between 0 and 1"),
            Error::SeasonalLengths => {
                f.write_str("seasonal_lengths must have the same length as periods")
            }
            Error::EmptyPeriods => f.write_str("periods must not be empty"),
        }
    }
}
