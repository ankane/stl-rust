use core::error;
use core::fmt;

/// An error.
#[derive(Debug, Eq, PartialEq)]
pub enum Error {
    Parameter(&'static str),
    Series(&'static str),
}

impl error::Error for Error {}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Error::Parameter(err) => f.write_str(err),
            Error::Series(err) => f.write_str(err),
        }
    }
}
