use std::fmt;

#[derive(Debug, PartialEq)]
pub enum Error {
    Parameter(String),
    Series(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Error::Parameter(ref err) => write!(f, "{}", err.as_str()),
            Error::Series(ref err) => write!(f, "{}", err.as_str()),
        }
    }
}
