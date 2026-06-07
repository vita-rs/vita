use std::fmt;
use std::io;

/// Position in a data source where a parse error was detected.
///
/// Line and column numbers are one-based. Text formats carry the exact line
/// and, when available, the column; binary formats carry the byte offset from
/// the start of the data.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Location {
    /// A position in a text source.
    Text {
        /// One-based line number.
        line: u32,
        /// One-based column number, when available.
        column: Option<u32>,
    },
    /// A byte position in a binary source.
    Binary {
        /// Zero-based byte offset from the start of the data.
        offset: u64,
    },
}

impl fmt::Display for Location {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Location::Text {
                line,
                column: Some(col),
            } => write!(f, "{}:{}", line, col),
            Location::Text { line, column: None } => write!(f, "{}", line),
            Location::Binary { offset } => write!(f, "offset {:#x}", offset),
        }
    }
}

/// A parse error at a known [`Location`] in the source.
#[derive(Clone, Debug, PartialEq)]
pub struct ParseError<K> {
    /// Where in the source the error was detected.
    pub location: Location,
    /// The format-specific error kind.
    pub kind: K,
}

impl<K> ParseError<K> {
    /// Constructs a `ParseError` at `location` with the given `kind`.
    #[inline]
    pub fn new(location: Location, kind: K) -> Self {
        Self { location, kind }
    }
}

impl<K: fmt::Display> fmt::Display for ParseError<K> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.location, self.kind)
    }
}

impl<K: fmt::Debug + fmt::Display> std::error::Error for ParseError<K> {}

/// The error type for Vita I/O operations, which may be either an underlying I/O error
/// or a parse error at a specific location in the source.
#[derive(Debug)]
pub enum Error<K> {
    /// An underlying I/O error.
    Io(io::Error),
    /// A parse error at a specific location in the source.
    Parse(ParseError<K>),
}

impl<K: fmt::Display> fmt::Display for Error<K> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::Io(e) => fmt::Display::fmt(e, f),
            Error::Parse(e) => fmt::Display::fmt(e, f),
        }
    }
}

impl<K: fmt::Debug + fmt::Display> std::error::Error for Error<K> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Error::Io(e) => Some(e),
            Error::Parse(_) => None,
        }
    }
}

impl<K> From<io::Error> for Error<K> {
    #[inline]
    fn from(e: io::Error) -> Self {
        Self::Io(e)
    }
}

impl<K> From<ParseError<K>> for Error<K> {
    #[inline]
    fn from(e: ParseError<K>) -> Self {
        Self::Parse(e)
    }
}
