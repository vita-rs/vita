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
        /// One-based byte column, when available.
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::error::Error as StdError;

    #[derive(Clone, Debug, PartialEq, Eq)]
    enum Kind {
        Unexpected,
        Invalid(String),
    }

    impl fmt::Display for Kind {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            match self {
                Kind::Unexpected => f.write_str("unexpected end of input"),
                Kind::Invalid(s) => write!(f, "invalid value: {s}"),
            }
        }
    }

    #[test]
    fn location_text_with_column() {
        assert_eq!(
            Location::Text {
                line: 3,
                column: Some(7)
            }
            .to_string(),
            "3:7",
        );
    }

    #[test]
    fn location_text_without_column() {
        assert_eq!(
            Location::Text {
                line: 3,
                column: None
            }
            .to_string(),
            "3",
        );
    }

    #[test]
    fn location_binary() {
        assert_eq!(
            Location::Binary { offset: 0x1a2b }.to_string(),
            "offset 0x1a2b",
        );
    }

    #[test]
    fn location_binary_zero() {
        assert_eq!(Location::Binary { offset: 0 }.to_string(), "offset 0x0",);
    }

    #[test]
    fn location_clone_and_eq() {
        let a = Location::Text {
            line: 1,
            column: Some(1),
        };
        assert_eq!(a.clone(), a);
        assert_ne!(
            a,
            Location::Text {
                line: 2,
                column: Some(1)
            }
        );
    }

    #[test]
    fn parse_error_display_text() {
        let e = ParseError::new(
            Location::Text {
                line: 5,
                column: Some(3),
            },
            Kind::Unexpected,
        );
        assert_eq!(e.to_string(), "5:3: unexpected end of input");
    }

    #[test]
    fn parse_error_display_binary() {
        let e = ParseError::new(
            Location::Binary { offset: 0xff },
            Kind::Invalid("bad".into()),
        );
        assert_eq!(e.to_string(), "offset 0xff: invalid value: bad");
    }

    #[test]
    fn parse_error_clone_and_eq() {
        let a = ParseError::new(
            Location::Text {
                line: 1,
                column: None,
            },
            Kind::Unexpected,
        );
        let b = ParseError::new(
            Location::Text {
                line: 2,
                column: None,
            },
            Kind::Unexpected,
        );
        assert_eq!(a.clone(), a);
        assert_ne!(a, b);
    }

    #[test]
    fn error_from_io() {
        let e: Error<Kind> = io::Error::new(io::ErrorKind::UnexpectedEof, "eof").into();
        assert!(matches!(e, Error::Io(_)));
    }

    #[test]
    fn error_from_parse() {
        let pe = ParseError::new(
            Location::Text {
                line: 1,
                column: None,
            },
            Kind::Unexpected,
        );
        let e: Error<Kind> = pe.into();
        assert!(matches!(e, Error::Parse(_)));
    }

    #[test]
    fn error_io_display() {
        let e: Error<Kind> = io::Error::new(io::ErrorKind::UnexpectedEof, "end of stream").into();
        assert!(e.to_string().contains("end of stream"));
    }

    #[test]
    fn error_parse_display() {
        let e: Error<Kind> = Error::Parse(ParseError::new(
            Location::Text {
                line: 2,
                column: Some(1),
            },
            Kind::Unexpected,
        ));
        assert_eq!(e.to_string(), "2:1: unexpected end of input");
    }

    #[test]
    fn error_source_io_is_some() {
        let e: Error<Kind> = io::Error::other("x").into();
        assert!(e.source().is_some());
    }

    #[test]
    fn error_source_parse_is_none() {
        let e: Error<Kind> = Error::Parse(ParseError::new(
            Location::Text {
                line: 1,
                column: None,
            },
            Kind::Unexpected,
        ));
        assert!(e.source().is_none());
    }
}
