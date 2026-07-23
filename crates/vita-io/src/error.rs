use std::{fmt, io};

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

    fn io_error() -> io::Error {
        io::Error::new(io::ErrorKind::UnexpectedEof, "boom")
    }

    fn parse_error() -> ParseError<&'static str> {
        ParseError::new(
            Location::Text {
                line: 3,
                column: Some(7),
            },
            "unexpected",
        )
    }

    #[test]
    fn text_location_without_a_column_shows_only_the_line() {
        let location = Location::Text {
            line: 12,
            column: None,
        };
        assert_eq!(location.to_string(), "12");
    }

    #[test]
    fn text_location_with_a_column_shows_line_and_column() {
        let location = Location::Text {
            line: 12,
            column: Some(5),
        };
        assert_eq!(location.to_string(), "12:5");
    }

    #[test]
    fn binary_location_shows_the_offset_in_hexadecimal() {
        let location = Location::Binary { offset: 255 };
        assert_eq!(location.to_string(), "offset 0xff");
    }

    #[test]
    fn new_stores_the_location_and_kind() {
        let error = ParseError::new(Location::Binary { offset: 8 }, "eof");
        assert_eq!(error.location, Location::Binary { offset: 8 });
        assert_eq!(error.kind, "eof");
    }

    #[test]
    fn parse_error_shows_the_location_then_the_kind() {
        assert_eq!(parse_error().to_string(), "3:7: unexpected");
    }

    #[test]
    fn an_io_error_converts_into_the_io_variant() {
        let error: Error<&str> = io_error().into();
        assert!(matches!(error, Error::Io(_)));
    }

    #[test]
    fn a_parse_error_converts_into_the_parse_variant() {
        let error: Error<&str> = parse_error().into();
        assert!(matches!(error, Error::Parse(_)));
    }

    #[test]
    fn io_variant_source_is_the_underlying_error() {
        let inner = io_error();
        let expected = inner.to_string();
        let error: Error<&str> = Error::Io(inner);
        let source = std::error::Error::source(&error);
        assert_eq!(source.map(|e| e.to_string()), Some(expected));
    }

    #[test]
    fn parse_variant_has_no_source() {
        let error: Error<&str> = Error::Parse(parse_error());
        assert!(std::error::Error::source(&error).is_none());
    }

    #[test]
    fn io_variant_displays_the_underlying_error() {
        let inner = io_error();
        let expected = inner.to_string();
        let error: Error<&str> = Error::Io(inner);
        assert_eq!(error.to_string(), expected);
    }

    #[test]
    fn parse_variant_displays_the_inner_parse_error() {
        let inner = parse_error();
        let expected = inner.to_string();
        let error: Error<&str> = Error::Parse(inner);
        assert_eq!(error.to_string(), expected);
    }
}
