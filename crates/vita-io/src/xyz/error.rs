use std::fmt;

/// What made an XYZ system invalid.
///
/// This is the format-specific error content; the surrounding [`ParseError`] pins it to
/// a [`Location`](crate::Location) in the input.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ErrorKind {
    /// The input ended before a complete system was read.
    UnexpectedEof,
    /// The atom-count line could not be parsed as a non-negative integer.
    AtomCount {
        /// The text found where a count was expected.
        found: Box<str>,
    },
    /// The atom count is a valid integer but exceeds the maximum representable (`u32::MAX`).
    AtomCountRange {
        /// The count that was declared.
        count: u64,
    },
    /// An atom line did not hold exactly four fields: a symbol and three coordinates.
    FieldCount {
        /// The number of whitespace-separated fields the line actually held.
        found: usize,
    },
    /// An atom's leading field was not a recognized element symbol.
    ElementSymbol {
        /// The unrecognized symbol.
        found: Box<str>,
    },
    /// An atom coordinate was not a finite real number.
    Coordinate {
        /// The text found where a coordinate was expected.
        found: Box<str>,
    },
}

impl fmt::Display for ErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ErrorKind::UnexpectedEof => f.write_str("unexpected end of input"),
            ErrorKind::AtomCount { found } => write!(f, "expected an atom count, found {found:?}"),
            ErrorKind::AtomCountRange { count } => write!(
                f,
                "atom count {count} exceeds the maximum ({max})",
                max = u32::MAX,
            ),
            ErrorKind::FieldCount { found } => {
                write!(f, "expected 4 fields (symbol x y z), found {found}")
            }
            ErrorKind::ElementSymbol { found } => write!(f, "unknown element symbol {found:?}"),
            ErrorKind::Coordinate { found } => write!(f, "invalid coordinate {found:?}"),
        }
    }
}

/// A parse error from reading XYZ: a [`Location`](crate::Location) paired with an
/// [`ErrorKind`].
pub type ParseError = crate::ParseError<ErrorKind>;

/// An error from reading or writing XYZ: either an I/O failure or a [`ParseError`].
pub type Error = crate::Error<ErrorKind>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unexpected_eof_reports_the_end_of_input() {
        assert_eq!(
            ErrorKind::UnexpectedEof.to_string(),
            "unexpected end of input",
        );
    }

    #[test]
    fn atom_count_reports_the_offending_text() {
        assert_eq!(
            ErrorKind::AtomCount {
                found: "two".into()
            }
            .to_string(),
            r#"expected an atom count, found "two""#,
        );
    }

    #[test]
    fn atom_count_range_reports_the_maximum() {
        assert_eq!(
            ErrorKind::AtomCountRange {
                count: 5_000_000_000,
            }
            .to_string(),
            "atom count 5000000000 exceeds the maximum (4294967295)",
        );
    }

    #[test]
    fn field_count_reports_the_number_of_fields_found() {
        assert_eq!(
            ErrorKind::FieldCount { found: 3 }.to_string(),
            "expected 4 fields (symbol x y z), found 3",
        );
    }

    #[test]
    fn element_symbol_reports_the_unknown_symbol() {
        assert_eq!(
            ErrorKind::ElementSymbol { found: "Xx".into() }.to_string(),
            r#"unknown element symbol "Xx""#,
        );
    }

    #[test]
    fn coordinate_reports_the_invalid_text() {
        assert_eq!(
            ErrorKind::Coordinate {
                found: "1.0.0".into()
            }
            .to_string(),
            r#"invalid coordinate "1.0.0""#,
        );
    }
}
