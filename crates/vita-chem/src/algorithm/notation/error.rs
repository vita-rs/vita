use std::fmt;

/// A parse failure: the kind of the first offending token, and where it
/// starts.
///
/// Every notation parser reports failures this way — the byte offset locates
/// the token in the input the caller already holds, and the kind, one of the
/// notation's own error enum, classifies it. Nothing is copied out, so the
/// error stays plain data.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct ParseError<K> {
    at: usize,
    kind: K,
}

impl<K> ParseError<K> {
    /// Constructs a parse error of `kind` at byte offset `at`.
    pub fn new(at: usize, kind: K) -> Self {
        ParseError { at, kind }
    }

    /// Byte offset of the offending token in the input.
    pub fn at(&self) -> usize {
        self.at
    }

    /// What went wrong there.
    pub fn kind(&self) -> K
    where
        K: Copy,
    {
        self.kind
    }
}

impl<K: fmt::Display> fmt::Display for ParseError<K> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} at byte {}", self.kind, self.at)
    }
}

impl<K: fmt::Debug + fmt::Display> std::error::Error for ParseError<K> {}

#[cfg(test)]
mod tests {
    use super::*;

    fn error() -> ParseError<&'static str> {
        ParseError::new(7, "bad token")
    }

    #[test]
    fn at_returns_the_offset() {
        assert_eq!(error().at(), 7);
    }

    #[test]
    fn kind_returns_the_kind() {
        assert_eq!(error().kind(), "bad token");
    }

    #[test]
    fn display_shows_the_kind_then_the_offset() {
        assert_eq!(error().to_string(), "bad token at byte 7");
    }
}
