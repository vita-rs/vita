//! Standard chemical notations: facts rendered as symbols.
//!
//! Each submodule is one notation — [`formula`] the Hill molecular
//! formula — writing a chemistry fact to its standard string and reading
//! the established vocabulary back. Writers are total and canonical, one
//! dialect deterministically; readers accept the notation's whole
//! vocabulary and repair nothing. A failed read locates itself as a
//! [`ParseError`], a byte offset paired with the notation's own error kind.

mod error;

pub mod formula;

pub use error::ParseError;
