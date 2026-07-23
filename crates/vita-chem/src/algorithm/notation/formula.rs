//! The molecular formula: a composition rendered in Hill order.
//!
//! [`write`](write()) gives a [`Composition`](crate::composition::Composition)
//! its one canonical formula — Hill order, bracketed nuclides, sign-first
//! charge. [`parse`] reads the notation's whole vocabulary back: units in
//! any order, repeated symbols accumulating, the `D`/`T` aliases, dotted
//! parts with coefficients. The two compose exactly —
//! `parse(&write(c)) == Ok(c)` — and a failed read classifies itself as an
//! [`ErrorKind`] located by a [`ParseError`].

// TODO: `smiles`, `cip`, `smarts`, ...

mod parse;
mod write;

pub use parse::{ErrorKind, ParseError, parse};
pub use write::write;
