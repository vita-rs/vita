//! I/O for the vita ecosystem: maps byte streams to vita capabilities and back.
//!
//! A format's reader produces a concrete type that implements exactly the
//! capabilities the format records — no more. Its writer accepts any type bounded
//! on those same capabilities, whatever its origin.
//!
//! # Formats
//!
//! | Module | Read output | Write bounds |
//! |--------|-------------|--------------|
//! | [`xyz`] | [`xyz::System`] | [`HasElements`](vita_core::HasElements) + [`HasPositions<V>`](vita_core::HasPositions) |
//!
//! # Errors
//!
//! Every format exposes [`Error<K>`](Error) — either a [`std::io::Error`] or a
//! [`ParseError<K>`](ParseError) pinning a format-specific kind to a [`Location`]
//! (text: line/column; binary: byte offset).

mod error;

pub mod xyz;

pub use error::{Error, Location, ParseError};
