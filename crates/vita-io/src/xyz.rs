//! The standard XYZ format.
//!
//! A text format for atomic symbols and Cartesian coordinates. A file contains one or
//! more frames, each:
//!
//! ```text
//! <count>
//! <comment>
//! <symbol> <x> <y> <z>
//!     ⋮  (count lines)
//! ```
//!
//! Fields are whitespace-separated; element symbols are case-insensitive.
//! Each frame yields a [`System`]: [`HasSites`](vita_core::HasSites),
//! [`HasElements`](vita_core::HasElements), [`HasPositions`](vita_core::HasPositions).
//!
//! # Reading
//!
//! [`read`] returns a [`Reader`]: one frame via [`Reader::system`], or a lazy frame
//! iterator via [`Reader::systems`].
//!
//! # Writing
//!
//! [`write`](write()) accepts any `S: HasElements + HasPositions<V>` — not only [`System`].

mod error;
mod read;
mod system;
mod write;

pub use error::{Error, ErrorKind, ParseError};
pub use read::{Reader, Systems, read};
pub use system::System;
pub use write::{Config, write};
