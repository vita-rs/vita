//! Reduction to the Bemis–Murcko framework.
//!
//! [`framework`] assigns every atom a structural [`Role`] — ring, linker,
//! or side chain — and keeps as the [`Framework`] the ring-and-linker core
//! left standing once the side chains fall away.

mod framework;

pub use framework::{Framework, Role, framework};
