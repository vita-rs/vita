//! The Hückel model of aromaticity, perceived and resolved.
//!
//! [`perceive`] finds the rings satisfying Hückel's rule as an
//! [`Aromaticity`], viewable over the molecule as [`WithAromaticity`];
//! [`kekulize`] goes the other way, resolving delocalized
//! [`Aromatic`](crate::BondOrder::Aromatic) bonds into one localized
//! single/double pattern — a [`Kekule`], viewable as [`WithKekule`]. Both
//! are overlays on an unchanged molecule.

mod kekulize;
mod perceive;

pub use kekulize::{Kekule, WithKekule, kekulize};
pub use perceive::{Aromaticity, WithAromaticity, perceive};
