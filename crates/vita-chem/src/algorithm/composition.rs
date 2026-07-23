//! The particle inventory: what a molecule is made of.
//!
//! A [`Composition`] counts atoms by [`Constituent`] — an element at its
//! natural isotopic mixture, or one specific nuclide — and carries the net
//! charge, the electron side of the inventory. Connectivity is forgotten by
//! definition: constitutional isomers share one composition. [`elemental`]
//! folds a molecule at element precision, [`isotopic`] at declared-nuclide
//! precision; compositions add as fragments do, and the molecular formula
//! ([`notation::formula`](crate::notation::formula)) is their standard
//! rendering.

mod constituent;
mod elemental;
mod isotopic;

pub use constituent::{Composition, Constituent};
pub use elemental::elemental;
pub use isotopic::isotopic;
