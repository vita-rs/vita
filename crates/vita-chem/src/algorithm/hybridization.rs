//! Hybridization perceived from the molecular graph.
//!
//! [`perceive`] fixes each site's label from its σ neighbors, its lone
//! pairs, and the conjugated systems it sits in — coordinates never enter —
//! and yields a [`Hybridizations`], viewable over the molecule as
//! [`WithHybridizations`]. Sites the model cannot count, foremost the d-
//! and f-block, carry no label.

mod perceive;

pub use perceive::{Hybridizations, WithHybridizations, perceive};
