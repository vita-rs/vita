//! Conjugated π systems, perceived from the Lewis graph alone.
//!
//! A site contributes a p-orbital where its octet arithmetic leaves one — a
//! multiple bond, a lone pair, an unpaired electron, or a vacancy — and
//! [`systems`] gathers the maximal interacting networks into a
//! [`ConjugatedSystems`] partition of [`ConjugatedSystem`]s. Orthogonality
//! the graph itself forces is honored — a cumulated site parts its two π
//! bonds into perpendicular systems — and every Lewis form of a molecule
//! yields the same partition.

mod systems;

pub use systems::{ConjugatedSystem, ConjugatedSystems, systems};
