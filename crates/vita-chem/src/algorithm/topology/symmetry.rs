//! The interchangeable sites of a molecule.
//!
//! [`orbits`] partitions the sites into [`Orbit`]s — the classes the bare
//! skeleton's automorphisms carry onto one another. Two sites in one orbit
//! are the same place in the graph, relabeled.

mod orbits;

pub use orbits::{Orbit, Orbits, orbits};
