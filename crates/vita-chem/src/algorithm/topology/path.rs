//! Walks over the bond skeleton, with distance counted in bond hops.
//!
//! [`bfs`] and [`dfs`] traverse lazily from a start site, nearest-first or
//! branch-first; [`path`] picks one shortest path between two sites,
//! [`paths`] keeps every one that ties, and [`distances`] answers all pairs
//! at once as a [`DistanceMatrix`]. Separation is a fact, not a failure: an
//! unreachable pair is an absent distance or an empty path set.

mod bfs;
mod dfs;
mod distances;
mod paths;
mod shortest;

pub use bfs::bfs;
pub use dfs::dfs;
pub use distances::{DistanceMatrix, distances};
pub use paths::{Paths, paths};
pub use shortest::path;
