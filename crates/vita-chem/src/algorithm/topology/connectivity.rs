//! Connected components and biconnected blocks: how a molecule holds
//! together.
//!
//! [`components`] partitions the sites into maximal mutually reachable
//! sets — the molecule's separate pieces. [`blocks`] partitions the bonds
//! one level finer, into maximal sets whose bonds pairwise share a cycle:
//! each [`Block`] is a lone bridge bond or a 2-connected subgraph.

mod blocks;
mod components;

pub use blocks::{Block, Blocks, blocks};
pub use components::{Component, Components, components};
