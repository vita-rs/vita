//! The molecular graph as pure structure: what sites and bonds alone
//! determine.
//!
//! [`connectivity`] partitions the graph into connected components and
//! biconnected blocks; [`path`] walks it, counting distance in bond hops;
//! [`ring`] perceives its cycles, from a membership bit up to the unique
//! ring families; [`scaffold`] strips it to the Bemis–Murcko framework;
//! [`symmetry`] groups the sites its automorphisms interchange. Elements
//! and bond orders never enter: every answer is a fact of the bare
//! skeleton.

pub mod connectivity;
pub mod path;
pub mod ring;
pub mod scaffold;
pub mod symmetry;
