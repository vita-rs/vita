//! Structure found inside structure.
//!
//! [`matches`](matches()) enumerates every [`Mapping`] of a pattern onto a
//! subgraph of a target — an injection of sites under which each pattern bond
//! has its counterpart; [`mcs`] finds the largest connected substructure two
//! molecules share, as a [`CommonSubgraph`] carrying the correspondence it
//! induces.

mod mcs;
mod subgraph;

pub use mcs::{CommonSubgraph, mcs};
pub use subgraph::{Mapping, matches};

use vita_core::SiteId;

use crate::algorithm::utils::{AdjacencyList, FxHashMap};
use crate::{BondId, HasBonds};

/// A molecule's sites and bonds in order, with the adjacency over their indices —
/// the dense `0..n` form the matching engine consumes.
struct Indexed {
    sites: Vec<SiteId>,
    bonds: Vec<BondId>,
    adjacency: AdjacencyList,
}

/// Indexes a molecule into contiguous vertices for the matching engine.
///
/// Site `sites[i]` becomes vertex `i` and bond `bonds[e]` becomes edge `e`, the
/// indices [`Indexed::adjacency`] is built over.
fn index<M: HasBonds>(mol: &M) -> Indexed {
    let sites: Vec<SiteId> = mol.sites().collect();
    let position: FxHashMap<SiteId, usize> =
        sites.iter().enumerate().map(|(i, &s)| (s, i)).collect();
    let bonds: Vec<BondId> = mol.bonds().collect();
    let adjacency = AdjacencyList::build(
        sites.len(),
        bonds.iter().enumerate().map(|(edge, &bond)| {
            let (a, b) = mol.bond_endpoints(bond);
            (edge, position[&a], position[&b])
        }),
    );
    Indexed {
        sites,
        bonds,
        adjacency,
    }
}
