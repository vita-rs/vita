use std::collections::HashMap;

use vita_core::SiteId;

use crate::{BondId, HasBonds};

/// A `(neighbour, bond index)` adjacency list over the sites `0..site_count`.
pub(super) type Adjacency = Vec<Vec<(usize, usize)>>;

/// A molecule's sites and bonds in order, with the adjacency over their indices —
/// the form the matching engine consumes.
pub(super) struct Indexed {
    pub(super) sites: Vec<SiteId>,
    pub(super) bonds: Vec<BondId>,
    pub(super) adjacency: Adjacency,
}

/// Indexes a molecule for the matching engine.
pub(super) fn index<M: HasBonds>(mol: &M) -> Indexed {
    let sites: Vec<SiteId> = mol.sites().collect();
    let position: HashMap<SiteId, usize> = sites.iter().enumerate().map(|(i, &s)| (s, i)).collect();
    let bonds: Vec<BondId> = mol.bonds().collect();
    let mut adjacency: Adjacency = vec![Vec::new(); sites.len()];
    for (edge, &bond) in bonds.iter().enumerate() {
        let (a, b) = mol.bond_endpoints(bond);
        adjacency[position[&a]].push((position[&b], edge));
        adjacency[position[&b]].push((position[&a], edge));
    }
    Indexed {
        sites,
        bonds,
        adjacency,
    }
}
