use std::collections::{HashMap, HashSet};

use vita_core::{HasSites, SiteId};

use crate::{BondId, HasBonds};

/// The ring membership of a molecule's sites and bonds.
///
/// A bond is a ring bond if and only if it is not a bridge; that is, removing
/// it would leave the molecule no less connected. A site is a ring site if and
/// only if it is incident to at least one ring bond.
///
/// Obtain via [`membership`].
pub struct RingMembership {
    sites: HashSet<SiteId>,
    bonds: HashSet<BondId>,
}

impl RingMembership {
    /// Returns `true` if `site` lies in at least one ring.
    ///
    /// Returns `false` if `site` is absent from the molecule or lies in no
    /// ring.
    pub fn site(&self, site: SiteId) -> bool {
        self.sites.contains(&site)
    }

    /// Returns `true` if `bond` lies in at least one ring.
    ///
    /// Returns `false` if `bond` is absent from the molecule or is a bridge.
    pub fn bond(&self, bond: BondId) -> bool {
        self.bonds.contains(&bond)
    }

    /// Iterates the sites that lie in at least one ring.
    pub fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
        self.sites.iter().copied()
    }

    /// Iterates the bonds that lie in at least one ring.
    pub fn bonds(&self) -> impl Iterator<Item = BondId> + '_ {
        self.bonds.iter().copied()
    }

    /// Returns `true` if the molecule contains no rings.
    pub fn is_acyclic(&self) -> bool {
        self.bonds.is_empty()
    }
}

/// Ring membership of every site and bond.
///
/// A bond is a ring bond exactly when it is not a bridge; a site is a ring site
/// exactly when it is incident to a ring bond. Computed with a single iterative
/// depth-first traversal using Tarjan's low-link bridge test, so no stack space
/// is consumed in proportion to traversal depth.
///
/// # Complexity
///
/// O(V + E) time and space.
pub fn membership<M: HasBonds + HasSites>(mol: &M) -> RingMembership {
    let sites: Vec<SiteId> = mol.sites().collect();
    let n = sites.len();

    let mut ring_sites: HashSet<SiteId> = HashSet::new();
    let mut ring_bonds: HashSet<BondId> = HashSet::new();

    if n == 0 {
        return RingMembership {
            sites: ring_sites,
            bonds: ring_bonds,
        };
    }

    let site_pos: HashMap<SiteId, usize> = sites.iter().enumerate().map(|(i, &s)| (s, i)).collect();

    let mut adj: Vec<Vec<(BondId, usize)>> = vec![vec![]; n];
    for bond in mol.bonds() {
        let (a, b) = mol.bond_endpoints(bond);
        let ai = site_pos[&a];
        let bi = site_pos[&b];
        adj[ai].push((bond, bi));
        adj[bi].push((bond, ai));
    }

    let mut disc = vec![usize::MAX; n];
    let mut low = vec![0usize; n];
    let mut timer = 0usize;

    let mut dfs_stack: Vec<(usize, Option<BondId>, usize)> = Vec::new();

    for start in 0..n {
        if disc[start] != usize::MAX {
            continue;
        }

        disc[start] = timer;
        low[start] = timer;
        timer += 1;
        dfs_stack.push((start, None, 0));

        while !dfs_stack.is_empty() {
            let (u, parent_bond, adj_pos) = *dfs_stack.last().unwrap();

            if adj_pos < adj[u].len() {
                let (bond, v) = adj[u][adj_pos];
                dfs_stack.last_mut().unwrap().2 += 1;

                if Some(bond) == parent_bond {
                    continue;
                }

                if disc[v] == usize::MAX {
                    disc[v] = timer;
                    low[v] = timer;
                    timer += 1;
                    dfs_stack.push((v, Some(bond), 0));
                } else if disc[v] < disc[u] {
                    if disc[v] < low[u] {
                        low[u] = disc[v];
                    }
                    ring_bonds.insert(bond);
                    ring_sites.insert(sites[u]);
                    ring_sites.insert(sites[v]);
                }
            } else {
                dfs_stack.pop();

                if let Some(&(pu, _, _)) = dfs_stack.last() {
                    if low[u] < low[pu] {
                        low[pu] = low[u];
                    }
                    if low[u] <= disc[pu] {
                        ring_bonds.insert(parent_bond.unwrap());
                        ring_sites.insert(sites[u]);
                        ring_sites.insert(sites[pu]);
                    }
                }
            }
        }
    }

    RingMembership {
        sites: ring_sites,
        bonds: ring_bonds,
    }
}
