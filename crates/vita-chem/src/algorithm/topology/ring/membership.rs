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

#[cfg(test)]
mod tests {
    use super::*;
    use vita_core::HasSites;

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn b(n: u32) -> BondId {
        BondId::new(n).unwrap()
    }

    struct Mol {
        sites: Vec<SiteId>,
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
    }

    impl HasSites for Mol {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }

    impl HasBonds for Mol {
        fn bonds(&self) -> impl Iterator<Item = BondId> + '_ {
            self.bonds.iter().copied()
        }

        fn bond_endpoints(&self, bond: BondId) -> (SiteId, SiteId) {
            let i = self.bonds.iter().position(|&x| x == bond).unwrap();
            self.endpoints[i]
        }
    }

    fn empty() -> Mol {
        Mol {
            sites: vec![],
            bonds: vec![],
            endpoints: vec![],
        }
    }

    fn single() -> Mol {
        Mol {
            sites: vec![s(1)],
            bonds: vec![],
            endpoints: vec![],
        }
    }

    fn chain() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            bonds: vec![b(1), b(2)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3))],
        }
    }

    fn triangle() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            bonds: vec![b(1), b(2), b(3)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3)), (s(1), s(3))],
        }
    }

    fn lollipop() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4)],
            bonds: vec![b(1), b(2), b(3), b(4)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3)), (s(1), s(3)), (s(1), s(4))],
        }
    }

    fn dumbbell() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6)],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6), b(7)],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(1), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(6)),
                (s(4), s(6)),
            ],
        }
    }

    fn two_triangles() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6)],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6)],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(1), s(3)),
                (s(4), s(5)),
                (s(5), s(6)),
                (s(4), s(6)),
            ],
        }
    }

    #[test]
    fn empty_molecule_is_acyclic() {
        assert!(membership(&empty()).is_acyclic());
    }

    #[test]
    fn single_site_is_acyclic() {
        assert!(membership(&single()).is_acyclic());
    }

    #[test]
    fn chain_is_acyclic() {
        assert!(membership(&chain()).is_acyclic());
    }

    #[test]
    fn triangle_is_not_acyclic() {
        assert!(!membership(&triangle()).is_acyclic());
    }

    #[test]
    fn chain_site_not_in_ring() {
        let m = membership(&chain());
        assert!(!m.site(s(1)));
        assert!(!m.site(s(2)));
        assert!(!m.site(s(3)));
    }

    #[test]
    fn chain_bond_not_in_ring() {
        let m = membership(&chain());
        assert!(!m.bond(b(1)));
        assert!(!m.bond(b(2)));
    }

    #[test]
    fn chain_has_no_ring_sites() {
        assert_eq!(membership(&chain()).sites().count(), 0);
    }

    #[test]
    fn triangle_every_site_in_ring() {
        let m = membership(&triangle());
        assert!([s(1), s(2), s(3)].iter().all(|&site| m.site(site)));
    }

    #[test]
    fn triangle_every_bond_in_ring() {
        let m = membership(&triangle());
        assert!([b(1), b(2), b(3)].iter().all(|&bond| m.bond(bond)));
    }

    #[test]
    fn triangle_ring_sites() {
        let mut sites: Vec<SiteId> = membership(&triangle()).sites().collect();
        sites.sort();
        assert_eq!(sites, vec![s(1), s(2), s(3)]);
    }

    #[test]
    fn triangle_ring_bonds() {
        let mut bonds: Vec<BondId> = membership(&triangle()).bonds().collect();
        bonds.sort();
        assert_eq!(bonds, vec![b(1), b(2), b(3)]);
    }

    #[test]
    fn lollipop_ring_sites() {
        let mut sites: Vec<SiteId> = membership(&lollipop()).sites().collect();
        sites.sort();
        assert_eq!(sites, vec![s(1), s(2), s(3)]);
    }

    #[test]
    fn lollipop_ring_bonds() {
        let mut bonds: Vec<BondId> = membership(&lollipop()).bonds().collect();
        bonds.sort();
        assert_eq!(bonds, vec![b(1), b(2), b(3)]);
    }

    #[test]
    fn lollipop_tail_site_not_in_ring() {
        assert!(!membership(&lollipop()).site(s(4)));
    }

    #[test]
    fn lollipop_tail_bond_not_in_ring() {
        assert!(!membership(&lollipop()).bond(b(4)));
    }

    #[test]
    fn dumbbell_all_sites_in_ring() {
        let m = membership(&dumbbell());
        assert!(
            [s(1), s(2), s(3), s(4), s(5), s(6)]
                .iter()
                .all(|&site| m.site(site))
        );
    }

    #[test]
    fn dumbbell_bridge_not_in_ring() {
        let m = membership(&dumbbell());
        assert!(!m.bond(b(4)));
        assert!(
            [b(1), b(2), b(3), b(5), b(6), b(7)]
                .iter()
                .all(|&bond| m.bond(bond))
        );
    }

    #[test]
    fn two_triangles_all_in_ring() {
        let m = membership(&two_triangles());
        assert_eq!(m.sites().count(), 6);
        assert_eq!(m.bonds().count(), 6);
    }

    #[test]
    fn unknown_site_not_in_ring() {
        assert!(!membership(&triangle()).site(s(99)));
    }

    #[test]
    fn unknown_bond_not_in_ring() {
        assert!(!membership(&triangle()).bond(b(99)));
    }
}
