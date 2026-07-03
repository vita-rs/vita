use vita_core::SiteId;

use crate::algorithm::utils::{AdjacencyList, FxHashMap};
use crate::{BondId, HasBonds};

/// Which of a molecule's sites and bonds lie in a ring.
///
/// A bond is a ring bond exactly when it is not a bridge — when its removal
/// would leave the number of connected components unchanged. A site is a ring
/// site exactly when it is incident to a ring bond. The ring sites and the ring
/// bonds are each held in ascending order.
///
/// Obtain via [`membership`].
pub struct RingMembership {
    sites: Vec<SiteId>,
    bonds: Vec<BondId>,
}

impl RingMembership {
    /// Returns `true` if `site` lies in a ring.
    ///
    /// Returns `false` if `site` is absent from the molecule or lies in no ring.
    pub fn contains_site(&self, site: SiteId) -> bool {
        self.sites.binary_search(&site).is_ok()
    }

    /// Returns `true` if `bond` lies in a ring.
    ///
    /// Returns `false` if `bond` is absent from the molecule or is a bridge.
    pub fn contains_bond(&self, bond: BondId) -> bool {
        self.bonds.binary_search(&bond).is_ok()
    }

    /// The ring sites, in ascending order.
    pub fn sites(&self) -> &[SiteId] {
        &self.sites
    }

    /// The ring bonds, in ascending order.
    pub fn bonds(&self) -> &[BondId] {
        &self.bonds
    }

    /// Returns `true` if the molecule has no rings.
    pub fn is_acyclic(&self) -> bool {
        self.bonds.is_empty()
    }

    /// Builds a membership from its ring sites and ring bonds.
    pub(super) fn from_sets(
        sites: impl IntoIterator<Item = SiteId>,
        bonds: impl IntoIterator<Item = BondId>,
    ) -> Self {
        let mut sites: Vec<SiteId> = sites.into_iter().collect();
        let mut bonds: Vec<BondId> = bonds.into_iter().collect();
        sites.sort_unstable();
        sites.dedup();
        bonds.sort_unstable();
        bonds.dedup();
        RingMembership { sites, bonds }
    }
}

/// The ring membership of a molecule's sites and bonds.
///
/// A bond is a ring bond exactly when it is not a bridge; a site is a ring site
/// exactly when it is incident to a ring bond. Bridges are found by Tarjan's
/// low-link test over an explicit stack, so recursion depth never bounds the
/// traversal.
///
/// # Complexity
///
/// O(V · log V + E · log E) time and O(V + E) space, over the molecule's `V`
/// sites and `E` bonds.
pub fn membership<M: HasBonds>(mol: &M) -> RingMembership {
    let sites: Vec<SiteId> = mol.sites().collect();
    let bonds: Vec<BondId> = mol.bonds().collect();
    let n = sites.len();
    let m = bonds.len();

    let index: FxHashMap<SiteId, usize> = sites.iter().enumerate().map(|(i, &s)| (s, i)).collect();
    let adjacency = AdjacencyList::build(
        n,
        bonds.iter().enumerate().map(|(e, &bond)| {
            let (a, b) = mol.bond_endpoints(bond);
            (e, index[&a], index[&b])
        }),
    );

    let is_bridge = bridges(n, m, &adjacency);

    let mut ring_bonds: Vec<BondId> = (0..m)
        .filter(|&e| !is_bridge[e])
        .map(|e| bonds[e])
        .collect();
    ring_bonds.sort_unstable();

    let mut ring_sites: Vec<SiteId> = (0..n)
        .filter(|&u| adjacency.neighbors(u).iter().any(|&(e, _)| !is_bridge[e]))
        .map(|u| sites[u])
        .collect();
    ring_sites.sort_unstable();

    RingMembership {
        sites: ring_sites,
        bonds: ring_bonds,
    }
}

/// Flags each edge that is a bridge: an edge on no cycle, whose removal would
/// disconnect its component.
///
/// Runs Tarjan's low-link test over an explicit stack, so recursion depth never
/// bounds the traversal. The flag at index `e` refers to the edge whose
/// identifier is `e` in `adjacency`.
fn bridges(n: usize, m: usize, adjacency: &AdjacencyList) -> Vec<bool> {
    let mut disc = vec![usize::MAX; n];
    let mut low = vec![0usize; n];
    let mut is_bridge = vec![false; m];
    let mut timer = 0;
    let mut stack: Vec<(usize, usize, usize)> = Vec::new();

    for start in 0..n {
        if disc[start] != usize::MAX {
            continue;
        }
        disc[start] = timer;
        low[start] = timer;
        timer += 1;
        stack.push((start, usize::MAX, 0));

        while let Some(&(u, parent_edge, cursor)) = stack.last() {
            if cursor < adjacency.neighbors(u).len() {
                stack.last_mut().unwrap().2 += 1;
                let (edge, v) = adjacency.neighbors(u)[cursor];
                if edge == parent_edge {
                    continue;
                }
                if disc[v] == usize::MAX {
                    disc[v] = timer;
                    low[v] = timer;
                    timer += 1;
                    stack.push((v, edge, 0));
                } else {
                    low[u] = low[u].min(disc[v]);
                }
            } else {
                stack.pop();
                if let Some(&(parent, _, _)) = stack.last() {
                    low[parent] = low[parent].min(low[u]);
                    if low[u] > disc[parent] {
                        is_bridge[parent_edge] = true;
                    }
                }
            }
        }
    }

    is_bridge
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::HasSites;

    use crate::BondId;

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
    fn every_bond_of_a_cycle_is_a_ring_bond() {
        let m = membership(&triangle());
        assert!([b(1), b(2), b(3)].iter().all(|&bond| m.contains_bond(bond)));
    }

    #[test]
    fn every_site_of_a_cycle_is_a_ring_site() {
        let m = membership(&triangle());
        assert!([s(1), s(2), s(3)].iter().all(|&site| m.contains_site(site)));
    }

    #[test]
    fn a_molecule_with_a_cycle_is_not_acyclic() {
        assert!(!membership(&triangle()).is_acyclic());
    }

    #[test]
    fn a_tree_is_acyclic() {
        assert!(membership(&chain()).is_acyclic());
    }

    #[test]
    fn a_bridge_is_not_a_ring_bond() {
        assert!(!membership(&lollipop()).contains_bond(b(4)));
    }

    #[test]
    fn a_site_incident_only_to_bridges_is_not_a_ring_site() {
        assert!(!membership(&lollipop()).contains_site(s(4)));
    }

    #[test]
    fn an_unknown_site_is_not_a_ring_site() {
        assert!(!membership(&triangle()).contains_site(s(99)));
    }

    #[test]
    fn an_unknown_bond_is_not_a_ring_bond() {
        assert!(!membership(&triangle()).contains_bond(b(99)));
    }

    #[test]
    fn a_site_shared_by_a_ring_and_a_bridge_is_a_ring_site() {
        assert!(membership(&dumbbell()).contains_site(s(3)));
    }

    #[test]
    fn a_bridge_joining_two_rings_is_not_a_ring_bond() {
        assert!(!membership(&dumbbell()).contains_bond(b(4)));
    }

    #[test]
    fn ring_sites_are_listed_in_ascending_order() {
        let m = membership(&dumbbell());
        assert_eq!(m.sites(), &[s(1), s(2), s(3), s(4), s(5), s(6)]);
    }

    #[test]
    fn ring_bonds_exclude_the_bridge_and_are_listed_in_ascending_order() {
        let m = membership(&dumbbell());
        assert_eq!(m.bonds(), &[b(1), b(2), b(3), b(5), b(6), b(7)]);
    }

    #[test]
    fn disjoint_cycles_are_each_perceived() {
        let m = membership(&two_triangles());
        assert_eq!(m.sites().len(), 6);
        assert_eq!(m.bonds().len(), 6);
    }

    #[test]
    fn membership_is_independent_of_input_order() {
        let shuffled = Mol {
            sites: vec![s(6), s(4), s(2), s(5), s(1), s(3)],
            bonds: vec![b(7), b(4), b(1), b(6), b(2), b(5), b(3)],
            endpoints: vec![
                (s(4), s(6)),
                (s(3), s(4)),
                (s(1), s(2)),
                (s(5), s(6)),
                (s(2), s(3)),
                (s(4), s(5)),
                (s(1), s(3)),
            ],
        };
        let parts = |mol: &Mol| {
            let m = membership(mol);
            (m.sites().to_vec(), m.bonds().to_vec())
        };
        assert_eq!(parts(&dumbbell()), parts(&shuffled));
    }
}
