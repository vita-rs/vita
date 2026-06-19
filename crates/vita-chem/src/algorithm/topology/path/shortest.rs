use std::collections::{HashMap, VecDeque};

use vita_core::SiteId;

use crate::HasBonds;

/// Shortest topological path between two sites, inclusive of both endpoints.
///
/// Uses BFS, guaranteeing the minimum number of bonds. The returned sequence
/// starts at `start` and ends at `end`. Returns `None` when `start` and `end`
/// lie in different connected components. Returns `Some` containing only
/// `start` when `start == end`.
///
/// # Complexity
///
/// O(V + E) time, O(V) auxiliary space.
pub fn path<M: HasBonds>(mol: &M, start: SiteId, end: SiteId) -> Option<Vec<SiteId>> {
    if start == end {
        return Some(vec![start]);
    }

    let mut parent: HashMap<SiteId, SiteId> = HashMap::new();
    let mut queue: VecDeque<SiteId> = VecDeque::new();

    parent.insert(start, start);
    queue.push_back(start);

    while let Some(site) = queue.pop_front() {
        for nb in mol.neighbors(site) {
            if parent.contains_key(&nb) {
                continue;
            }
            parent.insert(nb, site);
            if nb == end {
                let mut result = vec![end];
                let mut cur = end;
                while cur != start {
                    cur = parent[&cur];
                    result.push(cur);
                }
                result.reverse();
                return Some(result);
            }
            queue.push_back(nb);
        }
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BondId;
    use std::collections::HashSet;
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

    fn star() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4)],
            bonds: vec![b(1), b(2), b(3)],
            endpoints: vec![(s(1), s(2)), (s(1), s(3)), (s(1), s(4))],
        }
    }

    fn cyclobutane() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4)],
            bonds: vec![b(1), b(2), b(3), b(4)],
            endpoints: vec![(s(1), s(2)), (s(1), s(3)), (s(2), s(4)), (s(3), s(4))],
        }
    }

    fn pentane() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5)],
            bonds: vec![b(1), b(2), b(3), b(4)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3)), (s(3), s(4)), (s(4), s(5))],
        }
    }

    fn two_components() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            bonds: vec![b(1)],
            endpoints: vec![(s(1), s(2))],
        }
    }

    #[test]
    fn same_start_and_end() {
        assert_eq!(path(&chain(), s(1), s(1)), Some(vec![s(1)]));
        assert_eq!(path(&chain(), s(2), s(2)), Some(vec![s(2)]));
        assert_eq!(path(&chain(), s(3), s(3)), Some(vec![s(3)]));
    }

    #[test]
    fn single_site_to_itself() {
        let mol = Mol {
            sites: vec![s(1)],
            bonds: vec![],
            endpoints: vec![],
        };
        assert_eq!(path(&mol, s(1), s(1)), Some(vec![s(1)]));
    }

    #[test]
    fn adjacent_sites() {
        assert_eq!(path(&chain(), s(1), s(2)), Some(vec![s(1), s(2)]));
    }

    #[test]
    fn adjacent_sites_reversed() {
        assert_eq!(path(&chain(), s(2), s(1)), Some(vec![s(2), s(1)]));
    }

    #[test]
    fn chain_two_hops() {
        assert_eq!(path(&chain(), s(1), s(3)), Some(vec![s(1), s(2), s(3)]),);
    }

    #[test]
    fn chain_two_hops_reversed() {
        assert_eq!(path(&chain(), s(3), s(1)), Some(vec![s(3), s(2), s(1)]),);
    }

    #[test]
    fn disconnected_returns_none() {
        assert_eq!(path(&two_components(), s(1), s(3)), None);
        assert_eq!(path(&two_components(), s(3), s(1)), None);
    }

    #[test]
    fn isolated_site_returns_none() {
        assert_eq!(path(&two_components(), s(3), s(2)), None);
    }

    #[test]
    fn ring_uses_direct_bond() {
        assert_eq!(path(&triangle(), s(1), s(3)), Some(vec![s(1), s(3)]),);
    }

    #[test]
    fn ring_path_length_is_one_hop() {
        assert_eq!(path(&triangle(), s(1), s(3)).unwrap().len(), 2);
    }

    #[test]
    fn star_leaf_to_leaf_through_center() {
        assert_eq!(path(&star(), s(2), s(3)), Some(vec![s(2), s(1), s(3)]),);
    }

    #[test]
    fn star_leaf_to_center() {
        assert_eq!(path(&star(), s(4), s(1)), Some(vec![s(4), s(1)]));
    }

    #[test]
    fn cyclobutane_shortest_path() {
        assert_eq!(
            path(&cyclobutane(), s(1), s(4)),
            Some(vec![s(1), s(2), s(4)]),
        );
    }

    #[test]
    fn cyclobutane_path_length_is_minimal() {
        assert_eq!(path(&cyclobutane(), s(1), s(4)).unwrap().len(), 3);
    }

    #[test]
    fn pentane_full_path() {
        assert_eq!(
            path(&pentane(), s(1), s(5)),
            Some(vec![s(1), s(2), s(3), s(4), s(5)]),
        );
    }

    #[test]
    fn pentane_from_interior() {
        assert_eq!(path(&pentane(), s(3), s(5)), Some(vec![s(3), s(4), s(5)]),);
    }

    #[test]
    fn first_element_is_start() {
        assert_eq!(path(&chain(), s(1), s(3)).unwrap().first(), Some(&s(1)));
        assert_eq!(path(&star(), s(2), s(4)).unwrap().first(), Some(&s(2)));
        assert_eq!(path(&pentane(), s(3), s(5)).unwrap().first(), Some(&s(3)));
    }

    #[test]
    fn last_element_is_end() {
        assert_eq!(path(&chain(), s(1), s(3)).unwrap().last(), Some(&s(3)));
        assert_eq!(path(&star(), s(2), s(4)).unwrap().last(), Some(&s(4)));
        assert_eq!(path(&pentane(), s(3), s(5)).unwrap().last(), Some(&s(5)));
    }

    #[test]
    fn consecutive_sites_are_bonded() {
        let mol = cyclobutane();
        let p = path(&mol, s(1), s(4)).unwrap();
        for w in p.windows(2) {
            assert!(mol.bond_between(w[0], w[1]).is_some());
        }
    }

    #[test]
    fn path_contains_no_duplicate_sites() {
        let p = path(&cyclobutane(), s(1), s(4)).unwrap();
        let unique: HashSet<_> = p.iter().copied().collect();
        assert_eq!(unique.len(), p.len());
    }
}
