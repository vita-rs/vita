use std::collections::{VecDeque, hash_map::Entry};

use vita_core::SiteId;

use crate::HasBonds;
use crate::algorithm::utils::FxHashMap;

/// A shortest topological path from `start` to `end`, inclusive of both endpoints.
///
/// Returns the sites along a path with the fewest bonds, beginning at `start`
/// and ending at `end`; consecutive sites are bonded. When several shortest
/// paths exist, the one the breadth-first search reaches first is chosen. The
/// result is `Some([start])` when `start == end`, and `None` when `end` is
/// unreachable from `start`.
///
/// # Complexity
///
/// O(V + E) time and O(V) auxiliary space in the worst case, over the `V` sites
/// and `E` bonds reachable from `start`, assuming
/// [`neighbors`](HasBonds::neighbors) runs in O(degree).
pub fn path<M: HasBonds>(mol: &M, start: SiteId, end: SiteId) -> Option<Vec<SiteId>> {
    if start == end {
        return Some(vec![start]);
    }

    let mut parent: FxHashMap<SiteId, SiteId> = FxHashMap::default();
    parent.insert(start, start);
    let mut frontier: VecDeque<SiteId> = VecDeque::from([start]);

    while let Some(site) = frontier.pop_front() {
        for neighbor in mol.neighbors(site) {
            if let Entry::Vacant(slot) = parent.entry(neighbor) {
                slot.insert(site);
                if neighbor == end {
                    let mut route = vec![end];
                    let mut at = end;
                    while at != start {
                        at = parent[&at];
                        route.push(at);
                    }
                    route.reverse();
                    return Some(route);
                }
                frontier.push_back(neighbor);
            }
        }
    }

    None
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

    fn square() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4)],
            bonds: vec![b(1), b(2), b(3), b(4)],
            endpoints: vec![(s(1), s(4)), (s(1), s(2)), (s(2), s(3)), (s(3), s(4))],
        }
    }

    fn disconnected() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            bonds: vec![b(1)],
            endpoints: vec![(s(1), s(2))],
        }
    }

    #[test]
    fn start_equal_to_end_yields_a_single_site() {
        assert_eq!(path(&chain(), s(2), s(2)), Some(vec![s(2)]));
    }

    #[test]
    fn connects_two_sites_ordered_from_start_to_end() {
        assert_eq!(path(&chain(), s(1), s(3)), Some(vec![s(1), s(2), s(3)]));
        assert_eq!(path(&chain(), s(3), s(1)), Some(vec![s(3), s(2), s(1)]));
    }

    #[test]
    fn sites_in_different_components_have_no_path() {
        assert_eq!(path(&disconnected(), s(1), s(3)), None);
        assert_eq!(path(&disconnected(), s(3), s(1)), None);
    }

    #[test]
    fn path_to_an_absent_site_is_none() {
        assert_eq!(path(&chain(), s(1), s(99)), None);
    }

    #[test]
    fn adjacent_sites_yield_a_two_site_path() {
        assert_eq!(path(&chain(), s(1), s(2)), Some(vec![s(1), s(2)]));
    }

    #[test]
    fn takes_the_direct_bond_over_a_longer_route() {
        assert_eq!(path(&triangle(), s(1), s(3)), Some(vec![s(1), s(3)]));
    }

    #[test]
    fn chooses_among_equal_length_paths_by_search_order() {
        assert_eq!(path(&square(), s(1), s(3)), Some(vec![s(1), s(4), s(3)]));
    }

    #[test]
    fn path_sites_are_consecutively_bonded() {
        let mol = square();
        let route = path(&mol, s(1), s(3)).unwrap();
        for pair in route.windows(2) {
            assert!(mol.bond_between(pair[0], pair[1]).is_some());
        }
    }

    #[test]
    fn path_visits_no_site_twice() {
        let route = path(&square(), s(1), s(3)).unwrap();
        let mut unique = route.clone();
        unique.sort_unstable();
        unique.dedup();
        assert_eq!(unique.len(), route.len());
    }
}
