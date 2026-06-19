use std::collections::{HashSet, VecDeque};

use vita_core::SiteId;

use crate::HasBonds;

/// Breadth-first traversal order from `start`.
///
/// Yields every site reachable from `start` in breadth-first (level) order,
/// visiting sites closer to `start` before sites farther away. Sites in
/// disconnected components are never yielded. `start` is always the first
/// site yielded.
///
/// # Complexity
///
/// O(V + E) time, O(V) auxiliary space.
pub fn bfs<M: HasBonds>(mol: &M, start: SiteId) -> impl Iterator<Item = SiteId> + '_ {
    BfsIter::new(mol, start)
}

struct BfsIter<'a, M> {
    mol: &'a M,
    queue: VecDeque<SiteId>,
    visited: HashSet<SiteId>,
}

impl<'a, M: HasBonds> BfsIter<'a, M> {
    fn new(mol: &'a M, start: SiteId) -> Self {
        let mut visited = HashSet::new();
        visited.insert(start);
        let mut queue = VecDeque::new();
        queue.push_back(start);
        Self {
            mol,
            queue,
            visited,
        }
    }
}

impl<'a, M: HasBonds> Iterator for BfsIter<'a, M> {
    type Item = SiteId;

    fn next(&mut self) -> Option<SiteId> {
        let site = self.queue.pop_front()?;
        let mol = self.mol;
        for nb in mol.neighbors(site) {
            if self.visited.insert(nb) {
                self.queue.push_back(nb);
            }
        }
        Some(site)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.queue.len(), None)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BondId;
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
    fn start_is_first() {
        let mol = chain();
        assert_eq!(bfs(&mol, s(1)).next(), Some(s(1)));
        assert_eq!(bfs(&mol, s(2)).next(), Some(s(2)));
        assert_eq!(bfs(&mol, s(3)).next(), Some(s(3)));
    }

    #[test]
    fn chain_from_terminal() {
        assert_eq!(
            bfs(&chain(), s(1)).collect::<Vec<_>>(),
            vec![s(1), s(2), s(3)],
        );
    }

    #[test]
    fn chain_from_middle() {
        assert_eq!(
            bfs(&chain(), s(2)).collect::<Vec<_>>(),
            vec![s(2), s(1), s(3)],
        );
    }

    #[test]
    fn chain_from_other_terminal() {
        assert_eq!(
            bfs(&chain(), s(3)).collect::<Vec<_>>(),
            vec![s(3), s(2), s(1)],
        );
    }

    #[test]
    fn single_site() {
        let mol = Mol {
            sites: vec![s(1)],
            bonds: vec![],
            endpoints: vec![],
        };
        assert_eq!(bfs(&mol, s(1)).collect::<Vec<_>>(), vec![s(1)]);
    }

    #[test]
    fn disconnected_stays_in_component() {
        let mol = two_components();
        let from_1: Vec<_> = bfs(&mol, s(1)).collect();
        assert_eq!(from_1.len(), 2);
        assert!(!from_1.contains(&s(3)));
    }

    #[test]
    fn disconnected_isolated_site_yields_only_itself() {
        assert_eq!(bfs(&two_components(), s(3)).collect::<Vec<_>>(), vec![s(3)],);
    }

    #[test]
    fn ring_no_duplicates() {
        let order: Vec<_> = bfs(&triangle(), s(1)).collect();
        let unique: HashSet<_> = order.iter().copied().collect();
        assert_eq!(unique.len(), order.len(), "BFS produced duplicate sites");
    }

    #[test]
    fn ring_visits_all_sites() {
        let mut order: Vec<_> = bfs(&triangle(), s(1)).collect();
        order.sort();
        assert_eq!(order, vec![s(1), s(2), s(3)]);
    }

    #[test]
    fn ring_level_order() {
        let order: Vec<_> = bfs(&triangle(), s(1)).collect();
        assert_eq!(order[0], s(1));
        let mut rest = [order[1], order[2]];
        rest.sort();
        assert_eq!(rest, [s(2), s(3)]);
    }

    #[test]
    fn star_from_center_visits_leaves_at_depth_one() {
        let order: Vec<_> = bfs(&star(), s(1)).collect();
        assert_eq!(order[0], s(1));
        let mut leaves: Vec<_> = order[1..].to_vec();
        leaves.sort();
        assert_eq!(leaves, vec![s(2), s(3), s(4)]);
    }

    #[test]
    fn star_from_leaf_visits_center_then_other_leaves() {
        let order: Vec<_> = bfs(&star(), s(2)).collect();
        assert_eq!(order[0], s(2));
        assert_eq!(order[1], s(1));
        let mut rest: Vec<_> = order[2..].to_vec();
        rest.sort();
        assert_eq!(rest, vec![s(3), s(4)]);
    }

    #[test]
    fn cyclobutane_level_order_from_vertex() {
        assert_eq!(
            bfs(&cyclobutane(), s(1)).collect::<Vec<_>>(),
            vec![s(1), s(2), s(3), s(4)],
        );
    }

    #[test]
    fn cyclobutane_opposite_vertex_is_last() {
        let order: Vec<_> = bfs(&cyclobutane(), s(1)).collect();
        let pos = |x| order.iter().position(|&v| v == x).unwrap();
        assert!(pos(s(2)) < pos(s(4)));
        assert!(pos(s(3)) < pos(s(4)));
    }

    #[test]
    fn pentane_from_middle_alternates_outward() {
        assert_eq!(
            bfs(&pentane(), s(3)).collect::<Vec<_>>(),
            vec![s(3), s(2), s(4), s(1), s(5)],
        );
    }

    #[test]
    fn count_equals_component_size() {
        assert_eq!(bfs(&chain(), s(1)).count(), 3);
        assert_eq!(bfs(&chain(), s(2)).count(), 3);
        assert_eq!(bfs(&two_components(), s(1)).count(), 2);
        assert_eq!(bfs(&two_components(), s(3)).count(), 1);
    }

    #[test]
    fn size_hint_at_start_has_lower_bound_of_one() {
        let (lower, upper) = bfs(&chain(), s(1)).size_hint();
        assert!(lower >= 1, "frontier must hold at least the start site");
        assert_eq!(upper, None);
    }

    #[test]
    fn size_hint_lower_bound_never_exceeds_remaining() {
        let mol = chain();
        let mut iter = bfs(&mol, s(1));
        let mut remaining = 3usize;
        loop {
            let (lower, _) = iter.size_hint();
            assert!(
                lower <= remaining,
                "lower={lower} exceeded remaining={remaining}",
            );
            match iter.next() {
                Some(_) => remaining -= 1,
                None => break,
            }
        }
    }

    #[test]
    fn size_hint_exhausted_lower_bound_is_zero() {
        let mol = Mol {
            sites: vec![s(1)],
            bonds: vec![],
            endpoints: vec![],
        };
        let mut iter = bfs(&mol, s(1));
        iter.next();
        assert_eq!(iter.size_hint().0, 0);
    }

    #[test]
    fn iterator_adapters_work() {
        let id_sum: u32 = bfs(&chain(), s(1)).map(|site| site.get()).sum();
        assert_eq!(id_sum, 1 + 2 + 3);
    }

    #[test]
    fn take_stops_early() {
        let first_two: Vec<_> = bfs(&pentane(), s(1)).take(2).collect();
        assert_eq!(first_two, vec![s(1), s(2)]);
    }
}
