use std::collections::VecDeque;

use vita_core::SiteId;

use crate::HasBonds;
use crate::algorithm::utils::FxHashSet;

/// The sites reachable from `start`, in breadth-first order.
///
/// Yields `start` first, then every site reachable from it through bonds,
/// nearer sites before farther: every site at distance `d` precedes every site
/// at distance `d + 1`. Sites equidistant from `start` follow the order the
/// search reaches them. Each site is yielded once; sites in other connected
/// components never appear. The search is lazy — each step expands one site —
/// so adapters such as [`take`](Iterator::take) halt it early.
///
/// # Complexity
///
/// Exhausting the traversal is O(V + E) time and O(V) auxiliary space, over the
/// `V` sites and `E` bonds reachable from `start`, assuming
/// [`neighbors`](HasBonds::neighbors) runs in O(degree).
pub fn bfs<M: HasBonds>(mol: &M, start: SiteId) -> impl Iterator<Item = SiteId> + '_ {
    Bfs::new(mol, start)
}

/// The iterator returned by [`bfs`].
struct Bfs<'a, M> {
    mol: &'a M,
    frontier: VecDeque<SiteId>,
    visited: FxHashSet<SiteId>,
}

impl<'a, M: HasBonds> Bfs<'a, M> {
    fn new(mol: &'a M, start: SiteId) -> Self {
        let mut visited = FxHashSet::default();
        visited.insert(start);
        let frontier = VecDeque::from([start]);
        Bfs {
            mol,
            frontier,
            visited,
        }
    }
}

impl<M: HasBonds> Iterator for Bfs<'_, M> {
    type Item = SiteId;

    fn next(&mut self) -> Option<SiteId> {
        let site = self.frontier.pop_front()?;
        let mol = self.mol;
        for neighbor in mol.neighbors(site) {
            if self.visited.insert(neighbor) {
                self.frontier.push_back(neighbor);
            }
        }
        Some(site)
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.frontier.len(), None)
    }
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

    fn star() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4)],
            bonds: vec![b(1), b(2), b(3)],
            endpoints: vec![(s(1), s(3)), (s(1), s(4)), (s(1), s(2))],
        }
    }

    fn tree() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5)],
            bonds: vec![b(1), b(2), b(3), b(4)],
            endpoints: vec![(s(1), s(2)), (s(1), s(3)), (s(2), s(4)), (s(2), s(5))],
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
    fn single_site_yields_only_itself() {
        assert_eq!(bfs(&single(), s(1)).collect::<Vec<_>>(), vec![s(1)]);
    }

    #[test]
    fn start_is_yielded_first() {
        assert_eq!(bfs(&chain(), s(2)).next(), Some(s(2)));
        assert_eq!(bfs(&chain(), s(3)).next(), Some(s(3)));
    }

    #[test]
    fn reaches_every_site_in_a_connected_molecule() {
        let mut reached: Vec<SiteId> = bfs(&tree(), s(1)).collect();
        reached.sort_unstable();
        assert_eq!(reached, vec![s(1), s(2), s(3), s(4), s(5)]);
    }

    #[test]
    fn does_not_reach_other_components() {
        let reached: Vec<SiteId> = bfs(&disconnected(), s(1)).collect();
        assert!(!reached.contains(&s(3)));
    }

    #[test]
    fn isolated_start_yields_only_itself() {
        assert_eq!(bfs(&disconnected(), s(3)).collect::<Vec<_>>(), vec![s(3)]);
    }

    #[test]
    fn yields_nearer_sites_before_farther() {
        let distance = |site: SiteId| match site.get() {
            1 => 0,
            2 | 3 => 1,
            _ => 2,
        };
        let order: Vec<SiteId> = bfs(&tree(), s(1)).collect();
        let distances: Vec<u32> = order.iter().map(|&site| distance(site)).collect();
        assert!(distances.windows(2).all(|w| w[0] <= w[1]));
    }

    #[test]
    fn visits_a_cycle_without_revisiting_a_site() {
        let order: Vec<SiteId> = bfs(&triangle(), s(1)).collect();
        let unique: FxHashSet<SiteId> = order.iter().copied().collect();
        assert_eq!(order.len(), unique.len());
    }

    #[test]
    fn orders_a_path_outward_from_an_endpoint() {
        assert_eq!(
            bfs(&chain(), s(1)).collect::<Vec<_>>(),
            vec![s(1), s(2), s(3)],
        );
    }

    #[test]
    fn orders_equidistant_sites_by_the_order_their_bonds_are_reached() {
        assert_eq!(
            bfs(&star(), s(1)).collect::<Vec<_>>(),
            vec![s(1), s(3), s(4), s(2)],
        );
    }

    #[test]
    fn counts_the_sites_in_the_reachable_component() {
        assert_eq!(bfs(&chain(), s(1)).count(), 3);
        assert_eq!(bfs(&disconnected(), s(1)).count(), 2);
    }

    #[test]
    fn take_yields_only_the_requested_prefix() {
        assert_eq!(
            bfs(&chain(), s(1)).take(2).collect::<Vec<_>>(),
            vec![s(1), s(2)],
        );
    }

    #[test]
    fn size_hint_lower_bound_never_exceeds_the_sites_remaining() {
        let mol = tree();
        let mut iter = bfs(&mol, s(1));
        let mut remaining = 5;
        loop {
            let (lower, upper) = iter.size_hint();
            assert!(
                lower <= remaining,
                "lower {lower} exceeded remaining {remaining}"
            );
            assert_eq!(upper, None);
            if iter.next().is_none() {
                break;
            }
            remaining -= 1;
        }
        assert_eq!(iter.size_hint(), (0, None));
    }

    #[test]
    fn exhausted_iterator_keeps_yielding_none() {
        let mol = single();
        let mut iter = bfs(&mol, s(1));
        assert_eq!(iter.next(), Some(s(1)));
        assert_eq!(iter.next(), None);
        assert_eq!(iter.next(), None);
    }
}
