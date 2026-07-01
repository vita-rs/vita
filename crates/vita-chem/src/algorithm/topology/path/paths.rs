use std::collections::VecDeque;

use vita_core::SiteId;

use crate::HasBonds;
use crate::algorithm::utils::FxHashMap;

/// Every shortest path from `start` to `end`.
///
/// A shortest path uses the fewest bonds; when several tie, all are returned,
/// each listing its sites in order from `start` to `end` inclusive. The result
/// is empty when `end` is unreachable from `start`, and is the single path
/// `[start]` when `start == end`. Paths are returned in ascending order.
///
/// # Complexity
///
/// O(V + E + P) time and space, where `V` and `E` are the sites and bonds
/// reachable from `start` and `P` is the total length of the paths returned,
/// assuming [`neighbors`](HasBonds::neighbors) runs in O(degree).
pub fn paths<M: HasBonds>(mol: &M, start: SiteId, end: SiteId) -> Vec<Vec<SiteId>> {
    if start == end {
        return vec![vec![start]];
    }

    let mut dist: FxHashMap<SiteId, usize> = FxHashMap::default();
    let mut parents: FxHashMap<SiteId, Vec<SiteId>> = FxHashMap::default();
    let mut queue: VecDeque<SiteId> = VecDeque::new();
    dist.insert(start, 0);
    queue.push_back(start);

    while let Some(site) = queue.pop_front() {
        let d = dist[&site];
        for neighbor in mol.neighbors(site) {
            match dist.get(&neighbor).copied() {
                None => {
                    dist.insert(neighbor, d + 1);
                    parents.insert(neighbor, vec![site]);
                    queue.push_back(neighbor);
                }
                Some(nd) if nd == d + 1 => parents.get_mut(&neighbor).unwrap().push(site),
                Some(_) => {}
            }
        }
    }

    if !dist.contains_key(&end) {
        return Vec::new();
    }

    let mut result: Vec<Vec<SiteId>> = Vec::new();
    let mut path: Vec<SiteId> = Vec::new();
    collect(&parents, start, end, &mut path, &mut result);
    result.sort_unstable();
    result
}

/// Follows the predecessor lists `parents` from `site` back to `start`,
/// appending each complete path to `result`. `path` holds the sites chosen so
/// far, in reverse; a completed path is reversed before it is stored.
fn collect(
    parents: &FxHashMap<SiteId, Vec<SiteId>>,
    start: SiteId,
    site: SiteId,
    path: &mut Vec<SiteId>,
    result: &mut Vec<Vec<SiteId>>,
) {
    path.push(site);
    if site == start {
        let mut full = path.clone();
        full.reverse();
        result.push(full);
    } else {
        for &parent in &parents[&site] {
            collect(parents, start, parent, path, result);
        }
    }
    path.pop();
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
            endpoints: vec![(s(1), s(2)), (s(2), s(3)), (s(3), s(4)), (s(1), s(4))],
        }
    }

    fn disconnected() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4)],
            bonds: vec![b(1), b(2)],
            endpoints: vec![(s(1), s(2)), (s(3), s(4))],
        }
    }

    #[test]
    fn start_equal_to_end_yields_the_lone_site() {
        assert_eq!(paths(&chain(), s(2), s(2)), vec![vec![s(2)]]);
    }

    #[test]
    fn a_tree_yields_its_unique_shortest_path() {
        assert_eq!(paths(&chain(), s(1), s(3)), vec![vec![s(1), s(2), s(3)]]);
    }

    #[test]
    fn disconnected_sites_yield_no_path() {
        assert!(paths(&disconnected(), s(1), s(3)).is_empty());
    }

    #[test]
    fn excludes_paths_longer_than_the_shortest() {
        assert_eq!(paths(&triangle(), s(1), s(3)), vec![vec![s(1), s(3)]]);
    }

    #[test]
    fn a_cycle_yields_every_tied_shortest_path() {
        let result = paths(&square(), s(1), s(3));
        assert_eq!(result.len(), 2);
        assert!(result.contains(&vec![s(1), s(2), s(3)]));
        assert!(result.contains(&vec![s(1), s(4), s(3)]));
    }

    #[test]
    fn all_paths_share_the_minimal_length() {
        let result = paths(&square(), s(1), s(3));
        assert!(result.iter().all(|path| path.len() == 3));
    }

    #[test]
    fn every_path_runs_from_start_to_end() {
        for path in paths(&square(), s(3), s(1)) {
            assert_eq!(path.first(), Some(&s(3)));
            assert_eq!(path.last(), Some(&s(1)));
        }
    }

    #[test]
    fn consecutive_sites_in_each_path_are_bonded() {
        let mol = square();
        for path in paths(&mol, s(1), s(3)) {
            for window in path.windows(2) {
                assert!(mol.bond_between(window[0], window[1]).is_some());
            }
        }
    }

    #[test]
    fn paths_are_listed_in_ascending_order() {
        let result = paths(&square(), s(1), s(3));
        let mut sorted = result.clone();
        sorted.sort();
        assert_eq!(result, sorted);
    }

    #[test]
    fn output_is_independent_of_input_order() {
        let reordered = Mol {
            sites: vec![s(4), s(3), s(2), s(1)],
            bonds: vec![b(4), b(3), b(2), b(1)],
            endpoints: vec![(s(1), s(4)), (s(3), s(4)), (s(2), s(3)), (s(1), s(2))],
        };
        assert_eq!(paths(&square(), s(1), s(3)), paths(&reordered, s(1), s(3)));
    }
}
