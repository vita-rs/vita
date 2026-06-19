use std::collections::{HashMap, VecDeque};

use vita_core::SiteId;

use crate::HasBonds;

/// All shortest topological paths between two sites, each inclusive of both endpoints.
///
/// Uses BFS to determine the minimum distance, then collects every path of
/// that length. Returns an empty vector when `start` and `end` lie in
/// different connected components. Returns `vec![[start]]` when `start == end`.
/// The order of paths in the returned vector is unspecified.
///
/// # Complexity
///
/// O(V + E + P) time and O(V + P) auxiliary space, where P is the total
/// number of sites across all returned paths.
pub fn paths<M: HasBonds>(mol: &M, start: SiteId, end: SiteId) -> Vec<Vec<SiteId>> {
    if start == end {
        return vec![vec![start]];
    }

    let mut dist: HashMap<SiteId, usize> = HashMap::new();
    let mut queue: VecDeque<SiteId> = VecDeque::new();

    dist.insert(start, 0);
    queue.push_back(start);

    while let Some(site) = queue.pop_front() {
        let d = dist[&site];
        for nb in mol.neighbors(site) {
            if let std::collections::hash_map::Entry::Vacant(e) = dist.entry(nb) {
                e.insert(d + 1);
                queue.push_back(nb);
            }
        }
    }

    if !dist.contains_key(&end) {
        return vec![];
    }

    let mut result: Vec<Vec<SiteId>> = Vec::new();
    let mut stack: Vec<(SiteId, Vec<SiteId>)> = vec![(end, vec![end])];

    while let Some((site, partial)) = stack.pop() {
        if site == start {
            let mut path = partial;
            path.reverse();
            result.push(path);
            continue;
        }
        let d = dist[&site];
        for nb in mol.neighbors(site) {
            if dist.get(&nb) == Some(&(d - 1)) {
                let mut next = partial.clone();
                next.push(nb);
                stack.push((nb, next));
            }
        }
    }

    result
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
        assert_eq!(paths(&chain(), s(1), s(1)), vec![vec![s(1)]]);
        assert_eq!(paths(&chain(), s(2), s(2)), vec![vec![s(2)]]);
        assert_eq!(paths(&chain(), s(3), s(3)), vec![vec![s(3)]]);
    }

    #[test]
    fn single_site_to_itself() {
        let mol = Mol {
            sites: vec![s(1)],
            bonds: vec![],
            endpoints: vec![],
        };
        assert_eq!(paths(&mol, s(1), s(1)), vec![vec![s(1)]]);
    }

    #[test]
    fn adjacent_sites() {
        assert_eq!(paths(&chain(), s(1), s(2)), vec![vec![s(1), s(2)]]);
    }

    #[test]
    fn adjacent_sites_reversed() {
        assert_eq!(paths(&chain(), s(2), s(1)), vec![vec![s(2), s(1)]]);
    }

    #[test]
    fn chain_two_hops() {
        assert_eq!(paths(&chain(), s(1), s(3)), vec![vec![s(1), s(2), s(3)]],);
    }

    #[test]
    fn chain_two_hops_reversed() {
        assert_eq!(paths(&chain(), s(3), s(1)), vec![vec![s(3), s(2), s(1)]],);
    }

    #[test]
    fn disconnected_returns_empty() {
        assert!(paths(&two_components(), s(1), s(3)).is_empty());
        assert!(paths(&two_components(), s(3), s(1)).is_empty());
    }

    #[test]
    fn isolated_site_returns_empty() {
        assert!(paths(&two_components(), s(3), s(2)).is_empty());
    }

    #[test]
    fn ring_single_shortest_path() {
        assert_eq!(paths(&triangle(), s(1), s(3)), vec![vec![s(1), s(3)]],);
    }

    #[test]
    fn ring_path_length_is_one_hop() {
        assert_eq!(paths(&triangle(), s(1), s(3))[0].len(), 2);
    }

    #[test]
    fn star_leaf_to_leaf_through_center() {
        assert_eq!(paths(&star(), s(2), s(3)), vec![vec![s(2), s(1), s(3)]],);
    }

    #[test]
    fn star_leaf_to_center() {
        assert_eq!(paths(&star(), s(4), s(1)), vec![vec![s(4), s(1)]]);
    }

    #[test]
    fn cyclobutane_two_shortest_paths() {
        let mut ps = paths(&cyclobutane(), s(1), s(4));
        ps.sort();
        assert_eq!(ps, vec![vec![s(1), s(2), s(4)], vec![s(1), s(3), s(4)],]);
    }

    #[test]
    fn cyclobutane_all_paths_same_length() {
        let ps = paths(&cyclobutane(), s(1), s(4));
        let len = ps[0].len();
        assert!(ps.iter().all(|p| p.len() == len));
    }

    #[test]
    fn pentane_full_path() {
        assert_eq!(
            paths(&pentane(), s(1), s(5)),
            vec![vec![s(1), s(2), s(3), s(4), s(5)]],
        );
    }

    #[test]
    fn pentane_from_interior() {
        assert_eq!(paths(&pentane(), s(3), s(5)), vec![vec![s(3), s(4), s(5)]],);
    }

    #[test]
    fn all_paths_start_at_start() {
        for p in paths(&cyclobutane(), s(1), s(4)) {
            assert_eq!(p.first(), Some(&s(1)));
        }
        for p in paths(&star(), s(2), s(4)) {
            assert_eq!(p.first(), Some(&s(2)));
        }
        for p in paths(&pentane(), s(3), s(5)) {
            assert_eq!(p.first(), Some(&s(3)));
        }
    }

    #[test]
    fn all_paths_end_at_end() {
        for p in paths(&cyclobutane(), s(1), s(4)) {
            assert_eq!(p.last(), Some(&s(4)));
        }
        for p in paths(&star(), s(2), s(4)) {
            assert_eq!(p.last(), Some(&s(4)));
        }
        for p in paths(&pentane(), s(3), s(5)) {
            assert_eq!(p.last(), Some(&s(5)));
        }
    }

    #[test]
    fn all_consecutive_sites_bonded() {
        let mol = cyclobutane();
        for p in paths(&mol, s(1), s(4)) {
            for w in p.windows(2) {
                assert!(mol.bond_between(w[0], w[1]).is_some());
            }
        }
    }

    #[test]
    fn all_paths_contain_no_duplicate_sites() {
        for p in paths(&cyclobutane(), s(1), s(4)) {
            let unique: HashSet<_> = p.iter().copied().collect();
            assert_eq!(unique.len(), p.len());
        }
    }

    #[test]
    fn tree_yields_exactly_one_path() {
        assert_eq!(paths(&chain(), s(1), s(3)).len(), 1);
        assert_eq!(paths(&pentane(), s(1), s(5)).len(), 1);
        assert_eq!(paths(&star(), s(2), s(3)).len(), 1);
    }
}
