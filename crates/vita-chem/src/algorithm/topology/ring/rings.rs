use std::collections::VecDeque;

use vita_core::SiteId;

use super::{RingMembership, RingSystems};
use crate::algorithm::utils::{
    AdjacencyList, BitSet, DisjointSet, FxHashMap, Gf2Basis, SortedMultimap,
};
use crate::{BondId, HasBonds};

/// A single ring of a molecule.
///
/// The sites are ordered around the ring: `sites()[i]` is joined to
/// `sites()[i + 1]` (indices wrapping) by `bonds()[i]`. The order is canonical —
/// it starts at the ring's smallest [`SiteId`] and proceeds toward its smaller
/// neighbor — so the same ring always reads the same way.
///
/// Obtain via [`Rings`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Ring {
    sites: Vec<SiteId>,
    bonds: Vec<BondId>,
}

impl Ring {
    /// Sites in ring order.
    pub fn sites(&self) -> &[SiteId] {
        &self.sites
    }

    /// Bonds in ring order, with `bonds()[i]` joining `sites()[i]` to its
    /// successor.
    pub fn bonds(&self) -> &[BondId] {
        &self.bonds
    }
}

/// The minimum cycle basis (MCB) of a molecule.
///
/// Decomposes the cycle space into the fewest independent rings — the cycle
/// rank, which [`count`](fn@super::count) reports — chosen to have the least
/// total size. For nearly all molecules this is the chemically intended ring
/// perception; where several least-size bases exist (cage systems such as
/// cubane), the candidates' canonical order fixes one deterministically.
///
/// Obtain via [`rings`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Rings {
    rings: Vec<Ring>,
    site_index: SortedMultimap<SiteId, usize>,
    bond_index: SortedMultimap<BondId, usize>,
}

impl Rings {
    /// Number of rings in the basis (equals the cycle rank).
    pub fn len(&self) -> usize {
        self.rings.len()
    }

    /// Returns `true` if the molecule is acyclic.
    pub fn is_empty(&self) -> bool {
        self.rings.is_empty()
    }

    /// Iterates all rings of the basis, ordered by size then by their sites.
    pub fn iter(&self) -> impl Iterator<Item = &Ring> + '_ {
        self.rings.iter()
    }

    /// Size of the smallest ring (the graph girth).
    ///
    /// Returns `None` if the molecule is acyclic.
    pub fn girth(&self) -> Option<usize> {
        self.rings.iter().map(|ring| ring.sites.len()).min()
    }

    /// Iterates all rings that contain `site`, in the order of [`iter`](Rings::iter).
    ///
    /// Returns an empty iterator if `site` is absent from the molecule or lies
    /// in no ring.
    pub fn of_site(&self, site: SiteId) -> impl Iterator<Item = &Ring> + '_ {
        self.site_index
            .get(&site)
            .iter()
            .map(move |&i| &self.rings[i])
    }

    /// Iterates all rings that contain `bond`, in the order of [`iter`](Rings::iter).
    ///
    /// A bond shared by fused rings appears in more than one ring. Returns an
    /// empty iterator if `bond` is absent from the molecule or is a bridge.
    pub fn of_bond(&self, bond: BondId) -> impl Iterator<Item = &Ring> + '_ {
        self.bond_index
            .get(&bond)
            .iter()
            .map(move |&i| &self.rings[i])
    }

    /// Returns `true` if some ring contains both `a` and `b`.
    ///
    /// Returns `false` if either site is absent from the molecule or no ring
    /// holds them together.
    pub fn same(&self, a: SiteId, b: SiteId) -> bool {
        let a_rings = self.site_index.get(&a);
        let b_rings = self.site_index.get(&b);
        a_rings.iter().any(|i| b_rings.contains(i))
    }

    /// Derives ring membership from the basis.
    ///
    /// The result is independent of which basis was chosen: a site or bond lies
    /// in a ring exactly when it appears in any basis ring.
    pub fn membership(&self) -> RingMembership {
        RingMembership::from_sets(
            self.site_index.iter().map(|(&site, _)| site),
            self.bond_index.iter().map(|(&bond, _)| bond),
        )
    }

    /// The ring systems: maximal groups of rings connected through shared sites.
    ///
    /// Fused, bridged, and spiro rings all coalesce into one system. Systems are
    /// in ascending site order.
    pub fn systems(&self) -> RingSystems {
        let mut components = DisjointSet::new(self.rings.len());
        for (_, members) in self.site_index.iter() {
            for pair in members.windows(2) {
                components.union(pair[0], pair[1]);
            }
        }

        let mut systems: Vec<Vec<SiteId>> = components
            .groups()
            .into_iter()
            .map(|group| {
                let mut sites: Vec<SiteId> = group
                    .iter()
                    .flat_map(|&i| self.rings[i].sites.iter().copied())
                    .collect();
                sites.sort_unstable();
                sites.dedup();
                sites
            })
            .collect();
        systems.sort_unstable();
        RingSystems::new(systems)
    }
}

/// Minimum cycle basis of a molecule.
///
/// Builds the smallest set of independent rings spanning the cycle space using
/// Horton's algorithm: the fundamental cycles of every site's breadth-first tree
/// form the candidate pool, from which a least-total-size basis is drawn by
/// Gaussian elimination over GF(2). Candidates are ordered by size and then
/// canonically, so the chosen basis is deterministic across runs.
///
/// # Complexity
///
/// O(V · E³ / w) time and O(V · E² / w) space, over the molecule's `V` sites and
/// `E` bonds for word width `w` = 64. A breadth-first tree at each site yields
/// O(V · E) candidate cycles of `E` bits, which Gaussian elimination over GF(2)
/// reduces to the basis.
pub fn rings<M: HasBonds>(mol: &M) -> Rings {
    let mut sites: Vec<SiteId> = mol.sites().collect();
    sites.sort_unstable();
    let n = sites.len();

    let mut rows: Vec<(BondId, usize, usize)> = mol
        .bonds()
        .map(|bond| {
            let (a, b) = mol.bond_endpoints(bond);
            let i = sites.binary_search(&a).unwrap();
            let j = sites.binary_search(&b).unwrap();
            (bond, i.min(j), i.max(j))
        })
        .collect();
    rows.sort_unstable_by_key(|&(_, lo, hi)| (lo, hi));

    let graph = AdjacencyList::build(
        n,
        rows.iter()
            .enumerate()
            .map(|(edge, &(_, lo, hi))| (edge, lo, hi)),
    );

    let candidates = horton_candidates(&graph, &rows);
    let basis = minimum_basis(candidates, rows.len());

    let mut rings: Vec<Ring> = basis
        .iter()
        .map(|cycle| trace_ring(cycle, &rows, &sites))
        .collect();
    rings.sort_by(|a, b| {
        a.sites
            .len()
            .cmp(&b.sites.len())
            .then_with(|| a.sites.cmp(&b.sites))
            .then_with(|| a.bonds.cmp(&b.bonds))
    });

    let site_index = SortedMultimap::from_pairs(
        rings
            .iter()
            .enumerate()
            .flat_map(|(i, ring)| ring.sites.iter().map(move |&site| (site, i))),
    );
    let bond_index = SortedMultimap::from_pairs(
        rings
            .iter()
            .enumerate()
            .flat_map(|(i, ring)| ring.bonds.iter().map(move |&bond| (bond, i))),
    );

    Rings {
        rings,
        site_index,
        bond_index,
    }
}

/// Fundamental cycles of every site's breadth-first tree, as edge bit vectors,
/// sorted by size then canonically and deduplicated.
fn horton_candidates(graph: &AdjacencyList, rows: &[(BondId, usize, usize)]) -> Vec<BitSet> {
    let n = graph.len();
    let m = rows.len();
    let mut candidates: Vec<BitSet> = Vec::new();

    for root in 0..n {
        let mut distance = vec![usize::MAX; n];
        let mut tree_bond = vec![usize::MAX; n];
        let mut tree_parent = vec![usize::MAX; n];
        distance[root] = 0;
        let mut queue = VecDeque::new();
        queue.push_back(root);
        while let Some(u) = queue.pop_front() {
            for &(edge, v) in graph.neighbors(u) {
                if distance[v] == usize::MAX {
                    distance[v] = distance[u] + 1;
                    tree_bond[v] = edge;
                    tree_parent[v] = u;
                    queue.push_back(v);
                }
            }
        }

        for (edge, &(_, lo, hi)) in rows.iter().enumerate() {
            if distance[lo] == usize::MAX || distance[hi] == usize::MAX {
                continue;
            }
            if tree_bond[lo] == edge || tree_bond[hi] == edge {
                continue;
            }

            let mut cycle = BitSet::zeros(m);
            let mut x = lo;
            while x != root {
                cycle.toggle(tree_bond[x]);
                x = tree_parent[x];
            }
            let mut y = hi;
            while y != root {
                cycle.toggle(tree_bond[y]);
                y = tree_parent[y];
            }
            cycle.toggle(edge);

            candidates.push(cycle);
        }
    }

    candidates.sort_unstable_by(|a, b| a.count_ones().cmp(&b.count_ones()).then_with(|| a.cmp(b)));
    candidates.dedup();
    candidates
}

/// Greedy least-total-size basis of the candidate cycles over GF(2).
fn minimum_basis(candidates: Vec<BitSet>, dimension: usize) -> Vec<BitSet> {
    let mut basis = Gf2Basis::new(dimension);
    let mut chosen: Vec<BitSet> = Vec::new();
    for candidate in candidates {
        if basis.insert(candidate.clone()) {
            chosen.push(candidate);
        }
    }
    chosen
}

/// Walks an edge bit vector into a canonically ordered [`Ring`].
fn trace_ring(cycle: &BitSet, rows: &[(BondId, usize, usize)], sites: &[SiteId]) -> Ring {
    let mut local: FxHashMap<usize, Vec<(usize, usize)>> = FxHashMap::default();
    for edge in (0..cycle.len()).filter(|&e| cycle.test(e)) {
        let (_, lo, hi) = rows[edge];
        local.entry(lo).or_default().push((edge, hi));
        local.entry(hi).or_default().push((edge, lo));
    }

    debug_assert!(local.values().all(|neighbors| neighbors.len() == 2));

    let start = *local.keys().min().expect("a ring has at least one site");
    let mut ring_sites: Vec<usize> = vec![start];
    let mut ring_bonds: Vec<usize> = Vec::new();
    let mut previous = usize::MAX;
    let mut current = start;

    loop {
        let &(edge, next) = if previous == usize::MAX {
            local[&current].iter().min_by_key(|&&(_, nb)| nb).unwrap()
        } else {
            local[&current]
                .iter()
                .find(|&&(_, nb)| nb != previous)
                .unwrap()
        };
        ring_bonds.push(edge);
        if next == start {
            break;
        }
        ring_sites.push(next);
        previous = current;
        current = next;
    }

    debug_assert_eq!(
        ring_bonds.len(),
        local.len(),
        "basis vector is not a single cycle"
    );

    Ring {
        sites: ring_sites.into_iter().map(|i| sites[i]).collect(),
        bonds: ring_bonds.into_iter().map(|e| rows[e].0).collect(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::HasSites;

    use crate::topology::ring::membership;

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

    fn square() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4)],
            bonds: vec![b(1), b(2), b(3), b(4)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3)), (s(3), s(4)), (s(1), s(4))],
        }
    }

    fn tadpole() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4)],
            bonds: vec![b(1), b(2), b(3), b(4)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3)), (s(1), s(3)), (s(1), s(4))],
        }
    }

    fn triangle_and_square() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5)],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6)],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(1), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(1), s(5)),
            ],
        }
    }

    fn spiro() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5)],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6)],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(1), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(3), s(5)),
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

    fn cube() -> Mol {
        Mol {
            sites: (1..=8).map(s).collect(),
            bonds: (1..=12).map(b).collect(),
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(1), s(4)),
                (s(5), s(6)),
                (s(6), s(7)),
                (s(7), s(8)),
                (s(5), s(8)),
                (s(1), s(5)),
                (s(2), s(6)),
                (s(3), s(7)),
                (s(4), s(8)),
            ],
        }
    }

    #[test]
    fn empty_molecule_has_no_rings() {
        let r = rings(&empty());
        assert_eq!(r.len(), 0);
        assert!(r.is_empty());
    }

    #[test]
    fn single_site_has_no_rings() {
        assert!(rings(&single()).is_empty());
    }

    #[test]
    fn chain_has_no_rings() {
        assert!(rings(&chain()).is_empty());
    }

    #[test]
    fn triangle_has_one_ring() {
        assert_eq!(rings(&triangle()).len(), 1);
    }

    #[test]
    fn triangle_ring_spans_three_sites_and_three_bonds() {
        let r = rings(&triangle());
        let ring = r.iter().next().unwrap();
        assert_eq!(ring.sites().len(), 3);
        assert_eq!(ring.bonds().len(), 3);
    }

    #[test]
    fn square_has_one_ring_of_four_sites() {
        let r = rings(&square());
        assert_eq!(r.len(), 1);
        assert_eq!(r.iter().next().unwrap().sites().len(), 4);
    }

    #[test]
    fn a_ring_starts_at_its_smallest_site() {
        let r = rings(&square());
        assert_eq!(r.iter().next().unwrap().sites()[0], s(1));
    }

    #[test]
    fn ring_sites_are_consecutively_bonded() {
        let mol = square();
        let r = rings(&mol);
        let ring = r.iter().next().unwrap();
        let sites = ring.sites();
        let bonds = ring.bonds();
        for i in 0..sites.len() {
            let next = sites[(i + 1) % sites.len()];
            assert_eq!(mol.bond_between(sites[i], next), Some(bonds[i]));
        }
    }

    #[test]
    fn each_triangle_site_lies_in_one_ring() {
        let r = rings(&triangle());
        for site in [s(1), s(2), s(3)] {
            assert_eq!(r.of_site(site).count(), 1);
        }
    }

    #[test]
    fn of_site_of_an_absent_site_is_empty() {
        assert_eq!(rings(&triangle()).of_site(s(99)).count(), 0);
    }

    #[test]
    fn of_site_of_an_acyclic_site_is_empty() {
        assert_eq!(rings(&chain()).of_site(s(2)).count(), 0);
    }

    #[test]
    fn of_bond_of_a_bridge_is_empty() {
        assert_eq!(rings(&tadpole()).of_bond(b(4)).count(), 0);
    }

    #[test]
    fn of_bond_of_an_absent_bond_is_empty() {
        assert_eq!(rings(&triangle()).of_bond(b(99)).count(), 0);
    }

    #[test]
    fn same_is_false_for_an_absent_site() {
        assert!(!rings(&triangle()).same(s(1), s(99)));
    }

    #[test]
    fn same_is_false_across_distinct_rings() {
        assert!(!rings(&spiro()).same(s(1), s(4)));
    }

    #[test]
    fn ring_count_equals_the_cycle_rank() {
        assert_eq!(rings(&triangle()).len(), 1);
        assert_eq!(rings(&two_triangles()).len(), 2);
        assert_eq!(rings(&cube()).len(), 5);
    }

    #[test]
    fn girth_is_the_smallest_ring_size() {
        assert_eq!(rings(&triangle()).girth(), Some(3));
        assert_eq!(rings(&square()).girth(), Some(4));
        assert_eq!(rings(&triangle_and_square()).girth(), Some(3));
    }

    #[test]
    fn girth_of_an_acyclic_molecule_is_none() {
        assert_eq!(rings(&chain()).girth(), None);
    }

    #[test]
    fn a_bond_shared_by_fused_rings_lies_in_both() {
        assert_eq!(rings(&triangle_and_square()).of_bond(b(3)).count(), 2);
    }

    #[test]
    fn a_site_shared_by_spiro_rings_lies_in_both() {
        assert_eq!(rings(&spiro()).of_site(s(3)).count(), 2);
    }

    #[test]
    fn rings_are_ordered_by_size() {
        let r = rings(&triangle_and_square());
        let sizes: Vec<usize> = r.iter().map(|ring| ring.sites().len()).collect();
        assert_eq!(sizes, vec![3, 4]);
    }

    #[test]
    fn same_is_true_within_a_ring() {
        assert!(rings(&triangle()).same(s(1), s(2)));
        assert!(rings(&spiro()).same(s(1), s(3)));
    }

    #[test]
    fn a_minimum_basis_of_the_cube_is_all_squares() {
        let r = rings(&cube());
        assert!(r.iter().all(|ring| ring.sites().len() == 4));
    }

    #[test]
    fn membership_matches_the_free_function() {
        let mol = triangle_and_square();
        let derived = rings(&mol).membership();
        let direct = membership(&mol);
        assert!(derived.sites().eq(direct.sites()));
        assert!(derived.bonds().eq(direct.bonds()));
    }

    #[test]
    fn an_acyclic_molecule_has_acyclic_membership() {
        assert!(rings(&chain()).membership().is_acyclic());
    }

    #[test]
    fn systems_unite_fused_rings() {
        let r = rings(&triangle_and_square());
        let systems: Vec<Vec<SiteId>> =
            r.systems().iter().map(|sys| sys.iter().collect()).collect();
        assert_eq!(systems, vec![vec![s(1), s(2), s(3), s(4), s(5)]]);
    }

    #[test]
    fn systems_separate_disjoint_rings() {
        let r = rings(&two_triangles());
        let systems: Vec<Vec<SiteId>> =
            r.systems().iter().map(|sys| sys.iter().collect()).collect();
        assert_eq!(
            systems,
            vec![vec![s(1), s(2), s(3)], vec![s(4), s(5), s(6)]]
        );
    }

    #[test]
    fn a_ring_system_reports_and_counts_its_sites() {
        let r = rings(&triangle_and_square());
        let systems = r.systems();
        let system = systems.iter().next().unwrap();
        assert_eq!(system.len(), 5);
        assert!(system.contains(s(1)));
        assert!(!system.contains(s(99)));
    }

    #[test]
    fn an_acyclic_molecule_has_no_ring_systems() {
        assert!(rings(&chain()).systems().is_empty());
    }

    #[test]
    fn systems_partition_the_ring_sites() {
        let mol = spiro();
        let r = rings(&mol);
        let systems = r.systems();
        let mut from_systems: Vec<SiteId> = systems.iter().flat_map(|sys| sys.iter()).collect();
        from_systems.sort_unstable();
        let from_membership: Vec<SiteId> = r.membership().sites().collect();
        assert_eq!(from_systems, from_membership);
    }

    #[test]
    fn the_basis_is_deterministic() {
        let sites = |r: &Rings| -> Vec<Vec<SiteId>> {
            r.iter().map(|ring| ring.sites().to_vec()).collect()
        };
        assert_eq!(sites(&rings(&cube())), sites(&rings(&cube())));
    }

    #[test]
    fn the_basis_is_independent_of_input_order() {
        let shuffled = Mol {
            sites: vec![s(4), s(1), s(5), s(2), s(3)],
            bonds: vec![b(6), b(2), b(4), b(1), b(5), b(3)],
            endpoints: vec![
                (s(5), s(1)),
                (s(3), s(2)),
                (s(3), s(4)),
                (s(1), s(2)),
                (s(5), s(4)),
                (s(1), s(3)),
            ],
        };
        assert_eq!(rings(&triangle_and_square()), rings(&shuffled));
    }
}
