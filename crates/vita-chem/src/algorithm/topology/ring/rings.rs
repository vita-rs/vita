use std::collections::{HashMap, HashSet, VecDeque};

use vita_core::{HasSites, SiteId};

use super::bitset::Bits;
use super::membership::RingMembership;
use crate::{BondId, HasBonds};

/// A single ring of a molecule.
///
/// The sites are ordered around the ring: `sites()[i]` is joined to
/// `sites()[i + 1]` (indices wrapping) by `bonds()[i]`. The order is
/// canonical — it starts at the ring's smallest [`SiteId`] and proceeds
/// toward its smaller neighbour — so the same ring always reads the same way.
///
/// Obtain via [`Rings`].
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
/// perception; where several equivalent bases exist (cage systems such as
/// cubane), a canonical site-ordered tie-break makes the choice deterministic.
///
/// Obtain via [`rings`].
pub struct Rings {
    rings: Vec<Ring>,
    site_index: HashMap<SiteId, Vec<usize>>,
    bond_index: HashMap<BondId, Vec<usize>>,
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

    /// Iterates all rings of the basis.
    pub fn iter(&self) -> impl Iterator<Item = &Ring> + '_ {
        self.rings.iter()
    }

    /// Size of the smallest ring (the graph girth).
    ///
    /// Returns `None` if the molecule is acyclic.
    pub fn girth(&self) -> Option<usize> {
        self.rings.iter().map(|r| r.sites.len()).min()
    }

    /// Iterates all rings that contain `site`.
    ///
    /// Returns an empty iterator if `site` is absent from the molecule or lies
    /// in no ring.
    pub fn of_site(&self, site: SiteId) -> impl Iterator<Item = &Ring> + '_ {
        self.site_index
            .get(&site)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
            .iter()
            .map(move |&i| &self.rings[i])
    }

    /// Iterates all rings that contain `bond`.
    ///
    /// A bond shared by fused rings appears in more than one ring. Returns an
    /// empty iterator if `bond` is absent from the molecule or is a bridge.
    pub fn of_bond(&self, bond: BondId) -> impl Iterator<Item = &Ring> + '_ {
        self.bond_index
            .get(&bond)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
            .iter()
            .map(move |&i| &self.rings[i])
    }

    /// Returns `true` if some ring contains both `a` and `b`.
    ///
    /// Returns `false` if either site is absent from the molecule or no ring
    /// holds them together.
    pub fn same(&self, a: SiteId, b: SiteId) -> bool {
        match (self.site_index.get(&a), self.site_index.get(&b)) {
            (Some(ra), Some(rb)) => ra.iter().any(|i| rb.contains(i)),
            _ => false,
        }
    }

    /// Derives ring membership from the basis.
    ///
    /// The result is independent of which basis was chosen: a site or bond
    /// lies in a ring exactly when it appears in any basis ring.
    pub fn membership(&self) -> RingMembership {
        let mut sites: HashSet<SiteId> = HashSet::new();
        let mut bonds: HashSet<BondId> = HashSet::new();
        for ring in &self.rings {
            sites.extend(ring.sites.iter().copied());
            bonds.extend(ring.bonds.iter().copied());
        }
        RingMembership::from_sets(sites, bonds)
    }

    /// Iterates the ring systems, each as its sorted set of sites.
    ///
    /// A ring system is a maximal group of rings connected through shared
    /// sites; fused, bridged, and spiro rings all coalesce into one system.
    /// Systems are yielded in ascending site order.
    pub fn systems(&self) -> impl Iterator<Item = Vec<SiteId>> {
        let k = self.rings.len();
        let mut parent: Vec<usize> = (0..k).collect();

        fn find(parent: &mut [usize], mut x: usize) -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        }

        for members in self.site_index.values() {
            for w in members.windows(2) {
                let ra = find(&mut parent, w[0]);
                let rb = find(&mut parent, w[1]);
                if ra != rb {
                    parent[ra] = rb;
                }
            }
        }

        let mut systems: HashMap<usize, HashSet<SiteId>> = HashMap::new();
        for ri in 0..k {
            let root = find(&mut parent, ri);
            systems
                .entry(root)
                .or_default()
                .extend(self.rings[ri].sites.iter().copied());
        }

        let mut out: Vec<Vec<SiteId>> = systems
            .into_values()
            .map(|set| {
                let mut v: Vec<SiteId> = set.into_iter().collect();
                v.sort_unstable();
                v
            })
            .collect();
        out.sort_unstable();
        out.into_iter()
    }
}

/// Minimum cycle basis of a molecule.
///
/// Builds the smallest set of independent rings spanning the cycle space using
/// Horton's algorithm: the fundamental cycles of every vertex's breadth-first
/// tree form the candidate pool, from which a least-total-size basis is drawn
/// by Gaussian elimination over GF(2). Candidates are ordered by size and then
/// canonically, so the chosen basis is deterministic across runs.
///
/// # Complexity
///
/// O(V² · E) time and O(V · E / 64) space.
pub fn rings<M: HasBonds + HasSites>(mol: &M) -> Rings {
    let mut sites: Vec<SiteId> = mol.sites().collect();
    sites.sort_unstable();
    let n = sites.len();
    let pos: HashMap<SiteId, usize> = sites.iter().enumerate().map(|(i, &s)| (s, i)).collect();

    let mut rows: Vec<(BondId, usize, usize)> = mol
        .bonds()
        .map(|bond| {
            let (a, b) = mol.bond_endpoints(bond);
            let (i, j) = (pos[&a], pos[&b]);
            (bond, i.min(j), i.max(j))
        })
        .collect();
    rows.sort_unstable_by_key(|&(_, lo, hi)| (lo, hi));
    let m = rows.len();

    let mut adj: Vec<Vec<(usize, usize)>> = vec![vec![]; n];
    for (e, &(_, lo, hi)) in rows.iter().enumerate() {
        adj[lo].push((e, hi));
        adj[hi].push((e, lo));
    }
    for a in adj.iter_mut() {
        a.sort_unstable_by_key(|&(_, nb)| nb);
    }

    let candidates = horton_candidates(n, m, &adj, &rows);
    let basis = minimum_basis(candidates);

    let result: Vec<Ring> = basis
        .iter()
        .map(|bits| trace_ring(bits, m, &rows, &sites))
        .collect();

    let mut site_index: HashMap<SiteId, Vec<usize>> = HashMap::new();
    let mut bond_index: HashMap<BondId, Vec<usize>> = HashMap::new();
    for (i, ring) in result.iter().enumerate() {
        for &s in &ring.sites {
            site_index.entry(s).or_default().push(i);
        }
        for &b in &ring.bonds {
            bond_index.entry(b).or_default().push(i);
        }
    }

    Rings {
        rings: result,
        site_index,
        bond_index,
    }
}

/// Fundamental cycles of every vertex's BFS tree, as edge bit vectors.
fn horton_candidates(
    n: usize,
    m: usize,
    adj: &[Vec<(usize, usize)>],
    rows: &[(BondId, usize, usize)],
) -> Vec<Bits> {
    let mut candidates: Vec<Bits> = Vec::new();

    for root in 0..n {
        let mut dist = vec![usize::MAX; n];
        let mut pred_bond = vec![usize::MAX; n];
        let mut pred_node = vec![usize::MAX; n];
        dist[root] = 0;
        let mut queue = VecDeque::new();
        queue.push_back(root);
        while let Some(u) = queue.pop_front() {
            for &(e, v) in &adj[u] {
                if dist[v] == usize::MAX {
                    dist[v] = dist[u] + 1;
                    pred_bond[v] = e;
                    pred_node[v] = u;
                    queue.push_back(v);
                }
            }
        }

        for (e, &(_, lo, hi)) in rows.iter().enumerate() {
            if dist[lo] == usize::MAX || dist[hi] == usize::MAX {
                continue;
            }
            if pred_bond[lo] == e || pred_bond[hi] == e {
                continue;
            }

            let mut bits = Bits::zeros(m);
            let mut x = lo;
            while x != root {
                bits.toggle(pred_bond[x]);
                x = pred_node[x];
            }
            let mut y = hi;
            while y != root {
                bits.toggle(pred_bond[y]);
                y = pred_node[y];
            }
            bits.toggle(e);

            candidates.push(bits);
        }
    }

    candidates.sort_unstable_by(|a, b| a.count_ones().cmp(&b.count_ones()).then_with(|| a.cmp(b)));
    candidates.dedup();
    candidates
}

/// Greedy minimum-weight basis extraction over GF(2).
fn minimum_basis(candidates: Vec<Bits>) -> Vec<Bits> {
    let mut pivots: HashMap<usize, Bits> = HashMap::new();
    let mut basis: Vec<Bits> = Vec::new();

    for cand in candidates {
        let mut residue = cand.clone();
        while let Some(p) = residue.lowest_set() {
            match pivots.get(&p) {
                Some(reducer) => residue.xor(reducer),
                None => {
                    pivots.insert(p, residue);
                    basis.push(cand);
                    break;
                }
            }
        }
    }

    basis
}

/// Walks an edge bit vector into a canonically ordered [`Ring`].
fn trace_ring(bits: &Bits, m: usize, rows: &[(BondId, usize, usize)], sites: &[SiteId]) -> Ring {
    let mut local: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
    for e in (0..m).filter(|&e| bits.test(e)) {
        let (_, lo, hi) = rows[e];
        local.entry(lo).or_default().push((e, hi));
        local.entry(hi).or_default().push((e, lo));
    }

    debug_assert!(local.values().all(|nbrs| nbrs.len() == 2));

    let start = *local.keys().min().expect("a ring has at least one site");
    let mut ring_sites: Vec<usize> = vec![start];
    let mut ring_bonds: Vec<usize> = Vec::new();
    let mut prev = usize::MAX;
    let mut cur = start;

    loop {
        let nbrs = &local[&cur];
        let &(e, next) = if prev == usize::MAX {
            nbrs.iter().min_by_key(|&&(_, nb)| nb).unwrap()
        } else {
            nbrs.iter().find(|&&(_, nb)| nb != prev).unwrap()
        };
        ring_bonds.push(e);
        if next == start {
            break;
        }
        ring_sites.push(next);
        prev = cur;
        cur = next;
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

    fn fused() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6)],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6), b(7)],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(1), s(4)),
                (s(3), s(5)),
                (s(5), s(6)),
                (s(4), s(6)),
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
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6), s(7), s(8)],
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
    fn empty_has_no_rings() {
        let r = rings(&empty());
        assert_eq!(r.len(), 0);
        assert!(r.is_empty());
    }

    #[test]
    fn single_site_has_no_rings() {
        assert_eq!(rings(&single()).len(), 0);
    }

    #[test]
    fn chain_has_no_rings() {
        assert!(rings(&chain()).is_empty());
    }

    #[test]
    fn chain_girth_is_none() {
        assert_eq!(rings(&chain()).girth(), None);
    }

    #[test]
    fn triangle_has_one_ring() {
        assert_eq!(rings(&triangle()).len(), 1);
    }

    #[test]
    fn triangle_ring_size_is_three() {
        let r = rings(&triangle());
        let ring = r.iter().next().unwrap();
        assert_eq!(ring.sites().len(), 3);
        assert_eq!(ring.bonds().len(), 3);
    }

    #[test]
    fn triangle_girth_is_three() {
        assert_eq!(rings(&triangle()).girth(), Some(3));
    }

    #[test]
    fn square_has_one_ring() {
        assert_eq!(rings(&square()).len(), 1);
    }

    #[test]
    fn square_girth_is_four() {
        assert_eq!(rings(&square()).girth(), Some(4));
    }

    #[test]
    fn fused_has_two_rings() {
        assert_eq!(rings(&fused()).len(), 2);
    }

    #[test]
    fn spiro_has_two_rings() {
        assert_eq!(rings(&spiro()).len(), 2);
    }

    #[test]
    fn two_triangles_has_two_rings() {
        assert_eq!(rings(&two_triangles()).len(), 2);
    }

    #[test]
    fn cube_has_five_rings() {
        assert_eq!(rings(&cube()).len(), 5);
    }

    #[test]
    fn cube_rings_are_all_squares() {
        let r = rings(&cube());
        assert!(r.iter().all(|ring| ring.sites().len() == 4));
        assert_eq!(r.girth(), Some(4));
    }

    #[test]
    fn ring_starts_at_smallest_site() {
        let r = rings(&square());
        let ring = r.iter().next().unwrap();
        assert_eq!(ring.sites()[0], s(1));
    }

    #[test]
    fn ring_sites_are_consecutively_bonded() {
        let mol = square();
        let r = rings(&mol);
        let ring = r.iter().next().unwrap();
        let sites = ring.sites();
        let bonds = ring.bonds();
        assert_eq!(sites.len(), bonds.len());
        for i in 0..sites.len() {
            let a = sites[i];
            let b = sites[(i + 1) % sites.len()];
            assert_eq!(mol.bond_between(a, b), Some(bonds[i]));
        }
    }

    #[test]
    fn triangle_each_site_in_one_ring() {
        let r = rings(&triangle());
        for site in [s(1), s(2), s(3)] {
            assert_eq!(r.of_site(site).count(), 1);
        }
    }

    #[test]
    fn fused_shared_site_in_two_rings() {
        let r = rings(&fused());
        assert_eq!(r.of_site(s(3)).count(), 2);
        assert_eq!(r.of_site(s(4)).count(), 2);
    }

    #[test]
    fn fused_shared_bond_in_two_rings() {
        let r = rings(&fused());
        assert_eq!(r.of_bond(b(3)).count(), 2);
    }

    #[test]
    fn of_site_unknown_is_empty() {
        assert_eq!(rings(&triangle()).of_site(s(99)).count(), 0);
    }

    #[test]
    fn of_site_acyclic_is_empty() {
        assert_eq!(rings(&chain()).of_site(s(2)).count(), 0);
    }

    #[test]
    fn of_bond_bridge_is_empty() {
        assert_eq!(rings(&tadpole()).of_bond(b(4)).count(), 0);
    }

    #[test]
    fn of_bond_unknown_is_empty() {
        assert_eq!(rings(&triangle()).of_bond(b(99)).count(), 0);
    }

    #[test]
    fn triangle_same_is_true() {
        assert!(rings(&triangle()).same(s(1), s(2)));
    }

    #[test]
    fn spiro_same_within_ring() {
        assert!(rings(&spiro()).same(s(1), s(2)));
    }

    #[test]
    fn spiro_same_across_rings_is_false() {
        assert!(!rings(&spiro()).same(s(1), s(4)));
    }

    #[test]
    fn same_shared_site_is_false_across_distinct_rings() {
        let r = rings(&fused());
        assert!(!r.same(s(1), s(6)));
    }

    #[test]
    fn same_unknown_is_false() {
        assert!(!rings(&triangle()).same(s(1), s(99)));
    }

    #[test]
    fn membership_matches_free_function() {
        let mol = fused();
        let derived = rings(&mol).membership();
        let direct = super::super::membership(&mol);

        let mut a: Vec<SiteId> = derived.sites().collect();
        let mut c: Vec<SiteId> = direct.sites().collect();
        a.sort_unstable();
        c.sort_unstable();
        assert_eq!(a, c);

        let mut a: Vec<BondId> = derived.bonds().collect();
        let mut c: Vec<BondId> = direct.bonds().collect();
        a.sort_unstable();
        c.sort_unstable();
        assert_eq!(a, c);
    }

    #[test]
    fn acyclic_membership_is_acyclic() {
        assert!(rings(&chain()).membership().is_acyclic());
    }

    #[test]
    fn triangle_one_system() {
        assert_eq!(rings(&triangle()).systems().count(), 1);
    }

    #[test]
    fn fused_one_system() {
        let systems: Vec<Vec<SiteId>> = rings(&fused()).systems().collect();
        assert_eq!(systems.len(), 1);
        assert_eq!(systems[0], vec![s(1), s(2), s(3), s(4), s(5), s(6)]);
    }

    #[test]
    fn spiro_one_system() {
        assert_eq!(rings(&spiro()).systems().count(), 1);
    }

    #[test]
    fn two_triangles_two_systems() {
        let systems: Vec<Vec<SiteId>> = rings(&two_triangles()).systems().collect();
        assert_eq!(systems.len(), 2);
        assert_eq!(systems[0], vec![s(1), s(2), s(3)]);
        assert_eq!(systems[1], vec![s(4), s(5), s(6)]);
    }

    #[test]
    fn systems_partition_ring_sites() {
        let mol = fused();
        let r = rings(&mol);
        let from_systems: HashSet<SiteId> = r.systems().flatten().collect();
        let from_membership: HashSet<SiteId> = r.membership().sites().collect();
        assert_eq!(from_systems, from_membership);
    }

    #[test]
    fn cube_is_deterministic() {
        let first: Vec<Vec<SiteId>> = rings(&cube()).iter().map(|r| r.sites().to_vec()).collect();
        let second: Vec<Vec<SiteId>> = rings(&cube()).iter().map(|r| r.sites().to_vec()).collect();
        assert_eq!(first, second);
    }

    #[test]
    fn output_is_independent_of_input_order() {
        let shuffled = Mol {
            sites: vec![s(4), s(1), s(6), s(3), s(5), s(2)],
            bonds: vec![b(7), b(3), b(5), b(1), b(6), b(4), b(2)],
            endpoints: vec![
                (s(6), s(4)),
                (s(4), s(3)),
                (s(3), s(5)),
                (s(2), s(1)),
                (s(5), s(6)),
                (s(1), s(4)),
                (s(2), s(3)),
            ],
        };
        let canonical: Vec<Vec<SiteId>> =
            rings(&fused()).iter().map(|r| r.sites().to_vec()).collect();
        let reordered: Vec<Vec<SiteId>> = rings(&shuffled)
            .iter()
            .map(|r| r.sites().to_vec())
            .collect();
        assert_eq!(canonical, reordered);
    }
}
