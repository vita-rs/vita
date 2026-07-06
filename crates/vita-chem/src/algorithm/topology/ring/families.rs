use std::collections::VecDeque;

use vita_core::SiteId;

use super::{RingMembership, RingSystems};
use crate::algorithm::utils::{
    AdjacencyList, BitSet, DisjointSet, FxHashMap, FxHashSet, Gf2Basis, SortedMultimap,
};
use crate::topology::connectivity::blocks;
use crate::{BondId, HasBonds};

/// A unique ring family (URF) of a molecule.
///
/// An equivalence class of interchangeable relevant cycles: rings of one
/// [`size`](Self::size) that overlap in a bond and differ only by smaller
/// rings. Its [`sites`](Self::sites) and [`bonds`](Self::bonds) are the union
/// over the member cycles, so they may exceed `size` when the family holds
/// several interchangeable rings.
///
/// Obtain via [`RingFamilies`].
pub struct RingFamily {
    size: usize,
    sites: Vec<SiteId>,
    bonds: Vec<BondId>,
}

impl RingFamily {
    /// The number of bonds in each ring of the family.
    pub fn size(&self) -> usize {
        self.size
    }

    /// Returns `true` if `site` lies in some ring of the family.
    pub fn contains_site(&self, site: SiteId) -> bool {
        self.sites.binary_search(&site).is_ok()
    }

    /// Returns `true` if `bond` lies in some ring of the family.
    pub fn contains_bond(&self, bond: BondId) -> bool {
        self.bonds.binary_search(&bond).is_ok()
    }

    /// Number of sites lying in some ring of the family.
    pub fn site_count(&self) -> usize {
        self.sites.len()
    }

    /// Number of bonds lying in some ring of the family.
    pub fn bond_count(&self) -> usize {
        self.bonds.len()
    }

    /// Iterates the sites lying in some ring of the family, in ascending order.
    pub fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
        self.sites.iter().copied()
    }

    /// Iterates the bonds lying in some ring of the family, in ascending order.
    pub fn bonds(&self) -> impl Iterator<Item = BondId> + '_ {
        self.bonds.iter().copied()
    }
}

/// The unique ring families (URFs) of a molecule.
///
/// The canonical ring decomposition of Kolodzik, Urbaczek and Rarey
/// (J. Chem. Inf. Model. 2012, 52, 2013–2021). Unlike a minimum cycle basis
/// ([`Rings`](super::Rings)), it is fixed by the molecular graph alone, with no
/// tie-break: cubane yields six families, one per face, where a minimum basis
/// arbitrarily selects five. The family count is at least the cycle rank and may
/// exceed it.
///
/// A site or bond may lie in several families, so [`of_site`](Self::of_site) and
/// [`of_bond`](Self::of_bond) yield iterators.
///
/// Obtain via [`families`].
pub struct RingFamilies {
    families: Vec<RingFamily>,
    site_index: SortedMultimap<SiteId, usize>,
    bond_index: SortedMultimap<BondId, usize>,
}

impl RingFamilies {
    /// Number of unique ring families.
    pub fn len(&self) -> usize {
        self.families.len()
    }

    /// Returns `true` if the molecule has no rings.
    pub fn is_empty(&self) -> bool {
        self.families.is_empty()
    }

    /// Iterates the families, ordered by size then by their sites.
    pub fn iter(&self) -> impl Iterator<Item = &RingFamily> + '_ {
        self.families.iter()
    }

    /// Size of the smallest ring — the graph girth.
    ///
    /// Returns `None` if the molecule is acyclic.
    pub fn girth(&self) -> Option<usize> {
        self.families.iter().map(|f| f.size).min()
    }

    /// Iterates the families containing `site`, in the order of [`iter`](Self::iter).
    ///
    /// Empty if `site` is absent from the molecule or lies in no ring.
    pub fn of_site(&self, site: SiteId) -> impl Iterator<Item = &RingFamily> + '_ {
        self.site_index
            .get(&site)
            .iter()
            .map(|&i| &self.families[i])
    }

    /// Iterates the families containing `bond`, in the order of [`iter`](Self::iter).
    ///
    /// Empty if `bond` is absent from the molecule or is a bridge.
    pub fn of_bond(&self, bond: BondId) -> impl Iterator<Item = &RingFamily> + '_ {
        self.bond_index
            .get(&bond)
            .iter()
            .map(|&i| &self.families[i])
    }

    /// Returns `true` if some family holds both `a` and `b`.
    ///
    /// Returns `false` if either site is absent from the molecule or no family
    /// holds them together.
    pub fn same(&self, a: SiteId, b: SiteId) -> bool {
        let ra = self.site_index.get(&a);
        let rb = self.site_index.get(&b);
        ra.iter().any(|i| rb.contains(i))
    }

    /// Iterates the spiro atoms — sites where two rings meet at that site alone,
    /// sharing no bond — in ascending order.
    pub fn spiro_sites(&self) -> impl Iterator<Item = SiteId> + '_ {
        self.site_index
            .iter()
            .filter(move |&(_, fams)| self.meet_at_point(fams))
            .map(|(&site, _)| site)
    }

    /// Iterates the fusion bonds — bonds shared by two or more rings — in
    /// ascending order.
    pub fn fusion_bonds(&self) -> impl Iterator<Item = BondId> + '_ {
        self.bond_index
            .iter()
            .filter(|&(_, fams)| fams.len() >= 2)
            .map(|(&bond, _)| bond)
    }

    /// The ring membership implied by the families.
    ///
    /// A site or bond lies in a ring exactly when it belongs to some family.
    pub fn membership(&self) -> RingMembership {
        RingMembership::from_sets(
            self.families.iter().flat_map(|f| f.sites.iter().copied()),
            self.families.iter().flat_map(|f| f.bonds.iter().copied()),
        )
    }

    /// The ring systems: maximal sets of families joined through shared sites.
    ///
    /// Fused, bridged, and spiro rings coalesce into one system. Systems are
    /// ordered by their sites.
    pub fn systems(&self) -> RingSystems {
        let mut components = DisjointSet::new(self.families.len());
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
                    .flat_map(|&f| self.families[f].sites.iter().copied())
                    .collect();
                sites.sort_unstable();
                sites.dedup();
                sites
            })
            .collect();
        systems.sort_unstable();
        RingSystems::new(systems)
    }

    /// Returns `true` if two of the families share exactly one site.
    fn meet_at_point(&self, fams: &[usize]) -> bool {
        for (i, &a) in fams.iter().enumerate() {
            for &b in &fams[i + 1..] {
                let (a, b) = (&self.families[a].sites, &self.families[b].sites);
                if a.iter().filter(|&&s| b.contains(&s)).count() == 1 {
                    return true;
                }
            }
        }
        false
    }
}

/// Unique ring families of a molecule.
///
/// Decomposes each biconnected ring system on its own: Vismara's algorithm
/// enumerates the relevant-cycle prototypes, a GF(2) elimination by ascending
/// size keeps the relevant ones, and the transitive closure of the edge-overlap
/// relation among equal-size cosets fuses them into families.
///
/// # Complexity
///
/// O(V² · E) time and O(V · E) space, over the `V` sites and `E` bonds of each
/// biconnected ring system.
pub fn families<M: HasBonds>(mol: &M) -> RingFamilies {
    let mut families: Vec<RingFamily> = Vec::new();

    for block in blocks(mol).iter() {
        if !block.is_ring() {
            continue;
        }

        let mut sites: Vec<SiteId> = block.sites().collect();
        sites.sort_unstable();
        let n = sites.len();
        let pos: FxHashMap<SiteId, usize> =
            sites.iter().enumerate().map(|(i, &s)| (s, i)).collect();

        let bond_ids: Vec<BondId> = block.bonds().collect();
        let m = bond_ids.len();
        let endpoints: Vec<(usize, usize)> = bond_ids
            .iter()
            .map(|&bond| {
                let (a, b) = mol.bond_endpoints(bond);
                let (i, j) = (pos[&a], pos[&b]);
                (i.min(j), i.max(j))
            })
            .collect();

        let adjacency = AdjacencyList::build(
            n,
            endpoints
                .iter()
                .enumerate()
                .map(|(e, &(lo, hi))| (e, lo, hi)),
        );
        let edge_id: FxHashMap<(usize, usize), usize> = endpoints
            .iter()
            .enumerate()
            .map(|(e, &ends)| (ends, e))
            .collect();

        for (size, edges) in bcc_families(n, m, &adjacency, &edge_id) {
            let mut bonds: Vec<BondId> = edges.iter().map(|&e| bond_ids[e]).collect();
            let mut family_sites: FxHashSet<usize> = FxHashSet::default();
            for &e in &edges {
                let (lo, hi) = endpoints[e];
                family_sites.insert(lo);
                family_sites.insert(hi);
            }
            let mut family_sites: Vec<SiteId> =
                family_sites.into_iter().map(|i| sites[i]).collect();
            family_sites.sort_unstable();
            bonds.sort_unstable();
            families.push(RingFamily {
                size,
                sites: family_sites,
                bonds,
            });
        }
    }

    families.sort_by(|a, b| a.size.cmp(&b.size).then_with(|| a.sites.cmp(&b.sites)));

    let site_index = SortedMultimap::from_pairs(
        families
            .iter()
            .enumerate()
            .flat_map(|(i, fam)| fam.sites.iter().map(move |&s| (s, i))),
    );
    let bond_index = SortedMultimap::from_pairs(
        families
            .iter()
            .enumerate()
            .flat_map(|(i, fam)| fam.bonds.iter().map(move |&b| (b, i))),
    );

    RingFamilies {
        families,
        site_index,
        bond_index,
    }
}

/// A relevant-cycle prototype: shortest paths from `r` to `p` and to `q`, closed
/// by a direct edge (`x` is `None`) or through an apex vertex `x`.
struct Cf {
    r: usize,
    p: usize,
    q: usize,
    x: Option<usize>,
    weight: usize,
    prototype: BitSet,
}

/// All-pairs shortest-path tables under Vismara's degree-descending vertex order.
struct Apsp {
    /// Full-graph shortest distances.
    dist: Vec<Vec<usize>>,
    /// One shortest-path predecessor through following vertices only.
    pred: Vec<Vec<usize>>,
    /// Whether a vertex is reached by a shortest path of following vertices.
    reachable: Vec<Vec<bool>>,
}

impl Apsp {
    fn compute(n: usize, adjacency: &AdjacencyList, degree: &[usize]) -> Self {
        let mut dist = vec![vec![usize::MAX; n]; n];
        let mut pred = vec![vec![usize::MAX; n]; n];
        let mut reachable = vec![vec![false; n]; n];
        for r in 0..n {
            let dist_full = bfs_full(r, n, adjacency);
            let (dist_restr, pred_restr) = bfs_restricted(r, n, adjacency, degree);
            for v in 0..n {
                dist[r][v] = dist_full[v];
                pred[r][v] = pred_restr[v];
                reachable[r][v] = dist_restr[v] != usize::MAX && dist_restr[v] == dist_full[v];
            }
        }
        Apsp {
            dist,
            pred,
            reachable,
        }
    }
}

/// The families of one biconnected ring system, each as its ring size and the
/// set of edge indices it spans.
fn bcc_families(
    n: usize,
    m: usize,
    adjacency: &AdjacencyList,
    edge_id: &FxHashMap<(usize, usize), usize>,
) -> Vec<(usize, FxHashSet<usize>)> {
    let degree: Vec<usize> = (0..n).map(|u| adjacency.neighbors(u).len()).collect();
    let apsp = Apsp::compute(n, adjacency, &degree);

    let mut cfs = vismara(n, m, adjacency, edge_id, &degree, &apsp);
    cfs.sort_by_key(|c| c.weight);

    let mut basis = Gf2Basis::new(m);
    let mut urfs: Vec<(usize, FxHashSet<usize>)> = Vec::new();

    let mut i = 0;
    while i < cfs.len() {
        let weight = cfs[i].weight;
        let mut j = i;
        while j < cfs.len() && cfs[j].weight == weight {
            j += 1;
        }

        let mut relevant: Vec<(usize, BitSet)> = Vec::new();
        for (offset, cf) in cfs[i..j].iter().enumerate() {
            let reduced = basis.reduce(&cf.prototype);
            if !reduced.is_zero() {
                relevant.push((i + offset, reduced));
            }
        }

        let mut edges: FxHashMap<usize, FxHashSet<usize>> = FxHashMap::default();
        for (k, _) in &relevant {
            edges.insert(*k, family_edges(&cfs[*k], adjacency, edge_id, &apsp));
        }

        let mut groups: FxHashMap<BitSet, Vec<usize>> = FxHashMap::default();
        for (k, reduced) in &relevant {
            groups.entry(reduced.clone()).or_default().push(*k);
        }
        for members in groups.values() {
            for component in edge_overlap_components(members, &edges) {
                let mut union: FxHashSet<usize> = FxHashSet::default();
                for k in component {
                    union.extend(&edges[&k]);
                }
                urfs.push((weight, union));
            }
        }

        for (_, reduced) in relevant {
            basis.insert(reduced);
        }
        i = j;
    }

    urfs
}

/// Vismara's relevant-cycle prototypes for the biconnected system.
fn vismara(
    n: usize,
    m: usize,
    adjacency: &AdjacencyList,
    edge_id: &FxHashMap<(usize, usize), usize>,
    degree: &[usize],
    apsp: &Apsp,
) -> Vec<Cf> {
    let mut cfs: Vec<Cf> = Vec::new();

    for r in 0..n {
        for y in 0..n {
            if !apsp.reachable[r][y] {
                continue;
            }
            let mut even_cand: Vec<usize> = Vec::new();
            for &(_, z) in adjacency.neighbors(y) {
                if !apsp.reachable[r][z] {
                    continue;
                }
                if apsp.dist[r][z] + 1 == apsp.dist[r][y] {
                    even_cand.push(z);
                } else if apsp.dist[r][z] != apsp.dist[r][y] + 1
                    && (degree[z] < degree[y] || (degree[z] == degree[y] && z < y))
                    && paths_share_only_start(r, y, z, apsp)
                {
                    cfs.push(make_cf(r, y, z, None, m, edge_id, apsp));
                }
            }
            for a in 0..even_cand.len() {
                for b in (a + 1)..even_cand.len() {
                    let (p, q) = (even_cand[a], even_cand[b]);
                    if paths_share_only_start(r, p, q, apsp) {
                        cfs.push(make_cf(r, p, q, Some(y), m, edge_id, apsp));
                    }
                }
            }
        }
    }

    cfs
}

/// Full breadth-first distances from `root`.
fn bfs_full(root: usize, n: usize, adjacency: &AdjacencyList) -> Vec<usize> {
    let mut dist = vec![usize::MAX; n];
    dist[root] = 0;
    let mut queue = VecDeque::new();
    queue.push_back(root);
    while let Some(u) = queue.pop_front() {
        for &(_, v) in adjacency.neighbors(u) {
            if dist[v] == usize::MAX {
                dist[v] = dist[u] + 1;
                queue.push_back(v);
            }
        }
    }
    dist
}

/// Breadth-first distances and predecessors from `root`, passing only through
/// vertices that follow `root` in the degree-descending order.
fn bfs_restricted(
    root: usize,
    n: usize,
    adjacency: &AdjacencyList,
    degree: &[usize],
) -> (Vec<usize>, Vec<usize>) {
    let mut dist = vec![usize::MAX; n];
    let mut pred = vec![usize::MAX; n];
    dist[root] = 0;
    pred[root] = root;
    let mut queue = VecDeque::new();
    queue.push_back(root);
    while let Some(u) = queue.pop_front() {
        for &(_, v) in adjacency.neighbors(u) {
            if dist[v] == usize::MAX && follows(v, root, degree) {
                dist[v] = dist[u] + 1;
                pred[v] = u;
                queue.push_back(v);
            }
        }
    }
    (dist, pred)
}

/// Returns `true` if `v` follows `root` in the degree-descending vertex order.
fn follows(v: usize, root: usize, degree: &[usize]) -> bool {
    degree[v] < degree[root] || (degree[v] == degree[root] && v < root)
}

/// Returns `true` if the predecessor paths from `r` to `y` and to `z` meet only
/// at `r`.
fn paths_share_only_start(r: usize, y: usize, z: usize, apsp: &Apsp) -> bool {
    let mut on_ry: FxHashSet<usize> = FxHashSet::default();
    let mut v = y;
    on_ry.insert(v);
    while v != r {
        v = apsp.pred[r][v];
        on_ry.insert(v);
    }
    let mut shared = 0;
    let mut v = z;
    if on_ry.contains(&v) {
        shared += 1;
    }
    while v != r {
        v = apsp.pred[r][v];
        if on_ry.contains(&v) {
            shared += 1;
        }
    }
    shared == 1
}

/// Builds the prototype cycle of a family from the single predecessor paths.
fn make_cf(
    r: usize,
    p: usize,
    q: usize,
    x: Option<usize>,
    m: usize,
    edge_id: &FxHashMap<(usize, usize), usize>,
    apsp: &Apsp,
) -> Cf {
    let mut prototype = BitSet::zeros(m);
    trace_path(r, p, &mut prototype, edge_id, apsp);
    trace_path(r, q, &mut prototype, edge_id, apsp);
    let weight = match x {
        Some(apex) => {
            prototype.set(edge_id[&ends(p, apex)]);
            prototype.set(edge_id[&ends(q, apex)]);
            apsp.dist[r][p] + apsp.dist[r][q] + 2
        }
        None => {
            prototype.set(edge_id[&ends(p, q)]);
            apsp.dist[r][p] + apsp.dist[r][q] + 1
        }
    };
    Cf {
        r,
        p,
        q,
        x,
        weight,
        prototype,
    }
}

/// Sets the edges of the single predecessor path from `r` to `target`.
fn trace_path(
    r: usize,
    target: usize,
    bits: &mut BitSet,
    edge_id: &FxHashMap<(usize, usize), usize>,
    apsp: &Apsp,
) {
    let mut v = target;
    while v != r {
        let u = apsp.pred[r][v];
        bits.set(edge_id[&ends(v, u)]);
        v = u;
    }
}

/// All edges a family spans: every shortest-path edge from `r` to `p` and to
/// `q`, plus the closing edges.
fn family_edges(
    cf: &Cf,
    adjacency: &AdjacencyList,
    edge_id: &FxHashMap<(usize, usize), usize>,
    apsp: &Apsp,
) -> FxHashSet<usize> {
    let mut edges: FxHashSet<usize> = FxHashSet::default();
    collect_shortest_edges(cf.r, cf.p, adjacency, apsp, &mut edges);
    collect_shortest_edges(cf.r, cf.q, adjacency, apsp, &mut edges);
    match cf.x {
        Some(apex) => {
            edges.insert(edge_id[&ends(cf.p, apex)]);
            edges.insert(edge_id[&ends(cf.q, apex)]);
        }
        None => {
            edges.insert(edge_id[&ends(cf.p, cf.q)]);
        }
    }
    edges
}

/// Collects every edge lying on a shortest path from `r` to `target`.
fn collect_shortest_edges(
    r: usize,
    target: usize,
    adjacency: &AdjacencyList,
    apsp: &Apsp,
    edges: &mut FxHashSet<usize>,
) {
    let mut visited = vec![false; adjacency.len()];
    let mut stack = vec![target];
    visited[target] = true;
    while let Some(v) = stack.pop() {
        if v == r {
            continue;
        }
        for &(e, u) in adjacency.neighbors(v) {
            if apsp.reachable[r][u] && apsp.dist[r][u] + 1 == apsp.dist[r][v] {
                edges.insert(e);
                if !visited[u] {
                    visited[u] = true;
                    stack.push(u);
                }
            }
        }
    }
}

/// The connected components of a coset's families under edge overlap.
fn edge_overlap_components(
    members: &[usize],
    edges: &FxHashMap<usize, FxHashSet<usize>>,
) -> Vec<Vec<usize>> {
    let mut overlap = DisjointSet::new(members.len());
    for a in 0..members.len() {
        for b in (a + 1)..members.len() {
            if !edges[&members[a]].is_disjoint(&edges[&members[b]]) {
                overlap.union(a, b);
            }
        }
    }
    overlap
        .groups()
        .into_iter()
        .map(|group| group.into_iter().map(|i| members[i]).collect())
        .collect()
}

/// Orders a pair of vertices as `(min, max)` for edge lookup.
fn ends(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::HasSites;

    use crate::BondId;
    use crate::topology::ring::{membership, rings};

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

    fn mol(sites: &[u32], bonds: &[(u32, u32, u32)]) -> Mol {
        Mol {
            sites: sites.iter().map(|&n| s(n)).collect(),
            bonds: bonds.iter().map(|&(id, _, _)| b(id)).collect(),
            endpoints: bonds.iter().map(|&(_, u, v)| (s(u), s(v))).collect(),
        }
    }

    fn reversed(m: &Mol) -> Mol {
        Mol {
            sites: m.sites.iter().rev().copied().collect(),
            bonds: m.bonds.iter().rev().copied().collect(),
            endpoints: m.endpoints.iter().rev().copied().collect(),
        }
    }

    fn empty() -> Mol {
        mol(&[], &[])
    }

    fn single() -> Mol {
        mol(&[1], &[])
    }

    fn chain() -> Mol {
        mol(&[1, 2, 3], &[(1, 1, 2), (2, 2, 3)])
    }

    fn triangle() -> Mol {
        mol(&[1, 2, 3], &[(1, 1, 2), (2, 2, 3), (3, 1, 3)])
    }

    fn square() -> Mol {
        mol(&[1, 2, 3, 4], &[(1, 1, 2), (2, 2, 3), (3, 3, 4), (4, 1, 4)])
    }

    fn tadpole() -> Mol {
        mol(&[1, 2, 3, 4], &[(1, 1, 2), (2, 2, 3), (3, 1, 3), (4, 1, 4)])
    }

    fn fused() -> Mol {
        mol(
            &[1, 2, 3, 4, 5, 6],
            &[
                (1, 1, 2),
                (2, 2, 3),
                (3, 3, 4),
                (4, 1, 4),
                (5, 3, 5),
                (6, 5, 6),
                (7, 4, 6),
            ],
        )
    }

    fn spiro() -> Mol {
        mol(
            &[1, 2, 3, 4, 5],
            &[
                (1, 1, 2),
                (2, 2, 3),
                (3, 1, 3),
                (4, 3, 4),
                (5, 4, 5),
                (6, 3, 5),
            ],
        )
    }

    fn two_triangles() -> Mol {
        mol(
            &[1, 2, 3, 4, 5, 6],
            &[
                (1, 1, 2),
                (2, 2, 3),
                (3, 1, 3),
                (4, 4, 5),
                (5, 5, 6),
                (6, 4, 6),
            ],
        )
    }

    fn bridged_square() -> Mol {
        mol(
            &[1, 2, 3, 4, 5, 6],
            &[
                (1, 1, 2),
                (2, 2, 3),
                (3, 3, 4),
                (4, 1, 4),
                (5, 2, 5),
                (6, 5, 6),
                (7, 4, 6),
            ],
        )
    }

    fn k4() -> Mol {
        mol(
            &[1, 2, 3, 4],
            &[
                (1, 1, 2),
                (2, 1, 3),
                (3, 1, 4),
                (4, 2, 3),
                (5, 2, 4),
                (6, 3, 4),
            ],
        )
    }

    fn cube() -> Mol {
        mol(
            &[1, 2, 3, 4, 5, 6, 7, 8],
            &[
                (1, 1, 2),
                (2, 2, 3),
                (3, 3, 4),
                (4, 1, 4),
                (5, 5, 6),
                (6, 6, 7),
                (7, 7, 8),
                (8, 5, 8),
                (9, 1, 5),
                (10, 2, 6),
                (11, 3, 7),
                (12, 4, 8),
            ],
        )
    }

    fn bicyclo222() -> Mol {
        mol(
            &[1, 2, 3, 4, 5, 6, 7, 8],
            &[
                (1, 1, 3),
                (2, 3, 4),
                (3, 4, 2),
                (4, 1, 5),
                (5, 5, 6),
                (6, 6, 2),
                (7, 1, 7),
                (8, 7, 8),
                (9, 8, 2),
            ],
        )
    }

    #[test]
    fn empty_molecule_has_no_families() {
        let f = families(&empty());
        assert_eq!(f.len(), 0);
        assert!(f.is_empty());
    }

    #[test]
    fn single_site_has_no_families() {
        assert!(families(&single()).is_empty());
    }

    #[test]
    fn acyclic_molecule_has_no_families() {
        assert!(families(&chain()).is_empty());
    }

    #[test]
    fn a_triangle_is_one_family() {
        assert_eq!(families(&triangle()).len(), 1);
    }

    #[test]
    fn a_triangle_family_is_a_three_membered_ring() {
        let f = families(&triangle());
        let family = f.iter().next().unwrap();
        assert_eq!(family.size(), 3);
        assert_eq!(family.site_count(), 3);
        assert_eq!(family.bond_count(), 3);
    }

    #[test]
    fn a_family_reports_membership_of_its_sites_and_bonds() {
        let f = families(&triangle());
        let family = f.iter().next().unwrap();
        assert!(family.contains_site(s(1)));
        assert!(!family.contains_site(s(99)));
        assert!(family.contains_bond(b(1)));
        assert!(!family.contains_bond(b(99)));
    }

    #[test]
    fn a_square_is_one_family_of_size_four() {
        let f = families(&square());
        assert_eq!(f.len(), 1);
        assert_eq!(f.iter().next().unwrap().size(), 4);
    }

    #[test]
    fn fused_rings_yield_two_families() {
        assert_eq!(families(&fused()).len(), 2);
    }

    #[test]
    fn iter_orders_families_by_ascending_size() {
        let f = families(&bridged_square());
        let sizes: Vec<usize> = f.iter().map(|fam| fam.size()).collect();
        assert!(sizes.windows(2).all(|w| w[0] <= w[1]));
    }

    #[test]
    fn of_site_of_an_unknown_site_is_empty() {
        assert_eq!(families(&triangle()).of_site(s(99)).count(), 0);
    }

    #[test]
    fn of_site_in_an_acyclic_molecule_is_empty() {
        assert_eq!(families(&chain()).of_site(s(2)).count(), 0);
    }

    #[test]
    fn of_bond_of_a_bridge_is_empty() {
        assert_eq!(families(&tadpole()).of_bond(b(4)).count(), 0);
    }

    #[test]
    fn of_bond_of_an_unknown_bond_is_empty() {
        assert_eq!(families(&triangle()).of_bond(b(99)).count(), 0);
    }

    #[test]
    fn same_is_false_for_an_unknown_site() {
        assert!(!families(&triangle()).same(s(1), s(99)));
    }

    #[test]
    fn an_acyclic_molecule_has_no_girth() {
        assert_eq!(families(&chain()).girth(), None);
    }

    #[test]
    fn a_degenerate_ring_system_yields_a_family_per_relevant_cycle() {
        for (m, count, size) in [(cube(), 6, 4), (k4(), 4, 3), (bicyclo222(), 3, 6)] {
            let f = families(&m);
            assert_eq!(f.len(), count);
            assert!(f.iter().all(|fam| fam.size() == size));
        }
    }

    #[test]
    fn unique_ring_families_exceed_a_minimum_basis_when_degenerate() {
        for (m, rank) in [(cube(), 5), (k4(), 3), (bicyclo222(), 2)] {
            let basis = rings(&m).len();
            assert_eq!(basis, rank);
            assert!(families(&m).len() > basis);
        }
    }

    #[test]
    fn unique_ring_families_match_a_minimum_basis_when_unambiguous() {
        for m in [triangle(), square(), fused(), two_triangles()] {
            assert_eq!(families(&m).len(), rings(&m).len());
        }
    }

    #[test]
    fn a_bridged_system_fuses_interchangeable_rings() {
        let f = families(&bridged_square());
        assert_eq!(f.len(), 2);
        let mut sizes: Vec<usize> = f.iter().map(|fam| fam.size()).collect();
        sizes.sort_unstable();
        assert_eq!(sizes, vec![4, 5]);
        let five = f.iter().find(|fam| fam.size() == 5).unwrap();
        assert_eq!(five.site_count(), 6);
    }

    #[test]
    fn girth_is_the_smallest_family_size() {
        assert_eq!(families(&triangle()).girth(), Some(3));
        assert_eq!(families(&square()).girth(), Some(4));
        assert_eq!(families(&bridged_square()).girth(), Some(4));
    }

    #[test]
    fn a_site_of_one_ring_lies_in_one_family() {
        let f = families(&triangle());
        for site in [s(1), s(2), s(3)] {
            assert_eq!(f.of_site(site).count(), 1);
        }
    }

    #[test]
    fn a_shared_site_lies_in_two_families() {
        let f = families(&fused());
        assert_eq!(f.of_site(s(3)).count(), 2);
        assert_eq!(f.of_site(s(4)).count(), 2);
    }

    #[test]
    fn a_shared_bond_lies_in_two_families() {
        assert_eq!(families(&fused()).of_bond(b(3)).count(), 2);
    }

    #[test]
    fn two_sites_of_one_ring_are_the_same() {
        assert!(families(&triangle()).same(s(1), s(2)));
    }

    #[test]
    fn two_sites_in_distinct_rings_are_not_the_same() {
        assert!(!families(&fused()).same(s(1), s(6)));
    }

    #[test]
    fn a_spiro_atom_is_reported() {
        let f = families(&spiro());
        assert_eq!(f.spiro_sites().collect::<Vec<_>>(), vec![s(3)]);
    }

    #[test]
    fn fused_rings_have_no_spiro_atom() {
        assert_eq!(families(&fused()).spiro_sites().count(), 0);
    }

    #[test]
    fn a_fusion_bond_is_reported() {
        let f = families(&fused());
        assert_eq!(f.fusion_bonds().collect::<Vec<_>>(), vec![b(3)]);
    }

    #[test]
    fn spiro_rings_have_no_fusion_bond() {
        assert_eq!(families(&spiro()).fusion_bonds().count(), 0);
    }

    #[test]
    fn membership_matches_the_membership_function() {
        let m = fused();
        let derived = families(&m).membership();
        let direct = membership(&m);
        assert!(derived.sites().eq(direct.sites()));
        assert!(derived.bonds().eq(direct.bonds()));
    }

    #[test]
    fn fused_rings_form_one_system() {
        let f = families(&fused());
        let systems: Vec<Vec<SiteId>> =
            f.systems().iter().map(|sys| sys.iter().collect()).collect();
        assert_eq!(systems, vec![vec![s(1), s(2), s(3), s(4), s(5), s(6)]]);
    }

    #[test]
    fn spiro_rings_form_one_system() {
        assert_eq!(families(&spiro()).systems().len(), 1);
    }

    #[test]
    fn separate_rings_form_separate_systems() {
        let f = families(&two_triangles());
        let systems: Vec<Vec<SiteId>> =
            f.systems().iter().map(|sys| sys.iter().collect()).collect();
        assert_eq!(
            systems,
            vec![vec![s(1), s(2), s(3)], vec![s(4), s(5), s(6)]]
        );
    }

    #[test]
    fn systems_partition_the_ring_sites() {
        let m = fused();
        let f = families(&m);
        let systems = f.systems();
        let mut from_systems: Vec<SiteId> = systems.iter().flat_map(|sys| sys.iter()).collect();
        from_systems.sort_unstable();
        let from_membership: Vec<SiteId> = f.membership().sites().collect();
        assert_eq!(from_systems, from_membership);
    }

    #[test]
    fn families_are_independent_of_input_order() {
        let shape = |m: &Mol| -> Vec<(usize, Vec<SiteId>)> {
            families(m)
                .iter()
                .map(|f| (f.size(), f.sites().collect::<Vec<_>>()))
                .collect()
        };
        for m in [fused(), bridged_square(), cube()] {
            assert_eq!(shape(&m), shape(&reversed(&m)));
        }
    }

    #[test]
    fn derived_queries_are_independent_of_input_order() {
        let query = |m: &Mol| -> (Vec<SiteId>, Vec<BondId>, Vec<Vec<SiteId>>) {
            let f = families(m);
            let systems: Vec<Vec<SiteId>> =
                f.systems().iter().map(|sys| sys.iter().collect()).collect();
            (
                f.spiro_sites().collect(),
                f.fusion_bonds().collect(),
                systems,
            )
        };
        for m in [spiro(), fused(), two_triangles()] {
            assert_eq!(query(&m), query(&reversed(&m)));
        }
    }
}
