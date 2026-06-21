use std::collections::{HashMap, HashSet, VecDeque};

use vita_core::{HasSites, SiteId};

use super::bitset::Bits;
use super::membership::RingMembership;
use crate::topology::connectivity::blocks;
use crate::{BondId, HasBonds};

/// A unique ring family (URF) of a molecule.
///
/// A URF is an equivalence class of relevant cycles that are interchangeable:
/// they share the same size and differ only by smaller rings while overlapping
/// in at least one bond. Every cycle of a family has the same size, reported by
/// [`Self::size`]. The family's [`Self::sites`] and [`Self::bonds`] are the
/// union over all its member cycles, so they may exceed `size` when the family
/// holds several interchangeable rings.
///
/// Obtain via [`RingFamilies`].
pub struct RingFamily {
    size: usize,
    sites: Vec<SiteId>,
    bonds: Vec<BondId>,
}

impl RingFamily {
    /// The size of every ring in the family (number of bonds in each cycle).
    pub fn size(&self) -> usize {
        self.size
    }

    /// All sites that lie in some ring of the family, in ascending order.
    pub fn sites(&self) -> &[SiteId] {
        &self.sites
    }

    /// All bonds that lie in some ring of the family, in ascending order.
    pub fn bonds(&self) -> &[BondId] {
        &self.bonds
    }
}

/// The unique ring families (URFs) of a molecule.
///
/// Decomposes the ring topology into the canonical unique ring families of
/// Kolodzik, Urbaczek and Rarey (J. Chem. Inf. Model. 2012, 52, 2013–2021).
/// Unlike a minimum cycle basis ([`Rings`](super::Rings)), the decomposition is
/// determined solely by the molecular graph, independent of any tie-break. The
/// number of families is at least the cycle rank and may exceed it: cubane
/// yields six families, one per face, where a minimum cycle basis arbitrarily
/// selects five.
///
/// A site or bond may belong to several families, so [`Self::of_site`] and
/// [`Self::of_bond`] yield iterators. Counting the families on an atom is the
/// graph-unique form of the SMARTS ring-membership query.
///
/// Obtain via [`families`].
pub struct RingFamilies {
    families: Vec<RingFamily>,
    site_index: HashMap<SiteId, Vec<usize>>,
    bond_index: HashMap<BondId, Vec<usize>>,
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

    /// Iterates all ring families.
    pub fn iter(&self) -> impl Iterator<Item = &RingFamily> + '_ {
        self.families.iter()
    }

    /// Size of the smallest ring (the graph girth).
    ///
    /// Returns `None` if the molecule is acyclic.
    pub fn girth(&self) -> Option<usize> {
        self.families.iter().map(|f| f.size).min()
    }

    /// Iterates all families that contain `site`.
    ///
    /// Returns an empty iterator if `site` is absent from the molecule or lies
    /// in no ring.
    pub fn of_site(&self, site: SiteId) -> impl Iterator<Item = &RingFamily> + '_ {
        self.site_index
            .get(&site)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
            .iter()
            .map(move |&i| &self.families[i])
    }

    /// Iterates all families that contain `bond`.
    ///
    /// Returns an empty iterator if `bond` is absent from the molecule or is a
    /// bridge.
    pub fn of_bond(&self, bond: BondId) -> impl Iterator<Item = &RingFamily> + '_ {
        self.bond_index
            .get(&bond)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
            .iter()
            .map(move |&i| &self.families[i])
    }

    /// Returns `true` if some family contains both `a` and `b`.
    ///
    /// Returns `false` if either site is absent from the molecule or no family
    /// holds them together.
    pub fn same(&self, a: SiteId, b: SiteId) -> bool {
        match (self.site_index.get(&a), self.site_index.get(&b)) {
            (Some(ra), Some(rb)) => ra.iter().any(|i| rb.contains(i)),
            _ => false,
        }
    }

    /// Iterates the spiro atoms of the molecule.
    ///
    /// A site is a spiro atom when two rings meet at it alone, sharing that one
    /// site and no bond.
    pub fn spiro_sites(&self) -> impl Iterator<Item = SiteId> + '_ {
        self.site_index
            .iter()
            .filter(|&(_, fams)| self.meet_at_point(fams))
            .map(|(&s, _)| s)
    }

    /// Iterates the fusion bonds of the molecule.
    ///
    /// A bond is a fusion bond when two or more rings share it, so it lies in
    /// two or more families.
    pub fn fusion_bonds(&self) -> impl Iterator<Item = BondId> + '_ {
        self.bond_index
            .iter()
            .filter(|(_, v)| v.len() >= 2)
            .map(|(&b, _)| b)
    }

    /// Derives ring membership from the families.
    ///
    /// A site or bond lies in a ring exactly when it appears in some family.
    pub fn membership(&self) -> RingMembership {
        let mut sites: HashSet<SiteId> = HashSet::new();
        let mut bonds: HashSet<BondId> = HashSet::new();
        for family in &self.families {
            sites.extend(family.sites.iter().copied());
            bonds.extend(family.bonds.iter().copied());
        }
        RingMembership::from_sets(sites, bonds)
    }

    /// Iterates the ring systems, each as its sorted set of sites.
    ///
    /// A ring system is a maximal group of families connected through shared
    /// sites; fused, bridged, and spiro rings all coalesce into one system.
    /// Systems are yielded in ascending site order.
    pub fn systems(&self) -> impl Iterator<Item = Vec<SiteId>> {
        let k = self.families.len();
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
        for fi in 0..k {
            let root = find(&mut parent, fi);
            systems
                .entry(root)
                .or_default()
                .extend(self.families[fi].sites.iter().copied());
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
/// Computes the canonical URF decomposition. Each biconnected ring system is
/// processed independently: Vismara's algorithm enumerates relevant cycle
/// prototypes, a GF(2) elimination by ascending size marks the relevant ones
/// and the linearly dependent pairs, an edge-overlap test confirms the
/// URF-pair-relation, and its transitive closure yields the families.
///
/// # Complexity
///
/// Polynomial in the size of each ring system.
pub fn families<M: HasBonds + HasSites>(mol: &M) -> RingFamilies {
    let bccs = blocks(mol);
    let mut families: Vec<RingFamily> = Vec::new();

    for block in bccs.iter() {
        if !block.is_ring() {
            continue;
        }

        let mut sites: Vec<SiteId> = block.sites().to_vec();
        sites.sort_unstable();
        let n = sites.len();
        let pos: HashMap<SiteId, usize> = sites.iter().enumerate().map(|(i, &s)| (s, i)).collect();

        let mut rows: Vec<(BondId, usize, usize)> = block
            .bonds()
            .iter()
            .map(|&bond| {
                let (a, b) = mol.bond_endpoints(bond);
                let (i, j) = (pos[&a], pos[&b]);
                (bond, i.min(j), i.max(j))
            })
            .collect();
        rows.sort_unstable_by_key(|&(_, lo, hi)| (lo, hi));
        let m = rows.len();

        let mut adj: Vec<Vec<(usize, usize)>> = vec![vec![]; n];
        let mut edge_id: HashMap<(usize, usize), usize> = HashMap::new();
        for (e, &(_, lo, hi)) in rows.iter().enumerate() {
            adj[lo].push((e, hi));
            adj[hi].push((e, lo));
            edge_id.insert((lo, hi), e);
        }
        for a in adj.iter_mut() {
            a.sort_unstable_by_key(|&(_, nb)| nb);
        }

        for (size, edges) in bcc_families(n, m, &adj, &edge_id) {
            let mut bonds: Vec<BondId> = edges.iter().map(|&e| rows[e].0).collect();
            let mut family_sites: HashSet<usize> = HashSet::new();
            for &e in &edges {
                let (_, lo, hi) = rows[e];
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

    let mut site_index: HashMap<SiteId, Vec<usize>> = HashMap::new();
    let mut bond_index: HashMap<BondId, Vec<usize>> = HashMap::new();
    for (i, family) in families.iter().enumerate() {
        for &s in &family.sites {
            site_index.entry(s).or_default().push(i);
        }
        for &b in &family.bonds {
            bond_index.entry(b).or_default().push(i);
        }
    }

    RingFamilies {
        families,
        site_index,
        bond_index,
    }
}

/// A relevant cycle family prototype, defined by Vismara's vertices.
struct Cf {
    r: usize,
    p: usize,
    q: usize,
    x: Option<usize>,
    weight: usize,
    prototype: Bits,
}

/// All-pairs shortest paths of a component under Vismara's vertex ordering.
struct Apsp {
    /// Full-graph shortest distances.
    dist: Vec<Vec<usize>>,
    /// One shortest-path predecessor through following vertices.
    pred: Vec<Vec<usize>>,
    /// Whether a vertex is reached by a shortest path of following vertices.
    reachable: Vec<Vec<bool>>,
}

impl Apsp {
    fn compute(n: usize, adj: &[Vec<(usize, usize)>], degree: &[usize]) -> Self {
        let mut dist = vec![vec![usize::MAX; n]; n];
        let mut pred = vec![vec![usize::MAX; n]; n];
        let mut reachable = vec![vec![false; n]; n];
        for r in 0..n {
            let dist_full = bfs_full(r, n, adj);
            let (dist_restr, pred_restr) = bfs_restricted(r, n, adj, degree);
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

/// Unique ring families within a single biconnected component, each as its
/// ring size and the set of edge indices spanned by the family.
fn bcc_families(
    n: usize,
    m: usize,
    adj: &[Vec<(usize, usize)>],
    edge_id: &HashMap<(usize, usize), usize>,
) -> Vec<(usize, HashSet<usize>)> {
    let degree: Vec<usize> = adj.iter().map(|a| a.len()).collect();
    let apsp = Apsp::compute(n, adj, &degree);

    let mut cfs = vismara(n, m, adj, edge_id, &degree, &apsp);
    cfs.sort_by_key(|c| c.weight);

    let mut basis: Vec<(usize, Bits)> = Vec::new();
    let mut urfs: Vec<(usize, HashSet<usize>)> = Vec::new();

    let mut i = 0;
    while i < cfs.len() {
        let weight = cfs[i].weight;
        let mut j = i;
        while j < cfs.len() && cfs[j].weight == weight {
            j += 1;
        }

        let mut relevant: Vec<(usize, Bits)> = Vec::new();
        for (k, cf) in (i..j).zip(&cfs[i..j]) {
            let mut reduced = cf.prototype.clone();
            reduce(&mut reduced, &basis);
            if !reduced.is_zero() {
                relevant.push((k, reduced));
            }
        }

        let mut edges: HashMap<usize, HashSet<usize>> = HashMap::new();
        for (k, _) in &relevant {
            edges.insert(*k, family_edges(&cfs[*k], n, adj, edge_id, &apsp));
        }

        let mut groups: HashMap<Bits, Vec<usize>> = HashMap::new();
        for (k, reduced) in &relevant {
            groups.entry(reduced.clone()).or_default().push(*k);
        }
        for members in groups.values() {
            for component in edge_overlap_components(members, &edges) {
                let mut union: HashSet<usize> = HashSet::new();
                for k in component {
                    union.extend(&edges[&k]);
                }
                urfs.push((weight, union));
            }
        }

        for (_, reduced) in relevant {
            extend(&mut basis, reduced);
        }

        i = j;
    }

    urfs
}

/// Vismara's relevant cycle prototypes for the component.
fn vismara(
    n: usize,
    m: usize,
    adj: &[Vec<(usize, usize)>],
    edge_id: &HashMap<(usize, usize), usize>,
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
            for &(_, z) in &adj[y] {
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
fn bfs_full(root: usize, n: usize, adj: &[Vec<(usize, usize)>]) -> Vec<usize> {
    let mut dist = vec![usize::MAX; n];
    dist[root] = 0;
    let mut queue = VecDeque::new();
    queue.push_back(root);
    while let Some(u) = queue.pop_front() {
        for &(_, v) in &adj[u] {
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
    adj: &[Vec<(usize, usize)>],
    degree: &[usize],
) -> (Vec<usize>, Vec<usize>) {
    let mut dist = vec![usize::MAX; n];
    let mut pred = vec![usize::MAX; n];
    dist[root] = 0;
    pred[root] = root;
    let mut queue = VecDeque::new();
    queue.push_back(root);
    while let Some(u) = queue.pop_front() {
        for &(_, v) in &adj[u] {
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
    let mut on_ry: HashSet<usize> = HashSet::new();
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

/// Builds the prototype cycle of a cycle family from single shortest paths.
fn make_cf(
    r: usize,
    p: usize,
    q: usize,
    x: Option<usize>,
    m: usize,
    edge_id: &HashMap<(usize, usize), usize>,
    apsp: &Apsp,
) -> Cf {
    let mut prototype = Bits::zeros(m);
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
    bits: &mut Bits,
    edge_id: &HashMap<(usize, usize), usize>,
    apsp: &Apsp,
) {
    let mut v = target;
    while v != r {
        let u = apsp.pred[r][v];
        bits.set(edge_id[&ends(v, u)]);
        v = u;
    }
}

/// All edges spanned by a cycle family: every shortest path from `r` to `p` and
/// to `q`, plus the closing edges.
fn family_edges(
    cf: &Cf,
    n: usize,
    adj: &[Vec<(usize, usize)>],
    edge_id: &HashMap<(usize, usize), usize>,
    apsp: &Apsp,
) -> HashSet<usize> {
    let mut edges: HashSet<usize> = HashSet::new();
    collect_shortest_edges(cf.r, cf.p, n, adj, apsp, &mut edges);
    collect_shortest_edges(cf.r, cf.q, n, adj, apsp, &mut edges);
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
    n: usize,
    adj: &[Vec<(usize, usize)>],
    apsp: &Apsp,
    edges: &mut HashSet<usize>,
) {
    let mut visited = vec![false; n];
    let mut stack = vec![target];
    visited[target] = true;
    while let Some(v) = stack.pop() {
        if v == r {
            continue;
        }
        for &(e, u) in &adj[v] {
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

/// Reduces `cycle` modulo the span of a reduced-echelon `basis`.
fn reduce(cycle: &mut Bits, basis: &[(usize, Bits)]) {
    for (pivot, vector) in basis {
        if cycle.test(*pivot) {
            cycle.xor(vector);
        }
    }
}

/// Extends a reduced-echelon `basis` with `vector` if it is independent.
fn extend(basis: &mut Vec<(usize, Bits)>, mut vector: Bits) {
    reduce(&mut vector, basis);
    let Some(pivot) = vector.lowest_set() else {
        return;
    };
    for (_, existing) in basis.iter_mut() {
        if existing.test(pivot) {
            existing.xor(&vector);
        }
    }
    basis.push((pivot, vector));
}

/// Connected components of the cycle families under edge overlap.
fn edge_overlap_components(
    members: &[usize],
    edges: &HashMap<usize, HashSet<usize>>,
) -> Vec<Vec<usize>> {
    let k = members.len();
    let mut parent: Vec<usize> = (0..k).collect();

    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }

    for a in 0..k {
        for b in (a + 1)..k {
            if !edges[&members[a]].is_disjoint(&edges[&members[b]]) {
                let ra = find(&mut parent, a);
                let rb = find(&mut parent, b);
                if ra != rb {
                    parent[ra] = rb;
                }
            }
        }
    }

    let mut components: HashMap<usize, Vec<usize>> = HashMap::new();
    for (idx, &member) in members.iter().enumerate() {
        let root = find(&mut parent, idx);
        components.entry(root).or_default().push(member);
    }
    components.into_values().collect()
}

/// Orders a pair of vertices as `(min, max)` for edge lookup.
fn ends(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
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

    fn mol(sites: &[u32], bonds: &[(u32, u32, u32)]) -> Mol {
        Mol {
            sites: sites.iter().map(|&n| s(n)).collect(),
            bonds: bonds.iter().map(|&(id, _, _)| b(id)).collect(),
            endpoints: bonds.iter().map(|&(_, u, v)| (s(u), s(v))).collect(),
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

    #[test]
    fn empty_has_no_families() {
        let f = families(&empty());
        assert_eq!(f.len(), 0);
        assert!(f.is_empty());
    }

    #[test]
    fn single_site_has_no_families() {
        assert_eq!(families(&single()).len(), 0);
    }

    #[test]
    fn chain_has_no_families() {
        assert!(families(&chain()).is_empty());
    }

    #[test]
    fn triangle_has_one_family() {
        assert_eq!(families(&triangle()).len(), 1);
    }

    #[test]
    fn triangle_family_size_is_three() {
        let f = families(&triangle());
        assert_eq!(f.iter().next().unwrap().size(), 3);
    }

    #[test]
    fn square_has_one_family() {
        assert_eq!(families(&square()).len(), 1);
    }

    #[test]
    fn fused_has_two_families() {
        assert_eq!(families(&fused()).len(), 2);
    }

    #[test]
    fn spiro_has_two_families() {
        assert_eq!(families(&spiro()).len(), 2);
    }

    #[test]
    fn two_triangles_has_two_families() {
        assert_eq!(families(&two_triangles()).len(), 2);
    }

    #[test]
    fn bicyclo222_has_three_families() {
        let f = families(&bicyclo222());
        assert_eq!(f.len(), 3);
        assert!(f.iter().all(|fam| fam.size() == 6));
    }

    #[test]
    fn cube_has_six_families() {
        let f = families(&cube());
        assert_eq!(f.len(), 6);
        assert!(f.iter().all(|fam| fam.size() == 4));
    }

    #[test]
    fn k4_has_four_families() {
        let f = families(&k4());
        assert_eq!(f.len(), 4);
        assert!(f.iter().all(|fam| fam.size() == 3));
    }

    #[test]
    fn bridged_square_merges_interchangeable_rings() {
        let f = families(&bridged_square());
        assert_eq!(f.len(), 2);
        let mut sizes: Vec<usize> = f.iter().map(|fam| fam.size()).collect();
        sizes.sort_unstable();
        assert_eq!(sizes, vec![4, 5]);
    }

    #[test]
    fn bridged_square_size_five_family_spans_all_atoms() {
        let f = families(&bridged_square());
        let five = f.iter().find(|fam| fam.size() == 5).unwrap();
        assert_eq!(five.sites().len(), 6);
    }

    #[test]
    fn girth_is_smallest_family() {
        assert_eq!(families(&triangle()).girth(), Some(3));
        assert_eq!(families(&square()).girth(), Some(4));
        assert_eq!(families(&bridged_square()).girth(), Some(4));
    }

    #[test]
    fn acyclic_girth_is_none() {
        assert_eq!(families(&chain()).girth(), None);
    }

    #[test]
    fn triangle_each_site_in_one_family() {
        let f = families(&triangle());
        for site in [s(1), s(2), s(3)] {
            assert_eq!(f.of_site(site).count(), 1);
        }
    }

    #[test]
    fn fused_shared_bond_in_two_families() {
        let f = families(&fused());
        assert_eq!(f.of_bond(b(3)).count(), 2);
    }

    #[test]
    fn fused_shared_site_in_two_families() {
        let f = families(&fused());
        assert_eq!(f.of_site(s(3)).count(), 2);
        assert_eq!(f.of_site(s(4)).count(), 2);
    }

    #[test]
    fn of_site_unknown_is_empty() {
        assert_eq!(families(&triangle()).of_site(s(99)).count(), 0);
    }

    #[test]
    fn of_bond_bridge_is_empty() {
        let tadpole = mol(&[1, 2, 3, 4], &[(1, 1, 2), (2, 2, 3), (3, 1, 3), (4, 1, 4)]);
        assert_eq!(families(&tadpole).of_bond(b(4)).count(), 0);
    }

    #[test]
    fn of_bond_unknown_is_empty() {
        assert_eq!(families(&triangle()).of_bond(b(99)).count(), 0);
    }

    #[test]
    fn spiro_same_within_ring_only() {
        let f = families(&spiro());
        assert!(f.same(s(1), s(2)));
        assert!(!f.same(s(1), s(4)));
    }

    #[test]
    fn same_unknown_is_false() {
        assert!(!families(&triangle()).same(s(1), s(99)));
    }

    #[test]
    fn spiro_site_is_the_shared_atom() {
        let f = families(&spiro());
        assert_eq!(f.spiro_sites().collect::<Vec<_>>(), vec![s(3)]);
    }

    #[test]
    fn fused_has_no_spiro_sites() {
        assert_eq!(families(&fused()).spiro_sites().count(), 0);
    }

    #[test]
    fn bicyclo222_has_no_spiro_sites() {
        assert_eq!(families(&bicyclo222()).spiro_sites().count(), 0);
    }

    #[test]
    fn fusion_bond_is_the_shared_edge() {
        let f = families(&fused());
        assert_eq!(f.fusion_bonds().collect::<Vec<_>>(), vec![b(3)]);
    }

    #[test]
    fn spiro_has_no_fusion_bonds() {
        assert_eq!(families(&spiro()).fusion_bonds().count(), 0);
    }

    #[test]
    fn membership_matches_free_function() {
        let m = fused();
        let derived = families(&m).membership();
        let direct = super::super::membership(&m);

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
    fn fused_one_system() {
        assert_eq!(families(&fused()).systems().count(), 1);
    }

    #[test]
    fn spiro_one_system() {
        assert_eq!(families(&spiro()).systems().count(), 1);
    }

    #[test]
    fn two_triangles_two_systems() {
        let systems: Vec<Vec<SiteId>> = families(&two_triangles()).systems().collect();
        assert_eq!(systems.len(), 2);
        assert_eq!(systems[0], vec![s(1), s(2), s(3)]);
        assert_eq!(systems[1], vec![s(4), s(5), s(6)]);
    }

    #[test]
    fn systems_partition_ring_sites() {
        let m = fused();
        let f = families(&m);
        let from_systems: HashSet<SiteId> = f.systems().flatten().collect();
        let from_membership: HashSet<SiteId> = f.membership().sites().collect();
        assert_eq!(from_systems, from_membership);
    }

    #[test]
    fn matches_minimum_basis_when_non_degenerate() {
        for m in [triangle(), square(), fused(), spiro(), two_triangles()] {
            assert_eq!(families(&m).len(), super::super::rings(&m).len());
        }
    }

    #[test]
    fn exceeds_minimum_basis_when_degenerate() {
        assert!(families(&cube()).len() > super::super::rings(&cube()).len());
    }

    #[test]
    fn cube_is_deterministic() {
        let first: Vec<Vec<SiteId>> = families(&cube())
            .iter()
            .map(|f| f.sites().to_vec())
            .collect();
        let second: Vec<Vec<SiteId>> = families(&cube())
            .iter()
            .map(|f| f.sites().to_vec())
            .collect();
        assert_eq!(first, second);
    }

    #[test]
    fn output_is_independent_of_input_order() {
        let shuffled = mol(
            &[6, 4, 2, 5, 3, 1],
            &[
                (7, 6, 4),
                (3, 4, 3),
                (1, 2, 1),
                (5, 5, 3),
                (2, 3, 2),
                (6, 6, 5),
                (4, 4, 1),
            ],
        );
        let canonical: Vec<(usize, Vec<SiteId>)> = families(&fused())
            .iter()
            .map(|f| (f.size(), f.sites().to_vec()))
            .collect();
        let reordered: Vec<(usize, Vec<SiteId>)> = families(&shuffled)
            .iter()
            .map(|f| (f.size(), f.sites().to_vec()))
            .collect();
        assert_eq!(canonical, reordered);
    }
}
