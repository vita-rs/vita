use std::collections::{HashMap, HashSet, VecDeque};

use vita_core::{HasSites, SiteId};

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
/// Decomposes the cycle space into the fewest independent rings — exactly
/// `cycle_rank` of them — chosen to have the least total size. For nearly all
/// molecules this is the chemically intended ring perception; where several
/// equivalent bases exist (cage systems such as cubane), a canonical
/// site-ordered tie-break makes the choice deterministic.
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
    let words = m.div_ceil(64);

    let mut adj: Vec<Vec<(usize, usize)>> = vec![vec![]; n];
    for (e, &(_, lo, hi)) in rows.iter().enumerate() {
        adj[lo].push((e, hi));
        adj[hi].push((e, lo));
    }
    for a in adj.iter_mut() {
        a.sort_unstable_by_key(|&(_, nb)| nb);
    }

    let candidates = horton_candidates(n, words, &adj, &rows);
    let basis = minimum_basis(words, candidates);

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
    words: usize,
    adj: &[Vec<(usize, usize)>],
    rows: &[(BondId, usize, usize)],
) -> Vec<Vec<u64>> {
    let mut candidates: Vec<Vec<u64>> = Vec::new();

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

            let mut bits = vec![0u64; words];
            let mut x = lo;
            while x != root {
                let pb = pred_bond[x];
                bits[pb >> 6] ^= 1 << (pb & 63);
                x = pred_node[x];
            }
            let mut y = hi;
            while y != root {
                let pb = pred_bond[y];
                bits[pb >> 6] ^= 1 << (pb & 63);
                y = pred_node[y];
            }
            bits[e >> 6] ^= 1 << (e & 63);

            candidates.push(bits);
        }
    }

    candidates.sort_unstable_by(|a, b| {
        let pa: u32 = a.iter().map(|w| w.count_ones()).sum();
        let pb: u32 = b.iter().map(|w| w.count_ones()).sum();
        pa.cmp(&pb).then_with(|| a.cmp(b))
    });
    candidates.dedup();
    candidates
}

/// Greedy minimum-weight basis extraction over GF(2).
fn minimum_basis(words: usize, candidates: Vec<Vec<u64>>) -> Vec<Vec<u64>> {
    let mut pivots: HashMap<usize, Vec<u64>> = HashMap::new();
    let mut basis: Vec<Vec<u64>> = Vec::new();

    for cand in candidates {
        let mut residue = cand.clone();
        while let Some(p) = lowest_set_bit(&residue) {
            match pivots.get(&p) {
                Some(reducer) => {
                    for k in 0..words {
                        residue[k] ^= reducer[k];
                    }
                }
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

/// Index of the lowest set bit, or `None` if all bits are clear.
fn lowest_set_bit(bits: &[u64]) -> Option<usize> {
    bits.iter()
        .enumerate()
        .find(|&(_, &w)| w != 0)
        .map(|(k, &w)| k * 64 + w.trailing_zeros() as usize)
}

/// Walks an edge bit vector into a canonically ordered [`Ring`].
fn trace_ring(bits: &[u64], m: usize, rows: &[(BondId, usize, usize)], sites: &[SiteId]) -> Ring {
    let mut local: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
    for e in (0..m).filter(|&e| (bits[e >> 6] >> (e & 63)) & 1 == 1) {
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
