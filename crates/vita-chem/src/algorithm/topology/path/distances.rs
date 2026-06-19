use std::collections::{HashMap, VecDeque};

use vita_core::{HasSites, SiteId};

use crate::HasBonds;

/// All-pairs topological distance matrix for a molecule.
///
/// Distances are measured in bond hops (minimum number of bonds to traverse).
/// The distance from a site to itself is zero. Pairs in different connected
/// components have no finite distance and yield `None` from [`Self::get`].
///
/// Obtain via [`distances`].
pub struct DistanceMatrix {
    sites: Vec<SiteId>,
    index: HashMap<SiteId, usize>,
    mat: Vec<usize>,
}

impl DistanceMatrix {
    fn n(&self) -> usize {
        self.sites.len()
    }

    fn idx(&self, s: SiteId) -> Option<usize> {
        self.index.get(&s).copied()
    }

    fn raw(&self, i: usize, j: usize) -> usize {
        self.mat[i * self.n() + j]
    }

    /// Distance from `a` to `b` in bond hops.
    ///
    /// Returns `None` if either site is absent or if `a` and `b` lie in
    /// different connected components.
    pub fn get(&self, a: SiteId, b: SiteId) -> Option<usize> {
        let i = self.idx(a)?;
        let j = self.idx(b)?;
        let d = self.raw(i, j);
        if d == usize::MAX { None } else { Some(d) }
    }

    /// Greatest distance from `s` to any site.
    ///
    /// Returns `None` if `s` is absent or cannot reach every site in the
    /// molecule.
    pub fn eccentricity(&self, s: SiteId) -> Option<usize> {
        let i = self.idx(s)?;
        let n = self.n();
        let mut max_d = 0usize;
        for j in 0..n {
            let d = self.raw(i, j);
            if d == usize::MAX {
                return None;
            }
            if d > max_d {
                max_d = d;
            }
        }
        Some(max_d)
    }

    /// Greatest eccentricity over all sites.
    ///
    /// Returns `None` if the molecule is disconnected or contains no sites.
    pub fn diameter(&self) -> Option<usize> {
        let mut max_ecc = None::<usize>;
        for &s in &self.sites {
            match self.eccentricity(s) {
                None => return None,
                Some(e) => max_ecc = Some(max_ecc.map_or(e, |m| m.max(e))),
            }
        }
        max_ecc
    }

    /// Least eccentricity over all sites.
    ///
    /// Returns `None` if the molecule is disconnected or contains no sites.
    pub fn radius(&self) -> Option<usize> {
        let mut min_ecc = None::<usize>;
        for &s in &self.sites {
            match self.eccentricity(s) {
                None => return None,
                Some(e) => min_ecc = Some(min_ecc.map_or(e, |m| m.min(e))),
            }
        }
        min_ecc
    }

    /// Sites whose eccentricity equals the radius.
    ///
    /// Yields nothing if the molecule is disconnected.
    pub fn center(&self) -> impl Iterator<Item = SiteId> + '_ {
        let r = self.radius();
        self.sites
            .iter()
            .copied()
            .filter(move |&s| r.is_some_and(|d| self.eccentricity(s) == Some(d)))
    }

    /// Sites whose eccentricity equals the diameter.
    ///
    /// Yields nothing if the molecule is disconnected.
    pub fn peripheral(&self) -> impl Iterator<Item = SiteId> + '_ {
        let d = self.diameter();
        self.sites
            .iter()
            .copied()
            .filter(move |&s| d.is_some_and(|diam| self.eccentricity(s) == Some(diam)))
    }

    /// Sum of all finite pairwise distances (Wiener index).
    ///
    /// Each unordered pair is counted once. Disconnected pairs contribute
    /// nothing. Returns zero for a molecule with fewer than two sites.
    pub fn wiener(&self) -> u64 {
        let n = self.n();
        let mut sum = 0u64;
        for i in 0..n {
            for j in (i + 1)..n {
                let d = self.raw(i, j);
                if d != usize::MAX {
                    sum += d as u64;
                }
            }
        }
        sum
    }
}

/// All-pairs shortest topological distances.
///
/// Runs BFS from every site to build the complete distance matrix. Sites in
/// different connected components record no finite distance between them;
/// [`DistanceMatrix::get`] returns `None` for those pairs.
///
/// # Complexity
///
/// O(V(V + E)) time, O(V²) space.
pub fn distances<M: HasBonds + HasSites>(mol: &M) -> DistanceMatrix {
    let sites: Vec<SiteId> = mol.sites().collect();
    let n = sites.len();
    let index: HashMap<SiteId, usize> = sites.iter().enumerate().map(|(i, &s)| (s, i)).collect();

    let mut mat = vec![usize::MAX; n * n];

    for i in 0..n {
        mat[i * n + i] = 0;
    }

    for (i, &start) in sites.iter().enumerate() {
        let mut queue = VecDeque::new();
        queue.push_back(start);

        while let Some(site) = queue.pop_front() {
            let si = index[&site];
            let d = mat[i * n + si];
            for nb in mol.neighbors(site) {
                let ni = index[&nb];
                if mat[i * n + ni] == usize::MAX {
                    mat[i * n + ni] = d + 1;
                    queue.push_back(nb);
                }
            }
        }
    }

    DistanceMatrix { sites, index, mat }
}
