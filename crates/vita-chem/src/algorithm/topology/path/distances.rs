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
    mat: Vec<usize>,
}

impl DistanceMatrix {
    fn n(&self) -> usize {
        self.sites.len()
    }

    fn idx(&self, s: SiteId) -> Option<usize> {
        self.sites.binary_search(&s).ok()
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

    /// Sites whose eccentricity equals the radius, in ascending order.
    ///
    /// Yields nothing if the molecule is disconnected.
    pub fn center(&self) -> impl Iterator<Item = SiteId> + '_ {
        let r = self.radius();
        self.sites
            .iter()
            .copied()
            .filter(move |&s| r.is_some_and(|d| self.eccentricity(s) == Some(d)))
    }

    /// Sites whose eccentricity equals the diameter, in ascending order.
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
    let mut sites: Vec<SiteId> = mol.sites().collect();
    sites.sort_unstable();
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

    DistanceMatrix { sites, mat }
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
            endpoints: vec![(s(1), s(2)), (s(1), s(3)), (s(1), s(4))],
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
    fn self_distance_is_zero() {
        let dm = distances(&chain());
        assert_eq!(dm.get(s(1), s(1)), Some(0));
        assert_eq!(dm.get(s(2), s(2)), Some(0));
        assert_eq!(dm.get(s(3), s(3)), Some(0));
    }

    #[test]
    fn adjacent_distance_is_one() {
        assert_eq!(distances(&chain()).get(s(1), s(2)), Some(1));
        assert_eq!(distances(&star()).get(s(1), s(4)), Some(1));
    }

    #[test]
    fn non_adjacent_distance() {
        assert_eq!(distances(&chain()).get(s(1), s(3)), Some(2));
        assert_eq!(distances(&star()).get(s(2), s(3)), Some(2));
    }

    #[test]
    fn distance_is_symmetric() {
        let dm = distances(&chain());
        assert_eq!(dm.get(s(1), s(3)), dm.get(s(3), s(1)));
        let dm = distances(&star());
        assert_eq!(dm.get(s(2), s(4)), dm.get(s(4), s(2)));
    }

    #[test]
    fn disconnected_pair_returns_none() {
        let dm = distances(&two_components());
        assert_eq!(dm.get(s(1), s(3)), None);
        assert_eq!(dm.get(s(3), s(1)), None);
        assert_eq!(dm.get(s(3), s(2)), None);
    }

    #[test]
    fn unknown_site_returns_none() {
        let dm = distances(&chain());
        assert_eq!(dm.get(s(99), s(1)), None);
        assert_eq!(dm.get(s(1), s(99)), None);
        assert_eq!(dm.eccentricity(s(99)), None);
    }

    #[test]
    fn chain_eccentricities() {
        let dm = distances(&chain());
        assert_eq!(dm.eccentricity(s(1)), Some(2));
        assert_eq!(dm.eccentricity(s(2)), Some(1));
        assert_eq!(dm.eccentricity(s(3)), Some(2));
    }

    #[test]
    fn star_eccentricities() {
        let dm = distances(&star());
        assert_eq!(dm.eccentricity(s(1)), Some(1));
        assert_eq!(dm.eccentricity(s(2)), Some(2));
        assert_eq!(dm.eccentricity(s(3)), Some(2));
        assert_eq!(dm.eccentricity(s(4)), Some(2));
    }

    #[test]
    fn eccentricity_is_none_when_disconnected() {
        let dm = distances(&two_components());
        assert_eq!(dm.eccentricity(s(1)), None);
        assert_eq!(dm.eccentricity(s(3)), None);
    }

    #[test]
    fn chain_diameter() {
        assert_eq!(distances(&chain()).diameter(), Some(2));
    }

    #[test]
    fn pentane_diameter() {
        assert_eq!(distances(&pentane()).diameter(), Some(4));
    }

    #[test]
    fn disconnected_diameter_is_none() {
        assert_eq!(distances(&two_components()).diameter(), None);
    }

    #[test]
    fn chain_radius() {
        assert_eq!(distances(&chain()).radius(), Some(1));
    }

    #[test]
    fn triangle_radius() {
        assert_eq!(distances(&triangle()).radius(), Some(1));
    }

    #[test]
    fn disconnected_radius_is_none() {
        assert_eq!(distances(&two_components()).radius(), None);
    }

    #[test]
    fn chain_center() {
        let center: Vec<SiteId> = distances(&chain()).center().collect();
        assert_eq!(center, vec![s(2)]);
    }

    #[test]
    fn star_center() {
        let center: Vec<SiteId> = distances(&star()).center().collect();
        assert_eq!(center, vec![s(1)]);
    }

    #[test]
    fn triangle_all_sites_are_center() {
        let center: Vec<SiteId> = distances(&triangle()).center().collect();
        assert_eq!(center, vec![s(1), s(2), s(3)]);
    }

    #[test]
    fn disconnected_center_is_empty() {
        assert_eq!(distances(&two_components()).center().count(), 0);
    }

    #[test]
    fn chain_peripheral() {
        let peripheral: Vec<SiteId> = distances(&chain()).peripheral().collect();
        assert_eq!(peripheral, vec![s(1), s(3)]);
    }

    #[test]
    fn star_peripheral() {
        let peripheral: Vec<SiteId> = distances(&star()).peripheral().collect();
        assert_eq!(peripheral, vec![s(2), s(3), s(4)]);
    }

    #[test]
    fn metrics_are_independent_of_input_order() {
        let shuffled = Mol {
            sites: vec![s(4), s(3), s(2), s(1)],
            bonds: vec![b(3), b(2), b(1)],
            endpoints: vec![(s(1), s(4)), (s(1), s(3)), (s(1), s(2))],
        };
        let metrics = |m: &Mol| -> (Vec<SiteId>, Vec<SiteId>) {
            let d = distances(m);
            (d.center().collect(), d.peripheral().collect())
        };
        assert_eq!(metrics(&star()), metrics(&shuffled));
    }

    #[test]
    fn disconnected_peripheral_is_empty() {
        assert_eq!(distances(&two_components()).peripheral().count(), 0);
    }

    #[test]
    fn chain_wiener() {
        assert_eq!(distances(&chain()).wiener(), 4);
    }

    #[test]
    fn pentane_wiener() {
        assert_eq!(distances(&pentane()).wiener(), 20);
    }

    #[test]
    fn wiener_skips_disconnected_pairs() {
        assert_eq!(distances(&two_components()).wiener(), 1);
    }

    #[test]
    fn single_site_trivial_metrics() {
        let dm = distances(&single());
        assert_eq!(dm.get(s(1), s(1)), Some(0));
        assert_eq!(dm.eccentricity(s(1)), Some(0));
        assert_eq!(dm.diameter(), Some(0));
        assert_eq!(dm.radius(), Some(0));
        assert_eq!(dm.center().collect::<Vec<_>>(), vec![s(1)]);
        assert_eq!(dm.peripheral().collect::<Vec<_>>(), vec![s(1)]);
        assert_eq!(dm.wiener(), 0);
    }
}
