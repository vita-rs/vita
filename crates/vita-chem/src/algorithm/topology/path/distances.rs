use std::collections::VecDeque;

use vita_core::SiteId;

use crate::HasBonds;
use crate::algorithm::utils::{AdjacencyList, FxHashMap};

/// Matrix entry for a pair of sites with no connecting path.
const UNREACHABLE: u32 = u32::MAX;

/// The all-pairs topological distances of a molecule, in bond hops.
///
/// The distance between two sites is the least number of bonds on a path
/// between them, and zero from a site to itself. Sites in different connected
/// components have no path, so [`get`](Self::get) returns `None` for such a
/// pair. From the matrix follow the standard graph-distance metrics:
/// eccentricity, diameter, radius, center, peripheral sites, and the Wiener
/// index.
///
/// Obtain via [`distances`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DistanceMatrix {
    sites: Vec<SiteId>,
    matrix: Vec<u32>,
}

impl DistanceMatrix {
    /// The row and column index of `site`, or `None` if it is absent.
    fn position(&self, site: SiteId) -> Option<usize> {
        self.sites.binary_search(&site).ok()
    }

    /// The distance between the sites at indices `i` and `j`, or `None` if no
    /// path joins them.
    fn cell(&self, i: usize, j: usize) -> Option<usize> {
        match self.matrix[i * self.sites.len() + j] {
            UNREACHABLE => None,
            d => Some(d as usize),
        }
    }

    /// The eccentricity of every site in site order, or `None` if any site
    /// fails to reach the rest.
    fn eccentricities(&self) -> Option<Vec<usize>> {
        self.sites.iter().map(|&s| self.eccentricity(s)).collect()
    }

    /// Distance from `a` to `b` in bond hops.
    ///
    /// Returns `None` if either site is absent, or if the two lie in different
    /// connected components.
    pub fn get(&self, a: SiteId, b: SiteId) -> Option<usize> {
        let i = self.position(a)?;
        let j = self.position(b)?;
        self.cell(i, j)
    }

    /// Greatest distance from `s` to any site.
    ///
    /// Returns `None` if `s` is absent or cannot reach every site.
    pub fn eccentricity(&self, s: SiteId) -> Option<usize> {
        let i = self.position(s)?;
        let mut farthest = 0;
        for j in 0..self.sites.len() {
            farthest = farthest.max(self.cell(i, j)?);
        }
        Some(farthest)
    }

    /// Greatest eccentricity over all sites.
    ///
    /// Returns `None` if the molecule is empty or disconnected.
    pub fn diameter(&self) -> Option<usize> {
        self.eccentricities()?.into_iter().max()
    }

    /// Least eccentricity over all sites.
    ///
    /// Returns `None` if the molecule is empty or disconnected.
    pub fn radius(&self) -> Option<usize> {
        self.eccentricities()?.into_iter().min()
    }

    /// Sites whose eccentricity equals the radius, in ascending order.
    ///
    /// Empty when the molecule is disconnected.
    pub fn center(&self) -> impl Iterator<Item = SiteId> + '_ {
        let radius = self.radius();
        self.sites
            .iter()
            .copied()
            .filter(move |&s| radius.is_some_and(|r| self.eccentricity(s) == Some(r)))
    }

    /// Sites whose eccentricity equals the diameter, in ascending order.
    ///
    /// Empty when the molecule is disconnected.
    pub fn peripheral(&self) -> impl Iterator<Item = SiteId> + '_ {
        let diameter = self.diameter();
        self.sites
            .iter()
            .copied()
            .filter(move |&s| diameter.is_some_and(|d| self.eccentricity(s) == Some(d)))
    }

    /// The Wiener index: the sum of distances over every unordered pair of sites.
    ///
    /// Each pair counts once; pairs in different components are skipped. Zero
    /// when the molecule has fewer than two sites.
    pub fn wiener(&self) -> u64 {
        let n = self.sites.len();
        let mut sum = 0;
        for i in 0..n {
            for j in (i + 1)..n {
                if let Some(d) = self.cell(i, j) {
                    sum += d as u64;
                }
            }
        }
        sum
    }
}

/// All-pairs topological distances of a molecule.
///
/// Performs a breadth-first search from every site, filling one matrix row per
/// source. A pair in different components records no distance.
///
/// # Complexity
///
/// O(V · (V + E)) time and O(V²) space, over the molecule's `V` sites and `E`
/// bonds.
pub fn distances<M: HasBonds>(mol: &M) -> DistanceMatrix {
    let mut sites: Vec<SiteId> = mol.sites().collect();
    sites.sort_unstable();
    let n = sites.len();
    let index: FxHashMap<SiteId, usize> = sites.iter().enumerate().map(|(i, &s)| (s, i)).collect();

    let graph = AdjacencyList::build(
        n,
        mol.bonds().map(|bond| {
            let (a, b) = mol.bond_endpoints(bond);
            (0, index[&a], index[&b])
        }),
    );

    let mut matrix = vec![UNREACHABLE; n * n];
    let mut frontier = VecDeque::new();
    for source in 0..n {
        matrix[source * n + source] = 0;
        frontier.push_back(source);
        while let Some(u) = frontier.pop_front() {
            let d = matrix[source * n + u];
            for &(_, v) in graph.neighbors(u) {
                if matrix[source * n + v] == UNREACHABLE {
                    matrix[source * n + v] = d + 1;
                    frontier.push_back(v);
                }
            }
        }
    }

    DistanceMatrix { sites, matrix }
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

    fn reversed(m: &Mol) -> Mol {
        Mol {
            sites: m.sites.iter().rev().copied().collect(),
            bonds: m.bonds.iter().rev().copied().collect(),
            endpoints: m.endpoints.iter().rev().copied().collect(),
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

    fn path() -> Mol {
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

    fn disconnected() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            bonds: vec![b(1)],
            endpoints: vec![(s(1), s(2))],
        }
    }

    #[test]
    fn single_site_diameter_is_zero() {
        assert_eq!(distances(&single()).diameter(), Some(0));
    }

    #[test]
    fn single_site_radius_is_zero() {
        assert_eq!(distances(&single()).radius(), Some(0));
    }

    #[test]
    fn single_site_is_its_own_center() {
        assert_eq!(
            distances(&single()).center().collect::<Vec<_>>(),
            vec![s(1)]
        );
    }

    #[test]
    fn single_site_is_its_own_periphery() {
        assert_eq!(
            distances(&single()).peripheral().collect::<Vec<_>>(),
            vec![s(1)],
        );
    }

    #[test]
    fn single_site_has_a_zero_wiener_index() {
        assert_eq!(distances(&single()).wiener(), 0);
    }

    #[test]
    fn empty_molecule_has_no_diameter() {
        assert_eq!(distances(&empty()).diameter(), None);
    }

    #[test]
    fn empty_molecule_has_no_radius() {
        assert_eq!(distances(&empty()).radius(), None);
    }

    #[test]
    fn distance_from_a_site_to_itself_is_zero() {
        let d = distances(&path());
        assert_eq!(d.get(s(1), s(1)), Some(0));
        assert_eq!(d.get(s(2), s(2)), Some(0));
    }

    #[test]
    fn distance_between_adjacent_sites_is_one() {
        assert_eq!(distances(&path()).get(s(1), s(2)), Some(1));
    }

    #[test]
    fn distance_counts_the_bonds_along_a_path() {
        assert_eq!(distances(&path()).get(s(1), s(3)), Some(2));
    }

    #[test]
    fn eccentricity_is_the_distance_to_the_farthest_site() {
        let d = distances(&path());
        assert_eq!(d.eccentricity(s(1)), Some(2));
        assert_eq!(d.eccentricity(s(2)), Some(1));
        assert_eq!(d.eccentricity(s(3)), Some(2));
    }

    #[test]
    fn diameter_is_the_greatest_eccentricity() {
        assert_eq!(distances(&star()).diameter(), Some(2));
    }

    #[test]
    fn radius_is_the_least_eccentricity() {
        assert_eq!(distances(&star()).radius(), Some(1));
    }

    #[test]
    fn center_holds_the_least_eccentric_sites() {
        assert_eq!(distances(&star()).center().collect::<Vec<_>>(), vec![s(1)]);
    }

    #[test]
    fn periphery_holds_the_most_eccentric_sites() {
        assert_eq!(
            distances(&star()).peripheral().collect::<Vec<_>>(),
            vec![s(2), s(3), s(4)],
        );
    }

    #[test]
    fn wiener_index_sums_every_pairwise_distance() {
        assert_eq!(distances(&star()).wiener(), 9);
    }

    #[test]
    fn distance_between_separate_components_is_none() {
        assert_eq!(distances(&disconnected()).get(s(1), s(3)), None);
    }

    #[test]
    fn distance_to_an_absent_site_is_none() {
        let d = distances(&path());
        assert_eq!(d.get(s(1), s(99)), None);
        assert_eq!(d.get(s(99), s(1)), None);
    }

    #[test]
    fn eccentricity_of_a_site_that_cannot_reach_all_is_none() {
        assert_eq!(distances(&disconnected()).eccentricity(s(1)), None);
    }

    #[test]
    fn eccentricity_of_an_absent_site_is_none() {
        assert_eq!(distances(&path()).eccentricity(s(99)), None);
    }

    #[test]
    fn disconnected_molecule_has_no_diameter() {
        assert_eq!(distances(&disconnected()).diameter(), None);
    }

    #[test]
    fn disconnected_molecule_has_no_radius() {
        assert_eq!(distances(&disconnected()).radius(), None);
    }

    #[test]
    fn disconnected_molecule_has_an_empty_center() {
        assert_eq!(distances(&disconnected()).center().count(), 0);
    }

    #[test]
    fn disconnected_molecule_has_an_empty_periphery() {
        assert_eq!(distances(&disconnected()).peripheral().count(), 0);
    }

    #[test]
    fn wiener_index_skips_pairs_in_separate_components() {
        assert_eq!(distances(&disconnected()).wiener(), 1);
    }

    #[test]
    fn distance_around_a_cycle_takes_the_shorter_arc() {
        assert_eq!(distances(&triangle()).get(s(1), s(3)), Some(1));
    }

    #[test]
    fn every_site_of_a_cycle_is_central() {
        assert_eq!(
            distances(&triangle()).center().collect::<Vec<_>>(),
            vec![s(1), s(2), s(3)],
        );
    }

    #[test]
    fn center_and_periphery_split_a_path_into_middle_and_ends() {
        let d = distances(&path());
        assert_eq!(d.center().collect::<Vec<_>>(), vec![s(2)]);
        assert_eq!(d.peripheral().collect::<Vec<_>>(), vec![s(1), s(3)]);
    }

    #[test]
    fn distance_is_symmetric() {
        let connected = distances(&star());
        assert_eq!(connected.get(s(2), s(3)), connected.get(s(3), s(2)));
        assert_eq!(connected.get(s(1), s(4)), connected.get(s(4), s(1)));
        let split = distances(&disconnected());
        assert_eq!(split.get(s(1), s(3)), split.get(s(3), s(1)));
    }

    #[test]
    fn radius_never_exceeds_diameter() {
        for mol in [path(), triangle(), star()] {
            let d = distances(&mol);
            assert!(d.radius() <= d.diameter());
        }
    }

    #[test]
    fn metrics_are_independent_of_input_order() {
        let profile = |m: &Mol| {
            let d = distances(m);
            (
                d.diameter(),
                d.radius(),
                d.center().collect::<Vec<_>>(),
                d.peripheral().collect::<Vec<_>>(),
                d.wiener(),
            )
        };
        assert_eq!(profile(&star()), profile(&reversed(&star())));
        assert_eq!(profile(&path()), profile(&reversed(&path())));
    }
}
