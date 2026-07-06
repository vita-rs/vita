use vita_core::SiteId;

use crate::HasBonds;
use crate::algorithm::utils::{FxHashMap, SortedMap, labeling};

/// A set of sites interchangeable by a symmetry of the molecule.
///
/// A relabeling of the molecule onto itself carries any member of the orbit onto
/// any other. *Which* symmetry — of the bare bond skeleton, or of the coloured
/// molecule — is fixed by the function that produced the orbit, not by this type.
///
/// Obtain from [`Orbits`] (skeleton symmetry) or from a
/// [`Canonical`](crate::canonical::Canonical) labeling (colour-preserving
/// symmetry).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Orbit {
    sites: Vec<SiteId>,
}

impl Orbit {
    /// Wraps an ascending set of sites; the sites must be sorted so that
    /// [`contains`](Self::contains) can bisect.
    pub(crate) fn new(sites: Vec<SiteId>) -> Self {
        Orbit { sites }
    }

    /// Number of sites in the orbit.
    pub fn len(&self) -> usize {
        self.sites.len()
    }

    /// Returns `true` if the orbit contains no sites. Always `false` — an orbit
    /// is non-empty — but provided alongside [`len`](Self::len).
    pub fn is_empty(&self) -> bool {
        self.sites.is_empty()
    }

    /// Returns `true` if `site` lies in this orbit.
    pub fn contains(&self, site: SiteId) -> bool {
        self.sites.binary_search(&site).is_ok()
    }

    /// Iterates the orbit's sites in ascending order.
    pub fn iter(&self) -> impl Iterator<Item = SiteId> + '_ {
        self.sites.iter().copied()
    }
}

/// The symmetry classes of a molecule's sites.
///
/// Each class is an [`Orbit`] — a maximal set of sites a graph automorphism can
/// interchange. The classification is topological, blind to the elements at the
/// sites. An empty molecule has no classes.
///
/// Obtain via [`orbits`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Orbits {
    groups: Vec<Orbit>,
    index: SortedMap<SiteId, usize>,
}

impl Orbits {
    /// Number of symmetry classes.
    pub fn len(&self) -> usize {
        self.groups.len()
    }

    /// Returns `true` if the molecule contains no sites.
    pub fn is_empty(&self) -> bool {
        self.groups.is_empty()
    }

    /// Iterates the classes, ordered by their sites, ascending within each.
    pub fn iter(&self) -> impl Iterator<Item = &Orbit> + '_ {
        self.groups.iter()
    }

    /// Returns the class containing `site`.
    ///
    /// Returns `None` if `site` is absent from the molecule.
    pub fn get(&self, site: SiteId) -> Option<&Orbit> {
        self.index.get(&site).map(|&group| &self.groups[group])
    }

    /// Returns `true` if `a` and `b` are interchangeable by a symmetry of the
    /// molecule.
    ///
    /// Returns `false` if either site is absent from the molecule.
    pub fn same(&self, a: SiteId, b: SiteId) -> bool {
        match (self.index.get(&a), self.index.get(&b)) {
            (Some(a), Some(b)) => a == b,
            _ => false,
        }
    }
}

/// Symmetry classes of a molecule.
///
/// Groups the sites into orbits under the automorphisms of the bond skeleton:
/// two sites share a class exactly when some relabeling of the molecule onto
/// itself maps one to the other. Colors play no part — the classification is
/// purely topological. Color refinement settles the sites by their
/// neighborhoods, and individualization resolves any class symmetry leaves
/// unsplit. Classes are ordered by their sites, ascending within each.
///
/// # Complexity
///
/// O(V · (V + E) · log V) time per refinement and O(V + E) space, over the
/// molecule's `V` sites and `E` bonds; near-linear in practice, exponential in
/// the worst case under individualization backtracking.
pub fn orbits<M: HasBonds>(mol: &M) -> Orbits {
    let sites: Vec<SiteId> = mol.sites().collect();
    let n = sites.len();
    let pos: FxHashMap<SiteId, usize> = sites.iter().enumerate().map(|(i, &s)| (s, i)).collect();

    let mut adjacency: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n];
    for bond in mol.bonds() {
        let (a, b) = mol.bond_endpoints(bond);
        adjacency[pos[&a]].push((pos[&b], 0));
        adjacency[pos[&b]].push((pos[&a], 0));
    }

    let seed = vec![0; n];
    let classes = labeling(&adjacency, &seed).orbits().to_vec();

    let mut groups: Vec<Vec<SiteId>> = vec![Vec::new(); n];
    for (vertex, &site) in sites.iter().enumerate() {
        groups[classes[vertex]].push(site);
    }
    groups.retain(|group| !group.is_empty());
    for group in &mut groups {
        group.sort_unstable();
    }
    groups.sort_unstable();

    let index = SortedMap::from_pairs(
        groups
            .iter()
            .enumerate()
            .flat_map(|(g, group)| group.iter().map(move |&site| (site, g))),
    );
    let groups = groups.into_iter().map(Orbit::new).collect();

    Orbits { groups, index }
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

    fn isolated_pair() -> Mol {
        Mol {
            sites: vec![s(1), s(2)],
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

    fn star() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4)],
            bonds: vec![b(1), b(2), b(3)],
            endpoints: vec![(s(1), s(2)), (s(1), s(3)), (s(1), s(4))],
        }
    }

    fn triangle() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            bonds: vec![b(1), b(2), b(3)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3)), (s(1), s(3))],
        }
    }

    #[test]
    fn empty_molecule_has_no_classes() {
        let orbits = orbits(&empty());
        assert_eq!(orbits.len(), 0);
        assert!(orbits.is_empty());
    }

    #[test]
    fn a_single_site_is_one_class() {
        let orbits = orbits(&single());
        assert_eq!(orbits.len(), 1);
        assert!(!orbits.is_empty());
        let class = orbits.get(s(1)).unwrap();
        assert_eq!(class.len(), 1);
        assert!(class.contains(s(1)));
    }

    #[test]
    fn symmetric_sites_are_in_the_same_class() {
        assert!(orbits(&path()).same(s(1), s(3)));
    }

    #[test]
    fn get_returns_the_class_of_a_site() {
        let orbits = orbits(&path());
        let class = orbits.get(s(1)).unwrap();
        assert_eq!(class.len(), 2);
        assert!(class.contains(s(1)) && class.contains(s(3)));
    }

    #[test]
    fn iter_lists_the_classes_in_site_order() {
        let orbits = orbits(&star());
        let classes: Vec<Vec<SiteId>> = orbits.iter().map(|class| class.iter().collect()).collect();
        assert_eq!(classes, vec![vec![s(1)], vec![s(2), s(3), s(4)]]);
    }

    #[test]
    fn asymmetric_sites_are_in_different_classes() {
        assert!(!orbits(&path()).same(s(1), s(2)));
    }

    #[test]
    fn same_is_false_for_an_absent_site() {
        assert!(!orbits(&path()).same(s(1), s(99)));
    }

    #[test]
    fn get_is_none_for_an_absent_site() {
        assert!(orbits(&path()).get(s(99)).is_none());
    }

    #[test]
    fn a_symmetric_ring_is_a_single_class() {
        let orbits = orbits(&triangle());
        assert_eq!(orbits.len(), 1);
        assert!(orbits.same(s(1), s(2)));
        assert!(orbits.same(s(2), s(3)));
    }

    #[test]
    fn isolated_sites_are_interchangeable() {
        let orbits = orbits(&isolated_pair());
        assert_eq!(orbits.len(), 1);
        assert!(orbits.same(s(1), s(2)));
    }

    #[test]
    fn an_asymmetric_site_is_a_singleton_class() {
        let orbits = orbits(&path());
        let class = orbits.get(s(2)).unwrap();
        assert_eq!(class.len(), 1);
        assert!(class.contains(s(2)));
    }

    #[test]
    fn classes_partition_every_site() {
        let orbits = orbits(&star());
        assert_eq!(orbits.len(), 2);
        let mut sites: Vec<SiteId> = orbits.iter().flat_map(|class| class.iter()).collect();
        sites.sort_unstable();
        assert_eq!(sites, vec![s(1), s(2), s(3), s(4)]);
    }

    #[test]
    fn orbits_are_independent_of_input_order() {
        let reordered = Mol {
            sites: vec![s(3), s(1), s(2)],
            bonds: vec![b(2), b(1)],
            endpoints: vec![(s(2), s(3)), (s(1), s(2))],
        };
        let classes = |mol: &Mol| -> Vec<Vec<SiteId>> {
            orbits(mol)
                .iter()
                .map(|class| class.iter().collect())
                .collect()
        };
        assert_eq!(classes(&path()), classes(&reordered));
    }
}
