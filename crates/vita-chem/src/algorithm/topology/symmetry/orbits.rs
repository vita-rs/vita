use std::collections::HashMap;

use vita_core::{HasSites, SiteId};

use crate::HasBonds;

/// The symmetry classes of a molecule.
///
/// Each class is a maximal set of interchangeable sites: a symmetry of the
/// molecular graph — a relabelling that leaves it unchanged — maps any member
/// of a class onto any other. An empty molecule has no classes.
///
/// Obtain via [`orbits`].
pub struct Orbits {
    groups: Vec<Vec<SiteId>>,
    index: HashMap<SiteId, usize>,
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

    /// Iterates all classes, each as a slice of site identifiers.
    pub fn iter(&self) -> impl Iterator<Item = &[SiteId]> + '_ {
        self.groups.iter().map(|g| g.as_slice())
    }

    /// Returns the class containing `site`.
    ///
    /// Returns `None` if `site` is not present in the molecule.
    pub fn get(&self, site: SiteId) -> Option<&[SiteId]> {
        let &g = self.index.get(&site)?;
        Some(&self.groups[g])
    }

    /// Returns `true` if `a` and `b` belong to the same symmetry class.
    ///
    /// Returns `false` if either site is absent from the molecule.
    pub fn same(&self, a: SiteId, b: SiteId) -> bool {
        match (self.index.get(&a), self.index.get(&b)) {
            (Some(ia), Some(ib)) => ia == ib,
            _ => false,
        }
    }
}

/// Symmetry classes of a molecule.
///
/// Refines the sites by 1-dimensional Weisfeiler–Leman colouring until no class
/// splits further: two sites stay together only while their neighbours' classes
/// match. Classes are ordered by their sites, which within each class are
/// ascending.
///
/// # Complexity
///
/// O(V · (V + E)) time.
pub fn orbits<M: HasBonds + HasSites>(mol: &M) -> Orbits {
    let sites: Vec<SiteId> = mol.sites().collect();
    let n = sites.len();
    if n == 0 {
        return Orbits {
            groups: Vec::new(),
            index: HashMap::new(),
        };
    }

    let pos: HashMap<SiteId, usize> = sites.iter().enumerate().map(|(i, &s)| (s, i)).collect();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for bond in mol.bonds() {
        let (a, b) = mol.bond_endpoints(bond);
        adj[pos[&a]].push(pos[&b]);
        adj[pos[&b]].push(pos[&a]);
    }

    let mut class = vec![0usize; n];
    let mut count = 1;
    loop {
        let mut ids: HashMap<(usize, Vec<usize>), usize> = HashMap::new();
        let mut next = vec![0usize; n];
        for v in 0..n {
            let mut neighbours: Vec<usize> = adj[v].iter().map(|&u| class[u]).collect();
            neighbours.sort_unstable();
            let id = ids.len();
            next[v] = *ids.entry((class[v], neighbours)).or_insert(id);
        }
        if ids.len() == count {
            break;
        }
        count = ids.len();
        class = next;
    }

    let mut groups: Vec<Vec<SiteId>> = vec![Vec::new(); count];
    for (v, &site) in sites.iter().enumerate() {
        groups[class[v]].push(site);
    }
    for group in &mut groups {
        group.sort_unstable();
    }
    groups.sort_unstable();

    let mut index: HashMap<SiteId, usize> = HashMap::new();
    for (i, group) in groups.iter().enumerate() {
        for &site in group {
            index.insert(site, i);
        }
    }

    Orbits { groups, index }
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

    fn isopentane() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5)],
            bonds: vec![b(1), b(2), b(3), b(4)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3)), (s(3), s(4)), (s(2), s(5))],
        }
    }

    #[test]
    fn empty_molecule_has_no_classes() {
        let orb = orbits(&empty());
        assert_eq!(orb.len(), 0);
        assert!(orb.is_empty());
    }

    #[test]
    fn single_site_is_one_class() {
        assert_eq!(orbits(&single()).len(), 1);
    }

    #[test]
    fn chain_has_two_classes() {
        assert_eq!(orbits(&chain()).len(), 2);
    }

    #[test]
    fn cycle_is_one_class() {
        assert_eq!(orbits(&triangle()).len(), 1);
    }

    #[test]
    fn branched_distinguishes_by_environment() {
        assert_eq!(orbits(&isopentane()).len(), 4);
    }

    #[test]
    fn classes_partition_all_sites() {
        let mol = isopentane();
        let orb = orbits(&mol);
        let via_iter: HashSet<SiteId> = orb.iter().flat_map(|c| c.iter().copied()).collect();
        let via_mol: HashSet<SiteId> = mol.sites().collect();
        assert_eq!(via_iter, via_mol);
        let total: usize = orb.iter().map(|c| c.len()).sum();
        assert_eq!(total, via_mol.len());
    }

    #[test]
    fn site_in_own_class() {
        let orb = orbits(&isopentane());
        for site in [s(1), s(2), s(3), s(4), s(5)] {
            assert!(orb.get(site).unwrap().contains(&site));
        }
    }

    #[test]
    fn unknown_site_returns_none() {
        assert!(orbits(&chain()).get(s(99)).is_none());
    }

    #[test]
    fn chain_ends_are_equivalent() {
        assert!(orbits(&chain()).same(s(1), s(3)));
    }

    #[test]
    fn chain_end_and_centre_are_not_equivalent() {
        assert!(!orbits(&chain()).same(s(1), s(2)));
    }

    #[test]
    fn equivalent_methyls_are_same() {
        assert!(orbits(&isopentane()).same(s(1), s(5)));
    }

    #[test]
    fn inequivalent_methyls_are_not_same() {
        assert!(!orbits(&isopentane()).same(s(1), s(4)));
    }

    #[test]
    fn self_same_is_true() {
        assert!(orbits(&chain()).same(s(2), s(2)));
    }

    #[test]
    fn unknown_site_not_same_as_known() {
        assert!(!orbits(&chain()).same(s(99), s(1)));
        assert!(!orbits(&chain()).same(s(1), s(99)));
    }

    #[test]
    fn same_is_symmetric() {
        let orb = orbits(&isopentane());
        assert_eq!(orb.same(s(1), s(5)), orb.same(s(5), s(1)));
        assert_eq!(orb.same(s(1), s(4)), orb.same(s(4), s(1)));
    }

    #[test]
    fn classes_are_independent_of_input_order() {
        let shuffled = Mol {
            sites: vec![s(5), s(3), s(1), s(4), s(2)],
            bonds: vec![b(4), b(2), b(1), b(3)],
            endpoints: vec![(s(2), s(5)), (s(2), s(3)), (s(1), s(2)), (s(3), s(4))],
        };
        let canonical: Vec<Vec<SiteId>> =
            orbits(&isopentane()).iter().map(|c| c.to_vec()).collect();
        let reordered: Vec<Vec<SiteId>> = orbits(&shuffled).iter().map(|c| c.to_vec()).collect();
        assert_eq!(canonical, reordered);
    }
}
