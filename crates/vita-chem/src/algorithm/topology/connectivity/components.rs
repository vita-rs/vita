use std::collections::HashMap;

use vita_core::{HasSites, SiteId};

use crate::HasBonds;

/// The connected components of a molecule.
///
/// Each component is a maximal set of sites that are mutually reachable
/// through bonds. Sites with no bonds form singleton components. An empty
/// molecule has no components.
///
/// Obtain via [`components`].
pub struct Components {
    groups: Vec<Vec<SiteId>>,
    index: HashMap<SiteId, usize>,
}

impl Components {
    /// Number of connected components.
    pub fn len(&self) -> usize {
        self.groups.len()
    }

    /// Returns `true` if the molecule contains no sites.
    pub fn is_empty(&self) -> bool {
        self.groups.is_empty()
    }

    /// Returns `true` if the molecule is a single connected component.
    pub fn is_connected(&self) -> bool {
        self.groups.len() == 1
    }

    /// Iterates all components, each as a slice of site identifiers.
    pub fn iter(&self) -> impl Iterator<Item = &[SiteId]> + '_ {
        self.groups.iter().map(|g| g.as_slice())
    }

    /// Returns the component containing `site`.
    ///
    /// Returns `None` if `site` is not present in the molecule.
    pub fn get(&self, site: SiteId) -> Option<&[SiteId]> {
        let &g = self.index.get(&site)?;
        Some(&self.groups[g])
    }

    /// Returns `true` if `a` and `b` belong to the same connected component.
    ///
    /// Returns `false` if either site is absent from the molecule.
    pub fn same(&self, a: SiteId, b: SiteId) -> bool {
        match (self.index.get(&a), self.index.get(&b)) {
            (Some(ia), Some(ib)) => ia == ib,
            _ => false,
        }
    }
}

/// Connected components of a molecule.
///
/// Returns every maximal set of mutually reachable sites. Sites with no bonds
/// form singleton components. The order of components follows `mol.sites()`;
/// the order of sites within each component is DFS discovery order.
///
/// # Complexity
///
/// O(V + E) time and space.
pub fn components<M: HasBonds + HasSites>(mol: &M) -> Components {
    let mut index: HashMap<SiteId, usize> = HashMap::new();
    let mut groups: Vec<Vec<SiteId>> = Vec::new();

    for start in mol.sites() {
        if index.contains_key(&start) {
            continue;
        }

        let g = groups.len();
        let mut group = Vec::new();
        let mut stack = vec![start];
        index.insert(start, g);

        while let Some(site) = stack.pop() {
            group.push(site);
            for nb in mol.neighbors(site) {
                if let std::collections::hash_map::Entry::Vacant(e) = index.entry(nb) {
                    e.insert(g);
                    stack.push(nb);
                }
            }
        }

        groups.push(group);
    }

    Components { groups, index }
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

    fn two_components() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            bonds: vec![b(1)],
            endpoints: vec![(s(1), s(2))],
        }
    }

    fn three_isolated() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            bonds: vec![],
            endpoints: vec![],
        }
    }

    #[test]
    fn empty_molecule_has_no_components() {
        let cmps = components(&empty());
        assert_eq!(cmps.len(), 0);
        assert!(cmps.is_empty());
        assert!(!cmps.is_connected());
    }

    #[test]
    fn single_site_is_one_component() {
        assert_eq!(components(&single()).len(), 1);
    }

    #[test]
    fn chain_is_one_component() {
        assert_eq!(components(&chain()).len(), 1);
    }

    #[test]
    fn two_components_has_two() {
        assert_eq!(components(&two_components()).len(), 2);
    }

    #[test]
    fn three_isolated_each_own_component() {
        assert_eq!(components(&three_isolated()).len(), 3);
    }

    #[test]
    fn chain_is_connected() {
        assert!(components(&chain()).is_connected());
    }

    #[test]
    fn single_site_is_connected() {
        assert!(components(&single()).is_connected());
    }

    #[test]
    fn two_components_is_not_connected() {
        assert!(!components(&two_components()).is_connected());
    }

    #[test]
    fn components_partition_all_sites() {
        let mol = two_components();
        let cmps = components(&mol);
        let via_iter: HashSet<SiteId> = cmps.iter().flat_map(|c| c.iter().copied()).collect();
        let via_mol: HashSet<SiteId> = mol.sites().collect();
        assert_eq!(via_iter, via_mol);
        let total: usize = cmps.iter().map(|c| c.len()).sum();
        assert_eq!(total, via_mol.len());
    }

    #[test]
    fn two_components_sizes() {
        let cmps = components(&two_components());
        let mut sizes: Vec<usize> = cmps.iter().map(|c| c.len()).collect();
        sizes.sort();
        assert_eq!(sizes, vec![1, 2]);
        assert_eq!(cmps.iter().count(), cmps.len());
    }

    #[test]
    fn site_in_own_component() {
        let cmps = components(&chain());
        for site in [s(1), s(2), s(3)] {
            assert!(cmps.get(site).unwrap().contains(&site));
        }
    }

    #[test]
    fn chain_connected_sites_share_component() {
        let cmps = components(&chain());
        assert_eq!(cmps.get(s(1)), cmps.get(s(3)));
    }

    #[test]
    fn isolated_site_has_singleton_component() {
        let cmps = components(&two_components());
        assert_eq!(cmps.get(s(3)).unwrap(), &[s(3)]);
    }

    #[test]
    fn unknown_site_returns_none() {
        assert!(components(&chain()).get(s(99)).is_none());
    }

    #[test]
    fn site_component_is_enumerated() {
        let mol = two_components();
        let cmps = components(&mol);
        for site in mol.sites() {
            let via_get = cmps.get(site).unwrap();
            assert!(cmps.iter().any(|c| c == via_get));
        }
    }

    #[test]
    fn self_same_is_true() {
        let cmps = components(&chain());
        assert!(cmps.same(s(1), s(1)));
        assert!(cmps.same(s(2), s(2)));
    }

    #[test]
    fn directly_bonded_sites_are_same() {
        assert!(components(&chain()).same(s(1), s(2)));
    }

    #[test]
    fn transitively_connected_sites_are_same() {
        assert!(components(&chain()).same(s(1), s(3)));
    }

    #[test]
    fn disconnected_sites_are_not_same() {
        assert!(!components(&two_components()).same(s(1), s(3)));
    }

    #[test]
    fn unknown_site_not_same_as_known() {
        assert!(!components(&chain()).same(s(99), s(1)));
        assert!(!components(&chain()).same(s(1), s(99)));
    }

    #[test]
    fn unknown_sites_are_not_same_as_each_other() {
        assert!(!components(&chain()).same(s(99), s(100)));
    }

    #[test]
    fn same_is_symmetric() {
        let cmps = components(&two_components());
        assert_eq!(cmps.same(s(1), s(2)), cmps.same(s(2), s(1)));
        assert_eq!(cmps.same(s(1), s(3)), cmps.same(s(3), s(1)));
    }
}
