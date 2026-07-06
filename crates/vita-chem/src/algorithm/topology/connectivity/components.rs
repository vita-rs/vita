use vita_core::SiteId;

use crate::HasBonds;
use crate::algorithm::utils::{FxHashSet, SortedMap};

/// A connected component: a maximal set of sites mutually reachable through
/// bonds.
///
/// Obtain from [`Components::iter`] or [`Components::get`].
pub struct Component {
    sites: Vec<SiteId>,
}

impl Component {
    /// Number of sites in the component.
    pub fn len(&self) -> usize {
        self.sites.len()
    }

    /// Returns `true` if the component contains no sites. Always `false` — a
    /// component is non-empty — but provided alongside [`len`](Self::len).
    pub fn is_empty(&self) -> bool {
        self.sites.is_empty()
    }

    /// Returns `true` if `site` lies in this component.
    pub fn contains(&self, site: SiteId) -> bool {
        self.sites.binary_search(&site).is_ok()
    }

    /// Iterates the component's sites in ascending order.
    pub fn iter(&self) -> impl Iterator<Item = SiteId> + '_ {
        self.sites.iter().copied()
    }
}

/// The connected components of a molecule.
///
/// Each component is a maximal set of sites mutually reachable through bonds.
/// Sites with no bonds form singleton components; an empty molecule has none.
///
/// Obtain via [`components`].
pub struct Components {
    groups: Vec<Component>,
    index: SortedMap<SiteId, usize>,
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

    /// Iterates the components, ordered by their sites.
    pub fn iter(&self) -> impl Iterator<Item = &Component> + '_ {
        self.groups.iter()
    }

    /// Returns the component containing `site`.
    ///
    /// Returns `None` if `site` is absent from the molecule.
    pub fn get(&self, site: SiteId) -> Option<&Component> {
        let &group = self.index.get(&site)?;
        Some(&self.groups[group])
    }

    /// Returns `true` if `a` and `b` lie in the same connected component.
    ///
    /// Returns `false` if either site is absent from the molecule.
    pub fn same(&self, a: SiteId, b: SiteId) -> bool {
        match (self.index.get(&a), self.index.get(&b)) {
            (Some(ga), Some(gb)) => ga == gb,
            _ => false,
        }
    }
}

/// Connected components of a molecule.
///
/// Returns every maximal set of mutually reachable sites, ordered by their
/// sites and ascending within each.
///
/// # Complexity
///
/// O(V · log V + E) time and O(V) space, over the molecule's `V` sites and `E`
/// bonds, assuming [`neighbors`](HasBonds::neighbors) runs in O(degree); the
/// log factor orders the components and their sites canonically.
pub fn components<M: HasBonds>(mol: &M) -> Components {
    let mut visited = FxHashSet::default();
    let mut groups: Vec<Vec<SiteId>> = Vec::new();

    for start in mol.sites() {
        if !visited.insert(start) {
            continue;
        }
        let mut group = vec![start];
        let mut stack = vec![start];
        while let Some(site) = stack.pop() {
            for neighbor in mol.neighbors(site) {
                if visited.insert(neighbor) {
                    group.push(neighbor);
                    stack.push(neighbor);
                }
            }
        }
        group.sort_unstable();
        groups.push(group);
    }

    groups.sort_unstable();
    let index = SortedMap::from_pairs(
        groups
            .iter()
            .enumerate()
            .flat_map(|(g, group)| group.iter().map(move |&site| (site, g))),
    );
    let groups = groups
        .into_iter()
        .map(|sites| Component { sites })
        .collect();

    Components { groups, index }
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

    fn isolated() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
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

    #[test]
    fn empty_molecule_has_no_components() {
        let c = components(&empty());
        assert_eq!(c.len(), 0);
        assert!(c.is_empty());
        assert!(!c.is_connected());
    }

    #[test]
    fn single_site_is_one_component() {
        assert_eq!(components(&single()).len(), 1);
    }

    #[test]
    fn isolated_sites_are_each_their_own_component() {
        assert_eq!(components(&isolated()).len(), 3);
    }

    #[test]
    fn connected_molecule_is_a_single_component() {
        assert_eq!(components(&chain()).len(), 1);
    }

    #[test]
    fn adjacent_sites_share_a_component() {
        assert!(components(&chain()).same(s(1), s(2)));
    }

    #[test]
    fn transitively_connected_sites_share_a_component() {
        assert!(components(&chain()).same(s(1), s(3)));
    }

    #[test]
    fn a_site_shares_a_component_with_itself() {
        assert!(components(&chain()).same(s(1), s(1)));
    }

    #[test]
    fn disconnected_pieces_are_separate_components() {
        assert_eq!(components(&two_components()).len(), 2);
    }

    #[test]
    fn sites_in_different_components_are_not_the_same() {
        assert!(!components(&two_components()).same(s(1), s(3)));
    }

    #[test]
    fn same_is_false_when_a_site_is_absent() {
        let c = components(&chain());
        assert!(!c.same(s(1), s(99)));
        assert!(!c.same(s(99), s(1)));
    }

    #[test]
    fn single_component_molecule_is_connected() {
        assert!(components(&chain()).is_connected());
    }

    #[test]
    fn multi_component_molecule_is_not_connected() {
        assert!(!components(&two_components()).is_connected());
    }

    #[test]
    fn get_returns_the_component_containing_the_site() {
        let c = components(&two_components());
        let pair = c.get(s(1)).unwrap();
        assert_eq!(pair.len(), 2);
        assert!(pair.contains(s(1)) && pair.contains(s(2)));
        let lone = c.get(s(3)).unwrap();
        assert_eq!(lone.len(), 1);
        assert!(lone.contains(s(3)));
    }

    #[test]
    fn get_of_an_absent_site_is_none() {
        assert!(components(&chain()).get(s(99)).is_none());
    }

    #[test]
    fn iter_yields_components_ordered_by_their_sites() {
        let c = components(&two_components());
        let groups: Vec<Vec<SiteId>> = c.iter().map(|group| group.iter().collect()).collect();
        assert_eq!(groups, vec![vec![s(1), s(2)], vec![s(3)]]);
    }

    #[test]
    fn components_partition_every_site() {
        let mol = two_components();
        let c = components(&mol);
        let total: usize = c.iter().map(|group| group.len()).sum();
        assert_eq!(total, mol.sites().count());
        for site in mol.sites() {
            assert!(c.get(site).is_some());
        }
    }

    #[test]
    fn same_is_symmetric() {
        let c = components(&two_components());
        for (a, b) in [(s(1), s(2)), (s(1), s(3)), (s(2), s(3)), (s(1), s(99))] {
            assert_eq!(c.same(a, b), c.same(b, a));
        }
    }

    #[test]
    fn output_is_independent_of_input_order() {
        let canonical = Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5)],
            bonds: vec![b(1), b(2), b(3)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3)), (s(4), s(5))],
        };
        let shuffled = Mol {
            sites: vec![s(5), s(3), s(1), s(4), s(2)],
            bonds: vec![b(3), b(1), b(2)],
            endpoints: vec![(s(4), s(5)), (s(1), s(2)), (s(2), s(3))],
        };
        let groups = |m: &Mol| -> Vec<Vec<SiteId>> {
            components(m)
                .iter()
                .map(|group| group.iter().collect())
                .collect()
        };
        assert_eq!(groups(&canonical), groups(&shuffled));
    }
}
