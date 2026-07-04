use vita_core::SiteId;

use super::{Indexed, index};
use crate::algorithm::utils::{SortedMap, embeddings};
use crate::{BondId, HasBonds};

/// A subgraph match: an injection of a pattern's sites onto a target's under
/// which every pattern bond has a counterpart in the target.
///
/// Maps each pattern site to the target site it stands for.
///
/// Obtain via [`matches`](matches()).
pub struct Mapping {
    pairs: SortedMap<SiteId, SiteId>,
}

impl Mapping {
    /// Number of matched sites.
    pub fn len(&self) -> usize {
        self.pairs.len()
    }

    /// Returns `true` if the match maps no sites.
    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    /// Returns the target site `pattern_site` is matched to.
    ///
    /// Returns `None` if `pattern_site` is absent from the pattern.
    pub fn get(&self, pattern_site: SiteId) -> Option<SiteId> {
        self.pairs.get(&pattern_site).copied()
    }

    /// Iterates the matched `(pattern site, target site)` pairs, ordered by
    /// pattern site.
    pub fn iter(&self) -> impl Iterator<Item = (SiteId, SiteId)> + '_ {
        self.pairs
            .iter()
            .map(|(&pattern, &target)| (pattern, target))
    }
}

/// Finds every way `pattern` occurs as a subgraph of `target`.
///
/// The match is the caller's to define: `site_match` decides which pattern site
/// may stand for which target site, and `bond_match` which pattern bond for which
/// target bond — pass element and bond-order equality to match by constitution,
/// loosen either to honour a query's own rules, and the library imposes no
/// default. Each match injects every pattern site onto a distinct target site so
/// that every pattern bond meets a `bond_match`-ing target bond; the target may
/// bear further bonds among the matched sites, so the match is a subgraph, not an
/// induced one.
///
/// The matches stream lazily, so the first, a count, or all of them each cost
/// only the search they need.
///
/// # Complexity
///
/// O(T^P) time in the worst case and O(P · T) auxiliary space, over a pattern of
/// `P` sites and a target of `T` sites; near-linear in practice for the
/// connected patterns chemistry poses.
pub fn matches<P, T>(
    pattern: &P,
    target: &T,
    site_match: impl Fn(SiteId, SiteId) -> bool,
    bond_match: impl Fn(BondId, BondId) -> bool,
) -> impl Iterator<Item = Mapping>
where
    P: HasBonds,
    T: HasBonds,
{
    let Indexed {
        sites: pattern_sites,
        bonds: pattern_bonds,
        adjacency: pattern_adjacency,
    } = index(pattern);
    let Indexed {
        sites: target_sites,
        bonds: target_bonds,
        adjacency: target_adjacency,
    } = index(target);
    let compat_pattern_sites = pattern_sites.clone();
    let compat_target_sites = target_sites.clone();

    embeddings(
        &pattern_adjacency,
        target_adjacency,
        move |p, t| site_match(compat_pattern_sites[p], compat_target_sites[t]),
        move |pe, te| bond_match(pattern_bonds[pe], target_bonds[te]),
    )
    .map(move |mapping| Mapping {
        pairs: SortedMap::from_pairs(
            mapping
                .iter()
                .enumerate()
                .map(|(p, &t)| (pattern_sites[p], target_sites[t])),
        ),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::{Element, HasElements, HasSites};

    use crate::{BondOrder, HasBondOrders};

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn b(n: u32) -> BondId {
        BondId::new(n).unwrap()
    }

    fn elem(symbol: &str) -> Element {
        Element::from_symbol(symbol).unwrap()
    }

    struct Mol {
        sites: Vec<SiteId>,
        elements: Vec<Element>,
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
        orders: Vec<BondOrder>,
    }

    impl HasSites for Mol {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }

    impl HasElements for Mol {
        fn element(&self, site: SiteId) -> Element {
            let i = self.sites.iter().position(|&x| x == site).unwrap();
            self.elements[i]
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

    impl HasBondOrders for Mol {
        fn bond_order(&self, bond: BondId) -> BondOrder {
            let i = self.bonds.iter().position(|&x| x == bond).unwrap();
            self.orders[i]
        }
    }

    fn found(pattern: &Mol, target: &Mol) -> Vec<Mapping> {
        matches(
            pattern,
            target,
            |p, t| pattern.element(p) == target.element(t),
            |pe, te| pattern.bond_order(pe) == target.bond_order(te),
        )
        .collect()
    }

    fn empty() -> Mol {
        Mol {
            sites: vec![],
            elements: vec![],
            bonds: vec![],
            endpoints: vec![],
            orders: vec![],
        }
    }

    fn carbon() -> Mol {
        Mol {
            sites: vec![s(1)],
            elements: vec![elem("C")],
            bonds: vec![],
            endpoints: vec![],
            orders: vec![],
        }
    }

    fn carbon_pair() -> Mol {
        Mol {
            sites: vec![s(1), s(2)],
            elements: vec![elem("C"), elem("C")],
            bonds: vec![],
            endpoints: vec![],
            orders: vec![],
        }
    }

    fn carbon_carbon() -> Mol {
        Mol {
            sites: vec![s(1), s(2)],
            elements: vec![elem("C"), elem("C")],
            bonds: vec![b(1)],
            endpoints: vec![(s(1), s(2))],
            orders: vec![BondOrder::Single],
        }
    }

    fn carbon_oxygen() -> Mol {
        Mol {
            sites: vec![s(1), s(2)],
            elements: vec![elem("C"), elem("O")],
            bonds: vec![b(1)],
            endpoints: vec![(s(1), s(2))],
            orders: vec![BondOrder::Single],
        }
    }

    fn ethene() -> Mol {
        Mol {
            sites: vec![s(1), s(2)],
            elements: vec![elem("C"), elem("C")],
            bonds: vec![b(1)],
            endpoints: vec![(s(1), s(2))],
            orders: vec![BondOrder::Double],
        }
    }

    fn propane() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            elements: vec![elem("C"), elem("C"), elem("C")],
            bonds: vec![b(1), b(2)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3))],
            orders: vec![BondOrder::Single, BondOrder::Single],
        }
    }

    fn ethanol() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            elements: vec![elem("C"), elem("C"), elem("O")],
            bonds: vec![b(1), b(2)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3))],
            orders: vec![BondOrder::Single, BondOrder::Single],
        }
    }

    fn cyclohexane() -> Mol {
        Mol {
            sites: (1..=6).map(s).collect(),
            elements: vec![elem("C"); 6],
            bonds: (1..=6).map(b).collect(),
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(6)),
                (s(6), s(1)),
            ],
            orders: vec![BondOrder::Single; 6],
        }
    }

    #[test]
    fn empty_pattern_matches_once() {
        let found = found(&empty(), &ethanol());
        assert_eq!(found.len(), 1);
        assert!(found[0].is_empty());
    }

    #[test]
    fn a_single_atom_matches_each_like_atom() {
        assert_eq!(found(&carbon(), &ethanol()).len(), 2);
    }

    #[test]
    fn a_bond_matches_in_both_directions() {
        assert_eq!(found(&carbon_carbon(), &ethanol()).len(), 2);
    }

    #[test]
    fn a_path_matches_every_walk_of_a_ring() {
        assert_eq!(found(&propane(), &cyclohexane()).len(), 12);
    }

    #[test]
    fn a_pattern_larger_than_the_target_does_not_match() {
        assert!(found(&cyclohexane(), &ethanol()).is_empty());
    }

    #[test]
    fn an_element_constrains_the_match() {
        assert!(found(&carbon_oxygen(), &cyclohexane()).is_empty());
    }

    #[test]
    fn a_bond_order_constrains_the_match() {
        assert!(found(&ethene(), &cyclohexane()).is_empty());
    }

    #[test]
    fn an_absent_substructure_yields_nothing() {
        assert!(found(&carbon_oxygen(), &propane()).is_empty());
    }

    #[test]
    fn a_disjoint_pattern_matches_distinct_atoms() {
        let found = found(&carbon_pair(), &ethanol());
        assert_eq!(found.len(), 2);
        assert!(found.iter().all(|m| m.get(s(1)) != m.get(s(2))));
    }

    #[test]
    fn a_mapping_resolves_each_pattern_site() {
        let found = found(&carbon_oxygen(), &ethanol());
        assert_eq!(found.len(), 1);
        let mapping = &found[0];
        assert_eq!(mapping.len(), 2);
        assert_eq!(mapping.get(s(1)), Some(s(2)));
        assert_eq!(mapping.get(s(2)), Some(s(3)));
        assert_eq!(
            mapping.iter().collect::<Vec<_>>(),
            vec![(s(1), s(2)), (s(2), s(3))]
        );
    }

    #[test]
    fn a_mapping_has_no_target_for_an_absent_site() {
        let found = found(&carbon_oxygen(), &ethanol());
        assert_eq!(found[0].get(s(99)), None);
    }

    #[test]
    fn the_first_match_is_found_without_enumerating_all() {
        let carbon = carbon();
        let ring = cyclohexane();
        let mut search = matches(
            &carbon,
            &ring,
            |p, t| carbon.element(p) == ring.element(t),
            |pe, te| carbon.bond_order(pe) == ring.bond_order(te),
        );
        assert!(search.next().is_some());
    }
}
