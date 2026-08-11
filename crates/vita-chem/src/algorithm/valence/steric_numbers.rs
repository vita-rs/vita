use vita_core::{HasElements, SiteId};

use super::lone_pairs;
use crate::algorithm::conjugation::systems;
use crate::algorithm::utils::SortedMap;
use crate::{HasBondOrders, HasFormalCharges, HasRadicalElectrons};

/// How many electron domains lie about each of a molecule's sites.
///
/// The number counts a site's bonds, one apiece whatever their order, plus the
/// lone pairs it keeps localized. A site whose arithmetic settles no exact
/// count goes uncounted.
///
/// Obtain via [`steric_numbers`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StericNumbers {
    counts: SortedMap<SiteId, u32>,
}

impl StericNumbers {
    /// Number of counted sites.
    pub fn len(&self) -> usize {
        self.counts.len()
    }

    /// Returns `true` if no site is counted.
    pub fn is_empty(&self) -> bool {
        self.counts.is_empty()
    }

    /// Returns the number of electron domains about `site`.
    ///
    /// Returns `None` if `site` is absent from the molecule or goes uncounted.
    pub fn get(&self, site: SiteId) -> Option<u32> {
        self.counts.get(&site).copied()
    }

    /// Iterates the `(site, steric number)` pairs, ordered by site.
    pub fn iter(&self) -> impl Iterator<Item = (SiteId, u32)> + '_ {
        self.counts.iter().map(|(&site, &count)| (site, count))
    }
}

/// The steric number of each of a molecule's sites.
///
/// A site's domains are the regions its valence electrons occupy: one per
/// bond, whatever its order — a multiple bond's π component lies along the
/// σ-bond already counted — and one per lone pair the site keeps localized. A
/// pair donated into a conjugated π network is not localized; it spreads along
/// the network and adds no domain. Discounting it leaves the total independent
/// of the Lewis form drawn: an amide nitrogen has three domains whether the
/// neutral form donates from it or the zwitterionic form donates from the
/// oxygen.
///
/// A site goes uncounted where [`lone_pairs`] settles no exact count — a d- or
/// f-block element, an incident aromatic bond, or arithmetic describing an
/// impossible structure.
///
/// # Complexity
///
/// O((V + E) · log (V + E)) time and O(V + E) space, over the molecule's `V`
/// sites and `E` bonds, assuming [`bonds_of`](crate::HasBonds::bonds_of) and
/// [`degree`](crate::HasBonds::degree) run in O(degree); perceiving the
/// conjugated systems dominates.
pub fn steric_numbers<M: HasBondOrders + HasElements + HasFormalCharges + HasRadicalElectrons>(
    mol: &M,
) -> StericNumbers {
    let conjugation = systems(mol);
    StericNumbers {
        counts: SortedMap::from_pairs(mol.sites().filter_map(|site| {
            let pairs = lone_pairs(mol, site)?;
            let donated: u32 = conjugation
                .of_site(site)
                .map(|system| system.donated_pairs(site))
                .sum();
            let kept = pairs
                .checked_sub(donated)
                .expect("a system draws its donations from the donor's own lone pairs");
            Some((site, mol.degree(site) as u32 + kept))
        })),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::{Element, HasSites};

    use crate::{BondId, BondOrder, HasBonds};

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
        formal_charges: Vec<i8>,
        radicals: Vec<u8>,
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

    impl HasFormalCharges for Mol {
        fn formal_charge(&self, site: SiteId) -> i8 {
            let i = self.sites.iter().position(|&x| x == site).unwrap();
            self.formal_charges[i]
        }
    }

    impl HasRadicalElectrons for Mol {
        fn radical_electron(&self, site: SiteId) -> u8 {
            let i = self.sites.iter().position(|&x| x == site).unwrap();
            self.radicals[i]
        }
    }

    fn molecule(atoms: &[(u32, &str, i8, u8)], bonds: &[(u32, u32, u32, BondOrder)]) -> Mol {
        Mol {
            sites: atoms.iter().map(|&(id, ..)| s(id)).collect(),
            elements: atoms.iter().map(|&(_, symbol, ..)| elem(symbol)).collect(),
            formal_charges: atoms.iter().map(|&(_, _, charge, _)| charge).collect(),
            radicals: atoms.iter().map(|&(.., radicals)| radicals).collect(),
            bonds: bonds.iter().map(|&(id, ..)| b(id)).collect(),
            endpoints: bonds.iter().map(|&(_, u, v, _)| (s(u), s(v))).collect(),
            orders: bonds.iter().map(|&(.., order)| order).collect(),
        }
    }

    fn reversed(m: &Mol) -> Mol {
        Mol {
            sites: m.sites.iter().rev().copied().collect(),
            elements: m.elements.iter().rev().copied().collect(),
            bonds: m.bonds.iter().rev().copied().collect(),
            endpoints: m.endpoints.iter().rev().copied().collect(),
            orders: m.orders.iter().rev().copied().collect(),
            formal_charges: m.formal_charges.iter().rev().copied().collect(),
            radicals: m.radicals.iter().rev().copied().collect(),
        }
    }

    fn empty() -> Mol {
        molecule(&[], &[])
    }

    fn methane() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "H", 0, 0),
                (3, "H", 0, 0),
                (4, "H", 0, 0),
                (5, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Single),
                (2, 1, 3, BondOrder::Single),
                (3, 1, 4, BondOrder::Single),
                (4, 1, 5, BondOrder::Single),
            ],
        )
    }

    fn water() -> Mol {
        molecule(
            &[(1, "O", 0, 0), (2, "H", 0, 0), (3, "H", 0, 0)],
            &[(1, 1, 2, BondOrder::Single), (2, 1, 3, BondOrder::Single)],
        )
    }

    fn ammonia() -> Mol {
        molecule(
            &[
                (1, "N", 0, 0),
                (2, "H", 0, 0),
                (3, "H", 0, 0),
                (4, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Single),
                (2, 1, 3, BondOrder::Single),
                (3, 1, 4, BondOrder::Single),
            ],
        )
    }

    fn formaldehyde() -> Mol {
        molecule(
            &[
                (1, "O", 0, 0),
                (2, "C", 0, 0),
                (3, "H", 0, 0),
                (4, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 2, 3, BondOrder::Single),
                (3, 2, 4, BondOrder::Single),
            ],
        )
    }

    fn hydrogen_cyanide() -> Mol {
        molecule(
            &[(1, "H", 0, 0), (2, "C", 0, 0), (3, "N", 0, 0)],
            &[(1, 1, 2, BondOrder::Single), (2, 2, 3, BondOrder::Triple)],
        )
    }

    fn formamide() -> Mol {
        molecule(
            &[
                (1, "O", 0, 0),
                (2, "C", 0, 0),
                (3, "N", 0, 0),
                (4, "H", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Double),
                (2, 2, 3, BondOrder::Single),
                (3, 2, 4, BondOrder::Single),
                (4, 3, 5, BondOrder::Single),
                (5, 3, 6, BondOrder::Single),
            ],
        )
    }

    fn formamide_polar() -> Mol {
        molecule(
            &[
                (1, "O", -1, 0),
                (2, "C", 0, 0),
                (3, "N", 1, 0),
                (4, "H", 0, 0),
                (5, "H", 0, 0),
                (6, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Single),
                (2, 2, 3, BondOrder::Double),
                (3, 2, 4, BondOrder::Single),
                (4, 3, 5, BondOrder::Single),
                (5, 3, 6, BondOrder::Single),
            ],
        )
    }

    fn azide() -> Mol {
        molecule(
            &[(1, "N", -1, 0), (2, "N", 1, 0), (3, "N", -1, 0)],
            &[(1, 1, 2, BondOrder::Double), (2, 2, 3, BondOrder::Double)],
        )
    }

    fn azide_unsymmetric() -> Mol {
        molecule(
            &[(1, "N", 0, 0), (2, "N", 1, 0), (3, "N", -2, 0)],
            &[(1, 1, 2, BondOrder::Triple), (2, 2, 3, BondOrder::Single)],
        )
    }

    fn sulfur_hexafluoride() -> Mol {
        molecule(
            &[
                (1, "S", 0, 0),
                (2, "F", 0, 0),
                (3, "F", 0, 0),
                (4, "F", 0, 0),
                (5, "F", 0, 0),
                (6, "F", 0, 0),
                (7, "F", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Single),
                (2, 1, 3, BondOrder::Single),
                (3, 1, 4, BondOrder::Single),
                (4, 1, 5, BondOrder::Single),
                (5, 1, 6, BondOrder::Single),
                (6, 1, 7, BondOrder::Single),
            ],
        )
    }

    fn aromatic_benzene() -> Mol {
        molecule(
            &[
                (1, "C", 0, 0),
                (2, "C", 0, 0),
                (3, "C", 0, 0),
                (4, "C", 0, 0),
                (5, "C", 0, 0),
                (6, "C", 0, 0),
                (7, "H", 0, 0),
                (8, "H", 0, 0),
                (9, "H", 0, 0),
                (10, "H", 0, 0),
                (11, "H", 0, 0),
                (12, "H", 0, 0),
            ],
            &[
                (1, 1, 2, BondOrder::Aromatic),
                (2, 2, 3, BondOrder::Aromatic),
                (3, 3, 4, BondOrder::Aromatic),
                (4, 4, 5, BondOrder::Aromatic),
                (5, 5, 6, BondOrder::Aromatic),
                (6, 6, 1, BondOrder::Aromatic),
                (7, 1, 7, BondOrder::Single),
                (8, 2, 8, BondOrder::Single),
                (9, 3, 9, BondOrder::Single),
                (10, 4, 10, BondOrder::Single),
                (11, 5, 11, BondOrder::Single),
                (12, 6, 12, BondOrder::Single),
            ],
        )
    }

    #[test]
    fn an_empty_molecule_counts_nothing() {
        let counted = steric_numbers(&empty());
        assert_eq!(counted.len(), 0);
        assert!(counted.is_empty());
    }

    #[test]
    fn counts_bonds_and_localized_lone_pairs() {
        assert_eq!(steric_numbers(&methane()).get(s(1)), Some(4));
        assert_eq!(steric_numbers(&water()).get(s(1)), Some(4));
        assert_eq!(steric_numbers(&ammonia()).get(s(1)), Some(4));
    }

    #[test]
    fn a_bondless_site_counts_only_its_lone_pairs() {
        let atom = molecule(&[(1, "O", 0, 0)], &[]);
        assert_eq!(steric_numbers(&atom).get(s(1)), Some(3));
    }

    #[test]
    fn a_multiple_bond_contributes_one_domain() {
        assert_eq!(steric_numbers(&formaldehyde()).get(s(2)), Some(3));
        assert_eq!(steric_numbers(&hydrogen_cyanide()).get(s(2)), Some(2));
    }

    #[test]
    fn a_donated_lone_pair_is_no_domain_of_its_own() {
        assert_eq!(steric_numbers(&formamide()).get(s(3)), Some(3));
    }

    #[test]
    fn orthogonal_pi_planes_each_take_their_donation() {
        let counted = steric_numbers(&azide_unsymmetric());
        assert_eq!(counted.get(s(3)), Some(2));
    }

    #[test]
    fn formamide_counts_the_same_in_either_lewis_form() {
        assert_eq!(
            steric_numbers(&formamide()),
            steric_numbers(&formamide_polar())
        );
    }

    #[test]
    fn azide_counts_the_same_in_either_lewis_form() {
        assert_eq!(
            steric_numbers(&azide()),
            steric_numbers(&azide_unsymmetric())
        );
    }

    #[test]
    fn a_hypervalent_site_counts_every_bond() {
        assert_eq!(steric_numbers(&sulfur_hexafluoride()).get(s(1)), Some(6));
    }

    #[test]
    fn a_d_block_site_goes_uncounted() {
        let atom = molecule(&[(1, "Fe", 0, 0)], &[]);
        assert!(steric_numbers(&atom).is_empty());
    }

    #[test]
    fn an_aromatic_ring_leaves_its_carbons_uncounted() {
        let counted = steric_numbers(&aromatic_benzene());
        assert_eq!(counted.get(s(1)), None);
        assert_eq!(counted.get(s(7)), Some(1));
        assert_eq!(counted.len(), 6);
    }

    #[test]
    fn an_absent_site_has_no_count() {
        assert_eq!(steric_numbers(&methane()).get(s(99)), None);
    }

    #[test]
    fn iter_yields_the_counted_sites_in_ascending_order() {
        let counted: Vec<(SiteId, u32)> = steric_numbers(&water()).iter().collect();
        assert_eq!(counted, vec![(s(1), 4), (s(2), 1), (s(3), 1)]);
    }

    #[test]
    fn the_counts_are_independent_of_input_order() {
        assert_eq!(
            steric_numbers(&formamide()),
            steric_numbers(&reversed(&formamide()))
        );
    }
}
