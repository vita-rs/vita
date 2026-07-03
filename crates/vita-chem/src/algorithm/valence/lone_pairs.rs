use vita_core::{HasElements, SiteId};

use super::valence;
use crate::algorithm::utils::valence_electrons;
use crate::{HasBondOrders, HasFormalCharges, HasRadicalElectrons};

/// Number of lone (non-bonding) electron pairs on `site`.
///
/// The valence electrons left once the bonding and unpaired ones are taken away,
/// paired up: `(valence electrons − formal charge − bond-order sum − radicals)
/// / 2`. Water's oxygen has two; the ammonium nitrogen has none.
///
/// Returns `None` when no exact count exists:
/// - `site` holds a d- or f-block element, whose valence-electron count is not
///   fixed;
/// - an incident bond is aromatic, so the bonding electrons are not localised
///   (see [`valence`], which the count builds on — kekulise the ring first);
/// - the arithmetic goes negative, describing an impossible structure.
///
/// # Complexity
///
/// O(d) time and O(1) space, where `d` is the degree of `site`, assuming
/// [`bonds_of`](crate::HasBonds::bonds_of) runs in O(degree).
pub fn lone_pairs<M: HasBondOrders + HasElements + HasFormalCharges + HasRadicalElectrons>(
    mol: &M,
    site: SiteId,
) -> Option<u32> {
    let electrons = valence_electrons(mol.element(site))? as i32;
    let bonding = valence(mol, site)? as i32;
    let charge = mol.formal_charge(site) as i32;
    let radicals = mol.radical_electron(site) as i32;
    let free = electrons - charge - bonding - radicals;
    if free < 0 {
        return None;
    }
    Some(free as u32 / 2)
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

    fn atom(symbol: &str, charge: i8, radicals: u8, orders: &[BondOrder]) -> Mol {
        let n = orders.len() as u32;
        let mut sites = vec![s(1)];
        let mut elements = vec![elem(symbol)];
        let mut formal_charges = vec![charge];
        let mut radical_counts = vec![radicals];
        for i in 2..=n + 1 {
            sites.push(s(i));
            elements.push(elem("H"));
            formal_charges.push(0);
            radical_counts.push(0);
        }
        Mol {
            sites,
            elements,
            bonds: (1..=n).map(b).collect(),
            endpoints: (2..=n + 1).map(|i| (s(1), s(i))).collect(),
            orders: orders.to_vec(),
            formal_charges,
            radicals: radical_counts,
        }
    }

    #[test]
    fn a_bare_atom_pairs_all_its_valence_electrons() {
        assert_eq!(lone_pairs(&atom("O", 0, 0, &[]), s(1)), Some(3));
    }

    #[test]
    fn counts_lone_pairs_of_neutral_atoms() {
        let water = atom("O", 0, 0, &[BondOrder::Single, BondOrder::Single]);
        assert_eq!(lone_pairs(&water, s(1)), Some(2));

        let ammonia = atom("N", 0, 0, &[BondOrder::Single; 3]);
        assert_eq!(lone_pairs(&ammonia, s(1)), Some(1));

        let methane = atom("C", 0, 0, &[BondOrder::Single; 4]);
        assert_eq!(lone_pairs(&methane, s(1)), Some(0));
    }

    #[test]
    fn positive_formal_charge_removes_electrons() {
        let ammonium = atom("N", 1, 0, &[BondOrder::Single; 4]);
        assert_eq!(lone_pairs(&ammonium, s(1)), Some(0));
    }

    #[test]
    fn negative_formal_charge_adds_electrons() {
        let hydroxide = atom("O", -1, 0, &[BondOrder::Single]);
        assert_eq!(lone_pairs(&hydroxide, s(1)), Some(3));
    }

    #[test]
    fn radical_electrons_stay_unpaired() {
        let methyl = atom("C", 0, 1, &[BondOrder::Single; 3]);
        assert_eq!(lone_pairs(&methyl, s(1)), Some(0));
    }

    #[test]
    fn undefined_for_a_d_block_element() {
        assert_eq!(lone_pairs(&atom("Fe", 0, 0, &[]), s(1)), None);
    }

    #[test]
    fn undefined_beside_an_aromatic_bond() {
        assert_eq!(
            lone_pairs(&atom("C", 0, 0, &[BondOrder::Aromatic]), s(1)),
            None
        );
    }

    #[test]
    fn undefined_when_the_electron_count_goes_negative() {
        assert_eq!(
            lone_pairs(&atom("N", 0, 0, &[BondOrder::Hextuple]), s(1)),
            None
        );
    }

    #[test]
    fn an_odd_free_electron_count_rounds_down() {
        let aminyl = atom("N", 0, 0, &[BondOrder::Single, BondOrder::Single]);
        assert_eq!(lone_pairs(&aminyl, s(1)), Some(1));
    }
}
