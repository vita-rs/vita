use vita_core::{Element, HasElements, SiteId};

use super::explicit::valence;
use crate::{HasBondOrders, HasFormalCharges, HasRadicalElectrons};

/// Number of lone (non-bonding) electron pairs on `site`.
///
/// What is left of the element's valence electrons once the bonding and
/// unpaired ones are removed: `(valence electrons − formal charge − bond-order
/// sum − radicals) / 2`. Oxygen in water has two; the ammonium nitrogen has
/// none.
///
/// Returns `None` when no exact count exists:
/// - `site` holds a d- or f-block element, whose valence-electron count is
///   ambiguous;
/// - an incident bond is aromatic, so the bonding electrons are not localised
///   (see [`valence`] — kekulise the ring first);
/// - the bookkeeping is negative, describing an impossible structure.
///
/// # Complexity
///
/// O(degree) time.
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

/// Valence (outer-shell s and p) electron count of a main-group element.
///
/// Returns `None` for the d- and f-block, where the count is ambiguous.
fn valence_electrons(element: Element) -> Option<u8> {
    Some(match element.atomic_number() {
        1 | 3 | 11 | 19 | 37 | 55 | 87 => 1,
        2 | 4 | 12 | 20 | 38 | 56 | 88 => 2,
        5 | 13 | 31 | 49 | 81 | 113 => 3,
        6 | 14 | 32 | 50 | 82 | 114 => 4,
        7 | 15 | 33 | 51 | 83 | 115 => 5,
        8 | 16 | 34 | 52 | 84 | 116 => 6,
        9 | 17 | 35 | 53 | 85 | 117 => 7,
        10 | 18 | 36 | 54 | 86 | 118 => 8,
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BondId, BondOrder, HasBonds};
    use vita_core::HasSites;

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
        charges: Vec<i8>,
        radicals: Vec<u8>,
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
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.elements[i]
        }
    }

    impl HasFormalCharges for Mol {
        fn formal_charge(&self, site: SiteId) -> i8 {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.charges[i]
        }
    }

    impl HasRadicalElectrons for Mol {
        fn radical_electron(&self, site: SiteId) -> u8 {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.radicals[i]
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

    fn atom(symbol: &str, charge: i8, radicals: u8, bonds: &[BondOrder]) -> Mol {
        let mut mol = Mol {
            sites: vec![s(1)],
            elements: vec![elem(symbol)],
            charges: vec![charge],
            radicals: vec![radicals],
            bonds: Vec::new(),
            endpoints: Vec::new(),
            orders: Vec::new(),
        };
        for (i, &order) in bonds.iter().enumerate() {
            let neighbour = s(i as u32 + 2);
            mol.sites.push(neighbour);
            mol.elements.push(elem("H"));
            mol.charges.push(0);
            mol.radicals.push(0);
            mol.bonds.push(b(i as u32 + 1));
            mol.endpoints.push((s(1), neighbour));
            mol.orders.push(order);
        }
        mol
    }

    #[test]
    fn neutral_atoms() {
        assert_eq!(
            lone_pairs(&atom("O", 0, 0, &[BondOrder::Single; 2]), s(1)),
            Some(2)
        );
        assert_eq!(
            lone_pairs(&atom("N", 0, 0, &[BondOrder::Single; 3]), s(1)),
            Some(1)
        );
        assert_eq!(
            lone_pairs(&atom("C", 0, 0, &[BondOrder::Single; 4]), s(1)),
            Some(0)
        );
        assert_eq!(
            lone_pairs(&atom("F", 0, 0, &[BondOrder::Single]), s(1)),
            Some(3)
        );
    }

    #[test]
    fn helium_has_one_pair() {
        assert_eq!(lone_pairs(&atom("He", 0, 0, &[]), s(1)), Some(1));
    }

    #[test]
    fn bond_order_enters_through_valence() {
        assert_eq!(
            lone_pairs(&atom("O", 0, 0, &[BondOrder::Double]), s(1)),
            Some(2)
        );
        assert_eq!(
            lone_pairs(&atom("N", 0, 0, &[BondOrder::Triple]), s(1)),
            Some(1)
        );
    }

    #[test]
    fn formal_charge_shifts_the_count() {
        assert_eq!(
            lone_pairs(&atom("N", 1, 0, &[BondOrder::Single; 4]), s(1)),
            Some(0)
        );
        assert_eq!(
            lone_pairs(&atom("O", 1, 0, &[BondOrder::Single; 3]), s(1)),
            Some(1)
        );
        assert_eq!(
            lone_pairs(&atom("O", -1, 0, &[BondOrder::Single]), s(1)),
            Some(3)
        );
        assert_eq!(
            lone_pairs(&atom("N", -1, 0, &[BondOrder::Single; 2]), s(1)),
            Some(2)
        );
    }

    #[test]
    fn radicals_are_subtracted() {
        assert_eq!(
            lone_pairs(&atom("C", 0, 1, &[BondOrder::Single; 3]), s(1)),
            Some(0)
        );
        assert_eq!(
            lone_pairs(&atom("O", 0, 1, &[BondOrder::Single]), s(1)),
            Some(2)
        );
    }

    #[test]
    fn electron_deficient_atoms_have_no_pairs() {
        assert_eq!(
            lone_pairs(&atom("B", 0, 0, &[BondOrder::Single; 3]), s(1)),
            Some(0)
        );
        assert_eq!(
            lone_pairs(&atom("C", 1, 0, &[BondOrder::Single; 3]), s(1)),
            Some(0)
        );
    }

    #[test]
    fn hypervalent_atoms_have_no_pairs() {
        assert_eq!(
            lone_pairs(&atom("S", 0, 0, &[BondOrder::Single; 6]), s(1)),
            Some(0)
        );
        assert_eq!(
            lone_pairs(&atom("P", 0, 0, &[BondOrder::Single; 5]), s(1)),
            Some(0)
        );
    }

    #[test]
    fn odd_electron_count_rounds_down() {
        assert_eq!(
            lone_pairs(&atom("N", 0, 0, &[BondOrder::Double]), s(1)),
            Some(1)
        );
    }

    #[test]
    fn transition_and_inner_metals_are_undefined() {
        assert_eq!(
            lone_pairs(&atom("Fe", 0, 0, &[BondOrder::Single; 2]), s(1)),
            None
        );
        assert_eq!(
            lone_pairs(&atom("U", 0, 0, &[BondOrder::Single]), s(1)),
            None
        );
    }

    #[test]
    fn aromatic_bonds_are_undefined() {
        assert_eq!(
            lone_pairs(&atom("C", 0, 0, &[BondOrder::Aromatic; 2]), s(1)),
            None
        );
    }

    #[test]
    fn impossible_valence_is_undefined() {
        assert_eq!(
            lone_pairs(&atom("C", 0, 0, &[BondOrder::Single; 5]), s(1)),
            None
        );
    }

    #[test]
    fn valence_electron_counts() {
        for (symbol, expected) in [
            ("H", Some(1)),
            ("He", Some(2)),
            ("Li", Some(1)),
            ("B", Some(3)),
            ("C", Some(4)),
            ("N", Some(5)),
            ("O", Some(6)),
            ("F", Some(7)),
            ("Ne", Some(8)),
            ("Cl", Some(7)),
            ("Pb", Some(4)),
            ("Rn", Some(8)),
            ("Fe", None),
            ("U", None),
        ] {
            assert_eq!(valence_electrons(elem(symbol)), expected, "{symbol}");
        }
    }
}
