use vita_core::HasElements;

use super::{Composition, Constituent};
use crate::HasFormalCharges;

/// The elemental composition of a molecule: every site counted under its
/// element at natural isotopic precision, with the net formal charge.
///
/// The coarser of the two folds — isotopic declarations, if any, are not
/// consulted. For a composition that keeps them, use
/// [`isotopic`](super::isotopic).
///
/// # Complexity
///
/// O(V · log V) time and O(V) space, over the molecule's `V` sites; the log
/// factor orders the counts canonically.
pub fn elemental<M: HasElements + HasFormalCharges>(mol: &M) -> Composition {
    Composition::from_counts(
        mol.elements()
            .map(|element| (Constituent::Element(element), 1)),
        mol.formal_charges().map(i32::from).sum(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::{Element, HasSites, SiteId};

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn elem(symbol: &str) -> Element {
        Element::from_symbol(symbol).unwrap()
    }

    fn natural(symbol: &str) -> Constituent {
        Constituent::Element(elem(symbol))
    }

    struct Mol {
        sites: Vec<SiteId>,
        elements: Vec<Element>,
        formal_charges: Vec<i8>,
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

    impl HasFormalCharges for Mol {
        fn formal_charge(&self, site: SiteId) -> i8 {
            let i = self.sites.iter().position(|&x| x == site).unwrap();
            self.formal_charges[i]
        }
    }

    fn no_sites() -> Mol {
        Mol {
            sites: vec![],
            elements: vec![],
            formal_charges: vec![],
        }
    }

    fn lone_hydrogen() -> Mol {
        Mol {
            sites: vec![s(1)],
            elements: vec![elem("H")],
            formal_charges: vec![0],
        }
    }

    fn water_mol() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            elements: vec![elem("O"), elem("H"), elem("H")],
            formal_charges: vec![0, 0, 0],
        }
    }

    fn hydroxide() -> Mol {
        Mol {
            sites: vec![s(1), s(2)],
            elements: vec![elem("O"), elem("H")],
            formal_charges: vec![-1, 0],
        }
    }

    #[test]
    fn a_molecule_without_sites_has_an_empty_elemental_composition() {
        assert_eq!(elemental(&no_sites()), Composition::from_counts([], 0));
    }

    #[test]
    fn a_single_site_molecule_counts_one_atom() {
        assert_eq!(elemental(&lone_hydrogen()).atom_count(), 1);
    }

    #[test]
    fn elemental_counts_every_site_under_its_element() {
        let expected = Composition::from_counts([(natural("H"), 2), (natural("O"), 1)], 0);
        assert_eq!(elemental(&water_mol()), expected);
    }

    #[test]
    fn elemental_sums_the_formal_charges() {
        assert_eq!(elemental(&hydroxide()).charge(), -1);
    }

    #[test]
    fn the_elemental_composition_is_independent_of_input_order() {
        let reordered = Mol {
            sites: vec![s(2), s(3), s(1)],
            elements: vec![elem("H"), elem("H"), elem("O")],
            formal_charges: vec![0, 0, 0],
        };
        assert_eq!(elemental(&reordered), elemental(&water_mol()));
    }
}
