use vita_core::HasIsotopes;

use super::{Composition, Constituent};
use crate::HasFormalCharges;

/// The isotopic composition of a molecule: every site counted under its
/// declared nuclide, with the net formal charge.
///
/// The finer of the two folds — [`HasIsotopes`] declares one nuclide per
/// site, and every count keeps it. For natural-mixture counting, use
/// [`elemental`](super::elemental).
///
/// # Complexity
///
/// O(V · log V) time and O(V) space, over the molecule's `V` sites; the log
/// factor orders the counts canonically.
pub fn isotopic<M: HasIsotopes + HasFormalCharges>(mol: &M) -> Composition {
    Composition::from_counts(
        mol.isotopes()
            .map(|isotope| (Constituent::Nuclide(isotope), 1)),
        mol.formal_charges().map(i32::from).sum(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::{Element, HasElements, HasSites, Isotope, SiteId};

    use crate::algorithm::composition::elemental;

    fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn elem(symbol: &str) -> Element {
        Element::from_symbol(symbol).unwrap()
    }

    fn nuclide(symbol: &str, mass_number: u16) -> Constituent {
        Constituent::Nuclide(Isotope::new(elem(symbol), mass_number).unwrap())
    }

    struct Mol {
        sites: Vec<SiteId>,
        elements: Vec<Element>,
        mass_numbers: Vec<u16>,
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

    impl HasIsotopes for Mol {
        fn isotope(&self, site: SiteId) -> Isotope {
            let i = self.sites.iter().position(|&x| x == site).unwrap();
            Isotope::new(self.elements[i], self.mass_numbers[i]).unwrap()
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
            mass_numbers: vec![],
            formal_charges: vec![],
        }
    }

    fn hydroxide() -> Mol {
        Mol {
            sites: vec![s(1), s(2)],
            elements: vec![elem("O"), elem("H")],
            mass_numbers: vec![16, 1],
            formal_charges: vec![-1, 0],
        }
    }

    fn semiheavy_water() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            elements: vec![elem("O"), elem("H"), elem("H")],
            mass_numbers: vec![16, 1, 2],
            formal_charges: vec![0, 0, 0],
        }
    }

    #[test]
    fn a_molecule_without_sites_has_an_empty_isotopic_composition() {
        assert_eq!(isotopic(&no_sites()), Composition::from_counts([], 0));
    }

    #[test]
    fn isotopic_counts_every_site_under_its_nuclide() {
        let expected = Composition::from_counts(
            [
                (nuclide("O", 16), 1),
                (nuclide("H", 1), 1),
                (nuclide("H", 2), 1),
            ],
            0,
        );
        assert_eq!(isotopic(&semiheavy_water()), expected);
    }

    #[test]
    fn isotopic_sums_the_formal_charges() {
        assert_eq!(isotopic(&hydroxide()).charge(), -1);
    }

    #[test]
    fn elemental_and_isotopic_agree_on_totals_and_charge() {
        let coarse = elemental(&semiheavy_water());
        let fine = isotopic(&semiheavy_water());
        assert_eq!(coarse.atom_count(), fine.atom_count());
        assert_eq!(coarse.charge(), fine.charge());
    }

    #[test]
    fn the_isotopic_composition_is_independent_of_input_order() {
        let reordered = Mol {
            sites: vec![s(3), s(1), s(2)],
            elements: vec![elem("H"), elem("O"), elem("H")],
            mass_numbers: vec![2, 16, 1],
            formal_charges: vec![0, 0, 0],
        };
        assert_eq!(isotopic(&reordered), isotopic(&semiheavy_water()));
    }
}
