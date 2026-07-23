use std::cmp::Ordering;

use vita_core::{HasElements, SiteId};

use crate::algorithm::utils::electronegativity_rank;
use crate::{BondOrder, HasBondOrders, HasFormalCharges};

/// Oxidation state of `site`: the charge left once every bond's electrons are
/// awarded wholly to the bond's more electronegative end.
///
/// The ionic approximation of the IUPAC 2016 recommendation: relative to the
/// formal charge, each bond adds its integer order when the partner outranks
/// `site` on the Allen electronegativity order (Pauling across the f-block),
/// subtracts it when `site` outranks the partner, and moves nothing between
/// equal ranks — homonuclear bonds above all. Methane's carbon is −4, the
/// carbon of CO₂ is +4, and the d-block answers where
/// [`lone_pairs`](super::lone_pairs) cannot. The recommendation's caveat — a
/// reversibly bonded Lewis-acid ligand keeps its donated pair — rests on a
/// reversibility the bond graph does not record, and is not applied.
///
/// Returns `None` when no exact state exists:
/// - an incident bond is [`Aromatic`](BondOrder::Aromatic), its electrons not
///   localized into an integer order (see [`valence`](super::valence) —
///   kekulize the ring first);
/// - an incident bond reaches an element beyond lawrencium, past every
///   measured electronegativity;
/// - the total leaves `i8`, describing an impossible structure.
///
/// # Complexity
///
/// O(d) time and O(1) space, where `d` is the degree of `site`, assuming
/// [`bonds_of`](crate::HasBonds::bonds_of) runs in O(degree).
pub fn oxidation_state<M: HasBondOrders + HasElements + HasFormalCharges>(
    mol: &M,
    site: SiteId,
) -> Option<i8> {
    let own = electronegativity_rank(mol.element(site));
    let mut state = i32::from(mol.formal_charge(site));
    for (bond, partner) in mol.bonds_of(site) {
        let order: i32 = match mol.bond_order(bond) {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 4,
            BondOrder::Quintuple => 5,
            BondOrder::Hextuple => 6,
            BondOrder::Aromatic => return None,
        };
        state += order
            * match own?.cmp(&electronegativity_rank(mol.element(partner))?) {
                Ordering::Less => 1,
                Ordering::Equal => 0,
                Ordering::Greater => -1,
            };
    }
    i8::try_from(state).ok()
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::{Element, HasSites};

    use crate::{BondId, HasBonds};

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

    fn molecule(atoms: &[(&str, i8)], bonds: &[(u32, u32, BondOrder)]) -> Mol {
        Mol {
            sites: (1..=atoms.len() as u32).map(s).collect(),
            elements: atoms.iter().map(|&(symbol, _)| elem(symbol)).collect(),
            bonds: (1..=bonds.len() as u32).map(b).collect(),
            endpoints: bonds.iter().map(|&(u, v, _)| (s(u), s(v))).collect(),
            orders: bonds.iter().map(|&(_, _, order)| order).collect(),
            formal_charges: atoms.iter().map(|&(_, charge)| charge).collect(),
        }
    }

    #[test]
    fn a_bond_free_site_keeps_its_formal_charge() {
        let sodium = molecule(&[("Na", 1)], &[]);
        assert_eq!(oxidation_state(&sodium, s(1)), Some(1));
    }

    #[test]
    fn a_bond_free_site_needs_no_ranking() {
        let rutherfordium = molecule(&[("Rf", 0)], &[]);
        assert_eq!(oxidation_state(&rutherfordium, s(1)), Some(0));
    }

    #[test]
    fn a_homonuclear_bond_moves_nothing() {
        let dihydrogen = molecule(&[("H", 0), ("H", 0)], &[(1, 2, BondOrder::Single)]);
        assert_eq!(oxidation_state(&dihydrogen, s(1)), Some(0));
        assert_eq!(oxidation_state(&dihydrogen, s(2)), Some(0));
    }

    #[test]
    fn methane_carbon_takes_all_four_bonds() {
        let methane = molecule(
            &[("C", 0), ("H", 0), ("H", 0), ("H", 0), ("H", 0)],
            &[
                (1, 2, BondOrder::Single),
                (1, 3, BondOrder::Single),
                (1, 4, BondOrder::Single),
                (1, 5, BondOrder::Single),
            ],
        );
        assert_eq!(oxidation_state(&methane, s(1)), Some(-4));
        assert_eq!(oxidation_state(&methane, s(2)), Some(1));
    }

    #[test]
    fn a_double_bond_moves_twice_the_electrons() {
        let carbon_dioxide = molecule(
            &[("O", 0), ("C", 0), ("O", 0)],
            &[(1, 2, BondOrder::Double), (2, 3, BondOrder::Double)],
        );
        assert_eq!(oxidation_state(&carbon_dioxide, s(2)), Some(4));
        assert_eq!(oxidation_state(&carbon_dioxide, s(1)), Some(-2));
    }

    #[test]
    fn a_triple_bond_moves_three_times_the_electrons() {
        let hydrogen_cyanide = molecule(
            &[("H", 0), ("C", 0), ("N", 0)],
            &[(1, 2, BondOrder::Single), (2, 3, BondOrder::Triple)],
        );
        assert_eq!(oxidation_state(&hydrogen_cyanide, s(3)), Some(-3));
        assert_eq!(oxidation_state(&hydrogen_cyanide, s(2)), Some(2));
    }

    #[test]
    fn higher_bond_orders_move_their_full_multiple() {
        for (order, moved) in [
            (BondOrder::Quadruple, 4),
            (BondOrder::Quintuple, 5),
            (BondOrder::Hextuple, 6),
        ] {
            let carbide = molecule(&[("W", 0), ("C", 0)], &[(1, 2, order)]);
            assert_eq!(oxidation_state(&carbide, s(1)), Some(moved));
            assert_eq!(oxidation_state(&carbide, s(2)), Some(-moved));
        }
    }

    #[test]
    fn the_formal_charge_seeds_the_state() {
        let ammonium = molecule(
            &[("N", 1), ("H", 0), ("H", 0), ("H", 0), ("H", 0)],
            &[
                (1, 2, BondOrder::Single),
                (1, 3, BondOrder::Single),
                (1, 4, BondOrder::Single),
                (1, 5, BondOrder::Single),
            ],
        );
        assert_eq!(oxidation_state(&ammonium, s(1)), Some(-3));
    }

    #[test]
    fn a_rank_tie_moves_nothing() {
        let chromium_osmium = molecule(&[("Cr", 0), ("Os", 0)], &[(1, 2, BondOrder::Single)]);
        assert_eq!(oxidation_state(&chromium_osmium, s(1)), Some(0));
        assert_eq!(oxidation_state(&chromium_osmium, s(2)), Some(0));
    }

    #[test]
    fn iodomethane_iodine_is_plus_one() {
        let iodomethane = molecule(
            &[("C", 0), ("I", 0), ("H", 0), ("H", 0), ("H", 0)],
            &[
                (1, 2, BondOrder::Single),
                (1, 3, BondOrder::Single),
                (1, 4, BondOrder::Single),
                (1, 5, BondOrder::Single),
            ],
        );
        assert_eq!(oxidation_state(&iodomethane, s(2)), Some(1));
        assert_eq!(oxidation_state(&iodomethane, s(1)), Some(-4));
    }

    #[test]
    fn defined_for_a_d_block_element() {
        let titanium_tetrachloride = molecule(
            &[("Ti", 0), ("Cl", 0), ("Cl", 0), ("Cl", 0), ("Cl", 0)],
            &[
                (1, 2, BondOrder::Single),
                (1, 3, BondOrder::Single),
                (1, 4, BondOrder::Single),
                (1, 5, BondOrder::Single),
            ],
        );
        assert_eq!(oxidation_state(&titanium_tetrachloride, s(1)), Some(4));
        assert_eq!(oxidation_state(&titanium_tetrachloride, s(2)), Some(-1));
    }

    #[test]
    fn a_hypervalent_site_still_answers() {
        let sulfur_hexafluoride = molecule(
            &[
                ("S", 0),
                ("F", 0),
                ("F", 0),
                ("F", 0),
                ("F", 0),
                ("F", 0),
                ("F", 0),
            ],
            &[
                (1, 2, BondOrder::Single),
                (1, 3, BondOrder::Single),
                (1, 4, BondOrder::Single),
                (1, 5, BondOrder::Single),
                (1, 6, BondOrder::Single),
                (1, 7, BondOrder::Single),
            ],
        );
        assert_eq!(oxidation_state(&sulfur_hexafluoride, s(1)), Some(6));
    }

    #[test]
    fn the_states_sum_to_the_net_charge() {
        let ammonium = molecule(
            &[("N", 1), ("H", 0), ("H", 0), ("H", 0), ("H", 0)],
            &[
                (1, 2, BondOrder::Single),
                (1, 3, BondOrder::Single),
                (1, 4, BondOrder::Single),
                (1, 5, BondOrder::Single),
            ],
        );
        let total: i32 = ammonium
            .sites()
            .map(|site| i32::from(oxidation_state(&ammonium, site).unwrap()))
            .sum();
        assert_eq!(total, 1);
    }

    #[test]
    fn undefined_beside_an_aromatic_bond() {
        let aromatic_pair = molecule(&[("C", 0), ("C", 0)], &[(1, 2, BondOrder::Aromatic)]);
        assert_eq!(oxidation_state(&aromatic_pair, s(1)), None);
    }

    #[test]
    fn undefined_beside_an_element_beyond_lawrencium() {
        let oxide = molecule(&[("Rf", 0), ("O", 0)], &[(1, 2, BondOrder::Double)]);
        assert_eq!(oxidation_state(&oxide, s(1)), None);
        assert_eq!(oxidation_state(&oxide, s(2)), None);
    }

    #[test]
    fn undefined_when_the_total_leaves_the_charge_range() {
        let mut atoms = vec![("C", 0)];
        let mut bonds = Vec::new();
        for i in 2..=129 {
            atoms.push(("F", 0));
            bonds.push((1, i, BondOrder::Single));
        }
        let fluoride_swarm = molecule(&atoms, &bonds);
        assert_eq!(oxidation_state(&fluoride_swarm, s(1)), None);
    }
}
