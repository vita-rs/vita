use vita_core::SiteId;

use crate::{BondOrder, HasBondOrders};

/// Explicit valence of `site`: the sum of the integer orders of its bonds.
///
/// A double bond counts two, a triple three, and so on; a site with no bonds has
/// valence zero. Returns `None` when any incident bond is
/// [`Aromatic`](BondOrder::Aromatic): a delocalised bond has no localised integer
/// order, leaving the valence undefined until the ring is kekulised.
///
/// # Complexity
///
/// O(d) time and O(1) space, where `d` is the degree of `site`, assuming
/// [`bonds_of`](crate::HasBonds::bonds_of) runs in O(degree).
pub fn valence<M: HasBondOrders>(mol: &M, site: SiteId) -> Option<u32> {
    let mut sum = 0;
    for (bond, _) in mol.bonds_of(site) {
        sum += match mol.bond_order(bond) {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 4,
            BondOrder::Quintuple => 5,
            BondOrder::Hextuple => 6,
            BondOrder::Aromatic => return None,
        };
    }
    Some(sum)
}

#[cfg(test)]
mod tests {
    use super::*;

    use vita_core::HasSites;

    use crate::{BondId, HasBonds};

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
        orders: Vec<BondOrder>,
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

    impl HasBondOrders for Mol {
        fn bond_order(&self, bond: BondId) -> BondOrder {
            let i = self.bonds.iter().position(|&x| x == bond).unwrap();
            self.orders[i]
        }
    }

    fn hub(orders: &[BondOrder]) -> Mol {
        let n = orders.len() as u32;
        Mol {
            sites: (1..=n + 1).map(s).collect(),
            bonds: (1..=n).map(b).collect(),
            endpoints: (2..=n + 1).map(|i| (s(1), s(i))).collect(),
            orders: orders.to_vec(),
        }
    }

    #[test]
    fn a_bondless_site_has_zero_valence() {
        assert_eq!(valence(&hub(&[]), s(1)), Some(0));
    }

    #[test]
    fn each_bond_order_adds_its_multiplicity() {
        for (order, multiplicity) in [
            (BondOrder::Single, 1),
            (BondOrder::Double, 2),
            (BondOrder::Triple, 3),
            (BondOrder::Quadruple, 4),
            (BondOrder::Quintuple, 5),
            (BondOrder::Hextuple, 6),
        ] {
            assert_eq!(valence(&hub(&[order]), s(1)), Some(multiplicity));
        }
    }

    #[test]
    fn valence_sums_all_incident_bonds() {
        let mol = hub(&[BondOrder::Single, BondOrder::Double, BondOrder::Triple]);
        assert_eq!(valence(&mol, s(1)), Some(6));
    }

    #[test]
    fn an_aromatic_bond_leaves_the_valence_undefined() {
        assert_eq!(valence(&hub(&[BondOrder::Aromatic]), s(1)), None);
    }

    #[test]
    fn one_aromatic_among_localised_bonds_is_undefined() {
        let mol = hub(&[BondOrder::Single, BondOrder::Aromatic, BondOrder::Double]);
        assert_eq!(valence(&mol, s(1)), None);
    }

    #[test]
    fn site_valences_sum_to_twice_the_total_bond_order() {
        let mol = hub(&[BondOrder::Single, BondOrder::Double, BondOrder::Triple]);
        let total: u32 = mol.sites().map(|site| valence(&mol, site).unwrap()).sum();
        assert_eq!(total, 2 * (1 + 2 + 3));
    }
}
