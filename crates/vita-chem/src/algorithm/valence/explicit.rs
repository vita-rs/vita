use vita_core::SiteId;

use crate::{BondOrder, HasBondOrders};

/// Explicit valence of `site`: the sum of the orders of its bonds.
///
/// A double bond counts as two, a triple as three, and so on; a site with no
/// bonds has valence zero. Returns `None` when any incident bond is
/// [`Aromatic`](BondOrder::Aromatic): a delocalised bond has no definite
/// localised order, leaving the valence undefined until the ring is kekulised.
///
/// # Complexity
///
/// O(degree) time.
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
    use crate::{BondId, HasBonds};
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

    fn methane() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5)],
            bonds: vec![b(1), b(2), b(3), b(4)],
            endpoints: vec![(s(1), s(2)), (s(1), s(3)), (s(1), s(4)), (s(1), s(5))],
            orders: vec![BondOrder::Single; 4],
        }
    }

    fn carbon_dioxide() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            bonds: vec![b(1), b(2)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3))],
            orders: vec![BondOrder::Double, BondOrder::Double],
        }
    }

    fn hydrogen_cyanide() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3)],
            bonds: vec![b(1), b(2)],
            endpoints: vec![(s(1), s(2)), (s(2), s(3))],
            orders: vec![BondOrder::Single, BondOrder::Triple],
        }
    }

    fn benzene() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6)],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6)],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(6)),
                (s(6), s(1)),
            ],
            orders: vec![BondOrder::Aromatic; 6],
        }
    }

    fn methylbenzene() -> Mol {
        Mol {
            sites: vec![s(1), s(2), s(3), s(4), s(5), s(6), s(7)],
            bonds: vec![b(1), b(2), b(3), b(4), b(5), b(6), b(7)],
            endpoints: vec![
                (s(1), s(2)),
                (s(2), s(3)),
                (s(3), s(4)),
                (s(4), s(5)),
                (s(5), s(6)),
                (s(6), s(1)),
                (s(1), s(7)),
            ],
            orders: vec![
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Aromatic,
                BondOrder::Single,
            ],
        }
    }

    #[test]
    fn isolated_site_is_zero() {
        let mol = Mol {
            sites: vec![s(1)],
            bonds: vec![],
            endpoints: vec![],
            orders: vec![],
        };
        assert_eq!(valence(&mol, s(1)), Some(0));
    }

    #[test]
    fn single_bonds() {
        let mol = methane();
        assert_eq!(valence(&mol, s(1)), Some(4));
        assert_eq!(valence(&mol, s(2)), Some(1));
    }

    #[test]
    fn double_counts_as_two() {
        let mol = carbon_dioxide();
        assert_eq!(valence(&mol, s(2)), Some(4));
        assert_eq!(valence(&mol, s(1)), Some(2));
    }

    #[test]
    fn triple_counts_as_three() {
        let mol = hydrogen_cyanide();
        assert_eq!(valence(&mol, s(1)), Some(1));
        assert_eq!(valence(&mol, s(2)), Some(4));
        assert_eq!(valence(&mol, s(3)), Some(3));
    }

    #[test]
    fn high_orders() {
        for (order, expected) in [
            (BondOrder::Quadruple, 4),
            (BondOrder::Quintuple, 5),
            (BondOrder::Hextuple, 6),
        ] {
            let mol = Mol {
                sites: vec![s(1), s(2)],
                bonds: vec![b(1)],
                endpoints: vec![(s(1), s(2))],
                orders: vec![order],
            };
            assert_eq!(valence(&mol, s(1)), Some(expected));
        }
    }

    #[test]
    fn aromatic_is_undefined() {
        let mol = benzene();
        assert_eq!(valence(&mol, s(1)), None);
    }

    #[test]
    fn aromatic_affects_only_incident_sites() {
        let mol = methylbenzene();
        assert_eq!(valence(&mol, s(7)), Some(1));
        assert_eq!(valence(&mol, s(1)), None);
    }
}
