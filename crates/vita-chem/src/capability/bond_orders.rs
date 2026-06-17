use crate::{BondId, BondOrder, HasBonds};

/// Per-bond multiplicity: the [`BondOrder`] of each bond.
///
/// Access is by keyed lookup: [`bond_order`](HasBondOrders::bond_order) maps a
/// [`BondId`] to its order.
/// [`bond_orders`](HasBondOrders::bond_orders) yields one order per bond in
/// [`bonds`](HasBonds::bonds) order.
///
/// # Contract
///
/// [`bond_order`](HasBondOrders::bond_order) is total over
/// [`bonds`](HasBonds::bonds): every bond has exactly one order.
/// [`bond_orders`](HasBondOrders::bond_orders) yields values in the same
/// order as [`bonds`](HasBonds::bonds).
pub trait HasBondOrders: HasBonds {
    /// Returns the order of `bond`.
    ///
    /// # Panics
    ///
    /// Panics if `bond` is not in [`bonds`](HasBonds::bonds).
    fn bond_order(&self, bond: BondId) -> BondOrder;

    /// Yields one order per bond, in [`bonds`](HasBonds::bonds) order.
    ///
    /// The default implementation looks up
    /// [`bond_order`](HasBondOrders::bond_order) per bond; override it when
    /// the orders can be produced directly.
    #[inline]
    fn bond_orders(&self) -> impl Iterator<Item = BondOrder> + '_ {
        self.bonds().map(move |bond| self.bond_order(bond))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vita_core::{HasSites, SiteId};

    fn site(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn bond(n: u32) -> BondId {
        BondId::new(n).unwrap()
    }

    struct Bare {
        sites: Vec<SiteId>,
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
        orders: Vec<BondOrder>,
    }
    impl HasSites for Bare {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }
    impl HasBonds for Bare {
        fn bonds(&self) -> impl Iterator<Item = BondId> + '_ {
            self.bonds.iter().copied()
        }

        fn bond_endpoints(&self, b: BondId) -> (SiteId, SiteId) {
            let i = self.bonds.iter().position(|&x| x == b).unwrap();
            self.endpoints[i]
        }
    }
    impl HasBondOrders for Bare {
        fn bond_order(&self, b: BondId) -> BondOrder {
            let i = self.bonds.iter().position(|&x| x == b).unwrap();
            self.orders[i]
        }
    }

    struct Columnar {
        sites: Vec<SiteId>,
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
        orders: Vec<BondOrder>,
    }
    impl HasSites for Columnar {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }
    impl HasBonds for Columnar {
        fn bonds(&self) -> impl Iterator<Item = BondId> + '_ {
            self.bonds.iter().copied()
        }

        fn bond_endpoints(&self, b: BondId) -> (SiteId, SiteId) {
            let i = self.bonds.iter().position(|&x| x == b).unwrap();
            self.endpoints[i]
        }
    }
    impl HasBondOrders for Columnar {
        fn bond_order(&self, b: BondId) -> BondOrder {
            let i = self.bonds.iter().position(|&x| x == b).unwrap();
            self.orders[i]
        }

        fn bond_orders(&self) -> impl Iterator<Item = BondOrder> + '_ {
            self.orders.iter().copied()
        }
    }

    fn hcn() -> Bare {
        Bare {
            sites: vec![site(1), site(2), site(3)],
            bonds: vec![bond(1), bond(2)],
            endpoints: vec![(site(1), site(2)), (site(2), site(3))],
            orders: vec![BondOrder::Single, BondOrder::Triple],
        }
    }

    #[test]
    fn bond_order() {
        let mol = hcn();
        assert_eq!(mol.bond_order(bond(1)), BondOrder::Single);
        assert_eq!(mol.bond_order(bond(2)), BondOrder::Triple);
    }

    #[test]
    fn bond_orders() {
        let mol = hcn();
        assert_eq!(
            mol.bond_orders().collect::<Vec<_>>(),
            vec![BondOrder::Single, BondOrder::Triple],
        );
    }

    #[test]
    fn bond_orders_empty() {
        let mol = Bare {
            sites: vec![],
            bonds: vec![],
            endpoints: vec![],
            orders: vec![],
        };
        assert_eq!(mol.bond_orders().count(), 0);
    }

    #[test]
    fn all_order_variants() {
        let mol = Bare {
            sites: vec![
                site(1),
                site(2),
                site(3),
                site(4),
                site(5),
                site(6),
                site(7),
            ],
            bonds: vec![bond(1), bond(2), bond(3), bond(4), bond(5), bond(6)],
            endpoints: vec![
                (site(1), site(2)),
                (site(2), site(3)),
                (site(3), site(4)),
                (site(4), site(5)),
                (site(5), site(6)),
                (site(6), site(7)),
            ],
            orders: vec![
                BondOrder::Double,
                BondOrder::Aromatic,
                BondOrder::Triple,
                BondOrder::Quadruple,
                BondOrder::Quintuple,
                BondOrder::Hextuple,
            ],
        };
        assert_eq!(mol.bond_order(bond(1)), BondOrder::Double);
        assert_eq!(mol.bond_order(bond(2)), BondOrder::Aromatic);
        assert_eq!(mol.bond_order(bond(3)), BondOrder::Triple);
        assert_eq!(mol.bond_order(bond(4)), BondOrder::Quadruple);
        assert_eq!(mol.bond_order(bond(5)), BondOrder::Quintuple);
        assert_eq!(mol.bond_order(bond(6)), BondOrder::Hextuple);
    }

    #[test]
    fn override_matches_default() {
        let sites = vec![site(1), site(2), site(3)];
        let bonds_vec = vec![bond(1), bond(2)];
        let endpoints_vec = vec![(site(1), site(2)), (site(2), site(3))];
        let orders_vec = vec![BondOrder::Single, BondOrder::Triple];

        let bare = Bare {
            sites: sites.clone(),
            bonds: bonds_vec.clone(),
            endpoints: endpoints_vec.clone(),
            orders: orders_vec.clone(),
        };
        let col = Columnar {
            sites,
            bonds: bonds_vec,
            endpoints: endpoints_vec,
            orders: orders_vec,
        };

        assert_eq!(
            bare.bond_orders().collect::<Vec<_>>(),
            col.bond_orders().collect::<Vec<_>>(),
        );
    }
}
