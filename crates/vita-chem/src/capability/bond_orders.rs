use crate::{BondId, BondOrder, HasBonds};

/// Per-bond multiplicity: the [`BondOrder`] of each bond.
///
/// Access is by lookup: [`bond_order`](HasBondOrders::bond_order) maps a
/// [`BondId`] to its order. [`bond_orders`](HasBondOrders::bond_orders)
/// iterates every `(bond, order)` pair.
///
/// # Contract
///
/// [`bond_order`](HasBondOrders::bond_order) is total over
/// [`bonds`](HasBonds::bonds): every bond has exactly one order.
pub trait HasBondOrders: HasBonds {
    /// Returns the order of `bond`.
    ///
    /// # Panics
    ///
    /// Panics if `bond` is not in [`bonds`](HasBonds::bonds).
    fn bond_order(&self, bond: BondId) -> BondOrder;

    /// Returns an iterator over every `(bond, order)` pair.
    ///
    /// Each order is yielded with its [`BondId`]. The default implementation looks
    /// up [`bond_order`](HasBondOrders::bond_order) per bond; override it when
    /// the pairs can be produced directly.
    #[inline]
    fn bond_orders(&self) -> impl Iterator<Item = (BondId, BondOrder)> + '_ {
        self.bonds().map(move |bond| (bond, self.bond_order(bond)))
    }
}
