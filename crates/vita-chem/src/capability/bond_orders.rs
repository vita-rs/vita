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
