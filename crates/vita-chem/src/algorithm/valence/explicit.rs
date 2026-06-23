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
