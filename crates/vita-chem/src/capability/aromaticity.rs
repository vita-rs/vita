use vita_core::SiteId;

use crate::{BondId, HasBonds};

/// Per-bond aromaticity: whether each bond lies in an aromatic π system.
///
/// Access is by keyed lookup:
/// [`is_aromatic`](HasAromaticity::is_aromatic) maps a [`BondId`] to whether
/// that bond is aromatic.
/// [`is_aromatic_site`](HasAromaticity::is_aromatic_site) reports the same for
/// a [`SiteId`], derived from the bonds incident to it.
///
/// # Contract
///
/// [`is_aromatic`](HasAromaticity::is_aromatic) is total over
/// [`bonds`](HasBonds::bonds): every bond is either aromatic or not.
pub trait HasAromaticity: HasBonds {
    /// Returns whether `bond` is aromatic.
    ///
    /// # Panics
    ///
    /// Panics if `bond` is not in [`bonds`](HasBonds::bonds).
    fn is_aromatic(&self, bond: BondId) -> bool;

    /// Returns whether `site` lies in an aromatic system.
    ///
    /// A site is aromatic exactly when one of its bonds is aromatic; in
    /// biphenyl the two sites joined by the central bond are thus aromatic
    /// while that bond is not. The default implementation scans the bonds
    /// incident to `site`; override it when site aromaticity can be determined
    /// directly.
    #[inline]
    fn is_aromatic_site(&self, site: SiteId) -> bool {
        self.bonds_of(site).any(|(bond, _)| self.is_aromatic(bond))
    }
}
