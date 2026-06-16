use vita_core::units::charge::{Charge, ChargeUnit};
use vita_core::{HasSites, Scalar, SiteId};

/// Per-site partial charge: the [`Charge`] of each site.
///
/// Access is by lookup:
/// [`partial_charge`](HasPartialCharges::partial_charge) maps a [`SiteId`]
/// to its charge, in any requested [unit](ChargeUnit).
/// [`partial_charges`](HasPartialCharges::partial_charges) iterates every
/// `(site, charge)` pair.
///
/// # Contract
///
/// [`partial_charge`](HasPartialCharges::partial_charge) is total over
/// [`sites`](HasSites::sites): every site has exactly one partial charge.
pub trait HasPartialCharges<V: Scalar>: HasSites {
    /// Returns the partial charge of `site`, in unit `U`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn partial_charge<U: ChargeUnit>(&self, site: SiteId) -> Charge<V, U>;

    /// Returns an iterator over every `(site, partial_charge)` pair, each
    /// charge in unit `U`.
    ///
    /// Each charge is yielded with its [`SiteId`]. The default
    /// implementation looks up
    /// [`partial_charge`](HasPartialCharges::partial_charge) per site;
    /// override it when the pairs can be produced directly.
    #[inline]
    fn partial_charges<U: ChargeUnit>(&self) -> impl Iterator<Item = (SiteId, Charge<V, U>)> + '_ {
        self.sites()
            .map(move |site| (site, self.partial_charge::<U>(site)))
    }
}
