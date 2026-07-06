use vita_core::{HasSites, SiteId};

/// Per-site Lewis-structure bookkeeping: the formal charge of each site.
///
/// Access is by keyed lookup:
/// [`formal_charge`](HasFormalCharges::formal_charge) maps a [`SiteId`] to
/// its formal charge.
/// [`formal_charges`](HasFormalCharges::formal_charges) yields one charge per
/// site in [`sites`](HasSites::sites) order.
///
/// # Contract
///
/// [`formal_charge`](HasFormalCharges::formal_charge) is total over
/// [`sites`](HasSites::sites): every site has exactly one formal charge.
/// [`formal_charges`](HasFormalCharges::formal_charges) yields values in the
/// same order as [`sites`](HasSites::sites).
pub trait HasFormalCharges: HasSites {
    /// Returns the formal charge of `site`, in units of elementary charge.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn formal_charge(&self, site: SiteId) -> i8;

    /// Yields one formal charge per site, in [`sites`](HasSites::sites) order.
    ///
    /// The default implementation looks up
    /// [`formal_charge`](HasFormalCharges::formal_charge) per site; override
    /// it when the charges can be produced directly.
    #[inline]
    fn formal_charges(&self) -> impl Iterator<Item = i8> + '_ {
        self.sites().map(move |site| self.formal_charge(site))
    }
}
