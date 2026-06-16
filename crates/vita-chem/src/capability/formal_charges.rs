use vita_core::{HasSites, SiteId};

/// Per-site Lewis-structure bookkeeping: the formal charge of each site.
///
/// Access is by lookup:
/// [`formal_charge`](HasFormalCharges::formal_charge) maps a [`SiteId`] to
/// its formal charge.
/// [`formal_charges`](HasFormalCharges::formal_charges) iterates every
/// `(site, formal_charge)` pair.
///
/// # Contract
///
/// [`formal_charge`](HasFormalCharges::formal_charge) is total over
/// [`sites`](HasSites::sites): every site has exactly one formal charge.
pub trait HasFormalCharges: HasSites {
    /// Returns the formal charge of `site`, in units of elementary charge.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn formal_charge(&self, site: SiteId) -> i8;

    /// Returns an iterator over every `(site, formal_charge)` pair.
    ///
    /// Each charge is yielded with its [`SiteId`]. The default implementation
    /// looks up [`formal_charge`](HasFormalCharges::formal_charge) per site;
    /// override it when the pairs can be produced directly.
    #[inline]
    fn formal_charges(&self) -> impl Iterator<Item = (SiteId, i8)> + '_ {
        self.sites()
            .map(move |site| (site, self.formal_charge(site)))
    }
}
