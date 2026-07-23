use vita_core::{HasSites, SiteId};

use crate::Hybridization;

/// Per-site orbital geometry: the [`Hybridization`] of each site.
///
/// Access is by keyed lookup:
/// [`hybridization`](HasHybridizations::hybridization) maps a [`SiteId`] to
/// its hybridization state.
/// [`hybridizations`](HasHybridizations::hybridizations) yields one
/// hybridization per site in [`sites`](HasSites::sites) order.
///
/// # Contract
///
/// [`hybridization`](HasHybridizations::hybridization) is total over
/// [`sites`](HasSites::sites): every site has exactly one hybridization.
/// [`hybridizations`](HasHybridizations::hybridizations) yields values in the
/// same order as [`sites`](HasSites::sites).
pub trait HasHybridizations: HasSites {
    /// Returns the hybridization of `site`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn hybridization(&self, site: SiteId) -> Hybridization;

    /// Yields one hybridization per site, in [`sites`](HasSites::sites) order.
    ///
    /// The default implementation looks up
    /// [`hybridization`](HasHybridizations::hybridization) per site; override
    /// it when the hybridizations can be produced directly.
    #[inline]
    fn hybridizations(&self) -> impl Iterator<Item = Hybridization> + '_ {
        self.sites().map(move |site| self.hybridization(site))
    }
}
