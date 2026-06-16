use vita_core::SiteId;

use crate::{HasBonds, Hybridization};

/// Per-site orbital geometry: the [`Hybridization`] of each site.
///
/// Access is by lookup:
/// [`hybridization`](HasHybridizations::hybridization) maps a [`SiteId`] to
/// its hybridization state.
/// [`hybridizations`](HasHybridizations::hybridizations) iterates every
/// `(site, hybridization)` pair.
///
/// # Contract
///
/// [`hybridization`](HasHybridizations::hybridization) is total over
/// [`sites`](vita_core::HasSites::sites): every site has exactly one hybridization.
pub trait HasHybridizations: HasBonds {
    /// Returns the hybridization of `site`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](vita_core::HasSites::sites).
    fn hybridization(&self, site: SiteId) -> Hybridization;

    /// Returns an iterator over every `(site, hybridization)` pair.
    ///
    /// Each hybridization is yielded with its [`SiteId`]. The default
    /// implementation looks up
    /// [`hybridization`](HasHybridizations::hybridization) per site;
    /// override it when the pairs can be produced directly.
    #[inline]
    fn hybridizations(&self) -> impl Iterator<Item = (SiteId, Hybridization)> + '_ {
        self.sites()
            .map(move |site| (site, self.hybridization(site)))
    }
}
