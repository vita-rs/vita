use crate::{Element, HasSites, SiteId};

/// Per-site chemical identity: the [`Element`] occupying each site.
///
/// Access is by lookup: [`element`](HasElements::element) maps a [`SiteId`] to its
/// element. [`elements`](HasElements::elements) iterates every `(site, element)` pair;
/// a site together with its element constitutes an atom.
///
/// # Contract
///
/// [`element`](HasElements::element) is total over [`sites`](HasSites::sites): every site
/// has exactly one element.
pub trait HasElements: HasSites {
    /// Returns the element occupying `site`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not one of this configuration's [`sites`](HasSites::sites).
    fn element(&self, site: SiteId) -> Element;

    /// Returns an iterator over every `(site, element)` pair.
    ///
    /// Each element is yielded with its [`SiteId`]. The default implementation looks up
    /// [`element`](HasElements::element) per site; override it when the pairs can be
    /// produced directly.
    #[inline]
    fn elements(&self) -> impl Iterator<Item = (SiteId, Element)> + '_ {
        self.sites().map(move |site| (site, self.element(site)))
    }
}
