use crate::{Element, HasSites, SiteId};

/// Per-site chemical identity: the [`Element`] occupying each site.
///
/// Access is by keyed lookup: [`element`](HasElements::element) maps a [`SiteId`] to
/// its element. [`elements`](HasElements::elements) yields one element per site in
/// [`sites`](HasSites::sites) order; a site together with its element constitutes an
/// atom.
///
/// # Contract
///
/// [`element`](HasElements::element) is total over [`sites`](HasSites::sites): every site
/// has exactly one element.
/// [`elements`](HasElements::elements) yields values in the same order as
/// [`sites`](HasSites::sites).
pub trait HasElements: HasSites {
    /// Returns the element occupying `site`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn element(&self, site: SiteId) -> Element;

    /// Yields one element per site, in [`sites`](HasSites::sites) order.
    ///
    /// The default implementation looks up [`element`](HasElements::element) per site;
    /// override it when the elements can be produced directly.
    #[inline]
    fn elements(&self) -> impl Iterator<Item = Element> + '_ {
        self.sites().map(move |site| self.element(site))
    }
}
