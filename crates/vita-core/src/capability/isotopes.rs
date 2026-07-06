use crate::{HasElements, Isotope, SiteId};

/// Per-site nuclear identity: the [`Isotope`] occupying each site.
///
/// Access is by keyed lookup: [`isotope`](HasIsotopes::isotope) maps a [`SiteId`] to
/// its isotope. An isotope refines an [`Element`](crate::Element) with a mass number,
/// so this capability builds on [`HasElements`].
/// [`isotopes`](HasIsotopes::isotopes) yields one isotope per site in
/// [`sites`](crate::HasSites::sites) order.
///
/// # Contract
///
/// [`isotope`](HasIsotopes::isotope) is total over [`sites`](crate::HasSites::sites):
/// every site has exactly one isotope, whose [`element`](Isotope::element) equals that
/// site's [`element`](HasElements::element).
/// [`isotopes`](HasIsotopes::isotopes) yields values in the same order as
/// [`sites`](crate::HasSites::sites).
pub trait HasIsotopes: HasElements {
    /// Returns the isotope occupying `site`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](crate::HasSites::sites).
    fn isotope(&self, site: SiteId) -> Isotope;

    /// Yields one isotope per site, in [`sites`](crate::HasSites::sites) order.
    ///
    /// The default implementation looks up [`isotope`](HasIsotopes::isotope) per site;
    /// override it when the isotopes can be produced directly.
    #[inline]
    fn isotopes(&self) -> impl Iterator<Item = Isotope> + '_ {
        self.sites().map(move |site| self.isotope(site))
    }
}
