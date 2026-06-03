use crate::{HasElements, Isotope, SiteId};

/// Per-site nuclear identity: the [`Isotope`] occupying each site.
///
/// Access is by lookup: [`isotope`](HasIsotopes::isotope) maps a [`SiteId`] to its
/// isotope. An isotope refines an [`Element`](crate::Element) with a mass number, so this
/// capability builds on [`HasElements`].
///
/// # Contract
///
/// [`isotope`](HasIsotopes::isotope) is total over [`sites`](crate::HasSites::sites):
/// every site has exactly one isotope, whose [`element`](Isotope::element) equals that
/// site's [`element`](HasElements::element).
pub trait HasIsotopes: HasElements {
    /// Returns the isotope occupying `site`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not one of this configuration's
    /// [`sites`](crate::HasSites::sites).
    fn isotope(&self, site: SiteId) -> Isotope;

    /// Returns an iterator over every `(site, isotope)` pair.
    ///
    /// Each isotope is yielded with its [`SiteId`]. The default implementation looks up
    /// [`isotope`](HasIsotopes::isotope) per site; override it when the pairs can be
    /// produced directly.
    #[inline]
    fn isotopes(&self) -> impl Iterator<Item = (SiteId, Isotope)> + '_ {
        self.sites().map(move |site| (site, self.isotope(site)))
    }
}
