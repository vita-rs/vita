use crate::SiteId;

/// The identity skeleton of a configuration: the [`SiteId`]s it is built from.
///
/// Every per-site datum in the ecosystem — element, position, mass, velocity — is a
/// column keyed on [`SiteId`]. `HasSites` enumerates those keys, and is therefore the
/// supertrait of every per-site capability: a type cannot expose data *about* sites
/// without first declaring *which* sites exist.
///
/// # Contract
///
/// [`sites`](HasSites::sites) yields each identifier exactly once, with no duplicates.
pub trait HasSites {
    /// Returns an iterator over the identifier of every site in the configuration.
    fn sites(&self) -> impl Iterator<Item = SiteId> + '_;

    /// Returns the number of sites.
    ///
    /// The default implementation consumes [`sites`](HasSites::sites); override it when
    /// the count is known in `O(1)`.
    #[inline]
    fn site_count(&self) -> usize {
        self.sites().count()
    }

    /// Returns whether `site` belongs to this configuration.
    ///
    /// The default implementation scans [`sites`](HasSites::sites); override it when
    /// membership can be decided in better than `O(n)`.
    #[inline]
    fn contains_site(&self, site: SiteId) -> bool {
        self.sites().any(|s| s == site)
    }
}
