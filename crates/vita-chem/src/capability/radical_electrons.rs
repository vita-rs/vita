use vita_core::{HasSites, SiteId};

/// Per-site Lewis-structure bookkeeping: the unpaired-electron count of each
/// site.
///
/// Access is by keyed lookup:
/// [`radical_electron`](HasRadicalElectrons::radical_electron) maps a
/// [`SiteId`] to its unpaired-electron count.
/// [`radical_electrons`](HasRadicalElectrons::radical_electrons) yields one
/// count per site in [`sites`](HasSites::sites) order.
///
/// # Contract
///
/// [`radical_electron`](HasRadicalElectrons::radical_electron) is total over
/// [`sites`](HasSites::sites): every site has exactly one unpaired-electron
/// count.
/// [`radical_electrons`](HasRadicalElectrons::radical_electrons) yields values
/// in the same order as [`sites`](HasSites::sites).
pub trait HasRadicalElectrons: HasSites {
    /// Returns the number of unpaired electrons on `site`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn radical_electron(&self, site: SiteId) -> u8;

    /// Yields one unpaired-electron count per site, in
    /// [`sites`](HasSites::sites) order.
    ///
    /// The default implementation looks up
    /// [`radical_electron`](HasRadicalElectrons::radical_electron) per site;
    /// override it when the counts can be produced directly.
    #[inline]
    fn radical_electrons(&self) -> impl Iterator<Item = u8> + '_ {
        self.sites().map(move |site| self.radical_electron(site))
    }
}
