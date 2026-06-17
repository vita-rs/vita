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

#[cfg(test)]
mod tests {
    use super::*;

    fn site(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    struct Bare {
        sites: Vec<SiteId>,
        radical_electrons: Vec<u8>,
    }
    impl HasSites for Bare {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }
    impl HasRadicalElectrons for Bare {
        fn radical_electron(&self, site: SiteId) -> u8 {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.radical_electrons[i]
        }
    }

    struct Columnar {
        sites: Vec<SiteId>,
        radical_electrons: Vec<u8>,
    }
    impl HasSites for Columnar {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }
    impl HasRadicalElectrons for Columnar {
        fn radical_electron(&self, site: SiteId) -> u8 {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.radical_electrons[i]
        }

        fn radical_electrons(&self) -> impl Iterator<Item = u8> + '_ {
            self.radical_electrons.iter().copied()
        }
    }

    fn formyl_radical() -> Bare {
        Bare {
            sites: vec![site(1), site(2), site(3)],
            radical_electrons: vec![1, 0, 0],
        }
    }

    #[test]
    fn radical_electron() {
        let mol = formyl_radical();
        assert_eq!(mol.radical_electron(site(1)), 1);
        assert_eq!(mol.radical_electron(site(2)), 0);
        assert_eq!(mol.radical_electron(site(3)), 0);
    }

    #[test]
    fn radical_electrons() {
        let mol = formyl_radical();
        assert_eq!(mol.radical_electrons().collect::<Vec<_>>(), vec![1, 0, 0]);
    }

    #[test]
    fn radical_electrons_empty() {
        let mol = Bare {
            sites: vec![],
            radical_electrons: vec![],
        };
        assert_eq!(mol.radical_electrons().count(), 0);
    }

    #[test]
    fn override_matches_default() {
        let sites = vec![site(1), site(2), site(3)];
        let counts = vec![1u8, 0, 0];

        let bare = Bare {
            sites: sites.clone(),
            radical_electrons: counts.clone(),
        };
        let col = Columnar {
            sites,
            radical_electrons: counts,
        };

        assert_eq!(
            bare.radical_electrons().collect::<Vec<_>>(),
            col.radical_electrons().collect::<Vec<_>>(),
        );
    }
}
