use crate::SiteId;

/// The identity skeleton: the [`SiteId`]s every per-site datum is keyed on.
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
    /// Returns an iterator over the identifier of every site.
    fn sites(&self) -> impl Iterator<Item = SiteId> + '_;

    /// Returns the number of sites.
    ///
    /// The default implementation consumes [`sites`](HasSites::sites); override it when
    /// the count is known in `O(1)`.
    #[inline]
    fn site_count(&self) -> usize {
        self.sites().count()
    }

    /// Returns whether `site` is in [`sites`](HasSites::sites).
    ///
    /// The default implementation scans [`sites`](HasSites::sites); override it when
    /// membership can be decided in better than `O(n)`.
    #[inline]
    fn contains_site(&self, site: SiteId) -> bool {
        self.sites().any(|s| s == site)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn site(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    struct Bare(Vec<SiteId>);
    impl HasSites for Bare {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.0.iter().copied()
        }
    }

    struct Indexed(Vec<SiteId>);
    impl HasSites for Indexed {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.0.iter().copied()
        }

        fn site_count(&self) -> usize {
            self.0.len()
        }

        fn contains_site(&self, site: SiteId) -> bool {
            self.0.binary_search(&site).is_ok()
        }
    }

    #[test]
    fn sites() {
        let sys = Bare(vec![site(1), site(4), site(9)]);
        assert_eq!(
            sys.sites().collect::<Vec<_>>(),
            vec![site(1), site(4), site(9)]
        );
    }

    #[test]
    fn sites_empty() {
        let sys = Bare(vec![]);
        assert_eq!(sys.sites().count(), 0);
    }

    #[test]
    fn site_count() {
        let sys = Bare(vec![site(1), site(2), site(3)]);
        assert_eq!(sys.site_count(), 3);
    }

    #[test]
    fn site_count_empty_is_zero() {
        let sys = Bare(vec![]);
        assert_eq!(sys.site_count(), 0);
    }

    #[test]
    fn contains_site() {
        let sys = Bare(vec![site(1), site(5)]);
        assert!(sys.contains_site(site(5)));
        assert!(!sys.contains_site(site(2)));
    }

    #[test]
    fn override_matches_default() {
        let ids = vec![site(1), site(3), site(7)];
        let bare = Bare(ids.clone());
        let indexed = Indexed(ids);

        assert_eq!(bare.site_count(), indexed.site_count());
        for n in 0..10 {
            let id = site(n + 1);
            assert_eq!(bare.contains_site(id), indexed.contains_site(id));
        }
    }
}
