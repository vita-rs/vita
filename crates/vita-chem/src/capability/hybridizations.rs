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

#[cfg(test)]
mod tests {
    use super::*;

    fn site(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    struct Bare {
        sites: Vec<SiteId>,
        hybridizations: Vec<Hybridization>,
    }
    impl HasSites for Bare {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }
    impl HasHybridizations for Bare {
        fn hybridization(&self, site: SiteId) -> Hybridization {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.hybridizations[i]
        }
    }

    struct Columnar {
        sites: Vec<SiteId>,
        hybridizations: Vec<Hybridization>,
    }
    impl HasSites for Columnar {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }
    impl HasHybridizations for Columnar {
        fn hybridization(&self, site: SiteId) -> Hybridization {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.hybridizations[i]
        }

        fn hybridizations(&self) -> impl Iterator<Item = Hybridization> + '_ {
            self.hybridizations.iter().copied()
        }
    }

    fn hcn() -> Bare {
        Bare {
            sites: vec![site(1), site(2), site(3)],
            hybridizations: vec![Hybridization::S, Hybridization::Sp, Hybridization::Sp],
        }
    }

    #[test]
    fn hybridization() {
        let mol = hcn();
        assert_eq!(mol.hybridization(site(1)), Hybridization::S);
        assert_eq!(mol.hybridization(site(2)), Hybridization::Sp);
        assert_eq!(mol.hybridization(site(3)), Hybridization::Sp);
    }

    #[test]
    fn hybridizations() {
        let mol = hcn();
        assert_eq!(
            mol.hybridizations().collect::<Vec<_>>(),
            vec![Hybridization::S, Hybridization::Sp, Hybridization::Sp],
        );
    }

    #[test]
    fn hybridizations_empty() {
        let mol = Bare {
            sites: vec![],
            hybridizations: vec![],
        };
        assert_eq!(mol.hybridizations().count(), 0);
    }

    #[test]
    fn all_hybridization_variants() {
        let mol = Bare {
            sites: (1..=9).map(site).collect(),
            hybridizations: vec![
                Hybridization::S,
                Hybridization::Sp,
                Hybridization::Sp2,
                Hybridization::Sp3,
                Hybridization::Sp2d,
                Hybridization::Sp3d,
                Hybridization::Sp3d2,
                Hybridization::Sp3d3,
                Hybridization::Other,
            ],
        };
        assert_eq!(mol.hybridization(site(1)), Hybridization::S);
        assert_eq!(mol.hybridization(site(2)), Hybridization::Sp);
        assert_eq!(mol.hybridization(site(3)), Hybridization::Sp2);
        assert_eq!(mol.hybridization(site(4)), Hybridization::Sp3);
        assert_eq!(mol.hybridization(site(5)), Hybridization::Sp2d);
        assert_eq!(mol.hybridization(site(6)), Hybridization::Sp3d);
        assert_eq!(mol.hybridization(site(7)), Hybridization::Sp3d2);
        assert_eq!(mol.hybridization(site(8)), Hybridization::Sp3d3);
        assert_eq!(mol.hybridization(site(9)), Hybridization::Other);
    }

    #[test]
    fn override_matches_default() {
        let sites = vec![site(1), site(2), site(3)];
        let hybridizations_vec = vec![Hybridization::S, Hybridization::Sp, Hybridization::Sp];

        let bare = Bare {
            sites: sites.clone(),
            hybridizations: hybridizations_vec.clone(),
        };
        let col = Columnar {
            sites,
            hybridizations: hybridizations_vec,
        };

        assert_eq!(
            bare.hybridizations().collect::<Vec<_>>(),
            col.hybridizations().collect::<Vec<_>>(),
        );
    }
}
