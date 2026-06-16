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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BondId;
    use vita_core::HasSites;

    fn site(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn bond(n: u32) -> BondId {
        BondId::new(n).unwrap()
    }

    struct Bare {
        sites: Vec<SiteId>,
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
        hybridizations: Vec<Hybridization>,
    }

    impl HasSites for Bare {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }

    impl HasBonds for Bare {
        fn bonds(&self) -> impl Iterator<Item = BondId> + '_ {
            self.bonds.iter().copied()
        }

        fn bond_endpoints(&self, b: BondId) -> (SiteId, SiteId) {
            let i = self.bonds.iter().position(|&x| x == b).unwrap();
            self.endpoints[i]
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
        bonds: Vec<BondId>,
        endpoints: Vec<(SiteId, SiteId)>,
        hybridizations: Vec<Hybridization>,
    }

    impl HasSites for Columnar {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }

    impl HasBonds for Columnar {
        fn bonds(&self) -> impl Iterator<Item = BondId> + '_ {
            self.bonds.iter().copied()
        }

        fn bond_endpoints(&self, b: BondId) -> (SiteId, SiteId) {
            let i = self.bonds.iter().position(|&x| x == b).unwrap();
            self.endpoints[i]
        }
    }

    impl HasHybridizations for Columnar {
        fn hybridization(&self, site: SiteId) -> Hybridization {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.hybridizations[i]
        }

        fn hybridizations(&self) -> impl Iterator<Item = (SiteId, Hybridization)> + '_ {
            self.sites
                .iter()
                .copied()
                .zip(self.hybridizations.iter().copied())
        }
    }

    fn hcn() -> Bare {
        Bare {
            sites: vec![site(1), site(2), site(3)],
            bonds: vec![bond(1), bond(2)],
            endpoints: vec![(site(1), site(2)), (site(2), site(3))],
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
            vec![
                (site(1), Hybridization::S),
                (site(2), Hybridization::Sp),
                (site(3), Hybridization::Sp),
            ]
        );
    }

    #[test]
    fn hybridizations_empty() {
        let mol = Bare {
            sites: vec![],
            bonds: vec![],
            endpoints: vec![],
            hybridizations: vec![],
        };
        assert_eq!(mol.hybridizations().count(), 0);
    }

    #[test]
    fn all_hybridization_variants() {
        let mol = Bare {
            sites: vec![
                site(1),
                site(2),
                site(3),
                site(4),
                site(5),
                site(6),
                site(7),
                site(8),
                site(9),
            ],
            bonds: vec![
                bond(1),
                bond(2),
                bond(3),
                bond(4),
                bond(5),
                bond(6),
                bond(7),
                bond(8),
            ],
            endpoints: vec![
                (site(1), site(2)),
                (site(2), site(3)),
                (site(3), site(4)),
                (site(4), site(5)),
                (site(5), site(6)),
                (site(6), site(7)),
                (site(7), site(8)),
                (site(8), site(9)),
            ],
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
        use std::collections::BTreeMap;

        let sites = vec![site(1), site(2), site(3)];
        let bonds_vec = vec![bond(1), bond(2)];
        let endpoints_vec = vec![(site(1), site(2)), (site(2), site(3))];
        let hybridizations_vec = vec![Hybridization::S, Hybridization::Sp, Hybridization::Sp];

        let bare = Bare {
            sites: sites.clone(),
            bonds: bonds_vec.clone(),
            endpoints: endpoints_vec.clone(),
            hybridizations: hybridizations_vec.clone(),
        };
        let col = Columnar {
            sites,
            bonds: bonds_vec,
            endpoints: endpoints_vec,
            hybridizations: hybridizations_vec,
        };

        let bare_hybs: BTreeMap<_, _> = bare.hybridizations().collect();
        let col_hybs: BTreeMap<_, _> = col.hybridizations().collect();
        assert_eq!(bare_hybs, col_hybs);
    }
}
