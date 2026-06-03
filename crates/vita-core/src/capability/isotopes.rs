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
    /// Panics if `site` is not in [`sites`](crate::HasSites::sites).
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Element, HasSites};

    fn site(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn protium() -> Isotope {
        Isotope::new(Element::new(1).unwrap(), 1).unwrap()
    }

    fn oxygen_16() -> Isotope {
        Isotope::new(Element::new(8).unwrap(), 16).unwrap()
    }

    struct Bare {
        sites: Vec<SiteId>,
        isotopes: Vec<Isotope>,
    }
    impl HasSites for Bare {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }
    impl HasElements for Bare {
        fn element(&self, site: SiteId) -> Element {
            self.isotope(site).element()
        }
    }
    impl HasIsotopes for Bare {
        fn isotope(&self, site: SiteId) -> Isotope {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.isotopes[i]
        }
    }

    struct Columnar {
        sites: Vec<SiteId>,
        isotopes: Vec<Isotope>,
    }
    impl HasSites for Columnar {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }
    impl HasElements for Columnar {
        fn element(&self, site: SiteId) -> Element {
            self.isotope(site).element()
        }
    }
    impl HasIsotopes for Columnar {
        fn isotope(&self, site: SiteId) -> Isotope {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.isotopes[i]
        }

        fn isotopes(&self) -> impl Iterator<Item = (SiteId, Isotope)> + '_ {
            self.sites
                .iter()
                .copied()
                .zip(self.isotopes.iter().copied())
        }
    }

    fn water() -> Bare {
        Bare {
            sites: vec![site(1), site(2), site(3)],
            isotopes: vec![oxygen_16(), protium(), protium()],
        }
    }

    #[test]
    fn isotope() {
        let mol = water();
        assert_eq!(mol.isotope(site(1)), oxygen_16());
        assert_eq!(mol.isotope(site(2)), protium());
    }

    #[test]
    fn isotopes() {
        let mol = water();
        assert_eq!(
            mol.isotopes().collect::<Vec<_>>(),
            vec![
                (site(1), oxygen_16()),
                (site(2), protium()),
                (site(3), protium())
            ]
        );
    }

    #[test]
    fn isotopes_empty() {
        let mol = Bare {
            sites: vec![],
            isotopes: vec![],
        };
        assert_eq!(mol.isotopes().count(), 0);
    }

    #[test]
    fn override_matches_default() {
        use std::collections::BTreeMap;

        let sites = vec![site(1), site(2), site(3)];
        let isotopes = vec![oxygen_16(), protium(), protium()];
        let bare = Bare {
            sites: sites.clone(),
            isotopes: isotopes.clone(),
        };
        let columnar = Columnar { sites, isotopes };

        let bare_isos: BTreeMap<_, _> = bare.isotopes().collect();
        let columnar_isos: BTreeMap<_, _> = columnar.isotopes().collect();
        assert_eq!(bare_isos, columnar_isos);
    }
}
