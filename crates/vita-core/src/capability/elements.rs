use crate::{Element, HasSites, SiteId};

/// Per-site chemical identity: the [`Element`] occupying each site.
///
/// Access is by lookup: [`element`](HasElements::element) maps a [`SiteId`] to its
/// element. [`elements`](HasElements::elements) iterates every `(site, element)` pair;
/// a site together with its element constitutes an atom.
///
/// # Contract
///
/// [`element`](HasElements::element) is total over [`sites`](HasSites::sites): every site
/// has exactly one element.
pub trait HasElements: HasSites {
    /// Returns the element occupying `site`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not one of this configuration's [`sites`](HasSites::sites).
    fn element(&self, site: SiteId) -> Element;

    /// Returns an iterator over every `(site, element)` pair.
    ///
    /// Each element is yielded with its [`SiteId`]. The default implementation looks up
    /// [`element`](HasElements::element) per site; override it when the pairs can be
    /// produced directly.
    #[inline]
    fn elements(&self) -> impl Iterator<Item = (SiteId, Element)> + '_ {
        self.sites().map(move |site| (site, self.element(site)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn site(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn hydrogen() -> Element {
        Element::new(1).unwrap()
    }

    fn oxygen() -> Element {
        Element::new(8).unwrap()
    }

    struct Bare {
        sites: Vec<SiteId>,
        elements: Vec<Element>,
    }
    impl HasSites for Bare {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }
    impl HasElements for Bare {
        fn element(&self, site: SiteId) -> Element {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.elements[i]
        }
    }

    struct Columnar {
        sites: Vec<SiteId>,
        elements: Vec<Element>,
    }
    impl HasSites for Columnar {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }
    impl HasElements for Columnar {
        fn element(&self, site: SiteId) -> Element {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.elements[i]
        }

        fn elements(&self) -> impl Iterator<Item = (SiteId, Element)> + '_ {
            self.sites
                .iter()
                .copied()
                .zip(self.elements.iter().copied())
        }
    }

    fn water() -> Bare {
        Bare {
            sites: vec![site(1), site(2), site(3)],
            elements: vec![oxygen(), hydrogen(), hydrogen()],
        }
    }

    #[test]
    fn element() {
        let mol = water();
        assert_eq!(mol.element(site(1)), oxygen());
        assert_eq!(mol.element(site(2)), hydrogen());
    }

    #[test]
    fn elements() {
        let mol = water();
        assert_eq!(
            mol.elements().collect::<Vec<_>>(),
            vec![
                (site(1), oxygen()),
                (site(2), hydrogen()),
                (site(3), hydrogen())
            ]
        );
    }

    #[test]
    fn elements_empty() {
        let mol = Bare {
            sites: vec![],
            elements: vec![],
        };
        assert_eq!(mol.elements().count(), 0);
    }

    #[test]
    fn override_matches_default() {
        use std::collections::BTreeMap;

        let sites = vec![site(1), site(2), site(3)];
        let elements = vec![oxygen(), hydrogen(), hydrogen()];
        let bare = Bare {
            sites: sites.clone(),
            elements: elements.clone(),
        };
        let columnar = Columnar { sites, elements };

        let bare_elems: BTreeMap<_, _> = bare.elements().collect();
        let columnar_elems: BTreeMap<_, _> = columnar.elements().collect();
        assert_eq!(bare_elems, columnar_elems);
    }
}
