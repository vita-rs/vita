use crate::{Element, HasSites, SiteId};

/// Per-site chemical identity: the [`Element`] occupying each site.
///
/// Access is by keyed lookup: [`element`](HasElements::element) maps a [`SiteId`] to
/// its element. [`elements`](HasElements::elements) yields one element per site in
/// [`sites`](HasSites::sites) order; a site together with its element constitutes an
/// atom.
///
/// # Contract
///
/// [`element`](HasElements::element) is total over [`sites`](HasSites::sites): every site
/// has exactly one element.
/// [`elements`](HasElements::elements) yields values in the same order as
/// [`sites`](HasSites::sites).
pub trait HasElements: HasSites {
    /// Returns the element occupying `site`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn element(&self, site: SiteId) -> Element;

    /// Yields one element per site, in [`sites`](HasSites::sites) order.
    ///
    /// The default implementation looks up [`element`](HasElements::element) per site;
    /// override it when the elements can be produced directly.
    #[inline]
    fn elements(&self) -> impl Iterator<Item = Element> + '_ {
        self.sites().map(move |site| self.element(site))
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

        fn elements(&self) -> impl Iterator<Item = Element> + '_ {
            self.elements.iter().copied()
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
            vec![oxygen(), hydrogen(), hydrogen()],
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
        let sites = vec![site(1), site(2), site(3)];
        let elements = vec![oxygen(), hydrogen(), hydrogen()];
        let bare = Bare {
            sites: sites.clone(),
            elements: elements.clone(),
        };
        let columnar = Columnar { sites, elements };

        assert_eq!(
            bare.elements().collect::<Vec<_>>(),
            columnar.elements().collect::<Vec<_>>(),
        );
    }
}
