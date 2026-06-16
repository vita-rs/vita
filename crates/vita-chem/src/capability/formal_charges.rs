use vita_core::{HasSites, SiteId};

/// Per-site Lewis-structure bookkeeping: the formal charge of each site.
///
/// Access is by lookup:
/// [`formal_charge`](HasFormalCharges::formal_charge) maps a [`SiteId`] to
/// its formal charge.
/// [`formal_charges`](HasFormalCharges::formal_charges) iterates every
/// `(site, formal_charge)` pair.
///
/// # Contract
///
/// [`formal_charge`](HasFormalCharges::formal_charge) is total over
/// [`sites`](HasSites::sites): every site has exactly one formal charge.
pub trait HasFormalCharges: HasSites {
    /// Returns the formal charge of `site`, in units of elementary charge.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn formal_charge(&self, site: SiteId) -> i8;

    /// Returns an iterator over every `(site, formal_charge)` pair.
    ///
    /// Each charge is yielded with its [`SiteId`]. The default implementation
    /// looks up [`formal_charge`](HasFormalCharges::formal_charge) per site;
    /// override it when the pairs can be produced directly.
    #[inline]
    fn formal_charges(&self) -> impl Iterator<Item = (SiteId, i8)> + '_ {
        self.sites()
            .map(move |site| (site, self.formal_charge(site)))
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
        formal_charges: Vec<i8>,
    }

    impl HasSites for Bare {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }

    impl HasFormalCharges for Bare {
        fn formal_charge(&self, site: SiteId) -> i8 {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.formal_charges[i]
        }
    }

    struct Columnar {
        sites: Vec<SiteId>,
        formal_charges: Vec<i8>,
    }

    impl HasSites for Columnar {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }

    impl HasFormalCharges for Columnar {
        fn formal_charge(&self, site: SiteId) -> i8 {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.formal_charges[i]
        }

        fn formal_charges(&self) -> impl Iterator<Item = (SiteId, i8)> + '_ {
            self.sites
                .iter()
                .copied()
                .zip(self.formal_charges.iter().copied())
        }
    }

    fn diazomethane() -> Bare {
        Bare {
            sites: vec![site(1), site(2), site(3)],
            formal_charges: vec![0, 1, -1],
        }
    }

    #[test]
    fn formal_charge() {
        let mol = diazomethane();
        assert_eq!(mol.formal_charge(site(1)), 0);
        assert_eq!(mol.formal_charge(site(2)), 1);
        assert_eq!(mol.formal_charge(site(3)), -1);
    }

    #[test]
    fn formal_charges() {
        let mol = diazomethane();
        assert_eq!(
            mol.formal_charges().collect::<Vec<_>>(),
            vec![(site(1), 0), (site(2), 1), (site(3), -1)]
        );
    }

    #[test]
    fn formal_charges_empty() {
        let mol = Bare {
            sites: vec![],
            formal_charges: vec![],
        };
        assert_eq!(mol.formal_charges().count(), 0);
    }

    #[test]
    fn override_matches_default() {
        use std::collections::BTreeMap;

        let sites = vec![site(1), site(2), site(3)];
        let charges = vec![0i8, 1, -1];

        let bare = Bare {
            sites: sites.clone(),
            formal_charges: charges.clone(),
        };
        let col = Columnar {
            sites,
            formal_charges: charges,
        };

        let bare_charges: BTreeMap<_, _> = bare.formal_charges().collect();
        let col_charges: BTreeMap<_, _> = col.formal_charges().collect();
        assert_eq!(bare_charges, col_charges);
    }
}
