use vita_core::{HasSites, SiteId};

/// Per-site Lewis-structure bookkeeping: the formal charge of each site.
///
/// Access is by keyed lookup:
/// [`formal_charge`](HasFormalCharges::formal_charge) maps a [`SiteId`] to
/// its formal charge.
/// [`formal_charges`](HasFormalCharges::formal_charges) yields one charge per
/// site in [`sites`](HasSites::sites) order.
///
/// # Contract
///
/// [`formal_charge`](HasFormalCharges::formal_charge) is total over
/// [`sites`](HasSites::sites): every site has exactly one formal charge.
/// [`formal_charges`](HasFormalCharges::formal_charges) yields values in the
/// same order as [`sites`](HasSites::sites).
pub trait HasFormalCharges: HasSites {
    /// Returns the formal charge of `site`, in units of elementary charge.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn formal_charge(&self, site: SiteId) -> i8;

    /// Yields one formal charge per site, in [`sites`](HasSites::sites) order.
    ///
    /// The default implementation looks up
    /// [`formal_charge`](HasFormalCharges::formal_charge) per site; override
    /// it when the charges can be produced directly.
    #[inline]
    fn formal_charges(&self) -> impl Iterator<Item = i8> + '_ {
        self.sites().map(move |site| self.formal_charge(site))
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

        fn formal_charges(&self) -> impl Iterator<Item = i8> + '_ {
            self.formal_charges.iter().copied()
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
        assert_eq!(mol.formal_charges().collect::<Vec<_>>(), vec![0, 1, -1]);
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

        assert_eq!(
            bare.formal_charges().collect::<Vec<_>>(),
            col.formal_charges().collect::<Vec<_>>(),
        );
    }
}
