use vita_core::units::charge::{Charge, ChargeUnit};
use vita_core::{HasSites, Scalar, SiteId};

/// Per-site partial charge: the [`Charge`] of each site.
///
/// Access is by lookup:
/// [`partial_charge`](HasPartialCharges::partial_charge) maps a [`SiteId`]
/// to its charge, in any requested [unit](ChargeUnit).
/// [`partial_charges`](HasPartialCharges::partial_charges) iterates every
/// `(site, charge)` pair.
///
/// # Contract
///
/// [`partial_charge`](HasPartialCharges::partial_charge) is total over
/// [`sites`](HasSites::sites): every site has exactly one partial charge.
pub trait HasPartialCharges<V: Scalar>: HasSites {
    /// Returns the partial charge of `site`, in unit `U`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn partial_charge<U: ChargeUnit>(&self, site: SiteId) -> Charge<V, U>;

    /// Returns an iterator over every `(site, partial_charge)` pair, each
    /// charge in unit `U`.
    ///
    /// Each charge is yielded with its [`SiteId`]. The default
    /// implementation looks up
    /// [`partial_charge`](HasPartialCharges::partial_charge) per site;
    /// override it when the pairs can be produced directly.
    #[inline]
    fn partial_charges<U: ChargeUnit>(&self) -> impl Iterator<Item = (SiteId, Charge<V, U>)> + '_ {
        self.sites()
            .map(move |site| (site, self.partial_charge::<U>(site)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vita_core::units::charge::ElementaryCharge;

    fn site(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn e(value: f64) -> Charge<f64, ElementaryCharge> {
        Charge::new(value)
    }

    struct Bare {
        sites: Vec<SiteId>,
        partial_charges: Vec<Charge<f64, ElementaryCharge>>,
    }

    impl HasSites for Bare {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }

    impl HasPartialCharges<f64> for Bare {
        fn partial_charge<U: ChargeUnit>(&self, site: SiteId) -> Charge<f64, U> {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.partial_charges[i].to()
        }
    }

    struct Columnar {
        sites: Vec<SiteId>,
        partial_charges: Vec<Charge<f64, ElementaryCharge>>,
    }

    impl HasSites for Columnar {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }

    impl HasPartialCharges<f64> for Columnar {
        fn partial_charge<U: ChargeUnit>(&self, site: SiteId) -> Charge<f64, U> {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.partial_charges[i].to()
        }

        fn partial_charges<U: ChargeUnit>(
            &self,
        ) -> impl Iterator<Item = (SiteId, Charge<f64, U>)> + '_ {
            self.sites
                .iter()
                .copied()
                .zip(self.partial_charges.iter().copied().map(|q| q.to::<U>()))
        }
    }

    fn water() -> Bare {
        Bare {
            sites: vec![site(1), site(2), site(3)],
            partial_charges: vec![e(-0.8), e(0.4), e(0.4)],
        }
    }

    #[test]
    fn partial_charge() {
        let mol = water();
        assert_eq!(mol.partial_charge::<ElementaryCharge>(site(1)), e(-0.8));
        assert_eq!(mol.partial_charge::<ElementaryCharge>(site(2)), e(0.4));
        assert_eq!(mol.partial_charge::<ElementaryCharge>(site(3)), e(0.4));
    }

    #[test]
    fn partial_charges() {
        let mol = water();
        assert_eq!(
            mol.partial_charges::<ElementaryCharge>()
                .collect::<Vec<_>>(),
            vec![(site(1), e(-0.8)), (site(2), e(0.4)), (site(3), e(0.4))]
        );
    }

    #[test]
    fn partial_charges_empty() {
        let mol = Bare {
            sites: vec![],
            partial_charges: vec![],
        };
        assert_eq!(mol.partial_charges::<ElementaryCharge>().count(), 0);
    }

    #[test]
    fn override_matches_default() {
        use std::collections::BTreeMap;

        let sites = vec![site(1), site(2), site(3)];
        let charges = vec![e(-0.8), e(0.4), e(0.4)];

        let bare = Bare {
            sites: sites.clone(),
            partial_charges: charges.clone(),
        };
        let col = Columnar {
            sites,
            partial_charges: charges,
        };

        let bare_charges: BTreeMap<_, _> = bare.partial_charges::<ElementaryCharge>().collect();
        let col_charges: BTreeMap<_, _> = col.partial_charges::<ElementaryCharge>().collect();
        assert_eq!(bare_charges, col_charges);
    }
}
