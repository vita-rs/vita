use crate::units::mass::{Mass, MassUnit};
use crate::{HasSites, Scalar, SiteId};

/// Per-site mass: the [`Mass`] of the particle at each site.
///
/// Access is by keyed lookup: [`mass`](HasMasses::mass) maps a [`SiteId`] to its mass,
/// in any requested [unit](MassUnit). [`masses`](HasMasses::masses) yields one mass per
/// site in [`sites`](HasSites::sites) order.
///
/// # Contract
///
/// [`mass`](HasMasses::mass) is total over [`sites`](HasSites::sites): every site has
/// exactly one mass.
/// [`masses`](HasMasses::masses) yields values in the same order as
/// [`sites`](HasSites::sites).
pub trait HasMasses<V: Scalar>: HasSites {
    /// Returns the mass of `site`, in unit `U`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn mass<U: MassUnit>(&self, site: SiteId) -> Mass<V, U>;

    /// Yields one mass per site, in [`sites`](HasSites::sites) order.
    ///
    /// The default implementation looks up [`mass`](HasMasses::mass) per site; override
    /// it when the masses can be produced directly.
    #[inline]
    fn masses<U: MassUnit>(&self) -> impl Iterator<Item = Mass<V, U>> + '_ {
        self.sites().map(move |site| self.mass::<U>(site))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::mass::Dalton;

    fn site(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn dalton(value: f64) -> Mass<f64, Dalton> {
        Mass::new(value)
    }

    struct Bare {
        sites: Vec<SiteId>,
        masses: Vec<Mass<f64, Dalton>>,
    }
    impl HasSites for Bare {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }
    impl HasMasses<f64> for Bare {
        fn mass<U: MassUnit>(&self, site: SiteId) -> Mass<f64, U> {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.masses[i].to()
        }
    }

    struct Columnar {
        sites: Vec<SiteId>,
        masses: Vec<Mass<f64, Dalton>>,
    }
    impl HasSites for Columnar {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }
    impl HasMasses<f64> for Columnar {
        fn mass<U: MassUnit>(&self, site: SiteId) -> Mass<f64, U> {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.masses[i].to()
        }

        fn masses<U: MassUnit>(&self) -> impl Iterator<Item = Mass<f64, U>> + '_ {
            self.masses.iter().copied().map(|m| m.to::<U>())
        }
    }

    fn water() -> Bare {
        Bare {
            sites: vec![site(1), site(2), site(3)],
            masses: vec![dalton(15.999), dalton(1.008), dalton(1.008)],
        }
    }

    #[test]
    fn mass() {
        let mol = water();
        assert_eq!(mol.mass::<Dalton>(site(1)), dalton(15.999));
        assert_eq!(mol.mass::<Dalton>(site(2)), dalton(1.008));
    }

    #[test]
    fn masses() {
        let mol = water();
        assert_eq!(
            mol.masses::<Dalton>().collect::<Vec<_>>(),
            vec![dalton(15.999), dalton(1.008), dalton(1.008)],
        );
    }

    #[test]
    fn masses_empty() {
        let mol = Bare {
            sites: vec![],
            masses: vec![],
        };
        assert_eq!(mol.masses::<Dalton>().count(), 0);
    }

    #[test]
    fn override_matches_default() {
        let sites = vec![site(1), site(2), site(3)];
        let masses = vec![dalton(15.999), dalton(1.008), dalton(1.008)];
        let bare = Bare {
            sites: sites.clone(),
            masses: masses.clone(),
        };
        let columnar = Columnar { sites, masses };

        assert_eq!(
            bare.masses::<Dalton>().collect::<Vec<_>>(),
            columnar.masses::<Dalton>().collect::<Vec<_>>(),
        );
    }
}
