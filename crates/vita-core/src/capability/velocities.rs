use crate::tensor::Vector3;
use crate::units::velocity::{Velocity, VelocityUnit};
use crate::{HasSites, Scalar, SiteId};

/// Per-site velocity: the [`Vector3`] velocity of each site.
///
/// Access is by lookup: [`velocity`](HasVelocities::velocity) maps a [`SiteId`] to its
/// velocity, in any requested [unit](VelocityUnit). [`velocities`](HasVelocities::velocities)
/// iterates every `(site, velocity)` pair.
///
/// # Contract
///
/// [`velocity`](HasVelocities::velocity) is total over [`sites`](HasSites::sites): every
/// site has exactly one velocity.
pub trait HasVelocities<V: Scalar>: HasSites {
    /// Returns the velocity of `site`, in unit `U`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn velocity<U: VelocityUnit>(&self, site: SiteId) -> Vector3<Velocity<V, U>>;

    /// Returns an iterator over every `(site, velocity)` pair, each velocity in unit `U`.
    ///
    /// Each velocity is yielded with its [`SiteId`]. The default implementation looks up
    /// [`velocity`](HasVelocities::velocity) per site; override it when the pairs can be
    /// produced directly.
    #[inline]
    fn velocities<U: VelocityUnit>(
        &self,
    ) -> impl Iterator<Item = (SiteId, Vector3<Velocity<V, U>>)> + '_ {
        self.sites()
            .map(move |site| (site, self.velocity::<U>(site)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::velocity::AngstromPerPicosecond;

    fn site(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn angstrom_per_picosecond(
        x: f64,
        y: f64,
        z: f64,
    ) -> Vector3<Velocity<f64, AngstromPerPicosecond>> {
        Vector3::new(Velocity::new(x), Velocity::new(y), Velocity::new(z))
    }

    struct Bare {
        sites: Vec<SiteId>,
        velocities: Vec<Vector3<Velocity<f64, AngstromPerPicosecond>>>,
    }
    impl HasSites for Bare {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }
    impl HasVelocities<f64> for Bare {
        fn velocity<U: VelocityUnit>(&self, site: SiteId) -> Vector3<Velocity<f64, U>> {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.velocities[i].map(|v| v.to())
        }
    }

    struct Columnar {
        sites: Vec<SiteId>,
        velocities: Vec<Vector3<Velocity<f64, AngstromPerPicosecond>>>,
    }
    impl HasSites for Columnar {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }
    impl HasVelocities<f64> for Columnar {
        fn velocity<U: VelocityUnit>(&self, site: SiteId) -> Vector3<Velocity<f64, U>> {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.velocities[i].map(|v| v.to())
        }

        fn velocities<U: VelocityUnit>(
            &self,
        ) -> impl Iterator<Item = (SiteId, Vector3<Velocity<f64, U>>)> + '_ {
            self.sites.iter().copied().zip(
                self.velocities
                    .iter()
                    .copied()
                    .map(|v| v.map(|c| c.to::<U>())),
            )
        }
    }

    fn flow() -> Bare {
        Bare {
            sites: vec![site(1), site(2), site(3)],
            velocities: vec![
                angstrom_per_picosecond(0.0, 0.0, 0.0),
                angstrom_per_picosecond(1.5, -0.5, 0.0),
                angstrom_per_picosecond(-1.5, -0.5, 0.0),
            ],
        }
    }

    #[test]
    fn velocity() {
        let sys = flow();
        assert_eq!(
            sys.velocity::<AngstromPerPicosecond>(site(1)),
            angstrom_per_picosecond(0.0, 0.0, 0.0)
        );
        assert_eq!(
            sys.velocity::<AngstromPerPicosecond>(site(2)),
            angstrom_per_picosecond(1.5, -0.5, 0.0)
        );
    }

    #[test]
    fn velocities() {
        let sys = flow();
        assert_eq!(
            sys.velocities::<AngstromPerPicosecond>()
                .collect::<Vec<_>>(),
            vec![
                (site(1), angstrom_per_picosecond(0.0, 0.0, 0.0)),
                (site(2), angstrom_per_picosecond(1.5, -0.5, 0.0)),
                (site(3), angstrom_per_picosecond(-1.5, -0.5, 0.0))
            ]
        );
    }

    #[test]
    fn velocities_empty() {
        let sys = Bare {
            sites: vec![],
            velocities: vec![],
        };
        assert_eq!(sys.velocities::<AngstromPerPicosecond>().count(), 0);
    }

    #[test]
    fn override_matches_default() {
        use std::collections::BTreeMap;

        let sites = vec![site(1), site(2), site(3)];
        let velocities = vec![
            angstrom_per_picosecond(0.0, 0.0, 0.0),
            angstrom_per_picosecond(1.5, -0.5, 0.0),
            angstrom_per_picosecond(-1.5, -0.5, 0.0),
        ];
        let bare = Bare {
            sites: sites.clone(),
            velocities: velocities.clone(),
        };
        let columnar = Columnar { sites, velocities };

        let bare_velocities: BTreeMap<_, _> = bare.velocities::<AngstromPerPicosecond>().collect();
        let columnar_velocities: BTreeMap<_, _> =
            columnar.velocities::<AngstromPerPicosecond>().collect();
        assert_eq!(bare_velocities, columnar_velocities);
    }
}
