use crate::tensor::Vector3;
use crate::units::acceleration::{Acceleration, AccelerationUnit};
use crate::{HasSites, Scalar, SiteId};

/// Per-site acceleration: the [`Vector3`] acceleration of each site.
///
/// Access is by keyed lookup: [`acceleration`](HasAccelerations::acceleration) maps a
/// [`SiteId`] to its acceleration, in any requested [unit](AccelerationUnit).
/// [`accelerations`](HasAccelerations::accelerations) yields one acceleration per site
/// in [`sites`](HasSites::sites) order.
///
/// # Contract
///
/// [`acceleration`](HasAccelerations::acceleration) is total over [`sites`](HasSites::sites):
/// every site has exactly one acceleration.
/// [`accelerations`](HasAccelerations::accelerations) yields values in the same order as
/// [`sites`](HasSites::sites).
pub trait HasAccelerations<V: Scalar>: HasSites {
    /// Returns the acceleration of `site`, in unit `U`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn acceleration<U: AccelerationUnit>(&self, site: SiteId) -> Vector3<Acceleration<V, U>>;

    /// Yields one acceleration per site, in [`sites`](HasSites::sites) order.
    ///
    /// The default implementation looks up
    /// [`acceleration`](HasAccelerations::acceleration) per site; override it when the
    /// accelerations can be produced directly.
    #[inline]
    fn accelerations<U: AccelerationUnit>(
        &self,
    ) -> impl Iterator<Item = Vector3<Acceleration<V, U>>> + '_ {
        self.sites().map(move |site| self.acceleration::<U>(site))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::acceleration::AngstromPerSquarePicosecond;

    fn site(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn angstrom_per_square_picosecond(
        x: f64,
        y: f64,
        z: f64,
    ) -> Vector3<Acceleration<f64, AngstromPerSquarePicosecond>> {
        Vector3::new(
            Acceleration::new(x),
            Acceleration::new(y),
            Acceleration::new(z),
        )
    }

    struct Bare {
        sites: Vec<SiteId>,
        accelerations: Vec<Vector3<Acceleration<f64, AngstromPerSquarePicosecond>>>,
    }
    impl HasSites for Bare {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }
    impl HasAccelerations<f64> for Bare {
        fn acceleration<U: AccelerationUnit>(&self, site: SiteId) -> Vector3<Acceleration<f64, U>> {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.accelerations[i].map(|a| a.to())
        }
    }

    struct Columnar {
        sites: Vec<SiteId>,
        accelerations: Vec<Vector3<Acceleration<f64, AngstromPerSquarePicosecond>>>,
    }
    impl HasSites for Columnar {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }
    impl HasAccelerations<f64> for Columnar {
        fn acceleration<U: AccelerationUnit>(&self, site: SiteId) -> Vector3<Acceleration<f64, U>> {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.accelerations[i].map(|a| a.to())
        }

        fn accelerations<U: AccelerationUnit>(
            &self,
        ) -> impl Iterator<Item = Vector3<Acceleration<f64, U>>> + '_ {
            self.accelerations
                .iter()
                .copied()
                .map(|a| a.map(|c| c.to::<U>()))
        }
    }

    fn field() -> Bare {
        Bare {
            sites: vec![site(1), site(2), site(3)],
            accelerations: vec![
                angstrom_per_square_picosecond(0.0, 0.0, 0.0),
                angstrom_per_square_picosecond(2.0, -1.0, 0.0),
                angstrom_per_square_picosecond(-2.0, -1.0, 0.0),
            ],
        }
    }

    #[test]
    fn acceleration() {
        let sys = field();
        assert_eq!(
            sys.acceleration::<AngstromPerSquarePicosecond>(site(1)),
            angstrom_per_square_picosecond(0.0, 0.0, 0.0)
        );
        assert_eq!(
            sys.acceleration::<AngstromPerSquarePicosecond>(site(2)),
            angstrom_per_square_picosecond(2.0, -1.0, 0.0)
        );
    }

    #[test]
    fn accelerations() {
        let sys = field();
        assert_eq!(
            sys.accelerations::<AngstromPerSquarePicosecond>()
                .collect::<Vec<_>>(),
            vec![
                angstrom_per_square_picosecond(0.0, 0.0, 0.0),
                angstrom_per_square_picosecond(2.0, -1.0, 0.0),
                angstrom_per_square_picosecond(-2.0, -1.0, 0.0),
            ],
        );
    }

    #[test]
    fn accelerations_empty() {
        let sys = Bare {
            sites: vec![],
            accelerations: vec![],
        };
        assert_eq!(
            sys.accelerations::<AngstromPerSquarePicosecond>().count(),
            0
        );
    }

    #[test]
    fn override_matches_default() {
        let sites = vec![site(1), site(2), site(3)];
        let accelerations = vec![
            angstrom_per_square_picosecond(0.0, 0.0, 0.0),
            angstrom_per_square_picosecond(2.0, -1.0, 0.0),
            angstrom_per_square_picosecond(-2.0, -1.0, 0.0),
        ];
        let bare = Bare {
            sites: sites.clone(),
            accelerations: accelerations.clone(),
        };
        let columnar = Columnar {
            sites,
            accelerations,
        };

        assert_eq!(
            bare.accelerations::<AngstromPerSquarePicosecond>()
                .collect::<Vec<_>>(),
            columnar
                .accelerations::<AngstromPerSquarePicosecond>()
                .collect::<Vec<_>>(),
        );
    }
}
