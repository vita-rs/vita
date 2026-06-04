use crate::tensor::Point3;
use crate::units::length::{Length, LengthUnit};
use crate::{HasSites, Scalar, SiteId};

/// Per-site position: the [`Point3`] locating each site in space.
///
/// Access is by lookup: [`position`](HasPositions::position) maps a [`SiteId`] to its
/// position, in any requested [unit](LengthUnit). [`positions`](HasPositions::positions)
/// iterates every `(site, position)` pair.
///
/// # Contract
///
/// [`position`](HasPositions::position) is total over [`sites`](HasSites::sites): every
/// site has exactly one position.
pub trait HasPositions<V: Scalar>: HasSites {
    /// Returns the position of `site`, in unit `U`.
    ///
    /// # Panics
    ///
    /// Panics if `site` is not in [`sites`](HasSites::sites).
    fn position<U: LengthUnit>(&self, site: SiteId) -> Point3<Length<V, U>>;

    /// Returns an iterator over every `(site, position)` pair, each position in unit `U`.
    ///
    /// Each position is yielded with its [`SiteId`]. The default implementation looks up
    /// [`position`](HasPositions::position) per site; override it when the pairs can be
    /// produced directly.
    #[inline]
    fn positions<U: LengthUnit>(
        &self,
    ) -> impl Iterator<Item = (SiteId, Point3<Length<V, U>>)> + '_ {
        self.sites()
            .map(move |site| (site, self.position::<U>(site)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::units::length::Angstrom;

    fn site(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    fn angstrom(x: f64, y: f64, z: f64) -> Point3<Length<f64, Angstrom>> {
        Point3::new(Length::new(x), Length::new(y), Length::new(z))
    }

    struct Bare {
        sites: Vec<SiteId>,
        positions: Vec<Point3<Length<f64, Angstrom>>>,
    }
    impl HasSites for Bare {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }
    impl HasPositions<f64> for Bare {
        fn position<U: LengthUnit>(&self, site: SiteId) -> Point3<Length<f64, U>> {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.positions[i].map(|l| l.to())
        }
    }

    struct Columnar {
        sites: Vec<SiteId>,
        positions: Vec<Point3<Length<f64, Angstrom>>>,
    }
    impl HasSites for Columnar {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            self.sites.iter().copied()
        }
    }
    impl HasPositions<f64> for Columnar {
        fn position<U: LengthUnit>(&self, site: SiteId) -> Point3<Length<f64, U>> {
            let i = self.sites.iter().position(|&s| s == site).unwrap();
            self.positions[i].map(|l| l.to())
        }

        fn positions<U: LengthUnit>(
            &self,
        ) -> impl Iterator<Item = (SiteId, Point3<Length<f64, U>>)> + '_ {
            self.sites.iter().copied().zip(
                self.positions
                    .iter()
                    .copied()
                    .map(|p| p.map(|l| l.to::<U>())),
            )
        }
    }

    fn water() -> Bare {
        Bare {
            sites: vec![site(1), site(2), site(3)],
            positions: vec![
                angstrom(0.0, 0.0, 0.0),
                angstrom(0.757, 0.586, 0.0),
                angstrom(-0.757, 0.586, 0.0),
            ],
        }
    }

    #[test]
    fn position() {
        let mol = water();
        assert_eq!(mol.position::<Angstrom>(site(1)), angstrom(0.0, 0.0, 0.0));
        assert_eq!(
            mol.position::<Angstrom>(site(2)),
            angstrom(0.757, 0.586, 0.0)
        );
    }

    #[test]
    fn positions() {
        let mol = water();
        assert_eq!(
            mol.positions::<Angstrom>().collect::<Vec<_>>(),
            vec![
                (site(1), angstrom(0.0, 0.0, 0.0)),
                (site(2), angstrom(0.757, 0.586, 0.0)),
                (site(3), angstrom(-0.757, 0.586, 0.0))
            ]
        );
    }

    #[test]
    fn positions_empty() {
        let mol = Bare {
            sites: vec![],
            positions: vec![],
        };
        assert_eq!(mol.positions::<Angstrom>().count(), 0);
    }

    #[test]
    fn override_matches_default() {
        use std::collections::BTreeMap;

        let sites = vec![site(1), site(2), site(3)];
        let positions = vec![
            angstrom(0.0, 0.0, 0.0),
            angstrom(0.757, 0.586, 0.0),
            angstrom(-0.757, 0.586, 0.0),
        ];
        let bare = Bare {
            sites: sites.clone(),
            positions: positions.clone(),
        };
        let columnar = Columnar { sites, positions };

        let bare_positions: BTreeMap<_, _> = bare.positions::<Angstrom>().collect();
        let columnar_positions: BTreeMap<_, _> = columnar.positions::<Angstrom>().collect();
        assert_eq!(bare_positions, columnar_positions);
    }
}
