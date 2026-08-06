//! What lies near what.
//!
//! [`arrangement`] lays the sites out by where they lie, and every question of nearness
//! is put to that one [`Arrangement`]: which sites lie [`near`](Arrangement::near) a
//! place, which [`pairs`](Arrangement::pairs) of them fall within the cutoff, and — as
//! [`Neighbors`] — that same relation indexed by site, to be asked about rather than
//! read through once. Each answer carries the displacement it was decided by, so a
//! reading that weighs by separation measures nothing twice. [`periodic_arrangement`]
//! and [`periodic_neighbors`] answer on the torus a lattice defines, taking every
//! separation to its shortest image.
//!
//! A place and a site are asked after differently: an [`Arrangement`] is indexed by the
//! space the sites sit in, so any point at all can be put to it, while [`Neighbors`] is
//! indexed by the sites themselves and holds no geometry, only the relation. Cells no
//! narrower than the cutoff carry both, so the work follows the sites and the pairs
//! they make rather than the square of their number.

mod arrangement;
mod neighbors;

pub use arrangement::{Arrangement, arrangement, periodic_arrangement};
pub use neighbors::{Neighbors, neighbors, periodic_neighbors};

#[cfg(test)]
mod fixture {
    use crate::geometry::fixture::{System, configuration, s};
    use crate::tensor::{Point3, Vector3};
    use crate::units::length::{Angstrom, Length, LengthUnit};
    use crate::{HasLattice, HasPositions, HasSites, Lattice, SiteId};

    pub struct Crystal {
        pub points: Vec<[f64; 3]>,
        lattice: Lattice<f64>,
    }

    impl HasSites for Crystal {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            (1..=self.points.len() as u32).map(s)
        }
    }

    impl HasPositions<f64> for Crystal {
        fn position<U: LengthUnit>(&self, site: SiteId) -> Point3<Length<f64, U>> {
            Point3::from_array(self.points[site.get() as usize - 1])
                .map(|value| Length::<f64, Angstrom>::new(value).to::<U>())
        }
    }

    impl HasLattice<f64> for Crystal {
        fn lattice(&self) -> Lattice<f64> {
            self.lattice
        }
    }

    pub struct Descending(pub System);

    impl HasSites for Descending {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            (1..=self.0.site_count() as u32).rev().map(s)
        }
    }

    impl HasPositions<f64> for Descending {
        fn position<U: LengthUnit>(&self, site: SiteId) -> Point3<Length<f64, U>> {
            self.0.position(site)
        }
    }

    pub fn reach(value: f64) -> Length<f64, Angstrom> {
        Length::new(value)
    }

    pub fn spot(point: [f64; 3]) -> Point3<Length<f64, Angstrom>> {
        Point3::from_array(point).map(Length::new)
    }

    pub fn chain() -> System {
        configuration(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
        ])
    }

    pub fn strewn() -> Vec<[f64; 3]> {
        (0..60)
            .map(|step| {
                let count = f64::from(step);
                [
                    (count * 2.7) % 10.0,
                    (count * 4.3) % 10.0,
                    (count * 6.1) % 10.0,
                ]
            })
            .collect()
    }

    pub fn cube(points: &[[f64; 3]]) -> Crystal {
        Crystal {
            points: points.to_vec(),
            lattice: Lattice::cubic(reach(10.0)).unwrap(),
        }
    }

    pub fn sheared(points: &[[f64; 3]]) -> Crystal {
        let edge = |x, y, z| Vector3::new(reach(x), reach(y), reach(z));
        Crystal {
            points: points.to_vec(),
            lattice: Lattice::from_vectors(
                edge(10.0, 0.0, 0.0),
                edge(0.0, 10.0, 0.0),
                edge(0.0, 9.0, 10.0),
            )
            .unwrap(),
        }
    }

    pub fn sweep(system: &Crystal, cutoff: Length<f64, Angstrom>) -> Vec<(SiteId, SiteId)> {
        let count = system.points.len() as u32;
        let lattice = system.lattice();
        (1..=count)
            .flat_map(|a| (a + 1..=count).map(move |b| (s(a), s(b))))
            .filter(|&(a, b)| {
                let separation = system.position::<Angstrom>(b) - system.position::<Angstrom>(a);
                lattice.minimum_image(separation).norm() <= cutoff
            })
            .collect()
    }
}
