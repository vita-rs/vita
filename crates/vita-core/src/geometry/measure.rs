//! Readings taken from a fixed tuple of sites.
//!
//! [`distance`], [`angle`], and [`dihedral`] fix a configuration up to a rigid motion:
//! one, two, and three of them place each further site, so no further invariant of five
//! or more sites is independent of them. [`signed_volume`] reads the handedness of four
//! sites directly, with no chain to sight along, and [`displacement`] carries a
//! separation with its direction. An angle or a dihedral whose arms collapse is absent
//! rather than zero.

mod distance;

pub use distance::distance;

use crate::tensor::Point3;
use crate::units::length::Angstrom;
use crate::{HasPositions, Quantity, Scalar, SiteId};

#[cfg(test)]
mod fixture {
    use super::*;

    use crate::HasSites;
    use crate::units::length::{Length, LengthUnit};

    pub struct System(Vec<Point3<f64>>);

    impl HasSites for System {
        fn sites(&self) -> impl Iterator<Item = SiteId> + '_ {
            (1..=self.0.len() as u32).map(s)
        }
    }

    impl HasPositions<f64> for System {
        fn position<U: LengthUnit>(&self, site: SiteId) -> Point3<Length<f64, U>> {
            self.0[site.get() as usize - 1]
                .map(|value| Length::<f64, Angstrom>::new(value).to::<U>())
        }
    }

    pub fn s(n: u32) -> SiteId {
        SiteId::new(n).unwrap()
    }

    pub fn configuration(points: &[[f64; 3]]) -> System {
        System(points.iter().copied().map(Point3::from_array).collect())
    }

    pub fn close<Q: Quantity<Value = f64>>(a: Q, b: Q) -> bool {
        (a.value() - b.value()).abs() <= 1e-12
    }
}
