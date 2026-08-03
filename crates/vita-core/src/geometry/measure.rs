//! Readings taken from a fixed tuple of sites.
//!
//! [`distance`], [`angle`], and [`dihedral`] fix a configuration up to a rigid motion:
//! one, two, and three of them place each further site, so no further invariant of five
//! or more sites is independent of them. [`signed_volume`] reads the handedness of four
//! sites directly, with no chain to sight along, and [`displacement`] carries a
//! separation with its direction. An angle or a dihedral whose arms collapse is absent
//! rather than zero.

mod angle;
mod dihedral;
mod displacement;
mod distance;
mod signed_volume;

pub use angle::angle;
pub use dihedral::dihedral;
pub use displacement::displacement;
pub use distance::distance;
pub use signed_volume::signed_volume;

use crate::tensor::Point3;
use crate::units::length::Angstrom;
use crate::{HasPositions, Quantity, Scalar, SiteId};

/// The position of `site` in ångströms, stripped of its unit: the readings below
/// multiply coordinates together, and a product of lengths is not a length.
fn point<S: HasPositions<V>, V: Scalar>(system: &S, site: SiteId) -> Point3<V> {
    system.position::<Angstrom>(site).map(Quantity::value)
}
