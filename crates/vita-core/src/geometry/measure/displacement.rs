use crate::tensor::Vector3;
use crate::units::length::{Length, LengthUnit};
use crate::{HasPositions, Scalar, SiteId};

/// The displacement from `from` to `to`, in unit `U`.
///
/// The separation of two sites as a vector, direction and all; its norm is their
/// [`distance`](super::distance()). Under periodic boundaries the separation to take is
/// its shortest image, [`Lattice::minimum_image`](crate::Lattice::minimum_image).
///
/// # Panics
///
/// Panics if `from` or `to` is not in [`sites`](crate::HasSites::sites).
pub fn displacement<S, V, U>(system: &S, from: SiteId, to: SiteId) -> Vector3<Length<V, U>>
where
    S: HasPositions<V>,
    V: Scalar,
    U: LengthUnit,
{
    system.position::<U>(to) - system.position::<U>(from)
}
