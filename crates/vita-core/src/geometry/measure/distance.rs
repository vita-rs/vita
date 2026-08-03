use crate::units::length::{Length, LengthUnit};
use crate::{HasPositions, Scalar, SiteId};

/// The distance between `a` and `b`, in unit `U`.
///
/// The norm of their [`displacement`](super::displacement()), taken so that it stays
/// representable wherever the coordinates carrying it are.
///
/// # Panics
///
/// Panics if `a` or `b` is not in [`sites`](crate::HasSites::sites).
pub fn distance<S, V, U>(system: &S, a: SiteId, b: SiteId) -> Length<V, U>
where
    S: HasPositions<V>,
    V: Scalar,
    U: LengthUnit,
{
    system.position::<U>(a).distance(system.position::<U>(b))
}
