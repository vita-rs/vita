use super::point;
use crate::units::angle::{Angle, AngleUnit, Radian};
use crate::{HasPositions, Scalar, SiteId};

/// The angle at `vertex` between the arms reaching `a` and `b`, in unit `U`.
///
/// Unsigned and within a half turn: an angle at a vertex carries no sense of its own.
/// Returns `None` when either arm has no length — a site placed on the vertex leaves the
/// angle undefined rather than zero.
///
/// # Panics
///
/// Panics if `a`, `vertex`, or `b` is not in [`sites`](crate::HasSites::sites).
pub fn angle<S, V, U>(system: &S, a: SiteId, vertex: SiteId, b: SiteId) -> Option<Angle<V, U>>
where
    S: HasPositions<V>,
    V: Scalar,
    U: AngleUnit,
{
    let center = point(system, vertex);
    let first = (point(system, a) - center).try_normalize()?;
    let second = (point(system, b) - center).try_normalize()?;
    Some(Angle::<V, Radian>::new(first.angle_between(second)).to::<U>())
}
