use super::point;
use crate::units::angle::{Angle, AngleUnit, Radian};
use crate::{HasPositions, Scalar, SiteId};

/// The dihedral angle of the chain `a`–`b`–`c`–`d`, in unit `U`.
///
/// The angle between the plane through `a`, `b`, `c` and the plane through `b`, `c`,
/// `d`, signed and within a half turn either way: sighting from `b` toward `c`, a
/// positive angle turns the far arm clockwise from the near one, so a syn chain reads
/// zero and an anti chain a half turn. Returns `None` when either triple is collinear
/// and leaves its plane undefined.
///
/// # Panics
///
/// Panics if `a`, `b`, `c`, or `d` is not in [`sites`](crate::HasSites::sites).
pub fn dihedral<S, V, U>(
    system: &S,
    a: SiteId,
    b: SiteId,
    c: SiteId,
    d: SiteId,
) -> Option<Angle<V, U>>
where
    S: HasPositions<V>,
    V: Scalar,
    U: AngleUnit,
{
    let (first, second, third, fourth) = (
        point(system, a),
        point(system, b),
        point(system, c),
        point(system, d),
    );
    let bridge = third - second;
    let axis = bridge.try_normalize()?;
    let leading = (second - first).cross(bridge).try_normalize()?;
    let trailing = bridge.cross(fourth - third).try_normalize()?;
    let sine = leading.cross(trailing).dot(axis);
    let cosine = leading.dot(trailing);
    Some(Angle::<V, Radian>::new(sine.atan2(cosine)).to::<U>())
}
