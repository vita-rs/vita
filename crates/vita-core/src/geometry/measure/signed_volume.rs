use super::point;
use crate::units::volume::{CubicAngstrom, Volume, VolumeUnit};
use crate::{HasPositions, Scalar, SiteId};

/// The signed volume of the tetrahedron on `a`, `b`, `c`, and `d`, in unit `U`.
///
/// One sixth of the determinant of the three edges leaving `a`, positive when they form
/// a right-handed set. Its magnitude already follows from the six pairwise distances;
/// its sign does not, and is the whole of what separates four sites from their mirror
/// image. Coplanar sites span no volume.
///
/// # Panics
///
/// Panics if `a`, `b`, `c`, or `d` is not in [`sites`](crate::HasSites::sites).
pub fn signed_volume<S, V, U>(
    system: &S,
    a: SiteId,
    b: SiteId,
    c: SiteId,
    d: SiteId,
) -> Volume<V, U>
where
    S: HasPositions<V>,
    V: Scalar,
    U: VolumeUnit,
{
    let (first, second, third, fourth) = (
        point(system, a),
        point(system, b),
        point(system, c),
        point(system, d),
    );
    let determinant = (second - first).dot((third - first).cross(fourth - first));
    Volume::<V, CubicAngstrom>::new(determinant / V::from_f64(6.0)).to::<U>()
}
