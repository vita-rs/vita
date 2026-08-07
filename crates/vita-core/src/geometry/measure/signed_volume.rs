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

#[cfg(test)]
mod tests {
    use super::*;

    use crate::geometry::fixture::{System, close, configuration, s};
    use crate::units::volume::CubicNanometer;

    fn corner() -> System {
        configuration(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ])
    }

    #[test]
    fn coplanar_sites_span_no_volume() {
        let system = configuration(&[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ]);
        let span: Volume<f64, CubicAngstrom> = signed_volume(&system, s(1), s(2), s(3), s(4));
        assert_eq!(span, Volume::new(0.0));
    }

    #[test]
    fn a_unit_corner_spans_a_sixth_of_the_unit_cube() {
        let span: Volume<f64, CubicAngstrom> = signed_volume(&corner(), s(1), s(2), s(3), s(4));
        assert!(close(span, Volume::new(1.0 / 6.0)));
    }

    #[test]
    fn a_signed_volume_is_given_in_the_requested_unit() {
        let span: Volume<f64, CubicNanometer> = signed_volume(&corner(), s(1), s(2), s(3), s(4));
        assert!(close(span, Volume::new(1.0 / 6000.0)));
    }

    #[test]
    fn exchanging_two_sites_reverses_the_sign() {
        let system = corner();
        let span: Volume<f64, CubicAngstrom> = signed_volume(&system, s(1), s(2), s(3), s(4));
        let swapped: Volume<f64, CubicAngstrom> = signed_volume(&system, s(1), s(2), s(4), s(3));
        assert_eq!(swapped, -span);
    }

    #[test]
    fn a_cyclic_exchange_of_the_last_three_sites_preserves_the_volume() {
        let system = corner();
        let span: Volume<f64, CubicAngstrom> = signed_volume(&system, s(1), s(2), s(3), s(4));
        let cycled: Volume<f64, CubicAngstrom> = signed_volume(&system, s(1), s(3), s(4), s(2));
        assert_eq!(cycled, span);
    }

    #[test]
    fn a_signed_volume_is_unchanged_by_translation() {
        let moved = configuration(&[
            [5.0, -3.0, 2.0],
            [6.0, -3.0, 2.0],
            [5.0, -2.0, 2.0],
            [5.0, -3.0, 3.0],
        ]);
        let span: Volume<f64, CubicAngstrom> = signed_volume(&corner(), s(1), s(2), s(3), s(4));
        let translated: Volume<f64, CubicAngstrom> = signed_volume(&moved, s(1), s(2), s(3), s(4));
        assert_eq!(translated, span);
    }
}
