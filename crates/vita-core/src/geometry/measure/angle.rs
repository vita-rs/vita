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

#[cfg(test)]
mod tests {
    use super::*;

    use core::f64::consts::{FRAC_PI_2, PI};

    use crate::geometry::fixture::{close, configuration, s};
    use crate::units::angle::Degree;

    #[test]
    fn an_angle_between_arms_along_one_ray_is_zero() {
        let system = configuration(&[[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]);
        let opening: Angle<f64, Radian> = angle(&system, s(1), s(2), s(3)).unwrap();
        assert!(close(opening, Angle::new(0.0)));
    }

    #[test]
    fn an_angle_across_a_straight_chain_is_a_half_turn() {
        let system = configuration(&[[-1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]);
        let opening: Angle<f64, Radian> = angle(&system, s(1), s(2), s(3)).unwrap();
        assert!(close(opening, Angle::new(PI)));
    }

    #[test]
    fn perpendicular_arms_make_a_quarter_turn() {
        let system = configuration(&[[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        let opening: Angle<f64, Radian> = angle(&system, s(1), s(2), s(3)).unwrap();
        assert!(close(opening, Angle::new(FRAC_PI_2)));
    }

    #[test]
    fn an_angle_is_given_in_the_requested_unit() {
        let system = configuration(&[[1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        let opening: Angle<f64, Degree> = angle(&system, s(1), s(2), s(3)).unwrap();
        assert!(close(opening, Angle::new(90.0)));
    }

    #[test]
    fn an_arm_of_no_length_leaves_the_angle_absent() {
        let system = configuration(&[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]]);
        let collapsed: Option<Angle<f64, Radian>> = angle(&system, s(1), s(2), s(3));
        let swapped: Option<Angle<f64, Radian>> = angle(&system, s(3), s(2), s(1));
        assert!(collapsed.is_none());
        assert!(swapped.is_none());
    }

    #[test]
    fn an_angle_is_symmetric_in_its_arms() {
        let system = configuration(&[[1.0, 2.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 3.0]]);
        let forward: Angle<f64, Radian> = angle(&system, s(1), s(2), s(3)).unwrap();
        let reversed: Angle<f64, Radian> = angle(&system, s(3), s(2), s(1)).unwrap();
        assert_eq!(forward, reversed);
    }
}
