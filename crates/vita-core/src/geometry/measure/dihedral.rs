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

#[cfg(test)]
mod tests {
    use super::*;

    use core::f64::consts::{FRAC_PI_2, PI};

    use crate::geometry::measure::fixture::{System, close, configuration, s};
    use crate::units::angle::Degree;

    fn chain(far: [f64; 3]) -> System {
        configuration(&[[0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], far])
    }

    #[test]
    fn a_syn_planar_chain_has_a_dihedral_of_zero() {
        let system = chain([1.0, 1.0, 0.0]);
        let torsion: Angle<f64, Radian> = dihedral(&system, s(1), s(2), s(3), s(4)).unwrap();
        assert!(close(torsion, Angle::new(0.0)));
    }

    #[test]
    fn an_anti_planar_chain_has_a_dihedral_of_a_half_turn() {
        let system = chain([1.0, -1.0, 0.0]);
        let torsion: Angle<f64, Radian> = dihedral(&system, s(1), s(2), s(3), s(4)).unwrap();
        assert!(close(torsion, Angle::new(PI)));
    }

    #[test]
    fn a_quarter_turn_of_the_far_arm_measures_a_quarter_turn() {
        let system = chain([1.0, 0.0, 1.0]);
        let torsion: Angle<f64, Radian> = dihedral(&system, s(1), s(2), s(3), s(4)).unwrap();
        assert!(close(torsion, Angle::new(FRAC_PI_2)));
    }

    #[test]
    fn a_dihedral_reverses_sign_when_the_far_arm_turns_the_other_way() {
        let system = chain([1.0, 0.0, -1.0]);
        let torsion: Angle<f64, Radian> = dihedral(&system, s(1), s(2), s(3), s(4)).unwrap();
        assert!(close(torsion, Angle::new(-FRAC_PI_2)));
    }

    #[test]
    fn a_dihedral_is_given_in_the_requested_unit() {
        let system = chain([1.0, 0.0, 1.0]);
        let torsion: Angle<f64, Degree> = dihedral(&system, s(1), s(2), s(3), s(4)).unwrap();
        assert!(close(torsion, Angle::new(90.0)));
    }

    #[test]
    fn a_collinear_leading_triple_leaves_the_dihedral_absent() {
        let system = configuration(&[
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
        ]);
        let torsion: Option<Angle<f64, Radian>> = dihedral(&system, s(1), s(2), s(3), s(4));
        assert!(torsion.is_none());
    }

    #[test]
    fn a_collinear_trailing_triple_leaves_the_dihedral_absent() {
        let system = chain([2.0, 0.0, 0.0]);
        let torsion: Option<Angle<f64, Radian>> = dihedral(&system, s(1), s(2), s(3), s(4));
        assert!(torsion.is_none());
    }

    #[test]
    fn a_dihedral_is_unchanged_when_its_chain_is_reversed() {
        let system = chain([1.0, 0.0, 1.0]);
        let forward: Angle<f64, Radian> = dihedral(&system, s(1), s(2), s(3), s(4)).unwrap();
        let reversed: Angle<f64, Radian> = dihedral(&system, s(4), s(3), s(2), s(1)).unwrap();
        assert!(close(forward, reversed));
    }

    #[test]
    fn a_dihedral_is_unchanged_when_its_chain_is_scaled() {
        let system = chain([1.0, 1.0, 1.0]);
        let scaled = configuration(&[
            [0.0, 3.0, 0.0],
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [3.0, 3.0, 3.0],
        ]);
        let torsion: Angle<f64, Radian> = dihedral(&system, s(1), s(2), s(3), s(4)).unwrap();
        let enlarged: Angle<f64, Radian> = dihedral(&scaled, s(1), s(2), s(3), s(4)).unwrap();
        assert!(close(torsion, enlarged));
    }
}
