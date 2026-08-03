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

#[cfg(test)]
mod tests {
    use super::*;

    use crate::geometry::measure::fixture::{close, configuration, s};
    use crate::units::length::{Angstrom, Nanometer};

    #[test]
    fn a_displacement_from_a_site_to_itself_is_zero() {
        let system = configuration(&[[1.0, 2.0, 3.0]]);
        let separation: Vector3<Length<f64, Angstrom>> = displacement(&system, s(1), s(1));
        assert_eq!(separation, Vector3::ZERO);
    }

    #[test]
    fn a_displacement_runs_from_the_first_site_to_the_second() {
        let system = configuration(&[[1.0, 0.0, 0.0], [4.0, 2.0, -1.0]]);
        let separation: Vector3<Length<f64, Angstrom>> = displacement(&system, s(1), s(2));
        let expected = Vector3::new(Length::new(3.0), Length::new(2.0), Length::new(-1.0));
        assert_eq!(separation, expected);
    }

    #[test]
    fn a_displacement_is_given_in_the_requested_unit() {
        let system = configuration(&[[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]]);
        let separation: Vector3<Length<f64, Nanometer>> = displacement(&system, s(1), s(2));
        assert!(close(separation.x, Length::new(1.0)));
    }

    #[test]
    fn a_displacement_negates_when_its_endpoints_are_exchanged() {
        let system = configuration(&[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]);
        let forward: Vector3<Length<f64, Angstrom>> = displacement(&system, s(1), s(2));
        let backward: Vector3<Length<f64, Angstrom>> = displacement(&system, s(2), s(1));
        assert_eq!(forward, -backward);
    }
}
