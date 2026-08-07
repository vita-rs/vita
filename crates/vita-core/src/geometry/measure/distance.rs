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

#[cfg(test)]
mod tests {
    use super::*;

    use crate::geometry::fixture::{close, configuration, s};
    use crate::units::length::{Angstrom, Nanometer};

    #[test]
    fn the_distance_between_coincident_sites_is_zero() {
        let system = configuration(&[[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]]);
        let separation: Length<f64, Angstrom> = distance(&system, s(1), s(2));
        assert_eq!(separation, Length::new(0.0));
    }

    #[test]
    fn a_distance_is_the_norm_of_the_separation() {
        let system = configuration(&[[0.0, 0.0, 0.0], [3.0, 4.0, 0.0]]);
        let separation: Length<f64, Angstrom> = distance(&system, s(1), s(2));
        assert_eq!(separation, Length::new(5.0));
    }

    #[test]
    fn a_distance_is_given_in_the_requested_unit() {
        let system = configuration(&[[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]]);
        let separation: Length<f64, Nanometer> = distance(&system, s(1), s(2));
        assert!(close(separation, Length::new(1.0)));
    }

    #[test]
    fn a_distance_is_symmetric() {
        let system = configuration(&[[1.0, 2.0, 3.0], [-4.0, 0.5, 6.0]]);
        let forward: Length<f64, Angstrom> = distance(&system, s(1), s(2));
        let backward: Length<f64, Angstrom> = distance(&system, s(2), s(1));
        assert_eq!(forward, backward);
    }
}
