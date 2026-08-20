use super::{from_angstroms, mean, weighted_by_mass, weighted_evenly};
use crate::tensor::Point3;
use crate::units::length::{Length, LengthUnit};
use crate::{HasMasses, HasPositions, Scalar};

/// The centroid of the sites, in unit `U`.
///
/// Their mean position, counting each once. Returns `None` for a system with no
/// sites.
pub fn centroid<S, V, U>(system: &S) -> Option<Point3<Length<V, U>>>
where
    S: HasPositions<V>,
    V: Scalar,
    U: LengthUnit,
{
    mean(weighted_evenly(system)).map(|center| center.map(from_angstroms))
}

/// The center of mass of the sites, in unit `U`.
///
/// Their mean position, counting each by its mass. Returns `None` for a system with
/// no sites, or when no mass falls on any of them.
pub fn center_of_mass<S, V, U>(system: &S) -> Option<Point3<Length<V, U>>>
where
    S: HasPositions<V> + HasMasses<V>,
    V: Scalar,
    U: LengthUnit,
{
    mean(weighted_by_mass(system)).map(|center| center.map(from_angstroms))
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::geometry::fixture::{close, configuration, points_close, weighted};
    use crate::geometry::moment::moments;
    use crate::units::length::{Angstrom, Nanometer};

    #[test]
    fn the_centroid_of_an_empty_system_is_absent() {
        let center: Option<Point3<Length<f64, Angstrom>>> = centroid(&configuration(&[]));
        assert!(center.is_none());
    }

    #[test]
    fn the_centroid_of_one_site_is_that_site() {
        let system = configuration(&[[1.0, 2.0, 3.0]]);
        let center: Point3<Length<f64, Angstrom>> = centroid(&system).unwrap();
        let expected = Point3::new(Length::new(1.0), Length::new(2.0), Length::new(3.0));
        assert_eq!(center, expected);
    }

    #[test]
    fn the_centroid_is_the_mean_of_the_positions() {
        let system = configuration(&[[0.0, 0.0, 0.0], [2.0, 4.0, 6.0]]);
        let center: Point3<Length<f64, Angstrom>> = centroid(&system).unwrap();
        let expected = Point3::new(Length::new(1.0), Length::new(2.0), Length::new(3.0));
        assert_eq!(center, expected);
    }

    #[test]
    fn the_center_of_mass_leans_toward_the_heavier_site() {
        let system = weighted(&[([0.0, 0.0, 0.0], 1.0), ([10.0, 0.0, 0.0], 3.0)]);
        let center: Point3<Length<f64, Angstrom>> = center_of_mass(&system).unwrap();
        assert_eq!(center.x, Length::new(7.5));
    }

    #[test]
    fn the_centroid_is_given_in_the_requested_unit() {
        let system = configuration(&[[10.0, 0.0, 0.0]]);
        let center: Point3<Length<f64, Nanometer>> = centroid(&system).unwrap();
        assert!(close(center.x, Length::new(1.0)));
    }

    #[test]
    fn the_center_of_mass_is_given_in_the_requested_unit() {
        let system = weighted(&[([10.0, 0.0, 0.0], 2.0)]);
        let center: Point3<Length<f64, Nanometer>> = center_of_mass(&system).unwrap();
        assert!(close(center.x, Length::new(1.0)));
    }

    #[test]
    fn the_center_of_mass_of_weightless_sites_is_absent() {
        let system = weighted(&[([0.0, 0.0, 0.0], 0.0), ([1.0, 0.0, 0.0], 0.0)]);
        let center: Option<Point3<Length<f64, Angstrom>>> = center_of_mass(&system);
        assert!(center.is_none());
    }

    #[test]
    fn the_center_of_mass_of_equal_masses_is_the_centroid() {
        let system = weighted(&[([1.0, 0.0, 0.0], 2.0), ([3.0, 0.0, 0.0], 2.0)]);
        let weighed: Point3<Length<f64, Angstrom>> = center_of_mass(&system).unwrap();
        let counted: Point3<Length<f64, Angstrom>> = centroid(&system).unwrap();
        assert_eq!(weighed, counted);
    }

    #[test]
    fn the_centroid_agrees_with_the_center_the_moments_carry() {
        let system = configuration(&[[1.0, 2.0, 3.0], [-4.0, 0.5, 6.0], [7.0, -8.0, 0.25]]);
        let alone: Point3<Length<f64, Angstrom>> = centroid(&system).unwrap();
        let carried = moments(&system).unwrap().center::<Angstrom>();
        assert!(points_close(alone, carried));
    }
}
