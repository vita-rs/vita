//! The moments of a distribution of sites.
//!
//! [`centroid`] and [`center_of_mass`] give the first moment alone, counting each
//! site once and counting it by its mass. [`moments`] and [`mass_moments`] carry the
//! second alongside it as a [`Moments`]: the center the sites are spread about and
//! the covariance they are spread by, whose trace gives the
//! [`radius_of_gyration`](Moments::radius_of_gyration) and whose eigendecomposition
//! gives the shape. The second moment is the last a quadratic form on
//! three-dimensional space can hold; the named combinations of its eigenvalues that
//! each field keeps are the caller's.
//!
//! Every reading takes the coordinates as they stand. A moment is taken about a mean of
//! positions, and a torus carries their addition but not the division that follows it: a
//! periodic distribution has moments only once lifted back into space, and neither the
//! coordinates nor the cell settles the lift.

mod center;
mod moments;

pub use center::{center_of_mass, centroid};
pub use moments::{Moments, mass_moments, moments};

use crate::tensor::{Point3, Vector3};
use crate::units::length::{Angstrom, Length, LengthUnit};
use crate::units::mass::Dalton;
use crate::{HasMasses, HasPositions, Quantity, Scalar};

/// The positions in ångströms, each paired with a weight of one.
fn weighted_evenly<S, V>(system: &S) -> impl Iterator<Item = (Point3<V>, V)> + '_
where
    S: HasPositions<V>,
    V: Scalar,
{
    system
        .positions::<Angstrom>()
        .map(|position| (position.map(Quantity::value), V::ONE))
}

/// The positions in ångströms, each paired with the mass at its site.
fn weighted_by_mass<S, V>(system: &S) -> impl Iterator<Item = (Point3<V>, V)> + '_
where
    S: HasPositions<V> + HasMasses<V>,
    V: Scalar,
{
    system
        .positions::<Angstrom>()
        .zip(system.masses::<Dalton>())
        .map(|(position, mass)| (position.map(Quantity::value), mass.value()))
}

/// The weighted mean of the `weighted` positions, or `None` if no weight falls on
/// any of them.
///
/// The one place a center is settled, so that a center taken alone and a center
/// carried by [`Moments`] are the same reading down to the last bit.
fn mean<V: Scalar>(weighted: impl Iterator<Item = (Point3<V>, V)>) -> Option<Point3<V>> {
    let mut total = V::ZERO;
    let mut sum = Vector3::ZERO;
    for (position, weight) in weighted {
        total += weight;
        sum += position.to_vector() * weight;
    }
    (total > V::ZERO).then(|| Point3::from_vector(sum / total))
}

/// A bare ångström reading wrapped as a length in unit `U`.
fn from_angstroms<V: Scalar, U: LengthUnit>(value: V) -> Length<V, U> {
    Length::<V, Angstrom>::new(value).to::<U>()
}
