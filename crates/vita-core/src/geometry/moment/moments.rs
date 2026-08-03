use super::{from_angstroms, mean, weighted_by_mass, weighted_evenly};
use crate::tensor::{Matrix3, Point3, Vector3};
use crate::units::area::{Area, AreaUnit, SquareAngstrom};
use crate::units::length::{Length, LengthUnit};
use crate::{HasMasses, HasPositions, Scalar};

/// The first and second moments of a distribution of sites.
///
/// The [`center`](Self::center) the sites are spread about and the
/// [`covariance`](Self::covariance) they are spread by, each normalized by the
/// total weight. Every further reading is one of the covariance's own
/// invariants: its trace is the [`radius_of_gyration`](Self::radius_of_gyration),
/// its eigendecomposition the [`principal_axes`](Self::principal_axes) and
/// [`principal_moments`](Self::principal_moments).
///
/// Obtain via [`moments`] or [`mass_moments`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Moments<V> {
    center: Point3<V>,
    covariance: Matrix3<V>,
}

impl<V: Scalar> Moments<V> {
    /// The center the sites are spread about, in unit `U`.
    ///
    /// Their weighted mean position, and the point the
    /// [`covariance`](Self::covariance) is taken about.
    pub fn center<U: LengthUnit>(self) -> Point3<Length<V, U>> {
        self.center.map(from_angstroms)
    }

    /// The covariance of the sites about their center, in unit `U`.
    ///
    /// The second central moment, normalized by the total weight — the quadratic
    /// form whose level surface is the ellipsoid the sites fill.
    pub fn covariance<U: AreaUnit>(self) -> Matrix3<Area<V, U>> {
        self.covariance.map(from_square_angstroms)
    }

    /// The radius of gyration, in unit `U`.
    ///
    /// The quadratic mean of the distance from the center, and so the root of the
    /// [`covariance`](Self::covariance)'s trace.
    pub fn radius_of_gyration<U: LengthUnit>(self) -> Length<V, U> {
        from_angstroms(self.covariance.trace().sqrt())
    }

    /// The principal axes, as the columns of a right-handed orthonormal frame.
    ///
    /// The eigenvectors of the [`covariance`](Self::covariance), ordered
    /// alongside the [`principal_moments`](Self::principal_moments). Repeated
    /// moments leave their axes undetermined: any frame spanning the same
    /// eigenspace serves, and which one comes out is not fixed.
    pub fn principal_axes(self) -> Matrix3<V> {
        symmetric_eigendecomposition(self.covariance).0
    }

    /// The principal moments, ascending, in unit `U`.
    ///
    /// The spread along each [principal axis](Self::principal_axes). The smallest
    /// names the direction the sites depart from least, so its axis is the normal
    /// of the plane they fit best.
    pub fn principal_moments<U: AreaUnit>(self) -> Vector3<Area<V, U>> {
        symmetric_eigendecomposition(self.covariance)
            .1
            .map(from_square_angstroms)
    }
}

/// The [`Moments`] of the sites, counting each once.
///
/// Their center is the centroid and their covariance the gyration tensor.
/// Returns `None` for a system with no sites.
pub fn moments<S, V>(system: &S) -> Option<Moments<V>>
where
    S: HasPositions<V>,
    V: Scalar,
{
    accumulate(|| weighted_evenly(system))
}

/// The [`Moments`] of the sites, counting each by its mass.
///
/// Their center is the center of mass, and the covariance `C` carries the inertia
/// tensor as `M · (tr(C) · I − C)` for total mass `M`. Returns `None` for a system
/// with no sites, or when no mass falls on any of them.
pub fn mass_moments<S, V>(system: &S) -> Option<Moments<V>>
where
    S: HasPositions<V> + HasMasses<V>,
    V: Scalar,
{
    accumulate(|| weighted_by_mass(system))
}

/// Folds the `weighted` positions into their moments, or `None` if no weight falls
/// on any of them.
///
/// The center is settled first and the departures from it gathered after, which
/// keeps the sums from ever growing far enough apart for their difference to cancel
/// and makes every contribution an outer square — so the covariance comes out
/// symmetric by construction rather than by repair.
fn accumulate<V, I>(weighted: impl Fn() -> I) -> Option<Moments<V>>
where
    V: Scalar,
    I: Iterator<Item = (Point3<V>, V)>,
{
    let center = mean(weighted())?;
    let mut total = V::ZERO;
    let mut comoment = Matrix3::ZERO;
    for (position, weight) in weighted() {
        let departure = position - center;
        total += weight;
        comoment += Matrix3::outer_product(departure, departure) * weight;
    }
    Some(Moments {
        center,
        covariance: comoment / total,
    })
}

/// A bare square-ångström reading wrapped as an area in unit `U`.
fn from_square_angstroms<V: Scalar, U: AreaUnit>(value: V) -> Area<V, U> {
    Area::<V, SquareAngstrom>::new(value).to::<U>()
}

/// The eigendecomposition of the symmetric `matrix`: its eigenvectors as the
/// columns of a right-handed orthonormal frame, and its eigenvalues ascending
/// alongside them.
///
/// Cyclic Jacobi rotations, which leave the frame orthonormal by construction and
/// so stay defined where eigenvalues repeat — the ordinary case here, a collinear
/// or planar distribution having a rank-deficient covariance. Sweeps run only while
/// the off-diagonal keeps shrinking, so the rotations stop at the floor of whatever
/// precision `V` carries, on no threshold and no count: each sweep leaves a strictly
/// smaller non-negative value, of which there are finitely many. `matrix` must be
/// symmetric, as a covariance is by construction.
fn symmetric_eigendecomposition<V: Scalar>(matrix: Matrix3<V>) -> (Matrix3<V>, Vector3<V>) {
    let mut values = [
        matrix.row(0).to_array(),
        matrix.row(1).to_array(),
        matrix.row(2).to_array(),
    ];
    let mut frame = [
        [V::ONE, V::ZERO, V::ZERO],
        [V::ZERO, V::ONE, V::ZERO],
        [V::ZERO, V::ZERO, V::ONE],
    ];
    let off_diagonal = |values: &[[V; 3]; 3]| {
        values[0][1] * values[0][1] + values[0][2] * values[0][2] + values[1][2] * values[1][2]
    };
    let mut remaining = off_diagonal(&values);
    while remaining > V::ZERO {
        for (p, q) in [(0usize, 1usize), (0, 2), (1, 2)] {
            let off = values[p][q];
            if off == V::ZERO {
                continue;
            }
            let theta = (values[q][q] - values[p][p]) / (off + off);
            let tangent = theta.signum() / (theta.abs() + (theta * theta + V::ONE).sqrt());
            let cosine = (tangent * tangent + V::ONE).sqrt().recip();
            let sine = tangent * cosine;

            values[p][p] -= tangent * off;
            values[q][q] += tangent * off;
            values[p][q] = V::ZERO;
            values[q][p] = V::ZERO;

            let rest = 3 - p - q;
            let (near, far) = (values[p][rest], values[q][rest]);
            values[p][rest] = cosine * near - sine * far;
            values[rest][p] = values[p][rest];
            values[q][rest] = sine * near + cosine * far;
            values[rest][q] = values[q][rest];

            for row in &mut frame {
                let (near, far) = (row[p], row[q]);
                row[p] = cosine * near - sine * far;
                row[q] = sine * near + cosine * far;
            }
        }
        let shrunk = off_diagonal(&values);
        if shrunk >= remaining {
            break;
        }
        remaining = shrunk;
    }

    let mut eigenvalues = [values[0][0], values[1][1], values[2][2]];
    for (lower, upper) in [(0usize, 1usize), (0, 2), (1, 2)] {
        if eigenvalues[upper] < eigenvalues[lower] {
            eigenvalues.swap(lower, upper);
            for row in &mut frame {
                row.swap(lower, upper);
            }
        }
    }

    let axes = Matrix3::from_rows(
        Vector3::from_array(frame[0]),
        Vector3::from_array(frame[1]),
        Vector3::from_array(frame[2]),
    );
    let axes = if axes.determinant() < V::ZERO {
        Matrix3::from_cols(axes.col(0), axes.col(1), -axes.col(2))
    } else {
        axes
    };
    (axes, Vector3::from_array(eigenvalues))
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::geometry::fixture::{System, close, configuration, weighted};
    use crate::units::area::SquareNanometer;
    use crate::units::length::{Angstrom, Nanometer};

    fn pair() -> System {
        configuration(&[[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    }

    fn square() -> System {
        configuration(&[
            [1.0, 1.0, 0.0],
            [-1.0, 1.0, 0.0],
            [-1.0, -1.0, 0.0],
            [1.0, -1.0, 0.0],
        ])
    }

    fn rhombus() -> System {
        configuration(&[
            [2.0, 0.0, 0.0],
            [-2.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ])
    }

    fn diamond() -> System {
        configuration(&[
            [2.0, 2.0, 0.0],
            [-2.0, -2.0, 0.0],
            [1.0, -1.0, 0.0],
            [-1.0, 1.0, 0.0],
        ])
    }

    fn diagonal() -> Matrix3<f64> {
        Matrix3::from_diagonal(Vector3::new(2.0, 1.0, 3.0))
    }

    fn tilted() -> Matrix3<f64> {
        Matrix3::from_rows(
            Vector3::new(2.5, 1.5, 0.0),
            Vector3::new(1.5, 2.5, 0.0),
            Vector3::new(0.0, 0.0, 0.0),
        )
    }

    fn repeated() -> Matrix3<f64> {
        Matrix3::from_rows(
            Vector3::new(2.0, 1.0, 0.0),
            Vector3::new(1.0, 2.0, 0.0),
            Vector3::new(0.0, 0.0, 3.0),
        )
    }

    fn saddle() -> Matrix3<f64> {
        Matrix3::from_rows(
            Vector3::new(4.0, 1.0, 2.0),
            Vector3::new(1.0, 3.0, -1.0),
            Vector3::new(2.0, -1.0, -5.0),
        )
    }

    #[test]
    fn the_moments_of_an_empty_system_are_absent() {
        let taken: Option<Moments<f64>> = moments(&configuration(&[]));
        assert!(taken.is_none());
    }

    #[test]
    fn the_moments_of_one_site_center_on_it() {
        let moment = moments(&configuration(&[[1.0, 2.0, 3.0]])).unwrap();
        let expected = Point3::new(Length::new(1.0), Length::new(2.0), Length::new(3.0));
        assert_eq!(moment.center::<Angstrom>(), expected);
    }

    #[test]
    fn the_covariance_of_one_site_vanishes() {
        let moment = moments(&configuration(&[[1.0, 2.0, 3.0]])).unwrap();
        assert_eq!(moment.covariance::<SquareAngstrom>(), Matrix3::ZERO);
    }

    #[test]
    fn a_diagonal_matrix_keeps_its_diagonal_as_its_eigenvalues() {
        let (_, eigenvalues) = symmetric_eigendecomposition(diagonal());
        assert_eq!(eigenvalues, Vector3::new(1.0, 2.0, 3.0));
    }

    #[test]
    fn moments_center_on_the_mean_of_the_positions() {
        let moment = moments(&square()).unwrap();
        let center = moment.center::<Angstrom>();
        assert!(close(center.x, Length::new(0.0)) && close(center.y, Length::new(0.0)));
    }

    #[test]
    fn mass_moments_center_on_the_center_of_mass() {
        let system = weighted(&[([0.0, 0.0, 0.0], 1.0), ([10.0, 0.0, 0.0], 3.0)]);
        let moment = mass_moments(&system).unwrap();
        assert_eq!(moment.center::<Angstrom>().x, Length::new(7.5));
    }

    #[test]
    fn the_covariance_measures_the_spread_about_the_center() {
        let moment = moments(&pair()).unwrap();
        let expected =
            Matrix3::from_diagonal(Vector3::new(Area::new(1.0), Area::new(0.0), Area::new(0.0)));
        assert_eq!(moment.covariance::<SquareAngstrom>(), expected);
    }

    #[test]
    fn the_radius_of_gyration_is_the_root_of_the_covariance_trace() {
        let moment = moments(&square()).unwrap();
        let spread: Length<f64, Angstrom> = moment.radius_of_gyration();
        assert!(close(spread, Length::new(2.0_f64.sqrt())));
    }

    #[test]
    fn the_principal_moments_ascend() {
        let moment = moments(&rhombus()).unwrap();
        let spread = moment.principal_moments::<SquareAngstrom>();
        assert!(spread.x <= spread.y && spread.y <= spread.z);
    }

    #[test]
    fn the_principal_axes_follow_the_directions_of_greatest_spread() {
        let widest = moments(&rhombus()).unwrap().principal_axes().col(2);
        assert!(close(widest.dot(Vector3::X).abs(), 1.0));
    }

    #[test]
    fn the_principal_axes_turn_with_a_tilted_distribution() {
        let widest = moments(&diamond()).unwrap().principal_axes().col(2);
        let diagonal = Vector3::new(1.0, 1.0, 0.0).normalize();
        assert!(close(widest.dot(diagonal).abs(), 1.0));
    }

    #[test]
    fn the_eigenvalues_ascend() {
        let (_, eigenvalues) = symmetric_eigendecomposition(tilted());
        assert!(
            close(eigenvalues.x, 0.0) && close(eigenvalues.y, 1.0) && close(eigenvalues.z, 4.0)
        );
    }

    #[test]
    fn the_eigenvectors_are_orthonormal() {
        let (frame, _) = symmetric_eigendecomposition(saddle());
        assert!(
            (0..3).all(|index| close(frame.col(index).norm_squared(), 1.0))
                && close(frame.col(0).dot(frame.col(1)), 0.0)
                && close(frame.col(0).dot(frame.col(2)), 0.0)
                && close(frame.col(1).dot(frame.col(2)), 0.0)
        );
    }

    #[test]
    fn the_eigenvectors_form_a_right_handed_frame() {
        let (frame, _) = symmetric_eigendecomposition(diagonal());
        assert!(close(frame.determinant(), 1.0));
    }

    #[test]
    fn the_center_is_given_in_the_requested_unit() {
        let moment = moments(&configuration(&[[10.0, 0.0, 0.0]])).unwrap();
        assert!(close(moment.center::<Nanometer>().x, Length::new(1.0)));
    }

    #[test]
    fn the_covariance_is_given_in_the_requested_unit() {
        let moment = moments(&pair()).unwrap();
        assert!(close(
            moment.covariance::<SquareNanometer>().col(0).x,
            Area::new(0.01)
        ));
    }

    #[test]
    fn the_radius_of_gyration_is_given_in_the_requested_unit() {
        let moment = moments(&pair()).unwrap();
        assert!(close(
            moment.radius_of_gyration::<Nanometer>(),
            Length::new(0.1)
        ));
    }

    #[test]
    fn the_principal_moments_are_given_in_the_requested_unit() {
        let moment = moments(&pair()).unwrap();
        let spread = moment.principal_moments::<SquareNanometer>();
        assert!(close(spread.z, Area::new(0.01)));
    }

    #[test]
    fn the_mass_moments_of_weightless_sites_are_absent() {
        let system = weighted(&[([0.0, 0.0, 0.0], 0.0), ([1.0, 0.0, 0.0], 0.0)]);
        let taken: Option<Moments<f64>> = mass_moments(&system);
        assert!(taken.is_none());
    }

    #[test]
    fn collinear_sites_leave_two_principal_moments_at_zero() {
        let moment = moments(&pair()).unwrap();
        let spread = moment.principal_moments::<SquareAngstrom>();
        assert!(close(spread.x, Area::new(0.0)) && close(spread.y, Area::new(0.0)));
    }

    #[test]
    fn coplanar_sites_leave_the_narrowest_axis_along_the_plane_normal() {
        let narrowest = moments(&square()).unwrap().principal_axes().col(0);
        assert!(close(narrowest.dot(Vector3::Z).abs(), 1.0));
    }

    #[test]
    fn repeated_eigenvalues_leave_the_eigenvectors_orthonormal() {
        let (frame, eigenvalues) = symmetric_eigendecomposition(repeated());
        assert!(
            close(eigenvalues.y, eigenvalues.z)
                && (0..3).all(|index| close(frame.col(index).norm_squared(), 1.0))
                && close(frame.col(1).dot(frame.col(2)), 0.0)
        );
    }

    #[test]
    fn an_indefinite_matrix_decomposes_alike() {
        let (_, eigenvalues) = symmetric_eigendecomposition(saddle());
        assert!(eigenvalues.x < 0.0 && eigenvalues.z > 0.0);
    }

    #[test]
    fn the_covariance_is_symmetric() {
        let covariance = moments(&diamond()).unwrap().covariance::<SquareAngstrom>();
        assert_eq!(covariance, covariance.transpose());
    }

    #[test]
    fn the_decomposition_reconstructs_the_matrix() {
        let (frame, eigenvalues) = symmetric_eigendecomposition(saddle());
        let rebuilt = frame * Matrix3::from_diagonal(eigenvalues) * frame.transpose();
        assert!(
            rebuilt
                .to_cols_array()
                .iter()
                .zip(saddle().to_cols_array())
                .all(|(&taken, expected)| close(taken, expected))
        );
    }

    #[test]
    fn the_principal_moments_are_independent_of_the_order_of_the_sites() {
        let shuffled = configuration(&[
            [0.0, 1.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [-2.0, 0.0, 0.0],
        ]);
        let taken = moments(&rhombus()).unwrap();
        let again = moments(&shuffled).unwrap();
        let one = taken.principal_moments::<SquareAngstrom>();
        let other = again.principal_moments::<SquareAngstrom>();
        assert!(close(one.x, other.x) && close(one.y, other.y) && close(one.z, other.z));
    }
}
