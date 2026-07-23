use crate::Scalar;
use crate::tensor::{Matrix3, Point3, Vector3};
use crate::units::angle::{Angle, AngleUnit, Radian};
use crate::units::length::{Angstrom, Length, LengthUnit};
use crate::units::volume::{CubicAngstrom, Volume, VolumeUnit};

/// A periodic lattice: the three basis vectors **a**, **b**, **c** that
/// generate a system's translational symmetry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Lattice<V> {
    basis: Matrix3<V>,
}

impl<V: Scalar> Lattice<V> {
    /// Wraps `basis` (columns are lattice vectors in ångströms), returning
    /// `None` if the matrix is degenerate (zero or non-finite volume).
    #[inline]
    fn from_basis(basis: Matrix3<V>) -> Option<Self> {
        if basis.is_invertible() {
            Some(Self { basis })
        } else {
            None
        }
    }

    /// Constructs a lattice from its three basis vectors, returning `None` if
    /// they are coplanar (zero enclosed volume).
    #[inline]
    pub fn from_vectors<U: LengthUnit>(
        a: Vector3<Length<V, U>>,
        b: Vector3<Length<V, U>>,
        c: Vector3<Length<V, U>>,
    ) -> Option<Self> {
        Self::from_basis(Matrix3::from_cols(
            a.map(into_angstroms),
            b.map(into_angstroms),
            c.map(into_angstroms),
        ))
    }

    /// Constructs a cubic lattice with edge length `side`, returning `None` if
    /// `side` is not strictly positive and finite.
    #[inline]
    pub fn cubic<U: LengthUnit>(side: Length<V, U>) -> Option<Self> {
        let s = into_angstroms(side);
        if !(s > V::ZERO && s.is_finite()) {
            return None;
        }
        Self::from_basis(Matrix3::from_diagonal(Vector3::splat(s)))
    }

    /// Constructs an orthorhombic (axis-aligned) lattice from its three edge
    /// lengths, returning `None` if any length is not strictly positive and
    /// finite.
    #[inline]
    pub fn orthorhombic<U: LengthUnit>(
        a: Length<V, U>,
        b: Length<V, U>,
        c: Length<V, U>,
    ) -> Option<Self> {
        let (a, b, c) = (into_angstroms(a), into_angstroms(b), into_angstroms(c));
        let positive = |s: V| s > V::ZERO && s.is_finite();
        if !(positive(a) && positive(b) && positive(c)) {
            return None;
        }
        Self::from_basis(Matrix3::from_diagonal(Vector3::new(a, b, c)))
    }

    /// Constructs a lattice from the crystallographic parameters: edge lengths
    /// `a`, `b`, `c` and inter-edge angles `alpha` (∠**b**,**c**), `beta`
    /// (∠**a**,**c**), `gamma` (∠**a**,**b**).
    ///
    /// **a** is placed along the `x` axis and **b** in the `xy` plane
    /// (standard crystallographic orientation). Returns `None` if any edge
    /// length is not strictly positive and finite, or if the parameters
    /// describe no realizable cell (inconsistent angles or degenerate
    /// volume).
    pub fn from_parameters<L: LengthUnit, A: AngleUnit>(
        a: Length<V, L>,
        b: Length<V, L>,
        c: Length<V, L>,
        alpha: Angle<V, A>,
        beta: Angle<V, A>,
        gamma: Angle<V, A>,
    ) -> Option<Self> {
        let a = into_angstroms(a);
        let b = into_angstroms(b);
        let c = into_angstroms(c);
        let positive = |s: V| s > V::ZERO && s.is_finite();
        if !(positive(a) && positive(b) && positive(c)) {
            return None;
        }
        let (_, ca) = alpha.to::<Radian>().value().sin_cos();
        let (_, cb) = beta.to::<Radian>().value().sin_cos();
        let (sg, cg) = gamma.to::<Radian>().value().sin_cos();

        let col_a = Vector3::new(a, V::ZERO, V::ZERO);
        let col_b = Vector3::new(b * cg, b * sg, V::ZERO);
        let cx = c * cb;
        let cy = c * (ca - cb * cg) / sg;
        let cz = c
            * (V::ONE - ca * ca - cb * cb - cg * cg + V::from_f64(2.0) * ca * cb * cg).sqrt()
            / sg;

        Self::from_basis(Matrix3::from_cols(col_a, col_b, Vector3::new(cx, cy, cz)))
    }

    /// Returns the first lattice vector **a**, in unit `U`.
    #[inline]
    pub fn a<U: LengthUnit>(&self) -> Vector3<Length<V, U>> {
        self.basis.col(0).map(from_angstroms)
    }

    /// Returns the second lattice vector **b**, in unit `U`.
    #[inline]
    pub fn b<U: LengthUnit>(&self) -> Vector3<Length<V, U>> {
        self.basis.col(1).map(from_angstroms)
    }

    /// Returns the third lattice vector **c**, in unit `U`.
    #[inline]
    pub fn c<U: LengthUnit>(&self) -> Vector3<Length<V, U>> {
        self.basis.col(2).map(from_angstroms)
    }

    /// Returns the three edge lengths `(|a|, |b|, |c|)`, in unit `U`.
    #[inline]
    pub fn edge_lengths<U: LengthUnit>(&self) -> Vector3<Length<V, U>> {
        Vector3::new(
            from_angstroms(self.basis.col(0).norm()),
            from_angstroms(self.basis.col(1).norm()),
            from_angstroms(self.basis.col(2).norm()),
        )
    }

    /// Returns the three inter-edge angles `(α, β, γ)` — ∠(**b**,**c**),
    /// ∠(**a**,**c**), ∠(**a**,**b**) — in unit `U`.
    #[inline]
    pub fn angles<U: AngleUnit>(&self) -> Vector3<Angle<V, U>> {
        let a = self.basis.col(0);
        let b = self.basis.col(1);
        let c = self.basis.col(2);
        Vector3::new(
            Angle::<V, Radian>::new(b.angle_between(c)).to::<U>(),
            Angle::<V, Radian>::new(a.angle_between(c)).to::<U>(),
            Angle::<V, Radian>::new(a.angle_between(b)).to::<U>(),
        )
    }

    /// Returns the enclosed volume `|det(basis)|`, in unit `U`.
    #[inline]
    pub fn volume<U: VolumeUnit>(&self) -> Volume<V, U> {
        Volume::<V, CubicAngstrom>::new(self.basis.determinant().abs()).to::<U>()
    }

    /// Maps a fractional coordinate to its Cartesian position `r = basis · f`,
    /// in unit `U`.
    #[inline]
    pub fn to_cartesian<U: LengthUnit>(&self, fractional: Point3<V>) -> Point3<Length<V, U>> {
        Point3::from_vector((self.basis * fractional.to_vector()).map(from_angstroms))
    }

    /// Maps a Cartesian position to its (dimensionless) fractional coordinate
    /// `f = basis⁻¹ · r`.
    #[inline]
    pub fn to_fractional<U: LengthUnit>(&self, cartesian: Point3<Length<V, U>>) -> Point3<V> {
        Point3::from_vector(self.basis.inverse() * cartesian.to_vector().map(into_angstroms))
    }
}

/// Converts a length to its bare scalar value in ångströms.
#[inline]
fn into_angstroms<V: Scalar, U: LengthUnit>(length: Length<V, U>) -> V {
    length.to::<Angstrom>().value()
}

/// Wraps a bare ångström scalar as a length in unit `U`.
#[inline]
fn from_angstroms<V: Scalar, U: LengthUnit>(value: V) -> Length<V, U> {
    Length::<V, Angstrom>::new(value).to::<U>()
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::units::angle::Degree;
    use crate::units::length::Nanometer;

    const TOL: f64 = 1e-9;

    fn length(x: f64) -> Length<f64, Angstrom> {
        Length::new(x)
    }

    fn degrees(x: f64) -> Angle<f64, Degree> {
        Angle::new(x)
    }

    fn length_values<U: LengthUnit>(v: Vector3<Length<f64, U>>) -> [f64; 3] {
        v.to_array().map(|l| l.value())
    }

    fn angle_values(v: Vector3<Angle<f64, Degree>>) -> [f64; 3] {
        v.to_array().map(|a| a.value())
    }

    fn assert_close(a: f64, b: f64) {
        assert!((a - b).abs() < TOL, "{a} not close to {b}");
    }

    fn assert_all_close(actual: [f64; 3], expected: [f64; 3]) {
        for (a, b) in actual.iter().zip(expected.iter()) {
            assert!((a - b).abs() < TOL, "{a} not close to {b}");
        }
    }

    #[test]
    fn cubic_has_equal_edge_lengths() {
        let lat = Lattice::cubic(length(2.0)).unwrap();
        assert_all_close(
            length_values(lat.edge_lengths::<Angstrom>()),
            [2.0, 2.0, 2.0],
        );
    }

    #[test]
    fn cubic_is_orthogonal() {
        let lat = Lattice::cubic(length(2.0)).unwrap();
        assert_all_close(angle_values(lat.angles::<Degree>()), [90.0, 90.0, 90.0]);
    }

    #[test]
    fn cubic_volume_is_the_side_cubed() {
        let lat = Lattice::cubic(length(2.0)).unwrap();
        assert_close(lat.volume::<CubicAngstrom>().value(), 8.0);
    }

    #[test]
    fn cubic_basis_vectors_are_axis_aligned() {
        let lat = Lattice::cubic(length(2.0)).unwrap();
        assert_all_close(length_values(lat.a::<Angstrom>()), [2.0, 0.0, 0.0]);
        assert_all_close(length_values(lat.b::<Angstrom>()), [0.0, 2.0, 0.0]);
        assert_all_close(length_values(lat.c::<Angstrom>()), [0.0, 0.0, 2.0]);
    }

    #[test]
    fn orthorhombic_edge_lengths_match_its_arguments() {
        let lat = Lattice::orthorhombic(length(2.0), length(3.0), length(4.0)).unwrap();
        assert_all_close(
            length_values(lat.edge_lengths::<Angstrom>()),
            [2.0, 3.0, 4.0],
        );
    }

    #[test]
    fn orthorhombic_volume_is_the_product_of_edges() {
        let lat = Lattice::orthorhombic(length(2.0), length(3.0), length(4.0)).unwrap();
        assert_close(lat.volume::<CubicAngstrom>().value(), 24.0);
    }

    #[test]
    fn from_vectors_recovers_its_basis_vectors() {
        let a = Vector3::new(length(2.0), length(0.0), length(0.0));
        let b = Vector3::new(length(0.0), length(3.0), length(0.0));
        let c = Vector3::new(length(0.0), length(0.0), length(4.0));
        let lat = Lattice::from_vectors(a, b, c).unwrap();
        assert_all_close(length_values(lat.a::<Angstrom>()), [2.0, 0.0, 0.0]);
        assert_all_close(length_values(lat.b::<Angstrom>()), [0.0, 3.0, 0.0]);
        assert_all_close(length_values(lat.c::<Angstrom>()), [0.0, 0.0, 4.0]);
    }

    #[test]
    fn from_parameters_reproduces_a_cubic_cell() {
        let lat = Lattice::from_parameters(
            length(2.0),
            length(2.0),
            length(2.0),
            degrees(90.0),
            degrees(90.0),
            degrees(90.0),
        )
        .unwrap();
        assert_all_close(
            length_values(lat.edge_lengths::<Angstrom>()),
            [2.0, 2.0, 2.0],
        );
        assert_all_close(angle_values(lat.angles::<Degree>()), [90.0, 90.0, 90.0]);
        assert_close(lat.volume::<CubicAngstrom>().value(), 8.0);
    }

    #[test]
    fn to_cartesian_maps_fractional_coordinates_to_positions() {
        let lat = Lattice::cubic(length(2.0)).unwrap();
        let cart = lat.to_cartesian::<Angstrom>(Point3::new(0.5, 0.5, 0.5));
        assert_all_close(cart.to_array().map(|l| l.value()), [1.0, 1.0, 1.0]);
    }

    #[test]
    fn lattice_vectors_convert_to_the_requested_unit() {
        let lat = Lattice::cubic(Length::<f64, Nanometer>::new(1.0)).unwrap();
        assert_all_close(length_values(lat.a::<Angstrom>()), [10.0, 0.0, 0.0]);
    }

    #[test]
    fn volume_is_the_absolute_value_of_the_determinant() {
        let a = Vector3::new(length(2.0), length(0.0), length(0.0));
        let b = Vector3::new(length(0.0), length(0.0), length(3.0));
        let c = Vector3::new(length(0.0), length(4.0), length(0.0));
        let lat = Lattice::from_vectors(a, b, c).unwrap();
        assert_close(lat.volume::<CubicAngstrom>().value(), 24.0);
    }

    #[test]
    fn cubic_rejects_a_zero_side() {
        assert!(Lattice::cubic(length(0.0)).is_none());
    }

    #[test]
    fn cubic_rejects_a_nonfinite_side() {
        assert!(Lattice::cubic(length(f64::INFINITY)).is_none());
        assert!(Lattice::cubic(length(f64::NAN)).is_none());
    }

    #[test]
    fn orthorhombic_rejects_a_nonpositive_edge() {
        assert!(Lattice::orthorhombic(length(0.0), length(3.0), length(4.0)).is_none());
        assert!(Lattice::orthorhombic(length(2.0), length(-3.0), length(4.0)).is_none());
    }

    #[test]
    fn from_vectors_rejects_coplanar_vectors() {
        let a = Vector3::new(length(1.0), length(0.0), length(0.0));
        let b = Vector3::new(length(0.0), length(1.0), length(0.0));
        let c = Vector3::new(length(1.0), length(1.0), length(0.0));
        assert!(Lattice::from_vectors(a, b, c).is_none());
    }

    #[test]
    fn from_parameters_rejects_a_nonpositive_edge() {
        let lat = Lattice::from_parameters(
            length(0.0),
            length(2.0),
            length(2.0),
            degrees(90.0),
            degrees(90.0),
            degrees(90.0),
        );
        assert!(lat.is_none());
    }

    #[test]
    fn from_parameters_rejects_an_unrealizable_cell() {
        let lat = Lattice::from_parameters(
            length(2.0),
            length(2.0),
            length(2.0),
            degrees(170.0),
            degrees(170.0),
            degrees(170.0),
        );
        assert!(lat.is_none());
    }

    #[test]
    fn from_parameters_roundtrips_edge_lengths_and_angles() {
        let lat = Lattice::from_parameters(
            length(3.0),
            length(4.0),
            length(5.0),
            degrees(70.0),
            degrees(80.0),
            degrees(100.0),
        )
        .unwrap();
        assert_all_close(
            length_values(lat.edge_lengths::<Angstrom>()),
            [3.0, 4.0, 5.0],
        );
        assert_all_close(angle_values(lat.angles::<Degree>()), [70.0, 80.0, 100.0]);
    }

    #[test]
    fn to_fractional_inverts_to_cartesian() {
        let lat = Lattice::orthorhombic(length(2.0), length(3.0), length(4.0)).unwrap();
        let frac = Point3::new(0.1, 0.2, 0.7);
        let cart = lat.to_cartesian::<Angstrom>(frac);
        let back = lat.to_fractional::<Angstrom>(cart);
        assert_all_close(back.to_array(), frac.to_array());
    }
}
