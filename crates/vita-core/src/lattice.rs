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
    /// (standard crystallographic orientation). Returns `None` if the
    /// parameters describe no realizable cell (inconsistent angles or
    /// non-finite degenerate volume).
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
    use crate::units::length::{Bohr, Nanometer};
    use crate::units::volume::CubicNanometer;

    fn angstrom(v: f64) -> Length<f64, Angstrom> {
        Length::new(v)
    }

    fn cubic_two() -> Lattice<f64> {
        Lattice::cubic(angstrom(2.0)).unwrap()
    }

    #[test]
    fn from_vectors() {
        let l = Lattice::from_vectors(
            Vector3::new(angstrom(3.0), angstrom(0.0), angstrom(0.0)),
            Vector3::new(angstrom(0.0), angstrom(4.0), angstrom(0.0)),
            Vector3::new(angstrom(0.0), angstrom(0.0), angstrom(5.0)),
        )
        .unwrap();
        assert_eq!(l.volume::<CubicAngstrom>().value(), 60.0);
    }

    #[test]
    fn from_vectors_rejects_coplanar() {
        assert!(
            Lattice::from_vectors(
                Vector3::new(angstrom(1.0), angstrom(0.0), angstrom(0.0)),
                Vector3::new(angstrom(2.0), angstrom(0.0), angstrom(0.0)),
                Vector3::new(angstrom(0.0), angstrom(0.0), angstrom(1.0)),
            )
            .is_none()
        );
    }

    #[test]
    fn cubic() {
        assert_eq!(cubic_two().volume::<CubicAngstrom>().value(), 8.0);
    }

    #[test]
    fn cubic_rejects_zero_edge() {
        assert!(Lattice::cubic(angstrom(0.0)).is_none());
    }

    #[test]
    fn cubic_rejects_negative_edge() {
        assert!(Lattice::cubic(angstrom(-1.0)).is_none());
    }

    #[test]
    fn orthorhombic() {
        let l = Lattice::orthorhombic(angstrom(3.0), angstrom(4.0), angstrom(5.0)).unwrap();
        assert_eq!(
            l.edge_lengths::<Angstrom>(),
            Vector3::new(angstrom(3.0), angstrom(4.0), angstrom(5.0))
        );
    }

    #[test]
    fn orthorhombic_rejects_nonpositive_edge() {
        assert!(Lattice::orthorhombic(angstrom(1.0), angstrom(0.0), angstrom(1.0)).is_none());
        assert!(Lattice::orthorhombic(angstrom(1.0), angstrom(1.0), angstrom(-1.0)).is_none());
    }

    #[test]
    fn from_parameters_recovers_cubic() {
        let right = Angle::<f64, Degree>::new(90.0);
        let l = Lattice::from_parameters(
            angstrom(2.0),
            angstrom(2.0),
            angstrom(2.0),
            right,
            right,
            right,
        )
        .unwrap();
        assert!((l.volume::<CubicAngstrom>().value() - 8.0).abs() < 1e-12);
        assert!((l.b::<Angstrom>().y.value() - 2.0).abs() < 1e-12);
    }

    #[test]
    fn from_parameters_rejects_inconsistent_angles() {
        let wide = Angle::<f64, Degree>::new(170.0);
        assert!(
            Lattice::from_parameters(
                angstrom(1.0),
                angstrom(1.0),
                angstrom(1.0),
                wide,
                wide,
                wide
            )
            .is_none()
        );
    }

    #[test]
    fn a() {
        assert_eq!(
            cubic_two().a::<Angstrom>(),
            Vector3::new(angstrom(2.0), angstrom(0.0), angstrom(0.0))
        );
    }

    #[test]
    fn b() {
        assert_eq!(
            cubic_two().b::<Angstrom>(),
            Vector3::new(angstrom(0.0), angstrom(2.0), angstrom(0.0))
        );
    }

    #[test]
    fn c() {
        assert_eq!(
            cubic_two().c::<Angstrom>(),
            Vector3::new(angstrom(0.0), angstrom(0.0), angstrom(2.0))
        );
    }

    #[test]
    fn edge_lengths() {
        assert_eq!(
            cubic_two().edge_lengths::<Angstrom>(),
            Vector3::new(angstrom(2.0), angstrom(2.0), angstrom(2.0))
        );
    }

    #[test]
    fn angles() {
        let l = Lattice::orthorhombic(angstrom(3.0), angstrom(4.0), angstrom(5.0)).unwrap();
        let a = l.angles::<Degree>();
        assert!((a.x.value() - 90.0).abs() < 1e-12);
        assert!((a.y.value() - 90.0).abs() < 1e-12);
        assert!((a.z.value() - 90.0).abs() < 1e-12);
    }

    #[test]
    fn volume() {
        assert_eq!(cubic_two().volume::<CubicAngstrom>().value(), 8.0);
        assert!(
            (Lattice::cubic(angstrom(10.0))
                .unwrap()
                .volume::<CubicNanometer>()
                .value()
                - 1.0)
                .abs()
                < 1e-12
        );
    }

    #[test]
    fn to_cartesian() {
        assert_eq!(
            cubic_two().to_cartesian::<Angstrom>(Point3::new(0.5, 0.5, 0.5)),
            Point3::new(angstrom(1.0), angstrom(1.0), angstrom(1.0))
        );
    }

    #[test]
    fn to_fractional() {
        assert_eq!(
            cubic_two().to_fractional(Point3::new(angstrom(1.0), angstrom(1.0), angstrom(1.0))),
            Point3::new(0.5, 0.5, 0.5)
        );
    }

    #[test]
    fn cartesian_fractional_roundtrip() {
        let l = Lattice::from_vectors(
            Vector3::new(angstrom(4.0), angstrom(1.0), angstrom(0.0)),
            Vector3::new(angstrom(0.0), angstrom(5.0), angstrom(0.0)),
            Vector3::new(angstrom(1.0), angstrom(0.0), angstrom(6.0)),
        )
        .unwrap();
        let frac = Point3::new(0.2, 0.7, 0.4);
        let back = l.to_fractional(l.to_cartesian::<Angstrom>(frac));
        assert!((back.x - frac.x).abs() < 1e-12);
        assert!((back.y - frac.y).abs() < 1e-12);
        assert!((back.z - frac.z).abs() < 1e-12);
    }

    #[test]
    fn accessors_are_unit_generic() {
        let l = Lattice::cubic(Length::<f64, Nanometer>::new(0.5)).unwrap();
        assert!((l.a::<Angstrom>().x.value() - 5.0).abs() < 1e-12);
        assert!((l.a::<Bohr>().x.value() - 5.0 / 0.529_177_210_544).abs() < 1e-9);
    }

    #[test]
    fn supports_f32() {
        let l = Lattice::cubic(Length::<f32, Angstrom>::new(2.0)).unwrap();
        assert_eq!(l.volume::<CubicAngstrom>().value(), 8.0_f32);
        assert_eq!(
            l.to_cartesian::<Angstrom>(Point3::new(0.5_f32, 0.5, 0.5)),
            Point3::new(
                Length::<f32, Angstrom>::new(1.0),
                Length::new(1.0),
                Length::new(1.0)
            )
        );
    }

    #[test]
    fn copy_and_clone() {
        let a = cubic_two();
        let b = a;
        let c = ::core::clone::Clone::clone(&a);
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn eq() {
        assert_eq!(cubic_two(), cubic_two());
        assert_ne!(cubic_two(), Lattice::cubic(angstrom(3.0)).unwrap());
    }

    #[test]
    fn debug() {
        let s = format!("{:?}", cubic_two());
        assert!(s.starts_with("Lattice { basis: Matrix3 {"));
    }
}
