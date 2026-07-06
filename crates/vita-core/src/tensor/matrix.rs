use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use super::Vector3;
use crate::Scalar;

/// A 3×3 matrix stored as three column vectors.
pub struct Matrix3<T> {
    cols: [Vector3<T>; 3],
}

impl<T: ::core::marker::Copy> ::core::marker::Copy for Matrix3<T> {}

impl<T: ::core::clone::Clone> ::core::clone::Clone for Matrix3<T> {
    #[inline]
    fn clone(&self) -> Self {
        Self {
            cols: [
                self.cols[0].clone(),
                self.cols[1].clone(),
                self.cols[2].clone(),
            ],
        }
    }
}

impl<T: ::core::fmt::Debug> ::core::fmt::Debug for Matrix3<T> {
    fn fmt(&self, f: &mut ::core::fmt::Formatter<'_>) -> ::core::fmt::Result {
        f.debug_struct("Matrix3")
            .field("x_col", &self.cols[0])
            .field("y_col", &self.cols[1])
            .field("z_col", &self.cols[2])
            .finish()
    }
}

impl<T: ::core::cmp::PartialEq> ::core::cmp::PartialEq for Matrix3<T> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.cols[0] == other.cols[0]
            && self.cols[1] == other.cols[1]
            && self.cols[2] == other.cols[2]
    }
}

impl<T> Matrix3<T> {
    /// Constructs a matrix from its three column vectors.
    #[inline]
    pub const fn from_cols(x_col: Vector3<T>, y_col: Vector3<T>, z_col: Vector3<T>) -> Self {
        Self {
            cols: [x_col, y_col, z_col],
        }
    }

    /// Applies `f` to every element in column-major order, returning the
    /// resulting matrix.
    #[inline]
    pub fn map<U, F: FnMut(T) -> U>(self, mut f: F) -> Matrix3<U> {
        let [c0, c1, c2] = self.cols;
        Matrix3::from_cols(c0.map(&mut f), c1.map(&mut f), c2.map(&mut f))
    }

    /// Combines `self` and `rhs` element-wise in column-major order through
    /// `f`.
    #[inline]
    pub fn zip_map<U, R, F: FnMut(T, U) -> R>(self, rhs: Matrix3<U>, mut f: F) -> Matrix3<R> {
        let [c0, c1, c2] = self.cols;
        let [r0, r1, r2] = rhs.cols;
        Matrix3::from_cols(
            c0.zip_map(r0, &mut f),
            c1.zip_map(r1, &mut f),
            c2.zip_map(r2, &mut f),
        )
    }
}

impl<T: Copy> Matrix3<T> {
    /// Constructs a matrix from its three row vectors.
    #[inline]
    pub fn from_rows(x_row: Vector3<T>, y_row: Vector3<T>, z_row: Vector3<T>) -> Self {
        Self::from_cols(
            Vector3::new(x_row.x, y_row.x, z_row.x),
            Vector3::new(x_row.y, y_row.y, z_row.y),
            Vector3::new(x_row.z, y_row.z, z_row.z),
        )
    }

    /// Constructs a matrix from a column-major array of nine elements.
    #[inline]
    pub fn from_cols_array(m: &[T; 9]) -> Self {
        Self::from_cols(
            Vector3::new(m[0], m[1], m[2]),
            Vector3::new(m[3], m[4], m[5]),
            Vector3::new(m[6], m[7], m[8]),
        )
    }

    /// Returns the elements as a column-major array of nine elements.
    #[inline]
    pub fn to_cols_array(&self) -> [T; 9] {
        [
            self.cols[0].x,
            self.cols[0].y,
            self.cols[0].z,
            self.cols[1].x,
            self.cols[1].y,
            self.cols[1].z,
            self.cols[2].x,
            self.cols[2].y,
            self.cols[2].z,
        ]
    }

    /// Returns the column at `index`.
    ///
    /// # Panics
    ///
    /// Panics if `index` is greater than `2`.
    #[inline]
    pub fn col(&self, index: usize) -> Vector3<T> {
        self.cols[index]
    }

    /// Returns the row at `index`.
    ///
    /// # Panics
    ///
    /// Panics if `index` is greater than `2`.
    #[inline]
    pub fn row(&self, index: usize) -> Vector3<T> {
        Vector3::new(
            self.cols[0][index],
            self.cols[1][index],
            self.cols[2][index],
        )
    }

    /// Returns the transpose, exchanging rows and columns.
    #[inline]
    pub fn transpose(&self) -> Self {
        Self::from_cols(self.row(0), self.row(1), self.row(2))
    }
}

impl<T: Copy + Default> Matrix3<T> {
    /// Constructs a diagonal matrix whose diagonal is `diagonal` and whose
    /// off-diagonal elements are the default (zero) value of `T`.
    #[inline]
    pub fn from_diagonal(diagonal: Vector3<T>) -> Self {
        let zero = T::default();
        Self::from_cols(
            Vector3::new(diagonal.x, zero, zero),
            Vector3::new(zero, diagonal.y, zero),
            Vector3::new(zero, zero, diagonal.z),
        )
    }

    /// Returns the main diagonal as a vector.
    #[inline]
    pub fn diagonal(&self) -> Vector3<T> {
        Vector3::new(self.cols[0].x, self.cols[1].y, self.cols[2].z)
    }
}

impl<T: Default> Default for Matrix3<T> {
    /// Returns the zero matrix.
    #[inline]
    fn default() -> Self {
        Self::from_cols(Vector3::default(), Vector3::default(), Vector3::default())
    }
}

impl<T: Neg<Output = T>> Neg for Matrix3<T> {
    type Output = Self;
    /// Returns the component-wise negation of `self`.
    #[inline]
    fn neg(self) -> Self {
        let [c0, c1, c2] = self.cols;
        Self::from_cols(-c0, -c1, -c2)
    }
}

impl<T: Add<Output = T>> Add for Matrix3<T> {
    type Output = Self;
    /// Returns the component-wise sum of `self` and `rhs`.
    #[inline]
    fn add(self, rhs: Self) -> Self {
        let [a0, a1, a2] = self.cols;
        let [b0, b1, b2] = rhs.cols;
        Self::from_cols(a0 + b0, a1 + b1, a2 + b2)
    }
}

impl<T: AddAssign> AddAssign for Matrix3<T> {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        let [b0, b1, b2] = rhs.cols;
        self.cols[0] += b0;
        self.cols[1] += b1;
        self.cols[2] += b2;
    }
}

impl<T: Sub<Output = T>> Sub for Matrix3<T> {
    type Output = Self;
    /// Returns the component-wise difference of `self` and `rhs`.
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        let [a0, a1, a2] = self.cols;
        let [b0, b1, b2] = rhs.cols;
        Self::from_cols(a0 - b0, a1 - b1, a2 - b2)
    }
}

impl<T: SubAssign> SubAssign for Matrix3<T> {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        let [b0, b1, b2] = rhs.cols;
        self.cols[0] -= b0;
        self.cols[1] -= b1;
        self.cols[2] -= b2;
    }
}

impl<V: Scalar> Matrix3<V> {
    /// The zero matrix.
    pub const ZERO: Self = Self::from_cols(Vector3::ZERO, Vector3::ZERO, Vector3::ZERO);

    /// The identity matrix.
    pub const IDENTITY: Self = Self::from_cols(Vector3::X, Vector3::Y, Vector3::Z);

    /// Constructs a non-uniform scaling matrix from the per-axis factors in
    /// `scale`.
    #[inline]
    pub fn from_scale(scale: Vector3<V>) -> Self {
        Self::from_cols(
            Vector3::new(scale.x, V::ZERO, V::ZERO),
            Vector3::new(V::ZERO, scale.y, V::ZERO),
            Vector3::new(V::ZERO, V::ZERO, scale.z),
        )
    }

    /// Constructs the outer (dyadic) product `a ⊗ b`, the 3×3 matrix whose
    /// `(i, j)` entry is `aᵢ · bⱼ`.
    #[inline]
    pub fn outer_product(a: Vector3<V>, b: Vector3<V>) -> Self {
        Self::from_cols(a * b.x, a * b.y, a * b.z)
    }

    /// Constructs a right-handed rotation of `angle` radians about the `x`
    /// axis.
    #[inline]
    pub fn from_rotation_x(angle: V) -> Self {
        let (sin, cos) = angle.sin_cos();
        Self::from_cols(
            Vector3::new(V::ONE, V::ZERO, V::ZERO),
            Vector3::new(V::ZERO, cos, sin),
            Vector3::new(V::ZERO, -sin, cos),
        )
    }

    /// Constructs a right-handed rotation of `angle` radians about the `y`
    /// axis.
    #[inline]
    pub fn from_rotation_y(angle: V) -> Self {
        let (sin, cos) = angle.sin_cos();
        Self::from_cols(
            Vector3::new(cos, V::ZERO, -sin),
            Vector3::new(V::ZERO, V::ONE, V::ZERO),
            Vector3::new(sin, V::ZERO, cos),
        )
    }

    /// Constructs a right-handed rotation of `angle` radians about the `z`
    /// axis.
    #[inline]
    pub fn from_rotation_z(angle: V) -> Self {
        let (sin, cos) = angle.sin_cos();
        Self::from_cols(
            Vector3::new(cos, sin, V::ZERO),
            Vector3::new(-sin, cos, V::ZERO),
            Vector3::new(V::ZERO, V::ZERO, V::ONE),
        )
    }

    /// Constructs a right-handed rotation of `angle` radians about the unit
    /// vector `axis` (Rodrigues' rotation formula).
    ///
    /// `axis` is assumed to be normalized; a non-unit axis yields a matrix
    /// that also scales.
    #[inline]
    pub fn from_axis_angle(axis: Vector3<V>, angle: V) -> Self {
        let (sin, cos) = angle.sin_cos();
        let t = V::ONE - cos;
        let Vector3 { x, y, z } = axis;
        Self::from_cols(
            Vector3::new(t * x * x + cos, t * x * y + sin * z, t * x * z - sin * y),
            Vector3::new(t * x * y - sin * z, t * y * y + cos, t * y * z + sin * x),
            Vector3::new(t * x * z + sin * y, t * y * z - sin * x, t * z * z + cos),
        )
    }

    /// Returns the trace, the sum of the diagonal elements.
    #[inline]
    pub fn trace(&self) -> V {
        self.cols[0].x + self.cols[1].y + self.cols[2].z
    }

    /// Returns the determinant.
    #[inline]
    pub fn determinant(&self) -> V {
        self.cols[0].dot(self.cols[1].cross(self.cols[2]))
    }

    /// Returns `true` if the matrix has a finite, non-zero determinant and is
    /// therefore invertible.
    #[inline]
    pub fn is_invertible(&self) -> bool {
        let det = self.determinant();
        det != V::ZERO && det.is_finite()
    }

    /// Returns the inverse, or `None` if the matrix is not invertible.
    #[inline]
    pub fn try_inverse(&self) -> Option<Self> {
        let r0 = self.cols[1].cross(self.cols[2]);
        let r1 = self.cols[2].cross(self.cols[0]);
        let r2 = self.cols[0].cross(self.cols[1]);
        let det = self.cols[0].dot(r0);
        if det == V::ZERO || !det.is_finite() {
            return None;
        }
        let inv_det = det.recip();
        Some(Self::from_rows(r0 * inv_det, r1 * inv_det, r2 * inv_det))
    }

    /// Returns the inverse.
    ///
    /// # Panics
    ///
    /// Panics if the matrix is not invertible; use
    /// [`try_inverse`][Self::try_inverse] to handle singular matrices.
    #[inline]
    pub fn inverse(&self) -> Self {
        self.try_inverse().expect("matrix is not invertible")
    }
}

impl<T: Mul<S, Output = T> + Copy, S: Scalar> Mul<S> for Matrix3<T> {
    type Output = Self;
    /// Scales every element by the scalar `rhs`.
    #[inline]
    fn mul(self, rhs: S) -> Self {
        let [c0, c1, c2] = self.cols;
        Self::from_cols(c0 * rhs, c1 * rhs, c2 * rhs)
    }
}

impl<T: MulAssign<S> + Copy, S: Scalar> MulAssign<S> for Matrix3<T> {
    #[inline]
    fn mul_assign(&mut self, rhs: S) {
        self.cols[0] *= rhs;
        self.cols[1] *= rhs;
        self.cols[2] *= rhs;
    }
}

impl<T: Div<S, Output = T> + Copy, S: Scalar> Div<S> for Matrix3<T> {
    type Output = Self;
    /// Divides every element by the scalar `rhs`.
    #[inline]
    fn div(self, rhs: S) -> Self {
        let [c0, c1, c2] = self.cols;
        Self::from_cols(c0 / rhs, c1 / rhs, c2 / rhs)
    }
}

impl<T: DivAssign<S> + Copy, S: Scalar> DivAssign<S> for Matrix3<T> {
    #[inline]
    fn div_assign(&mut self, rhs: S) {
        self.cols[0] /= rhs;
        self.cols[1] /= rhs;
        self.cols[2] /= rhs;
    }
}

impl<V: Scalar> Mul<Vector3<V>> for Matrix3<V> {
    type Output = Vector3<V>;
    /// Applies the matrix to the vector, returning `self * rhs`.
    #[inline]
    fn mul(self, rhs: Vector3<V>) -> Vector3<V> {
        self.cols[0] * rhs.x + self.cols[1] * rhs.y + self.cols[2] * rhs.z
    }
}

impl<V: Scalar> Mul for Matrix3<V> {
    type Output = Self;
    /// Returns the matrix product `self * rhs`.
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        Self::from_cols(self * rhs.cols[0], self * rhs.cols[1], self * rhs.cols[2])
    }
}

impl<V: Scalar> MulAssign for Matrix3<V> {
    #[inline]
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::FRAC_PI_2;

    fn columns() -> Matrix3<f64> {
        Matrix3::from_cols(
            Vector3::new(1.0, 2.0, 3.0),
            Vector3::new(4.0, 5.0, 6.0),
            Vector3::new(7.0, 8.0, 9.0),
        )
    }

    #[test]
    fn from_cols() {
        let m = columns();
        assert_eq!(m.col(0), Vector3::new(1.0, 2.0, 3.0));
        assert_eq!(m.col(1), Vector3::new(4.0, 5.0, 6.0));
        assert_eq!(m.col(2), Vector3::new(7.0, 8.0, 9.0));
    }

    #[test]
    fn from_rows() {
        let m = Matrix3::from_rows(
            Vector3::new(1.0, 2.0, 3.0),
            Vector3::new(4.0, 5.0, 6.0),
            Vector3::new(7.0, 8.0, 9.0),
        );
        assert_eq!(m.row(0), Vector3::new(1.0, 2.0, 3.0));
        assert_eq!(m.row(1), Vector3::new(4.0, 5.0, 6.0));
        assert_eq!(m.row(2), Vector3::new(7.0, 8.0, 9.0));
    }

    #[test]
    fn from_cols_array() {
        let m = Matrix3::from_cols_array(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        assert_eq!(m, columns());
    }

    #[test]
    fn to_cols_array() {
        assert_eq!(
            columns().to_cols_array(),
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
        );
    }

    #[test]
    fn from_diagonal() {
        let m = Matrix3::from_diagonal(Vector3::new(1.0, 2.0, 3.0));
        assert_eq!(
            m.to_cols_array(),
            [1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0]
        );
    }

    #[test]
    fn diagonal() {
        assert_eq!(columns().diagonal(), Vector3::new(1.0, 5.0, 9.0));
    }

    #[test]
    fn map() {
        assert_eq!(
            columns().map(|e| e * 2.0).to_cols_array(),
            [2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0]
        );
    }

    #[test]
    fn zip_map() {
        let sum = columns().zip_map(columns(), |a, b| a + b);
        assert_eq!(sum, columns().map(|e| e * 2.0));
    }

    #[test]
    fn col() {
        assert_eq!(columns().col(2), Vector3::new(7.0, 8.0, 9.0));
    }

    #[test]
    #[should_panic]
    fn col_panics_when_out_of_bounds() {
        let _ = columns().col(3);
    }

    #[test]
    fn row() {
        assert_eq!(columns().row(1), Vector3::new(2.0, 5.0, 8.0));
    }

    #[test]
    fn transpose() {
        assert_eq!(
            columns().transpose().to_cols_array(),
            [1.0, 4.0, 7.0, 2.0, 5.0, 8.0, 3.0, 6.0, 9.0]
        );
        assert_eq!(columns().transpose().transpose(), columns());
    }

    #[test]
    fn default_is_zero() {
        assert_eq!(Matrix3::<f64>::default(), Matrix3::ZERO);
    }

    #[test]
    fn copy_and_clone() {
        let a = columns();
        let b = a;
        let c = ::core::clone::Clone::clone(&a);
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn eq() {
        assert_eq!(columns(), columns());
        assert_ne!(columns(), Matrix3::<f64>::IDENTITY);
    }

    #[test]
    fn debug() {
        assert_eq!(
            format!("{:?}", Matrix3::<f64>::IDENTITY),
            concat!(
                "Matrix3 { ",
                "x_col: Vector3 { x: 1.0, y: 0.0, z: 0.0 }, ",
                "y_col: Vector3 { x: 0.0, y: 1.0, z: 0.0 }, ",
                "z_col: Vector3 { x: 0.0, y: 0.0, z: 1.0 } }"
            )
        );
    }

    #[test]
    fn zero_constant() {
        assert_eq!(Matrix3::<f64>::ZERO.to_cols_array(), [0.0; 9]);
    }

    #[test]
    fn identity_constant() {
        assert_eq!(
            Matrix3::<f64>::IDENTITY.to_cols_array(),
            [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
        );
    }

    #[test]
    fn neg() {
        assert_eq!((-columns()), columns().map(|e| -e));
    }

    #[test]
    fn add() {
        assert_eq!(columns() + columns(), columns().map(|e| e * 2.0));
    }

    #[test]
    fn add_assign() {
        let mut m = columns();
        m += columns();
        assert_eq!(m, columns().map(|e| e * 2.0));
    }

    #[test]
    fn sub() {
        assert_eq!(columns() - columns(), Matrix3::ZERO);
    }

    #[test]
    fn sub_assign() {
        let mut m = columns();
        m -= columns();
        assert_eq!(m, Matrix3::ZERO);
    }

    #[test]
    fn mul_scalar() {
        assert_eq!(columns() * 2.0, columns().map(|e| e * 2.0));
    }

    #[test]
    fn mul_assign_scalar() {
        let mut m = columns();
        m *= 2.0;
        assert_eq!(m, columns().map(|e| e * 2.0));
    }

    #[test]
    fn div_scalar() {
        assert_eq!(columns().map(|e| e * 2.0) / 2.0, columns());
    }

    #[test]
    fn div_assign_scalar() {
        let mut m = columns().map(|e| e * 2.0);
        m /= 2.0;
        assert_eq!(m, columns());
    }

    #[test]
    fn mul_vector() {
        assert_eq!(
            Matrix3::<f64>::IDENTITY * Vector3::new(1.0, 2.0, 3.0),
            Vector3::new(1.0, 2.0, 3.0)
        );
        assert_eq!(
            Matrix3::from_scale(Vector3::new(2.0, 3.0, 4.0)) * Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(2.0, 3.0, 4.0)
        );
    }

    #[test]
    fn mul_matrix() {
        assert_eq!(Matrix3::<f64>::IDENTITY * columns(), columns());
        assert_eq!(
            Matrix3::from_scale(Vector3::new(2.0, 2.0, 2.0))
                * Matrix3::from_scale(Vector3::new(3.0, 3.0, 3.0)),
            Matrix3::from_scale(Vector3::new(6.0, 6.0, 6.0))
        );
    }

    #[test]
    fn mul_assign_matrix() {
        let mut m = columns();
        m *= Matrix3::IDENTITY;
        assert_eq!(m, columns());
    }

    #[test]
    fn from_scale() {
        assert_eq!(
            Matrix3::from_scale(Vector3::new(2.0, 3.0, 4.0)).to_cols_array(),
            [2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 4.0]
        );
    }

    #[test]
    fn outer_product() {
        assert_eq!(
            Matrix3::outer_product(Vector3::new(1.0, 2.0, 3.0), Vector3::new(4.0, 5.0, 6.0)),
            Matrix3::from_cols(
                Vector3::new(4.0, 8.0, 12.0),
                Vector3::new(5.0, 10.0, 15.0),
                Vector3::new(6.0, 12.0, 18.0),
            )
        );
    }

    #[test]
    fn from_rotation_x() {
        let m = Matrix3::from_rotation_x(FRAC_PI_2);
        assert!((m * Vector3::Y - Vector3::Z).norm() < 1e-12);
        assert!((m * Vector3::Z - (-Vector3::<f64>::Y)).norm() < 1e-12);
    }

    #[test]
    fn from_rotation_y() {
        let m = Matrix3::from_rotation_y(FRAC_PI_2);
        assert!((m * Vector3::Z - Vector3::X).norm() < 1e-12);
        assert!((m * Vector3::X - (-Vector3::<f64>::Z)).norm() < 1e-12);
    }

    #[test]
    fn from_rotation_z() {
        let m = Matrix3::from_rotation_z(FRAC_PI_2);
        assert!((m * Vector3::X - Vector3::Y).norm() < 1e-12);
        assert!((m * Vector3::Y - (-Vector3::<f64>::X)).norm() < 1e-12);
    }

    #[test]
    fn from_axis_angle() {
        let a = Matrix3::from_axis_angle(Vector3::Z, 0.7).to_cols_array();
        let b = Matrix3::from_rotation_z(0.7).to_cols_array();
        for i in 0..9 {
            assert!((a[i] - b[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn trace() {
        assert_eq!(columns().trace(), 15.0);
    }

    #[test]
    fn determinant() {
        assert_eq!(Matrix3::<f64>::IDENTITY.determinant(), 1.0);
        assert_eq!(
            Matrix3::from_diagonal(Vector3::new(2.0, 3.0, 4.0)).determinant(),
            24.0
        );
    }

    #[test]
    fn is_invertible() {
        assert!(Matrix3::<f64>::IDENTITY.is_invertible());
        assert!(!Matrix3::from_scale(Vector3::new(1.0, 0.0, 1.0)).is_invertible());
    }

    #[test]
    fn try_inverse() {
        assert_eq!(
            Matrix3::from_scale(Vector3::new(2.0, 4.0, 8.0)).try_inverse(),
            Some(Matrix3::from_scale(Vector3::new(0.5, 0.25, 0.125)))
        );
    }

    #[test]
    fn try_inverse_singular_is_none() {
        assert_eq!(
            Matrix3::from_scale(Vector3::new(1.0, 0.0, 1.0)).try_inverse(),
            None
        );
    }

    #[test]
    fn inverse() {
        assert_eq!(
            Matrix3::from_scale(Vector3::new(2.0, 4.0, 8.0)).inverse(),
            Matrix3::from_scale(Vector3::new(0.5, 0.25, 0.125))
        );
    }

    #[test]
    #[should_panic]
    fn inverse_panics_when_singular() {
        Matrix3::from_scale(Vector3::new(1.0, 0.0, 1.0)).inverse();
    }

    #[test]
    fn inverse_roundtrip() {
        let m = Matrix3::from_cols(
            Vector3::new(2.0, 1.0, 0.0),
            Vector3::new(1.0, 2.0, 1.0),
            Vector3::new(0.0, 1.0, 2.0),
        );
        let product = (m * m.inverse()).to_cols_array();
        let identity = Matrix3::<f64>::IDENTITY.to_cols_array();
        for i in 0..9 {
            assert!((product[i] - identity[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn f32_mul_vector() {
        assert_eq!(
            Matrix3::from_scale(Vector3::<f32>::new(2.0, 3.0, 4.0)) * Vector3::new(1.0, 1.0, 1.0),
            Vector3::new(2.0, 3.0, 4.0)
        );
    }
}
