use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use crate::Scalar;
use crate::tensor::Vector3;

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
